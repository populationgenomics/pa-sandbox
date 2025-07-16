from typing import TYPE_CHECKING

import hail as hl
from cpg_flow.targets import Cohort
from cpg_flow.utils import to_path  # type: ignore[ReportUnknownVariableType]
from cpg_utils.config import config_retrieve, genome_build, get_driver_image
from cpg_utils.hail_batch import get_batch, init_batch, output_path  # type: ignore[ReportUnknownVariableType]
from hailtop.batch.job import PythonJob, PythonResult
from loguru import logger

from ..gnomad_methods.bi_allelic_sites_inbreeding import bi_allelic_site_inbreeding_expr

if TYPE_CHECKING:
    from hail.vds.variant_dataset import VariantDataset


def _initalise_sites_table_job(cohort: Cohort, name: str, job_memory: str, job_cpus: int) -> PythonJob:
    job: PythonJob = get_batch().new_python_job(
        name=name,
        attributes=cohort.get_job_attrs() or {} | {'tool': 'Hail:LD_prune'},  # type: ignore[ReportUnknownVariableType]
    )
    job.image(image=get_driver_image())
    job.memory(job_memory)
    job.cpu(job_cpus)
    job.n_max_attempts(2)
    return job


def _initalise_sites_table_merge_job(cohort: Cohort, job_memory: str, job_cpus: int) -> PythonJob:
    job: PythonJob = get_batch().new_python_job(
        name=f'Merging per chromosome sites tables for {cohort.name}',
        attributes=cohort.get_job_attrs() or {} | {'tool': 'Hail:MergeSitesTables'},  # type: ignore[ReportUnknownVariableType]
    )
    job.image(image=get_driver_image())
    job.memory(job_memory)
    job.cpu(job_cpus)
    job.n_max_attempts(2)
    return job

def get_filtering_intervals(chromosome: str) -> list[hl.Interval]:

    exomes: bool = config_retrieve(['generate_sites_table', 'exomes'])
    intersected_bed_file: str = config_retrieve(['generate_sites_table', 'intersected_bed_file'])

    # If 'autosomes' is passed, keep intervals from all 22 chromosomes.
    if chromosome == 'autosomes':
        target = {f'chr{i}' for i in range(1, 22)}
    else:
        target = {chromosome}

    # Optionally filter to the exome regions.
    if exomes:
        if not intersected_bed_file:
            raise ValueError('If --exomes is set, you must provide at least one --capture-region-bed-files')

        # Read capture regions and filter to target chromosomes in one step
        capture_interval_ht = hl.import_bed(
            str(intersected_bed_file),
            reference_genome=genome_build()
        )

        # Filter intervals to target chromosomes and collect
        intervals = capture_interval_ht.filter(
            hl.literal(target).contains(capture_interval_ht.interval.start.contig)
        ).interval.collect()

        if not intervals:
            # No capture regions found for target chromosomes
            raise ValueError(f"No capture regions found for chromosomes: {target}")

    # Otherwise just subset to chromosomes.
    else:
        intervals = [
            hl.eval(hl.parse_locus_interval(chrom, reference_genome=genome_build()))
            for chrom in target
        ]

    return intervals

def generate_sites_table(cohort: Cohort, sites_table_outpath: str) -> PythonJob:
    job_memory: str = config_retrieve(['workflow', 'job_memory'])
    job_cpus: int = config_retrieve(['workflow', 'job_cpus'])
    cohort_name: str = cohort.name
    sites_jobs: list[PythonResult] = []
    chromosomes: list[str] = config_retrieve(['generate_sites_table', 'chromosome_list'])
    for chromosome in chromosomes:
        sites_jobs.append(
            _initalise_sites_table_job(
                cohort=cohort,
                name=f'Generate sites table for {chromosome} with {cohort.name}',
                job_memory=job_memory,
                job_cpus=job_cpus,
            ).call(
                _run_sites_per_chromosome,
                cohort_name=cohort_name,
                chromosome=chromosome,
            )
        )
    merge_job: PythonJob = _initalise_sites_table_merge_job(cohort=cohort, job_memory=job_memory, job_cpus=job_cpus)
    merge_job.call(
        _run_merge_sites_table, filtered_chromosome_tables=sites_jobs, sites_table_outpath=sites_table_outpath
    )
    return merge_job


def _run_sites_per_chromosome(cohort_name: str, chromosome: str) -> str:

    # Paths to the input VDS or dense MatrixTable.
    vds_path: str = config_retrieve(['generate_sites_table', 'vds_path'], None)
    dense_mt_path: str = config_retrieve(['generate_sites_table', 'dense_mt_path'], None)

    if not vds_path and not dense_mt_path:
            raise ValueError('One of vds_path or dense_mt_path must be provided')

    # Optional filtering to perform before variant QC (samples, variants, and intervals).
    exomes: bool = config_retrieve(['generate_sites_table', 'exomes'])
    external_sites_filter_table_path: str = config_retrieve(
        ['generate_sites_table', 'external_sites_filter_table_path'], None
    )
    subsample: bool = config_retrieve(['generate_sites_table', 'subsample'])
    subsample_n: int = config_retrieve(['generate_sites_table', 'subsample_n'])

    samples_to_drop = config_retrieve(['generate_sites_table', 'samples_to_drop'], None)

    # How many partitions to create after filtering.
    n_partitions = config_retrieve(['generate_sites_table', 'n_partitions'])

    # Variant QC thresholds.
    allele_frequency_min: float = config_retrieve(['generate_sites_table', 'allele_frequency_min'])
    call_rate_min: float = config_retrieve(['generate_sites_table', 'call_rate_min'])
    f_stat: float = config_retrieve(['generate_sites_table', 'f_stat'])
    p_value_hwe: float = config_retrieve(['generate_sites_table', 'p_value_hwe'])

    # LD pruning parameters.
    r2_value: float = config_retrieve(['generate_sites_table', 'r2_value'])
    bp_window_size: int = config_retrieve(['generate_sites_table', 'bp_window_size'])

    init_batch()

    # Determine paths for intermediate outputs.
    pre_ld_prune_path: str = output_path(
        f'cohort{cohort_name}_{chromosome}_dense_mt_{"exome_" if exomes else ""}pre_pruning.mt', 'tmp'
    )
    post_ld_prune_outpath: str = output_path(
        f'cohort{cohort_name}_{chromosome}_dense_mt_{"exome_" if exomes else ""}pruned.mt', 'tmp'
    )

    if not to_path(post_ld_prune_outpath).exists():

        if not to_path(pre_ld_prune_path).exists():

            # Fetch the filtering intervals.
            filtering_intervals: list[hl.Interval] = get_filtering_intervals(chromosome)

            # Fetch any additional site filtering tables.
            if external_sites_filter_table_path:
                external_sites_table: hl.Table = hl.read_table(external_sites_filter_table_path)

                # LC pipeline VQSR has AS_FilterStatus in info field. We need to annotate
                # sites in the external sites table with AS_FilterStatus if it does not exist
                # based on the `filters` field.
                if 'AS_FilterStatus' not in list(external_sites_table.info.keys()):
                    # if 'AS_FilterStatus' not in [k for k in external_sites_table.info.keys()]:
                    external_sites_table = external_sites_table.annotate(
                        info=external_sites_table.info.annotate(
                            AS_FilterStatus=hl.if_else(hl.len(external_sites_table.filters) == 0, 'PASS', 'FAIL'),
                        ),
                    )
                passed_variants = external_sites_table.filter(external_sites_table.info.AS_FilterStatus == 'PASS')

            # Fetch any samples that should be dropped.
            if samples_to_drop:
                sample_hts = [hl.read_table(path) for path in samples_to_drop]
                all_samples_to_drop = sample_hts[0]
                for ht in sample_hts[1:]:
                    all_samples_to_drop = all_samples_to_drop.union(ht)

            if not dense_mt_path:

                # Read VDS then filter to intervals, to avoid spanning ref blocks being dropped silently.
                vds: VariantDataset = hl.vds.read_vds(str(vds_path))
                vds = hl.vds.filter_intervals(vds, filtering_intervals, split_reference_blocks=False)

                # Filter to variant sites that pass QC.
                if external_sites_filter_table_path:
                    vds = hl.vds.filter_variants(vds, passed_variants)

                # Remove all multiallelic sites prior to densification.
                vds = hl.vds.filter_variants(
                    vds,
                    vds.variant_data.filter_rows(hl.len(vds.variant_data.alleles) == 2).rows(),
                    keep=True,
                )

                # Remove samples that are present in the samples_to_drop list.
                if samples_to_drop:
                    vds = hl.vds.filter_samples(vds, all_samples_to_drop, keep=False)

                # Recode the local genotypes as GT entries if required.
                if 'GT' not in vds.variant_data.entry:
                    vds.variant_data = vds.variant_data.annotate_entries(
                        GT=hl.vds.lgt_to_gt(vds.variant_data.LGT, vds.variant_data.LA)
                    )

                logger.info('Densifying VDS')
                cohort_dense_mt: hl.MatrixTable = hl.vds.to_dense_mt(vds)

            else:

                # Read the dense MatrixTable.
                cohort_dense_mt = hl.read_matrix_table(dense_mt_path)

                # Filter to intervals.
                cohort_dense_mt = hl.filter_intervals(cohort_dense_mt, filtering_intervals)

                # Filter to variant sites that pass QC.
                if external_sites_filter_table_path:
                    cohort_dense_mt = cohort_dense_mt.semi_join_rows(passed_variants)

                # Remove samples that are present in the samples_to_drop list.
                if samples_to_drop:
                    cohort_dense_mt = cohort_dense_mt.anti_join_cols(all_samples_to_drop)

            # Run variant QC
            logger.info('Running variant QC')
            cohort_dense_mt = hl.variant_qc(cohort_dense_mt)
            cohort_dense_mt = cohort_dense_mt.annotate_rows(
                IB=bi_allelic_site_inbreeding_expr(
                    cohort_dense_mt.GT,
                ),
            )

            logger.info('Performing variant filtering')
            cohort_dense_mt = cohort_dense_mt.filter_rows(
                (hl.is_snp(cohort_dense_mt.alleles[0], cohort_dense_mt.alleles[1]))
                & (cohort_dense_mt.locus.in_autosome())
                & (hl.min(cohort_dense_mt.variant_qc.AF[0], cohort_dense_mt.variant_qc.AF[1]) > allele_frequency_min)
                & (cohort_dense_mt.variant_qc.call_rate > call_rate_min)
                & (f_stat < cohort_dense_mt.IB)
                & (cohort_dense_mt.variant_qc.p_value_hwe > p_value_hwe)
            )
            logger.info('Filtering complete')

            # subsample
            if subsample:
                if not subsample_n:
                    raise ValueError('If --subsample is set, you must provide a value for --subsample-n')
                logger.info('Sub-sampling sites table before LD pruning')
                nrows = cohort_dense_mt.count_rows()
                logger.info(f'There are {nrows} before sub-sampling')
                cohort_dense_mt = cohort_dense_mt.sample_rows(
                    subsample_n / nrows,
                    seed=12345,
                )

            # Repartition the matrix table to account for variant filtering.
            logger.info('Repartioning the sites table before LD pruning')
            cohort_dense_mt = cohort_dense_mt.repartition(n_partitions)

            logger.info('Checkpointing pre-LD pruning matrix table')
            cohort_dense_mt = cohort_dense_mt.checkpoint(pre_ld_prune_path, overwrite=True)
            logger.info('Checkpointing of the pre-LD pruning matrix table is complete')

        else:
            cohort_dense_mt = hl.read_matrix_table(pre_ld_prune_path)

        # as per gnomAD, LD-prune variants with a cutoff of r2 = 0.1
        logger.info('Pruning sites table')
        pruned_variant_table = hl.ld_prune(
            cohort_dense_mt.GT,
            r2=r2_value,
            bp_window_size=bp_window_size,
        )

        logger.info(
            f'Pruning complete. Number of variants in pruned_variant_table: {pruned_variant_table.count()}'
        )

        post_ld_prune_outpath = output_path(
            f'cohort{cohort_name}_{chromosome}_dense_mt_{"exome_" if exomes else ""}pruned.mt', 'tmp'
        )
        pruned_variant_table.write(post_ld_prune_outpath, overwrite=True)
        logger.info('Done writing sites table')

    return post_ld_prune_outpath


def _run_merge_sites_table(filtered_chromosome_tables: list[str], sites_table_outpath: str) -> None:
    logger.info('Merging per chromosome sites tables into one')
    init_batch()
    merged_sites_tables_list: list[hl.Table] = [hl.read_table(table) for table in filtered_chromosome_tables]
    merged_sites_table: hl.Table = hl.Table.union(*merged_sites_tables_list)
    merged_sites_table = merged_sites_table.repartition(config_retrieve(['generate_sites_table', 'n_partitions']))
    merged_sites_table.write(sites_table_outpath, overwrite=True)
