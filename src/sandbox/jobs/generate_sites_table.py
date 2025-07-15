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


def _run_sites_per_chromosome(cohort_name: str, chromosome: str) -> str:  # noqa: PLR0915
    vds_path: str = config_retrieve(['generate_sites_table', 'vds_path'])
    exomes: bool = config_retrieve(['generate_sites_table', 'exomes'])
    intersected_bed_file: str = config_retrieve(['generate_sites_table', 'intersected_bed_file'])
    external_sites_filter_table_path: str = config_retrieve(
        ['generate_sites_table', 'external_sites_filter_table_path']
    )
    subsample: bool = config_retrieve(['generate_sites_table', 'subsample'])
    subsample_n: int = config_retrieve(['generate_sites_table', 'subsample_n'])

    samples_to_drop = config_retrieve(['generate_sites_table', 'samples_to_drop'])

    n_partitions = config_retrieve(['generate_sites_table', 'n_partitions'])

    allele_frequency_min: float = config_retrieve(['generate_sites_table', 'allele_frequency_min'])
    call_rate_min: float = config_retrieve(['generate_sites_table', 'call_rate_min'])
    f_stat: float = config_retrieve(['generate_sites_table', 'f_stat'])
    p_value_hwe: float = config_retrieve(['generate_sites_table', 'p_value_hwe'])

    r2_value: float = config_retrieve(['generate_sites_table', 'r2_value'])
    bp_window_size: int = config_retrieve(['generate_sites_table', 'bp_window_size'])

    intervals: list[hl.Interval] = []

    init_batch()

    pre_ld_prune_path: str = output_path(
        f'cohort{cohort_name}_{chromosome}_dense_mt_{"exome_" if exomes else ""}pre_pruning.mt', 'tmp'
    )
    post_ld_prune_outpath: str = output_path(
        f'cohort{cohort_name}_{chromosome}_dense_mt_{"exome_" if exomes else ""}pruned.mt', 'tmp'
    )

    if not to_path(post_ld_prune_outpath).exists():
        if not to_path(pre_ld_prune_path).exists():
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

            # exomes
            if exomes:
                if not intersected_bed_file:
                    raise ValueError('If --exomes is set, you must provide at least one --capture-region-bed-files')
                # Read in capture region bed files
                capture_interval_ht: hl.Table = hl.import_bed(
                    str(intersected_bed_file), reference_genome=genome_build()
                )
                # Generate list of intervals
                intervals = capture_interval_ht.interval.collect()

            tmp_intervals: list[hl.Interval] = []
            if not intervals:
                tmp_intervals = [hl.eval(hl.parse_locus_interval(chromosome, reference_genome=genome_build()))]
            else:
                for interval in intervals:
                    if interval.start.contig == chromosome:
                        tmp_intervals.append(interval)

            filtering_intervals: list[hl.Interval] = hl.eval(hl.array(tmp_intervals))

            # Read VDS then filter, to avoid ref blocks that span intervals being dropped silently
            vds: VariantDataset = hl.vds.read_vds(str(vds_path))
            vds = hl.vds.filter_intervals(vds, filtering_intervals, split_reference_blocks=False)

            # Filter to variant sites that pass VQSR
            passed_variants = external_sites_table.filter(external_sites_table.info.AS_FilterStatus == 'PASS')
            vds = hl.vds.filter_variants(vds, passed_variants)

            # Remove all multiallelic sites
            vds = hl.vds.filter_variants(
                vds,
                vds.variant_data.filter_rows(hl.len(vds.variant_data.alleles) == 2).rows(),  # noqa: PLR2004
                keep=True,
            )

            # Remove samples that are present in the samples_to_drop list
            sample_hts = [hl.read_table(path) for path in samples_to_drop]
            all_samples_to_drop = sample_hts[0]
            for ht in sample_hts[1:]:
                all_samples_to_drop = all_samples_to_drop.union(ht)
            vds = hl.vds.filter_samples(vds, all_samples_to_drop, keep=False)

            if 'GT' not in vds.variant_data.entry:
                vds.variant_data = vds.variant_data.annotate_entries(
                    GT=hl.vds.lgt_to_gt(vds.variant_data.LGT, vds.variant_data.LA)
                )

            logger.info('Densifying VDS')
            cohort_dense_mt: hl.MatrixTable = hl.vds.to_dense_mt(vds)
            logger.info('Done densifying VDS. Now running variant QC')

            # Run variant QC
            # choose variants based off of gnomAD v3 parameters
            # Inbreeding coefficient > -0.80 (no excess of heterozygotes)
            # Must be single nucleotide variants that are autosomal (i.e., no sex), and bi-allelic
            # Have an allele frequency above 1% (note deviation from gnomAD, which is 0.1%)
            # Have a call rate above 99%
            cohort_dense_mt = hl.variant_qc(cohort_dense_mt)

            logger.info('Done running variant QC. Now generating sites table and filtering using gnomAD v3 parameters')
            cohort_dense_mt = cohort_dense_mt.annotate_rows(
                IB=bi_allelic_site_inbreeding_expr(
                    cohort_dense_mt.GT,
                ),
            )
            cohort_dense_mt = cohort_dense_mt.filter_rows(
                (hl.is_snp(cohort_dense_mt.alleles[0], cohort_dense_mt.alleles[1]))
                & (cohort_dense_mt.locus.in_autosome())
                & (hl.min(cohort_dense_mt.variant_qc.AF[0], cohort_dense_mt.variant_qc.AF[1]) > allele_frequency_min)
                & (cohort_dense_mt.variant_qc.call_rate > call_rate_min)
                & (f_stat < cohort_dense_mt.IB)
                & (cohort_dense_mt.variant_qc.p_value_hwe > p_value_hwe)
            )
            logger.info('Done filtering using gnomAD v3 parameters')

            # downsize input variants for ld_prune
            # otherwise, persisting the pruned_variant_table will cause
            # script to fail. See https://github.com/populationgenomics/ancestry/pull/79

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

            # Repartition the matrix table to account for our aggressive variant filtering.
            logger.info('Repartioning the sites table pre-LD pruning')
            cohort_dense_mt = cohort_dense_mt.repartition(n_partitions)

            logger.info('Writing sites table pre-LD pruning')
            cohort_dense_mt = cohort_dense_mt.checkpoint(pre_ld_prune_path, overwrite=True)
            logger.info('Done writing sites table pre-LD pruning')

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
            f'Done pruning sites table. Number of variants in pruned_variant_table: {pruned_variant_table.count()}'
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
