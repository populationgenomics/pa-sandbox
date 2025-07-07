from typing import TYPE_CHECKING

import hail as hl
from cpg_flow.targets import Cohort
from cpg_utils.config import config_retrieve, genome_build, get_driver_image
from cpg_utils.hail_batch import get_batch, output_path
from hailtop.batch.job import PythonJob, PythonResult
from loguru import logger

if TYPE_CHECKING:
    from hail.vds.variant_dataset import VariantDataset


def _initalise_sites_table_job(cohort: Cohort, name: str) -> PythonJob:
    job: PythonJob = get_batch().new_python_job(
        name=name,
        attributes=cohort.get_job_attrs() or {} | {'tool': 'Hail:LD_prune'},  # type: ignore[ReportUnknownVariableType]
    )
    job.image(image=get_driver_image())
    job.memory('highmem')
    job.cpu(4)
    return job


def _initalise_sites_table_merge_job(cohort: Cohort) -> PythonJob:
    job: PythonJob = get_batch().new_python_job(
        name=f'Merging per chromosome sites tables for {cohort.name}',
        attributes=cohort.get_job_attrs() or {} | {'tool': 'Hail:MergeSitesTables'},  # type: ignore[ReportUnknownVariableType]
    )
    job.image(image=get_driver_image())
    job.memory('highmem')
    job.cpu(4)
    return job


def generate_sites_table(cohort: Cohort, sites_table_outpath: str) -> PythonJob:
    cohort_name: str = cohort.name
    sites_jobs: list[PythonResult] = []
    chromosomes: list[str] = [f'chr{x}' for x in [*list(range(1, 23)), 'X', 'Y', 'M']]
    for chromosome in chromosomes:
        sites_jobs.append(
            _initalise_sites_table_job(
                cohort=cohort, name=f'Generate sites table for chr{chromosome} with {cohort.name}'
            ).call(
                _run_sites_per_chromosome,
                cohort_name=cohort_name,
                chromosome=chromosome,
            )
        )
    merge_job: PythonJob = _initalise_sites_table_merge_job(cohort=cohort)
    merge_job.call(
        _run_merge_sites_table, filtered_chromosome_tables=sites_jobs, sites_table_outpath=sites_table_outpath
    )
    return merge_job


def _run_sites_per_chromosome(cohort_name: str, chromosome: str) -> str:
    vds_path: str = config_retrieve(['generate_sites_table', 'vds_path'])
    exomes: bool = config_retrieve(['generate_sites_table', 'exomes'])
    intersected_bed_file: str = config_retrieve(['generate_sites_table', 'intersected_bed_file'])
    external_sites_filter_table_path: str = config_retrieve(
        ['generate_sites_table', 'external_sites_filter_table_path']
    )
    subsample: bool = config_retrieve(['generate_sites_table', 'subsample'])
    subsample_n: int = config_retrieve(['generate_sites_table', 'subsample_n'])

    intervals: list[hl.Interval] = []
    external_sites_table = hl.read_table(external_sites_filter_table_path)

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
        capture_interval_ht: hl.Table = hl.import_bed(str(intersected_bed_file), reference_genome=genome_build())
        # Generate list of intervals
        intervals = capture_interval_ht.interval.collect()

    intervals = [
        interval
        if intervals and interval.start.contig == chromosome
        else hl.Interval(
            hl.eval(hl.locus(chromosome, reference_genome=genome_build())), include_start=True, include_end=True
        )
        for interval in intervals
    ]

    # Read VDS then filter, to avoid ref blocks that span intervals being dropped silently
    vds: VariantDataset = hl.vds.read_vds(str(vds_path))
    vds = hl.vds.filter_intervals(vds, intervals, split_reference_blocks=False)

    # Filter to variant sites that pass VQSR
    passed_variants = external_sites_table.filter(external_sites_table.info.AS_FilterStatus == 'PASS')
    vds = hl.vds.filter_variants(vds, passed_variants)

    # Remove all multiallelic sites
    vds = hl.vds.filter_variants(vds, hl.len(vds.variant_data.alleles) > 2)  # noqa: PLR2004

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
        IB=hl.agg.inbreeding(
            cohort_dense_mt.GT,
            cohort_dense_mt.variant_qc.AF[1],
        ),
    )
    cohort_dense_mt = cohort_dense_mt.filter_rows(
        (hl.is_snp(cohort_dense_mt.alleles[0], cohort_dense_mt.alleles[1]))
        & (cohort_dense_mt.locus.in_autosome())
        & (cohort_dense_mt.variant_qc.AF[1] > 0.01)  # noqa: PLR2004
        & (cohort_dense_mt.variant_qc.call_rate > 0.99)  # noqa: PLR2004
        & (cohort_dense_mt.IB.f_stat > -0.80)  # noqa: PLR2004
        & (cohort_dense_mt.variant_qc.p_value_hwe > 1e-8)  # noqa: PLR2004
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

    logger.info('Writing sites table pre-LD pruning')
    checkpoint_path: str = output_path(
        f'cohort{cohort_name}_{chromosome}_dense_mt_{"exome_" if exomes else ""}pre_pruning.mt', 'tmp'
    )
    cohort_dense_mt = cohort_dense_mt.checkpoint(checkpoint_path, overwrite=True)
    logger.info('Done writing sites table pre-LD pruning')

    # as per gnomAD, LD-prune variants with a cutoff of r2 = 0.1
    logger.info('Pruning sites table')
    pruned_variant_table = hl.ld_prune(
        cohort_dense_mt.GT,
        r2=0.1,
        bp_window_size=500000,
    )

    logger.info(f'Done pruning sites table. Number of variants in pruned_variant_table: {pruned_variant_table.count()}')

    tmp_chr_outpath: str = output_path(
        f'cohort{cohort_name}_{chromosome}_dense_mt_{"exome_" if exomes else ""}pruned.mt', 'tmp'
    )
    pruned_variant_table.write(tmp_chr_outpath, overwrite=True)
    logger.info('Done writing sites table')
    return tmp_chr_outpath


def _run_merge_sites_table(filtered_chromosome_tables: list[str], sites_table_outpath: str) -> None:
    logger.info('Merging per chromosome sites tables into one')
    merged_sites_tables_list: list[hl.MatrixTable] = [
        hl.read_table(table).to_matrix_table() for table in filtered_chromosome_tables
    ]
    merged_sites_table: hl.Table = hl.MatrixTable.union_rows(*merged_sites_tables_list).rows()
    merged_sites_table = merged_sites_table.repartition(100, shuffle=True)
    merged_sites_table.write(sites_table_outpath, overwrite=True)
