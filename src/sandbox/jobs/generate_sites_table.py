import hail as hl
from cpg_flow.targets import Cohort
from cpg_utils.config import config_retrieve, genome_build, get_driver_image
from cpg_utils.hail_batch import get_batch, init_batch, output_path
from hailtop.batch.job import PythonJob


def _initalise_sites_table_job(cohort: Cohort) -> PythonJob:
    job: PythonJob = get_batch().new_python_job(
        name=f'Generate sites table for {cohort.name}',
        attributes=cohort.get_job_attrs() or {} | {'tool': 'Hail:LD_prune'},  # type: ignore[ReportUnknownVariableType]
    )
    job.image(image=get_driver_image())
    return job


def generate_sites_table(cohort: Cohort, chromosome: str) -> PythonJob:
    job: PythonJob = _initalise_sites_table_job(cohort=cohort)
    job.call(_run, cohort, chromosome=chromosome)
    return job


def _run(cohort: Cohort, chromosome: str) -> None:
    vds_path: str = config_retrieve(['generate_sites_table', 'vds_path'])
    exomes: bool = config_retrieve(['generate_sites_table', 'exomes'])
    intersected_bed_file: str = config_retrieve(['generate_sites_table', 'intersected_bed_file'])
    external_sites_filter_table_path: str = config_retrieve(
        ['generate_sites_table', 'external_sites_filter_table_path']
    )
    subsample: bool = config_retrieve(['generate_sites_table', 'subsample'])
    subsample_n: int = config_retrieve(['generate_sites_table', 'subsample_n'])
    sites_table_outpath: str = config_retrieve(['generate_sites_table', 'sites_table_outpath'])

    print(f'Input vds_path: {vds_path}')

    print('Will be writing to pruned_variant_table_path:', sites_table_outpath)

    # Initialise batch
    init_batch(
        worker_memory='highmem',
        driver_memory='highmem',
        driver_cores=4,
    )
    external_sites_table = hl.read_table(external_sites_filter_table_path)

    # LC pipeline VQSR has AS_FilterStatus in info field. We need to annotate
    # sites in the external sites table with AS_FilterStatus if it does not exist
    # based on the `filters` field.
    if 'AS_FilterStatus' not in [k for k in external_sites_table.info.keys()]:
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
        intervals: list[hl.Interval] = capture_interval_ht.interval.collect()

    intervals = [
        interval
        if intervals and interval.start.contig == chromosome
        else hl.Interval(
            hl.eval(hl.locus(chromosome, reference_genome=genome_build())), include_start=True, include_end=True
        )
        for interval in intervals
    ]

    vds = hl.vds.read_vds(str(vds_path), intervals=intervals)

    # Filter to variant sites that pass VQSR
    passed_variants = external_sites_table.filter(external_sites_table.info.AS_FilterStatus == 'PASS')
    vds = hl.vds.filter_variants(vds, passed_variants)

    print('Splitting multi-allelic sites')
    vds = hl.vds.split_multi(vds, filter_changed_loci=True)
    print('Done splitting multi-allelic sites')

    print('Densifying VDS')
    cohort_dense_mt = hl.vds.to_dense_mt(vds)
    print('Done densifying VDS')

    # Run variant QC
    print('Running variant QC')
    # choose variants based off of gnomAD v3 parameters
    # Inbreeding coefficient > -0.80 (no excess of heterozygotes)
    # Must be single nucleotide variants that are autosomal (i.e., no sex), and bi-allelic
    # Have an allele frequency above 1% (note deviation from gnomAD, which is 0.1%)
    # Have a call rate above 99%
    cohort_dense_mt = hl.variant_qc(cohort_dense_mt)
    print('Done running variant QC')

    print('Generating sites table')
    print('Filtering using gnomAD v3 parameters')
    cohort_dense_mt = cohort_dense_mt.annotate_rows(
        IB=hl.agg.inbreeding(
            cohort_dense_mt.GT,
            cohort_dense_mt.variant_qc.AF[1],
        ),
    )
    cohort_dense_mt = cohort_dense_mt.filter_rows(
        (hl.len(cohort_dense_mt.alleles) == 2)
        & (cohort_dense_mt.locus.in_autosome())
        & (cohort_dense_mt.variant_qc.AF[1] > 0.01)
        & (cohort_dense_mt.variant_qc.call_rate > 0.99)
        & (cohort_dense_mt.IB.f_stat > -0.80),
    )
    print('Done filtering using gnomAD v3 parameters')

    # downsize input variants for ld_prune
    # otherwise, persisting the pruned_variant_table will cause
    # script to fail. See https://github.com/populationgenomics/ancestry/pull/79

    # subsample
    if subsample:
        if not subsample_n:
            raise ValueError('If --subsample is set, you must provide a value for --subsample-n')
        print('Sub-sampling sites table before LD pruning')
        nrows = cohort_dense_mt.count_rows()
        print(f'pre sub-sample rows = {nrows}')
        cohort_dense_mt = cohort_dense_mt.sample_rows(
            subsample_n / nrows,
            seed=12345,
        )

    print('Writing sites table pre-LD pruning')
    checkpoint_path = output_path(f'cohort{cohort.name}_dense_mt_{"exome_" if exomes else ""}pre_pruning.mt', 'tmp')
    cohort_dense_mt = cohort_dense_mt.checkpoint(checkpoint_path, overwrite=True)
    print('Done writing sites table pre-LD pruning')

    # as per gnomAD, LD-prune variants with a cutoff of r2 = 0.1
    print('Pruning sites table')
    pruned_variant_table = hl.ld_prune(
        cohort_dense_mt.GT,
        r2=0.1,
        bp_window_size=500000,
    )
    print('Done pruning sites table')

    print('Repartitioning sites table')
    # repartition table after pruning
    pruned_variant_table = pruned_variant_table.repartition(100, shuffle=False)
    print('Number of variants in pruned_variant_table:', pruned_variant_table.count())
    print('Done repartitioning sites table')

    print(f'Writing sites table to {sites_table_outpath}')
    pruned_variant_table.write(sites_table_outpath, overwrite=True)
    print('Done writing sites table')
