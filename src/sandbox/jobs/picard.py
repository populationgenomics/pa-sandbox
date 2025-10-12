"""
Create Hail Batch jobs to run Picard tools (marking duplicates, QC).
"""

import hailtop.batch as hb
from cpg_flow.resources import (
    STANDARD,
    storage_for_joint_vcf,
)
from cpg_flow.utils import can_reuse, exists
from cpg_utils import Path
from cpg_utils.config import get_config, image_path, reference_path
from cpg_utils.hail_batch import command, fasta_res_group
from hailtop.batch.job import Job


def get_intervals(
    b: hb.Batch,
    scatter_count: int,
    source_intervals_path: Path | None = None,
    exclude_intervals_path: Path | None = None,
    job_attrs: dict[str, str] | None = None,
    output_prefix: Path | None = None,
) -> tuple[Job | None, list[hb.ResourceFile]]:
    """
    Add a job that splits genome/exome intervals into sub-intervals to be used to
    parallelize variant calling.

    @param b: Hail Batch object,
    @param scatter_count: number of target sub-intervals,
    @param source_intervals_path: path to source intervals to split. Would check for
        config if not provided.
    @param exclude_intervals_path: path to file with intervals to exclude.
        Would check for config if not provided.
    @param job_attrs: attributes for Hail Batch job,
    @param output_prefix: path optionally to save split subintervals.

    The job calls picard IntervalListTools to scatter the input interval list
    into scatter_count sub-interval lists, inspired by this WARP task :
    https://github.com/broadinstitute/warp/blob/bc90b0db0138747685b459c83ce52c8576ce03cd/tasks/broad/Utilities.wdl

    Note that we use the mode INTERVAL_SUBDIVISION instead of
    BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW. Modes other than
    INTERVAL_SUBDIVISION produce an unpredictable number of intervals. WDL can
    handle that, but Hail Batch is not dynamic and expects a certain number
    of output files.
    """
    assert scatter_count > 0, scatter_count
    sequencing_type = get_config()['workflow']['sequencing_type']
    source_intervals_path = str(
        source_intervals_path or reference_path(f'broad/{sequencing_type}_calling_interval_lists'),
    )
    exclude_intervals_path = (
        exclude_intervals_path or reference_path('hg38_telomeres_and_centromeres_intervals/interval_list') or None
    )

    if scatter_count == 1:
        # Special case when we don't need to split
        return None, [b.read_input(source_intervals_path)]

    if output_prefix:
        interval_lists_exist = all(
            exists(output_prefix / f'{idx}.interval_list') for idx in range(1, scatter_count + 1)
        )
        if interval_lists_exist:
            return None, [b.read_input(str(output_prefix / f'{idx + 1}.interval_list')) for idx in range(scatter_count)]

    j = b.new_job(
        f'Make {scatter_count} intervals for {sequencing_type}',
        attributes=(job_attrs or {}) | dict(tool='picard IntervalListTools'),
    )
    j.image(image_path('picard'))
    STANDARD.set_resources(j, storage_gb=16, mem_gb=2)

    break_bands_at_multiples_of = {
        'genome': 100000,
        'exome': 0,
    }.get(sequencing_type, 0)

    extra_cmd = ''
    if exclude_intervals_path:
        # If there are intervals to exclude, subtract them from the source intervals
        extra_cmd = f"""-ACTION SUBTRACT \
        -SI {b.read_input(str(exclude_intervals_path))} \
        """

    cmd = f"""
    mkdir $BATCH_TMPDIR/out

    picard -Xms1000m -Xmx1500m \
    IntervalListTools \
    -SCATTER_COUNT {scatter_count} \
    -SUBDIVISION_MODE INTERVAL_SUBDIVISION \
    -UNIQUE true \
    -SORT true \
    -BREAK_BANDS_AT_MULTIPLES_OF {break_bands_at_multiples_of} \
    -I {b.read_input(source_intervals_path)} \
    {extra_cmd} \
    -OUTPUT $BATCH_TMPDIR/out
    ls $BATCH_TMPDIR/out
    ls $BATCH_TMPDIR/out/*
    """
    for idx in range(scatter_count):
        name = f'temp_{str(idx + 1).zfill(4)}_of_{scatter_count}'
        cmd += f"""
        ln $BATCH_TMPDIR/out/{name}/scattered.interval_list {j[f'{idx + 1}.interval_list']}
        """

    j.command(command(cmd))
    if output_prefix:
        for idx in range(scatter_count):
            b.write_output(
                j[f'{idx + 1}.interval_list'],
                str(output_prefix / f'{idx + 1}.interval_list'),
            )

    intervals: list[hb.ResourceFile] = []
    for idx in range(scatter_count):
        interval = j[f'{idx + 1}.interval_list']
        assert isinstance(interval, hb.ResourceFile)
        intervals.append(interval)
    return j, intervals


def vcf_qc(
    b: hb.Batch,
    vcf_or_gvcf: hb.ResourceGroup,
    is_gvcf: bool,
    sequencing_group_count: int | None = None,
    job_attrs: dict | None = None,
    output_summary_path: Path | None = None,
    output_detail_path: Path | None = None,
    overwrite: bool = False,
) -> Job | None:
    """
    Make job that runs Picard CollectVariantCallingMetrics.
    """
    if output_summary_path and can_reuse(output_detail_path, overwrite):
        return None

    job_attrs = (job_attrs or {}) | {'tool': 'picard CollectVariantCallingMetrics'}
    j = b.new_job('CollectVariantCallingMetrics', job_attrs)
    j.image(image_path('picard'))
    storage_gb = 20 if is_gvcf else storage_for_joint_vcf(sequencing_group_count, site_only=False)
    res = STANDARD.set_resources(j, storage_gb=storage_gb, mem_gb=3)
    reference = fasta_res_group(b)
    dbsnp_vcf = b.read_input_group(
        base=reference_path('broad/dbsnp_vcf'),
        index=reference_path('broad/dbsnp_vcf_index'),
    )
    sequencing_type = get_config()['workflow']['sequencing_type']
    intervals_file = b.read_input(reference_path(f'broad/{sequencing_type}_evaluation_interval_lists'))

    input_file = vcf_or_gvcf['g.vcf.gz'] if is_gvcf else vcf_or_gvcf['vcf.gz']

    cmd = f"""\
    picard {res.java_mem_options()} \
    CollectVariantCallingMetrics \
    INPUT={input_file} \
    OUTPUT=$BATCH_TMPDIR/prefix \
    DBSNP={dbsnp_vcf['base']} \
    SEQUENCE_DICTIONARY={reference['dict']} \
    TARGET_INTERVALS={intervals_file} \
    GVCF_INPUT={"true" if is_gvcf else "false"}

    cp $BATCH_TMPDIR/prefix.variant_calling_summary_metrics {j.summary}
    cp $BATCH_TMPDIR/prefix.variant_calling_detail_metrics {j.detail}
    """

    j.command(command(cmd))

    if output_summary_path:
        b.write_output(j.summary, str(output_summary_path))
    if output_detail_path:
        b.write_output(j.detail, str(output_detail_path))
    return j
