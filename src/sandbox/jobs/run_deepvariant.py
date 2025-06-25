"""
A job should contain the logic for a single Stage
"""

from typing import TYPE_CHECKING

from cpg_flow.filetypes import CramPath
from cpg_utils.config import image_path
from cpg_utils.hail_batch import fasta_res_group, get_batch

if TYPE_CHECKING:
    from hailtop.batch.job import Job


def run(output_vcf: str,
        output_gvcf: str,
        sequencing_group_name: str,
        cram_path: CramPath) -> 'Job':
    """
    This is a simple example of a job that writes a statement to a file.

    Args:
        statement (str): the intended file contents
        output_file (str): the path to write the file to

    Returns:
        the resulting job
    """
    b = get_batch()
    # create a job
    j = b.new_job('Pangenome Aware DeepVariant')

    # choose an image to run this job in (default is bare ubuntu)
    j.image(image_path('deepvariant_pangenome_aware'))
    j.memory('96Gi')
    j.storage('50Gi')
    j.declare_resource_group(
        output_gvcf={
            'g.vcf.gz': '{root}-' + sequencing_group_name + '.g.vcf.gz',
            'g.vcf.gz.tbi': '{root}-' + sequencing_group_name + '.g.vcf.gz.tbi',
        },
    )
    j.declare_resource_group(
        output_vcf={
            'vcf.gz': '{root}-' + sequencing_group_name + '.vcf.gz',
            'vcf.gz.tbi': '{root}-' + sequencing_group_name + '.vcf.gz.tbi',
        },
    )

    reference = fasta_res_group(b)
    # copy test data
    j.command(
        f"""
        set -ex
        CRAM=$BATCH_TMPDIR/{sequencing_group_name}.cram
        CRAI=$BATCH_TMPDIR/{sequencing_group_name}.cram.crai

        # Retrying copying to avoid google bandwidth limits
        retry_gs_cp {cram_path.path} $CRAM
        retry_gs_cp {cram_path.index_path} $CRAI
        ls -l $CRAM
        ls -l $CRAI

        ls -l ${{BATCH_TMPDIR}}

        HTTPDIR=https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38
        curl -L "${{HTTPDIR}}/hprc-v1.1-mc-grch38.gbz" -o ${{BATCH_TMPDIR}}/hprc-v1.1-mc-grch38.gbz

        # Run pangenome-aware DeepVariant
        mkdir -p ${{BATCH_TMPDIR}}/intermediate_results_dir
        pwd
        ls -l
        ls -l /opt
        ls -l /opt/deepvariant/
        ls -l /opt/deepvariant/bin/
        /opt/deepvariant/bin/run_pangenome_aware_deepvariant \
        --model_type WGS \
        --ref {reference.base} \
        --reads $BATCH_TMPDIR/{sequencing_group_name}.cram  \
        --pangenome ${{BATCH_TMPDIR}}/hprc-v1.1-mc-grch38.gbz \
        --output_vcf {j.output_vcf['vcf.gz']} \
        --output_gvcf {j.output_gvcf['g.vcf.gz']} \
        --num_shards $(nproc) \
        --intermediate_results_dir ${{BATCH_TMPDIR}}/intermediate_results_dir
        """,
    )

    # write the output to the expected location
    get_batch().write_output(j.output_vcf, output_vcf.replace('.vcf.gz', ''))
    get_batch().write_output(j.output_gvcf, output_gvcf.replace('.g.vcf.gz', ''))

    # return the job
    return j
