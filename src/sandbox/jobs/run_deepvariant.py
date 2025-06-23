"""
A job should contain the logic for a single Stage
"""

from typing import TYPE_CHECKING

from cpg_utils.config import image_path
from cpg_utils.hail_batch import get_batch

if TYPE_CHECKING:
    from hailtop.batch.job import Job


def run(output_vcf: str, output_gvcf: str) -> 'Job':
    """
    This is a simple example of a job that writes a statement to a file.

    Args:
        statement (str): the intended file contents
        output_file (str): the path to write the file to

    Returns:
        the resulting job
    """

    # create a job
    j = get_batch().new_job('DeepVariant')

    # choose an image to run this job in (default is bare ubuntu)
    j.image(image_path(['deepvariant']))

    # copy test data
    j.command(
         f"""
        INPUT_DIR="$BATCH_TMPDIR/quickstart-testdata"
        CHECKPOINT_DIR="$BATCH_TMPDIR/checkpoint"
        DATA_HTTP_DIR="https://storage.googleapis.com/deepvariant/quickstart-testdata"

        mkdir -p ${{INPUT_DIR}}
        mkdir -p ${{CHECKPOINT_DIR}}
        wget -P ${{INPUT_DIR}} "${{DATA_HTTP_DIR}}/NA12878_S1.chr20.10_10p1mb.bam"
        wget -P ${{INPUT_DIR}} "${{DATA_HTTP_DIR}}/NA12878_S1.chr20.10_10p1mb.bam.bai"
        wget -P ${{INPUT_DIR}} "${{DATA_HTTP_DIR}}/test_nist.b37_chr20_100kbp_at_10mb.bed"
        wget -P ${{INPUT_DIR}} "${{DATA_HTTP_DIR}}/test_nist.b37_chr20_100kbp_at_10mb.vcf.gz"
        wget -P ${{INPUT_DIR}} "${{DATA_HTTP_DIR}}/test_nist.b37_chr20_100kbp_at_10mb.vcf.gz.tbi"
        wget -P ${{INPUT_DIR}} "${{DATA_HTTP_DIR}}/ucsc.hg19.chr20.unittest.fasta"
        wget -P ${{INPUT_DIR}} "${{DATA_HTTP_DIR}}/ucsc.hg19.chr20.unittest.fasta.fai"
        wget -P ${{INPUT_DIR}} "${{DATA_HTTP_DIR}}/ucsc.hg19.chr20.unittest.fasta.gz"
        wget -P ${{INPUT_DIR}} "${{DATA_HTTP_DIR}}/ucsc.hg19.chr20.unittest.fasta.gz.fai"
        wget -P ${{INPUT_DIR}} "${{DATA_HTTP_DIR}}/ucsc.hg19.chr20.unittest.fasta.gz.gzi"

        # Run make_examples
        /opt/deepvariant/bin/run_deepvariant \
        --model_type=WGS \
        --vcf_stats_report=true \
        --ref=${{INPUT_DIR}}/ucsc.hg19.chr20.unittest.fasta \
        --reads=${{INPUT_DIR}}/NA12878_S1.chr20.10_10p1mb.bam \
        --regions "chr20:10,000,000-10,010,000" \
        --output_vcf={j.outVcf} \
        --output_gvcf={j.outGvcf} \
        --intermediate_results_dir /output/intermediate_results_dir \
        --num_shards=1
        """,
    )

    # write the output to the expected location
    get_batch().write_output(j.outVcf, output_vcf)
    get_batch().write_output(j.outGvcf, output_gvcf)

    # return the job
    return j
