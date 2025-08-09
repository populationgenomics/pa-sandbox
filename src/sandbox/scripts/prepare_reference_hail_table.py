"""Create a reference hail table containing every locus in the genome, excluding centromeres and telomeres.

This script takes a reference hail table, selects distinct loci only, removes centromeres and telomeres,
and repartitions it.

Arguments:
    --n-partitions: int. Number of partitions for the resulting reference hail table.
    --ref-table-outpath: str. The output location for the filtered reference hail table.
    --chromosome: str, optional. An optional chromosome to subset to.

Example usage:

analysis-runner --dataset agdd \
    --access-level full \
    --output-dir generate-reference-table \
    --description 'Generate a reference hail table' \
    python3 prepare_reference_hail_table.py \
    --n-partitions 5000 \
    --ref-table-outpath gs://cpg-common-main/references/large_cohort/grch38_no_centromeres_or_telomeres.ht
"""

import logging

import hail as hl

from argparse import ArgumentParser, Namespace
from cpg_utils.config import reference_path
from cpg_utils.hail_batch import init_batch

logging.basicConfig(
    format='%(asctime)s (%(name)s %(lineno)s): %(message)s',
    datefmt='%m/%d/%Y %I:%M:%S %p',
)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

def main(
    n_partitions: int,
    ref_table_outpath: str,
) -> None:
    init_batch()

    # Load the full seqr reference data table.
    ref_ht = hl.read_table(reference_path('seqr_combined_reference_data'))

    # Retain only the locus field, and subset to unique values.
    logger.info('Rekeying reference HT to distinct loci.')
    ref_ht = ref_ht.key_by('locus').select().distinct()

    # Filter out telomeres and centromeres.
    logger.info('Removing centromeres and telomeres.')
    tel_cent_ht = hl.read_table(reference_path('gnomad/tel_and_cent_ht'))
    ref_ht = hl.filter_intervals(
        ref_ht,
        tel_cent_ht.interval.collect(),
        keep=False,
    )

    # Repartition the reference table to reduce job overhead.
    logger.info('Repartitioning the reference table.')
    ref_ht = ref_ht.repartition(n_partitions)

    # Write.
    logger.info(f'Writing the reference to {ref_table_outpath}.')
    ref_ht.write(ref_table_outpath)


if __name__ == '__main__':
    parser = ArgumentParser(
        description='Generates a reference hail table, with centromeres and telomeres removed.'
    )
    parser.add_argument(
        '--n-partitions',
        help='Number of partitions for the resulting reference hail table.',
        required=True,
        default=10,
        type=int,
    )
    parser.add_argument(
        '--ref-table-outpath',
        help='The output location for the filtered reference hail table.',
        required=True,
        type=str,
    )
    args: Namespace = parser.parse_args()
    main(
        n_partitions=args.n_partitions,
        ref_table_outpath=args.ref_table_outpath,
    )
