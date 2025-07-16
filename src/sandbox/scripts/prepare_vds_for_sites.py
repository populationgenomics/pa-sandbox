"""Create a dense MatrixTable from a VDS, along with filtering and repartitioning.

This script subsets a given VDS by either a number of samples, a region or set of regions, or both.

Arguments:
    --vds-path: str. Path to the VDS in GCS.
    --samples-to-drop: str, optional. Comma separated list of hail tables with samples to drop
    --exome-intervals: str. The path to the BED file with the exome intervals that you want to use for filtering.
    --variant-table: str. Path to the hail table that contains variants that you want to filter on
    --n-partitions: int. The number of partitions to repartition the end MatrixTable into
    --sites-table-outpath: str. The path to write the resultant MatrixTable
    --chromosome: str, optional. An optional chromosome to subset to.

Example usage:

analysis-runner --dataset hgdp-1kg \
    --access-level test \
    --output-dir test-subset \
    --description 'Test VDS subsetting script' \
    python3 prepare_vds_for_sites.py --vds-path gs://cpg-hgdp-1kg-main/vds/v1-0.vds/ \
    --exome-intervals gs://cpg-common-main/references/exome-probesets/hg38/mackenzie_intersect_exome_regions.bed \
    --variant-table gs://cpg-common-main/references/gnomad/v4.1/ht/gnomad.genomes.v4.1.sites.ht \
    --n-partitions 100 \
    --sites-table-outpath gs://cpg-common-main/references/ancestry/hgdp-1kg-wgs-pruned_variants.mt
"""

from argparse import ArgumentParser, Namespace

import hail as hl
from cpg_utils.hail_batch import init_batch


def main(
    vds_path: str,
    samples_to_drop: str | None,
    exome_intervals: str,
    variant_table: str,
    n_partitions: int,
    sites_table_outpath: str,
    chromosome: str | None,
) -> None:
    init_batch()
    input_vds: hl.vds.VariantDataset = hl.vds.read_vds(vds_path)

    if input_vds.ref_block_max_length_field not in input_vds.reference_data.globals:
        hl.vds.store_ref_block_max_length(vds_path)

    # Always subset by interval first, if possible
    # https://discuss.hail.is/t/filtering-samples-from-vds-in-google-cloud/3718/6
    exome_interval_table: hl.Table = hl.import_bed(exome_intervals, reference_genome='GRCh38')
    if chromosome:
        exome_interval_table = exome_interval_table.filter(exome_interval_table.interval.start.contig == chromosome)
    input_vds = hl.vds.filter_intervals(input_vds, exome_interval_table, keep=True, split_reference_blocks=False)

    variant_qc_table: hl.Table = hl.read_table(variant_table)
    if 'AS_FilterStatus' not in list(variant_qc_table.info.keys()):
        variant_qc_table = variant_qc_table.annotate(
            info=variant_qc_table.info.annotate(
                AS_FilterStatus=hl.if_else(hl.len(variant_qc_table.filters) == 0, 'PASS', 'FAIL'),
            ),
        )

    # Filter to variant sites that pass VQSR
    passed_variants: hl.Table = variant_qc_table.filter(variant_qc_table.info.AS_FilterStatus == 'PASS')
    input_vds = hl.vds.filter_variants(input_vds, passed_variants)

    # Remove all multiallelic sites
    input_vds = hl.vds.filter_variants(
        input_vds,
        input_vds.variant_data.filter_rows(hl.len(input_vds.variant_data.alleles) == 2).rows(),  # noqa: PLR2004
        keep=True,
    )

    if 'GT' not in input_vds.variant_data.entry:
        input_vds.variant_data = input_vds.variant_data.annotate_entries(
            GT=hl.vds.lgt_to_gt(input_vds.variant_data.LGT, input_vds.variant_data.LA)
        )

    # Filter out any samples
    if samples_to_drop:
        sample_hts: list[hl.Table] = [hl.read_table(path) for path in samples_to_drop]
        all_samples_to_drop: hl.Table = sample_hts[0]
        for ht in sample_hts[1:]:
            all_samples_to_drop = all_samples_to_drop.union(ht)
        input_vds = hl.vds.filter_samples(input_vds, all_samples_to_drop, keep=False)

    dense_mt: hl.MatrixTable = hl.vds.to_dense_mt(input_vds)
    dense_mt = dense_mt.repartition(n_partitions)
    dense_mt.write(sites_table_outpath)


if __name__ == '__main__':
    parser = ArgumentParser(
        description='VDS patch and subsest script for sites generation. Generates a dense MatrixTable.'
    )
    parser.add_argument('--vds-path', help='Path to VDS in GCP.', required=True)
    parser.add_argument(
        '--samples-to-drop',
        help='Samples(s) to drop from the VDS provided either on the command line as a comma-separated list of IDs, or in a text file with one per line.',
        required=False,
    )
    parser.add_argument(
        '--exome-intervals',
        help='BED file of exome regions to subset the VDS to.',
        required=True,
    )
    parser.add_argument(
        '--variant-table',
        help='Hail table of variants to use for filtering the VDS',
        required=True,
    )
    parser.add_argument(
        '--n-partitions',
        help='Number of partitions for the resulting MatrixTable',
        required=True,
        default=10,
        type=int,
    )
    parser.add_argument(
        '--sites-table-outpath',
        help='The output location for the densified MatrixTable',
        required=True,
        type=str,
    )
    parser.add_argument(
        '--chromosome',
        help="A chromosome, for running smaller scale tests. GRCh38 format e.g. 'chr21'",
        required=False,
        type=str,
    )
    args: Namespace = parser.parse_args()
    main(
        vds_path=args.vds_path,
        samples_to_drop=args.samples_to_drop,
        exome_intervals=args.exome_intervals,
        variant_table=args.variant_table,
        n_partitions=args.n_partitions,
        sites_table_outpath=args.sites_table_outpath,
        chromosome=args.chromosome,
    )
