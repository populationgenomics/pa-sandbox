from typing import TYPE_CHECKING

import hail as hl
from cpg_flow.targets import Cohort
from cpg_utils.config import config_retrieve, get_driver_image
from cpg_utils.hail_batch import get_batch, init_batch, output_path  # type: ignore[ReportUnknownVariableType]
from gnomad.utils.sparse_mt import default_compute_info
from hailtop.batch.job import PythonJob
from loguru import logger

if TYPE_CHECKING:
    from hail.vds.variant_dataset import VariantDataset


INFO_VCF_AS_PIPE_DELIMITED_FIELDS = [
    "AS_QUALapprox",
    "AS_VarDP",
    "AS_MQ_DP",
    "AS_RAW_MQ",
    "AS_SB_TABLE",
]

FORMAT_DICT = {
    'GQ': {
    'Description': (
        'Phred-scaled confidence that the genotype assignment is correct. '
        'Value is the difference between the second lowest PL and the lowest PL '
        '(always normalized to 0).'
    ),
    'Type': 'Integer'
},
    'GT': {
        'Description': 'Genotype',
        'Type': 'String'
    },
    'DP': {
        'Description': 'Approximate read depth (reads with MQ=255 or with bad mates are filtered)',
        'Type': 'Integer'
    },
    'AD': {
        'Description': 'Allelic depths for the ref and alt alleles in the order listed',
        'Number': 'R',
        'Type': 'Integer'
    },
    'MIN_DP': {
        'Description': 'Minimum DP observed within the GVCF block',
        'Type': 'Integer'
    },
    'GP': {
        'Description': 'Phred-scaled posterior probabilities for genotypes as defined in the VCF specification'
    },
    'PGT': {
    'Description': (
        'Physical phasing haplotype information, describing how the alternate '
        'alleles are phased in relation to one another'
    ),
    'Type': 'String'
    },
    'RGQ': {
        'Description': (
            'Unconditional reference genotype confidence, encoded as a phred quality '
            '-10*log10 p(genotype call is wrong)'
        )
    },
    'PID': {
        'Description': (
            'Physical phasing ID information, where each unique ID within a given '
            'sample (but not across samples) connects records within a phasing group'
        ),
        'Type': 'String'
    },
    'PG': {
        'Description': 'genotype priors in Phred Scale'
    },
    'PL': {
        'Description': 'Normalized, phred-scaled likelihoods for genotypes as defined in the VCF specification',
        'Number': 'G',
        'Type': 'Integer'
    },
    'PS': {
        'Description': (
            'Physical phasing set information, identifying the position of the first '
            'variant in the phasing group'
        )
    },
    'SB': {
        'Description': (
            "Per-sample component statistics which comprise the Fisher's exact test to "
            "detect strand bias. Values are: depth of reference allele on forward strand, "
            "depth of reference allele on reverse strand, depth of alternate allele on "
            "forward strand, depth of alternate allele on reverse strand."
        ),
        'Type': 'Integer'
    }
}

INFO_DICT = {
    'AS_ReadPosRankSum': {
        'Description': 'AS Z-score from Wilcoxon rank sum test of alternate vs. reference read position bias'
    },
    'AS_MQRankSum': {
        'Description': 'AS Z-score from Wilcoxon rank sum test of alternate vs. reference read mapping qualities'
    },
    'AS_QUALapprox': {
        'Description': 'AS Sum of PL[0] values; used to approximate the QUAL score'
    },
    'AS_VarDP': {
        'Description': 'AS Depth over variant genotypes (does not include depth of reference samples)'
    },
    'AS_SB_TABLE': {
        'Number': '.',
        'Description': 'Allele-specific forward/reverse read counts for strand bias tests'
    },
    'AS_MQ': {
        'Description': 'AS Root mean square of the mapping quality of reads across all samples'
    },
    'AS_QD': {
        'Description': 'AS Variant call confidence normalized by depth of sample reads supporting a variant'
    },
    'AS_FS': {
        'Description': "AS Phred-scaled p-value of Fisher's exact test for strand bias"
    },
    'AS_SOR': {
        'Description': 'AS Strand bias estimated by the symmetric odds ratio test'
    },
    'AS_pab_max': {
        'Number': 'A',
        'Description': (
            'Maximum p-value over callset for binomial test of observed allele balance '
            'for a heterozygous genotype, given expectation of 0.5'
        )
    },
    'ReadPosRankSum': {
        'Description': 'Z-score from Wilcoxon rank sum test of alternate vs. reference read position bias'
    },
    'MQRankSum': {
        'Description': 'Z-score from Wilcoxon rank sum test of alternate vs. reference read mapping qualities'
    },
    'QUALapprox': {
        'Number': '1',
        'Description': 'Sum of PL[0] values; used to approximate the QUAL score'
    },
    'VarDP': {
        'Description': 'Depth over variant genotypes (does not include depth of reference samples)'
    },
    'MQ': {
        'Description': 'Root mean square of the mapping quality of reads across all samples'
    },
    'QD': {
        'Description': 'Variant call confidence normalized by depth of sample reads supporting a variant'
    },
    'FS': {
        'Description': "Phred-scaled p-value of Fisher's exact test for strand bias"
    },
    'SOR': {
        'Description': 'Strand bias estimated by the symmetric odds ratio test'
    }
}

def adjust_vcf_incompatible_types(
    ht: hl.Table,
    pipe_delimited_annotations: list[str] = INFO_VCF_AS_PIPE_DELIMITED_FIELDS,
) -> hl.Table:
    """
    Create a Table ready for vcf export.

    In particular, the following conversions are done:
        - All int64 are coerced to int32
        - Fields specified by `pipe_delimited_annotations` are converted from arrays to pipe-delimited strings

    :param ht: Input Table.
    :param pipe_delimited_annotations: List of info fields (they must be fields of the ht.info Struct).
    :return: Table ready for VCF export.
    """

    def get_pipe_expr(array_expr: hl.expr.ArrayExpression) -> hl.expr.StringExpression:
        return hl.delimit(array_expr.map(lambda x: hl.or_else(hl.str(x), "")), "|")

    # Make sure the HT is keyed by locus, alleles
    ht = ht.key_by("locus", "alleles")

    info_type_convert_expr = {}
    # Convert int64 fields to int32 (int64 isn't supported by VCF)
    for f, ft in ht.info.dtype.items():
        if ft == hl.dtype("int64"):
            logger.warning(
                "Coercing field info.%s from int64 to int32 for VCF output. Value"
                " will be capped at int32 max value.",
                f,
            )
            info_type_convert_expr.update(
                {f: hl.int32(hl.min(2**31 - 1, ht.info[f]))}
            )
        elif ft == hl.dtype("array<int64>"):
            logger.warning(
                "Coercing field info.%s from array<int64> to array<int32> for VCF"
                " output. Array values will be capped at int32 max value.",
                f,
            )
            info_type_convert_expr.update(
                {f: ht.info[f].map(lambda x: hl.int32(hl.min(2**31 - 1, x)))}
            )

    ht = ht.annotate(info=ht.info.annotate(**info_type_convert_expr))

    info_expr = {}

    # Make sure to pipe-delimit fields that need to.
    # Note: the expr needs to be prefixed by "|" because GATK expect one value for the ref (always empty)
    # Note2: this doesn't produce the correct annotation for AS_SB_TABLE, it
    # is handled below
    for f in pipe_delimited_annotations:
        if f in ht.info and f != "AS_SB_TABLE":
            info_expr[f] = "|" + get_pipe_expr(ht.info[f])

    # Flatten SB if it is an array of arrays
    if "SB" in ht.info and not isinstance(ht.info.SB, hl.expr.ArrayNumericExpression):
        info_expr["SB"] = ht.info.SB[0].extend(ht.info.SB[1])

    if "AS_SB_TABLE" in ht.info:
        info_expr["AS_SB_TABLE"] = get_pipe_expr(
            ht.info.AS_SB_TABLE.map(lambda x: hl.delimit(x, ","))
        )

    # Annotate with new expression
    return ht.annotate(info=ht.info.annotate(**info_expr))

def extract_as_pls(
    lpl_expr: hl.expr.ArrayExpression,
    allele_idx: hl.expr.Int32Expression,
) -> hl.expr.ArrayExpression:
    """
    Function pulled from `gnomad_qc/v4/annotations/generate_variant_qc_annotations.py`.

    Extract PLs for a specific allele from an LPL array expression.

    PL/LPL represents the normalized Phred-scaled likelihoods of the possible
    genotypes from all considered alleles (or local alleles).

    If three alleles are considered, LPL genotype indexes are:
    [0/0, 0/1, 1/1, 0/2, 1/2, 2/2].

    If we want to extract the PLs for each alternate allele, we need to extract:
        - allele 1: [0/0, 0/1, 1/1]
        - allele 2: [0/0, 0/2, 2/2]

    Example:
        - LPL: [138, 98, 154, 26, 0, 14]
        - Extract allele 1 PLs: [0/0, 0/1, 1/1] -> [138, 98, 154]
        - Extract allele 2 PLs: [0/0, 0/2, 2/2] -> [138, 26, 14]

    :param lpl_expr: LPL ArrayExpression.
    :param allele_idx: The index of the alternate allele to extract PLs for.
    :return: ArrayExpression of PLs for the specified allele.
    """
    calls_to_keep = hl.array(
        [hl.call(0, 0), hl.call(0, allele_idx), hl.call(allele_idx, allele_idx)],
    )
    return calls_to_keep.map(lambda c: lpl_expr[c.unphased_diploid_gt_index()])


def recompute_as_qualapprox_from_lpl(mt: hl.MatrixTable) -> hl.expr.ArrayExpression:
    """
    Function pulled from `gnomad_qc/v4/annotations/generate_variant_qc_annotations.py`.

    Recompute AS_QUALapprox from LPL.

    QUALapprox is the (Phred-scaled) probability that all reads at the site are hom-ref,
    so QUALapprox is PL[0]. To get the QUALapprox for just one allele, pull out the
    PLs for just that allele, then normalize by subtracting the smallest element from
    all the entries (so the best genotype is 0) and then use the normalized PL[0]
    value for that allele's QUALapprox.

    .. note::

        - The first element of AS_QUALapprox is always None.
        - If the allele is a star allele, we set QUALapprox for that allele to 0.
        - If GQ == 0 and PL[0] for the allele == 1, we set QUALapprox for the allele
          to 0.

    Example:
        Starting Values:
            - alleles: [‘G’, ‘*’, ‘A’, ‘C’, ‘GCTT’, ‘GT’, ‘T’]
            - LGT: 1/2
            - LA: [0, 1, 6]
            - LPL: [138, 98, 154, 26, 0, 14]
            - QUALapprox: 138

        Use `extract_as_pls` to get PLs for each allele:
            - allele 1: [138, 98, 154]
            - allele 2: [138, 26, 14]

        Normalize PLs by subtracting the smallest element from all the PLs:
            - allele 1: [138-98, 98-98, 154-98] -> [40, 0, 56]
            - allele 2: [138-14, 26-14, 14-14] -> [124, 12, 0]

        Use the first element of the allele specific PLs to generate AS_QUALapprox:
        [None, 40, 124]

        Set QUALapprox to 0 for the star allele: [None, 0, 124]

    :param mt: Input MatrixTable.
    :return: AS_QUALapprox ArrayExpression recomputed from LPL.
    """
    return hl.enumerate(mt.LA).map(
        lambda i: (
            hl.case()
            .when(mt.alleles[i[1]] == "*", 0)
            .when(
                i[0] > 0,
                hl.bind(
                    lambda pl_0: hl.if_else((mt.GQ == 0) & (pl_0 == 1), 0, pl_0),
                    hl.bind(lambda x: x[0] - hl.min(x), extract_as_pls(mt.LPL, i[0])),
                ),
            )
            .or_missing()
        ),
    )


def correct_as_annotations(
    mt: hl.MatrixTable,
    set_to_missing: bool = False,
) -> hl.expr.StructExpression:
    """
    Function pulled from `gnomad_qc/v4/annotations/generate_variant_qc_annotations.py`.

    Correct allele specific annotations that are longer than the length of LA.

    For some entries in the MatrixTable, the following annotations are longer than LA,
    when they should be the same length as LA:

        - AS_SB_TABLE
        - AS_RAW_MQ
        - AS_RAW_ReadPosRankSum
        - AS_RAW_MQRankSum

    This function corrects these annotations by either dropping the alternate allele
    with the index corresponding to the min value of AS_RAW_MQ, or setting them to
    missing if `set_to_missing` is True.

    :param mt: Input MatrixTable.
    :param set_to_missing: Whether to set the annotations to missing instead of
        correcting them.
    :return: StructExpression with corrected allele specific annotations.
    """
    annotations_to_correct = [
        "AS_SB_TABLE",
        "AS_RAW_MQ",
        "AS_RAW_ReadPosRankSum",
        "AS_RAW_MQRankSum",
    ]
    annotations_to_correct_dict = {a: mt.gvcf_info[a] for a in annotations_to_correct}

    # Identify index corresponding to min of AS_RAW_MQ, skipping the reference allele.
    as_raw_mq_no_ref = mt.gvcf_info.AS_RAW_MQ[1:]
    idx_remove = as_raw_mq_no_ref.index(hl.min(as_raw_mq_no_ref)) + 1

    corrected_annotations = {
        a: hl.if_else(
            hl.len(expr) > hl.len(mt.LA),
            hl.or_missing(
                not set_to_missing,
                expr[:idx_remove].extend(expr[idx_remove + 1 :]),
            ),
            expr,
        )
        for a, expr in annotations_to_correct_dict.items()
    }

    return hl.struct(**corrected_annotations)


def _filter_rows_and_add_tags(mt: hl.MatrixTable) -> hl.MatrixTable:
    # Filter to only non-reference sites.
    # An example of a variant with hl.len(mt.alleles) > 1 BUT NOT
    # hl.agg.any(mt.LGT.is_non_ref()) is a variant that spans a deletion,
    # which was however filtered out, so the LGT was set to NA, however the site
    # was preserved to account for the presence of that spanning deletion.
    # locus   alleles    LGT
    # chr1:1 ["GCT","G"] 0/1
    # chr1:3 ["T","*"]   NA
    mt = mt.filter_rows((hl.len(mt.alleles) > 1) & (hl.agg.any(mt.LGT.is_non_ref())))

    # annotate site level DP as site_dp onto the mt rows to avoid name collision
    mt = mt.annotate_rows(site_dp=hl.agg.sum(mt.DP))

    # Add AN tag as ANS
    return mt.annotate_rows(ANS=hl.agg.count_where(hl.is_defined(mt.LGT)) * 2)


def _create_info_ht(mt: hl.MatrixTable, n_partitions: int) -> hl.Table:
    """Create info table from vcf matrix table"""

    info_ht: hl.Table = default_compute_info(
        mt,
        as_annotations=True,
        site_annotations=True,
        n_partitions=n_partitions,
    )

    return info_ht.annotate(info=info_ht.info.annotate(DP=mt.rows()[info_ht.key].site_dp))


def _initalise_vds_export_to_vcf_job(cohort: Cohort, name: str, job_memory: str, job_cpus: int) -> PythonJob:
    job: PythonJob = get_batch().new_python_job(
        name=name,
        attributes=cohort.get_job_attrs() or {} | {'tool': 'Hail'},  # type: ignore[ReportUnknownVariableType]
    )
    job.image(image=get_driver_image())
    job.memory(job_memory)
    job.cpu(job_cpus)
    job.n_max_attempts(2)
    return job


def vds_to_vcf(cohort: Cohort, vds_path: str, vcf_outpath: str, chrom: str) -> PythonJob:
    job_memory: str = config_retrieve(['workflow', 'job_memory'])
    job_cpus: int = config_retrieve(['workflow', 'job_cpus'])

    vds_to_vcf_job: PythonJob = _initalise_vds_export_to_vcf_job(
        cohort=cohort,
        name=f'ExportVdsToVcf-{chrom}',
        job_memory=job_memory,
        job_cpus=job_cpus,
    ).call(
        _run_vds_to_vcf,
        vds_path=vds_path,
        vcf_outpath=vcf_outpath,
        chrom=chrom
    )
    return vds_to_vcf_job

def globalise_entries(vds: 'VariantDataset') -> 'VariantDataset':
    vds.variant_data = vds.variant_data.annotate_entries(
            GT=hl.vds.lgt_to_gt(vds.variant_data.LGT, vds.variant_data.LA),
        )

    vds.variant_data = vds.variant_data.annotate_entries(
            AD=hl.vds.local_to_global(
                vds.variant_data.LAD,
                vds.variant_data.LA,
                n_alleles=hl.len(vds.variant_data.alleles),
                fill_value=0,
                number='R',
            ),
        )

    vds.variant_data = vds.variant_data.annotate_entries(
            PL=hl.vds.local_to_global(
                vds.variant_data.LPL,
                vds.variant_data.LA,
                n_alleles=hl.len(vds.variant_data.alleles),
                fill_value=0,
                number='G',
            ),
        )

    vds.variant_data = vds.variant_data.annotate_entries(
        PGT=hl.vds.lgt_to_gt(vds.variant_data.LPGT, vds.variant_data.LA)
    )

    return vds

def generate_info_ht(mt: 'hl.MatrixTable', chrom: str) -> hl.Table:
    """
    Code taken from MakeSiteOnlyVcf stage in large_cohorts.py
    """
    # Compute and checkpoint the allele specific info annotations after recomputing
    # AS_QUALapprox from LPL, and fixing the length of AS_SB_TABLE, AS_RAW_MQ,
    # AS_RAW_ReadPosRankSum and AS_RAW_MQRankSum.
    mt = mt.annotate_rows(alt_alleles_range_array=hl.range(1, hl.len(mt.alleles)))

    correct_mt = mt.annotate_entries(
        gvcf_info=mt.gvcf_info.annotate(
            AS_QUALapprox=recompute_as_qualapprox_from_lpl(mt),
            **correct_as_annotations(mt),
        ),
    )

    correct_mt = _filter_rows_and_add_tags(correct_mt)

    info_ht = _create_info_ht(correct_mt, n_partitions=mt.n_partitions())

    info_ht = info_ht.drop('quasi_info', 'AS_lowqual', 'lowqual')

    info_ht = adjust_vcf_incompatible_types(
        info_ht,
        # With default INFO_VCF_AS_PIPE_DELIMITED_FIELDS, AS_VarDP will be converted
        # into a pipe-delimited value e.g.: VarDP=|132.1|140.2
        # which breaks VQSR parser (it doesn't recognise the delimiter and treats
        # it as an array with a single string value "|132.1|140.2", leading to
        # an IndexOutOfBound exception when trying to access value for second allele)
        pipe_delimited_annotations=[],
    )

    return info_ht.checkpoint(output_path(f'{chrom}_info_ht.ht', category='tmp'), overwrite=True)

def _run_vds_to_vcf(vds_path: str, vcf_outpath: str, chrom: str) -> str:
    init_batch()

    vds: VariantDataset = hl.vds.read_vds(vds_path)

    # vds = hl.vds.filter_chromosomes(
    #     vds,
    #     keep=chrom,
    # )

    logger.info('Filtering to spanning deletions...')
    spanning_deletions_mt = vds.variant_data.filter_rows(
        vds.variant_data.alleles.contains('*')
    )
    logger.info(
        f'Checkpointing spanning deletions MT to {output_path("all_chrom_spanning_deletions.mt", category="tmp", test=True)}'
    )
    spanning_deletions_mt = spanning_deletions_mt.repartition(2500)
    spanning_deletions_mt = spanning_deletions_mt.checkpoint(
        output_path('all_chrom_spanning_deletions.mt', category='tmp', test=True),
        overwrite=True
    )

    # logger.info('Globalising entries')
    # vds = globalise_entries(vds)

    # logger.info('Densifying...')
    # mt = hl.vds.to_dense_mt(vds)

    # logger.info(f'Checkpointing dense MT to {output_path(f"{chrom}_dense_mt.mt", category="tmp")}')
    # mt = mt.checkpoint(output_path(f'{chrom}_dense_mt.mt', category='tmp'), overwrite=True)


    # logger.info('Generating info_ht...')
    # info_ht = generate_info_ht(mt, chrom=chrom)

    # # Annotate back the MT with the info HT
    # logger.info('Annotating MT with info_ht')
    # mt = mt.annotate_rows(info=info_ht[mt.row_key].info)

    # # Drop local and unnecessary fields
    # logger.info('Dropping local and unnecessary fields')
    # mt = mt.drop('gvcf_info')
    # mt = mt.drop('LAD', 'LGT', 'LA', 'LPL', 'LPGT')

    # logger.info(f'Exporting {chrom} VCF to {vcf_outpath}')
    # metadata = {'info': INFO_DICT, 'format': FORMAT_DICT}
    # hl.export_vcf(mt, vcf_outpath, tabix=True, metadata=metadata)

    return vcf_outpath
