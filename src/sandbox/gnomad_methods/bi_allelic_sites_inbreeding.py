import hail as hl


def bi_allelic_site_inbreeding_expr(
    call: hl.expr.CallExpression | None = None,
    callstats_expr: hl.expr.StructExpression | None = None,
) -> hl.expr.Float32Expression:
    """
    Return the site inbreeding coefficient as an expression to be computed on a MatrixTable.

    This is implemented based on the GATK InbreedingCoeff metric:
    https://software.broadinstitute.org/gatk/documentation/article.php?id=8032

    .. note::

        The computation is run based on the counts of alternate alleles and thus should only be run on bi-allelic sites.

    :param call: Expression giving the calls in the MT
    :param callstats_expr: StructExpression containing only alternate allele AC, AN, and homozygote_count as integers. If passed, used to create expression in place of GT calls.
    :return: Site inbreeding coefficient expression
    """
    if call is None and callstats_expr is None:
        raise ValueError('One of `call` or `callstats_expr` must be passed.')

    def inbreeding_coeff(
        gt_counts: hl.expr.DictExpression,
    ) -> hl.expr.Float32Expression:
        n = gt_counts.get(0, 0) + gt_counts.get(1, 0) + gt_counts.get(2, 0)
        p = (2 * gt_counts.get(0, 0) + gt_counts.get(1, 0)) / (2 * n)
        q = (2 * gt_counts.get(2, 0) + gt_counts.get(1, 0)) / (2 * n)
        return 1 - (gt_counts.get(1, 0) / (2 * p * q * n))

    if callstats_expr is not None:
        # Check that AC, AN, and homozygote count are all ints
        if not (
            ((callstats_expr.AC.dtype == hl.tint32) | (callstats_expr.AC.dtype == hl.tint64))
            & ((callstats_expr.AN.dtype == hl.tint32) | (callstats_expr.AN.dtype == hl.tint64))
            & (
                (callstats_expr.homozygote_count.dtype == hl.tint32)
                | (callstats_expr.homozygote_count.dtype == hl.tint64)
            )
        ):
            raise ValueError(
                "callstats_expr must be a StructExpression containing fields 'AC',"
                " 'AN', and 'homozygote_count' of types int32 or int64."
            )
        n = callstats_expr.AN / 2
        q = callstats_expr.AC / callstats_expr.AN
        p = 1 - q
        return 1 - (callstats_expr.AC - (2 * callstats_expr.homozygote_count)) / (2 * p * q * n)
    return hl.bind(inbreeding_coeff, hl.agg.counter(call.n_alt_alleles()))
