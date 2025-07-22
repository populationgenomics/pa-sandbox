from loguru import logger
import itertools
import hail as hl
from typing import List, Optional, Dict

def generate_freq_group_membership_array(
    ht: hl.Table,
    strata_expr: List[Dict[str, hl.expr.StringExpression]],
    downsamplings: Optional[List[int]] = None,
    ds_pop_counts: Optional[Dict[str, int]] = None,
    remove_zero_sample_groups: bool = False,
    no_raw_group: bool = False,
) -> hl.Table:
    """
    Generate a Table with a 'group_membership' array for each sample indicating whether the sample belongs to specific stratification groups.

    .. note::
        This function is primarily used through `annotate_freq` but can be used
        independently if desired. Please see the `annotate_freq` function for more
        complete documentation.

    The following global annotations are added to the returned Table:
        - freq_meta: Each element of the list contains metadata on a stratification
          group.
        - freq_meta_sample_count: sample count per grouping defined in `freq_meta`.
        - If downsamplings or ds_pop_counts are specified, they are also added as
          global annotations on the returned Table.

    Each sample is annotated with a 'group_membership' array indicating whether the
    sample belongs to specific stratification groups. All possible value combinations
    are determined for each stratification grouping in the `strata_expr` list.

    :param ht: Input Table that contains Expressions specified by `strata_expr`.
    :param strata_expr: List of dictionaries specifying stratification groups where
        the keys of each dictionary are strings and the values are corresponding
        expressions that define the values to stratify frequency calculations by.
    :param downsamplings: List of downsampling values to include in the stratifications.
    :param ds_pop_counts: Dictionary of population counts for each downsampling value.
    :param remove_zero_sample_groups: Whether to remove groups with a sample count of 0.
        Default is False.
    :param no_raw_group: Whether to remove the raw group from the 'group_membership'
        annotation and the 'freq_meta' and 'freq_meta_sample_count' global annotations.
        Default is False.
    :return: Table with the 'group_membership' array annotation.
    """
    errors = []
    ds_in_strata = any("downsampling" in s for s in strata_expr)
    global_idx_in_ds_expr = any(
        "global_idx" in s["downsampling"] for s in strata_expr if "downsampling" in s
    )
    pop_in_strata = any("pop" in s for s in strata_expr)
    pop_idx_in_ds_expr = any(
        "pop_idx" in s["downsampling"]
        for s in strata_expr
        if "downsampling" in s and ds_pop_counts is not None
    )

    if downsamplings is not None and not ds_in_strata:
        errors.append(
            "Strata must contain a downsampling expression when downsamplings"
            "are provided."
        )
    if downsamplings is not None and not global_idx_in_ds_expr:
        errors.append(
            "Strata must contain a downsampling expression with 'global_idx' when "
            "downsamplings are provided."
        )
    if ds_pop_counts is not None and not pop_in_strata:
        errors.append(
            "Strata must contain a population expression 'pop' when ds_pop_counts "
            " are provided."
        )
    if ds_pop_counts is not None and not pop_idx_in_ds_expr:
        errors.append(
            "Strata must contain a downsampling expression with 'pop_idx' when "
            "ds_pop_counts are provided."
        )

    if errors:
        raise ValueError("The following errors were found: \n" + "\n".join(errors))

    # Get counters for all strata.
    strata_counts = ht.aggregate(
        hl.struct(
            **{
                k: hl.agg.filter(hl.is_defined(v), hl.agg.counter({k: v}))
                for strata in strata_expr
                for k, v in strata.items()
            }
        )
    )

    # Add all desired strata to sample group filters.
    sample_group_filters = [({}, True)]
    for strata in strata_expr:
        downsampling_expr = strata.get("downsampling")
        strata_values = []
        # Add to all downsampling groups, both global and population-specific, to
        # strata.
        for s in strata:
            if s == "downsampling":
                v = [("downsampling", d) for d in downsamplings]
            else:
                v = [(s, k[s]) for k in strata_counts.get(s, {})]
                if s == "pop" and downsampling_expr is not None:
                    v.append(("pop", "global"))
            strata_values.append(v)

        # Get all combinations of strata values.
        strata_combinations = itertools.product(*strata_values)
        # Create sample group filters that are evaluated on each sample for each strata
        # combination. Strata combinations are evaluated as a logical AND, e.g.
        # {"pop":nfe, "downsampling":1000} or "nfe-10000" creates the filter expression
        # pop == nfe AND downsampling pop_idx < 10000.
        for combo in strata_combinations:
            combo = dict(combo)
            ds = combo.get("downsampling")
            pop = combo.get("pop")
            # If combo contains downsampling, determine the downsampling index
            # annotation to use.
            downsampling_idx = "global_idx"
            if ds is not None:
                if pop is not None and pop != "global":
                    # Don't include population downsamplings where the downsampling is
                    # larger than the number of samples in the population.
                    if ds > ds_pop_counts[pop]:
                        continue
                    downsampling_idx = "pop_idx"

            # If combo contains downsampling, add downsampling filter expression.
            combo_filter_exprs = []
            for s, v in combo.items():
                if s == "downsampling":
                    combo_filter_exprs.append(downsampling_expr[downsampling_idx] < v)
                else:
                    if s != "pop" or v != "global":
                        combo_filter_exprs.append(strata[s] == v)
            combo = {k: str(v) for k, v in combo.items()}
            sample_group_filters.append((combo, hl.all(combo_filter_exprs)))

    n_groups = len(sample_group_filters)
    logger.info("number of filters: %i", n_groups)

    # Get sample count per strata group.
    freq_meta_sample_count = ht.aggregate(
        [hl.agg.count_where(x[1]) for x in sample_group_filters]
    )

    if remove_zero_sample_groups:
        filter_freq = hl.enumerate(freq_meta_sample_count).filter(lambda x: x[1] > 0)
        freq_meta_sample_count = filter_freq.map(lambda x: x[1])
        idx_keep = hl.eval(filter_freq.map(lambda x: x[0]))
        sample_group_filters = [sample_group_filters[i] for i in idx_keep]

    # Annotate columns with group_membership.
    ht = ht.select(group_membership=[x[1] for x in sample_group_filters])

    # Create and annotate global expression with meta and sample count information.
    freq_meta = [
        dict(**sample_group[0], group="adj") for sample_group in sample_group_filters
    ]

    if not no_raw_group:
        # Sample group membership for the "raw" group, representing all samples, is
        # the same as the first group in the group_membership array.
        ht = ht.annotate(
            group_membership=hl.array([ht.group_membership[0]]).extend(
                ht.group_membership
            )
        )
        # Add the "raw" group, representing all samples, to the freq_meta_expr list.
        freq_meta.insert(1, {"group": "raw"})
        freq_meta_sample_count = hl.array([freq_meta_sample_count[0]]).extend(
            freq_meta_sample_count
        )

    global_expr = {
        "freq_meta": freq_meta,
        "freq_meta_sample_count": freq_meta_sample_count,
    }

    if downsamplings is not None:
        global_expr["downsamplings"] = downsamplings
    if ds_pop_counts is not None:
        global_expr["ds_pop_counts"] = ds_pop_counts

    ht = ht.select_globals(**global_expr)
    ht = ht.checkpoint(hl.utils.new_temp_file("group_membership", "ht"))

    return ht
