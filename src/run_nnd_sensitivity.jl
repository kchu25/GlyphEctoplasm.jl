"""
NND sensitivity analysis: sweep k values across all motifs and apply FDR correction.

Uses `nnd_sensitivity_batch_1d` to evaluate all k values per motif in a single
pass (shared sort, shared distance table, shared permutations), giving a
~length(ks)× speedup over calling `nnd_permutation_test_1d` separately.
"""

"""
    run_nnd_sensitivity_analysis(contributions_df, dfs, all_indices, pts, motif_sizes;
                                 ks=SENSITIVITY_KS, save_path="tmp", page_title="n/a")

Run NND permutation tests across multiple k values for every motif (singletons + multi-motifs).
Applies Benjamini–Hochberg FDR correction across all tests and saves results to CSV.

# Arguments
- `contributions_df`: DataFrame with singleton contribution data (already sorted/remapped)
- `dfs`: Vector of DataFrames for multi-motif data (already remapped)
- `all_indices`: Training data point indices
- `pts`: Training data NamedTuple with `.labels`
- `motif_sizes`: Vector of multi-motif sizes (e.g., [2, 3])

# Keyword Arguments
- `ks::Vector{Int}`: k values to sweep (default: `SENSITIVITY_KS` = [1,3,5,7,10,15,20])
- `save_path::String`: Directory for saving CSV output
- `page_title::String`: Dataset identifier written into the CSV

# Output
Saves `nnd_sensitivity.csv` with columns:
  dataset, motif_id, motif_type, k, n_subpop, n_background, obs_mNND, raw_pvalue, adjusted_pvalue

# Returns
DataFrame with all sensitivity analysis results (same as saved CSV content)
"""
function run_nnd_sensitivity_analysis(
        contributions_df, dfs, all_indices, pts, motif_sizes;
        ks::Vector{Int} = SENSITIVITY_KS,
        save_path::String = "tmp",
        page_title::String = "n/a"
    )
    @info "Running NND sensitivity analysis across k ∈ $ks …"
    background = pts.labels
    N_bg = length(background)
    
    # Accumulate rows: (dataset, motif_id, motif_type, k, n_subpop, n_background, obs_mNND, raw_pvalue)
    RowType = NamedTuple{(:dataset,:motif_id,:motif_type,:k,:n_subpop,:n_background,:obs_mNND,:raw_pvalue),
                          Tuple{String,String,String,Int,Int,Int,Float64,Float64}}
    rows = RowType[]

    # Helper: run batched NND and append rows for one motif
    function _process_motif!(subset_labels, motif_id::String, motif_type::String)
        m = length(subset_labels)
        m < 2 && return
        # Filter ks that are valid for this pool size
        valid_ks = filter(kv -> kv < m + N_bg, ks)
        isempty(valid_ks) && return

        batch_results = nnd_sensitivity_batch_1d(subset_labels, background; ks=valid_ks)
        for res in batch_results
            push!(rows, (dataset=page_title, motif_id=motif_id, motif_type=motif_type,
                         k=res.k, n_subpop=m, n_background=N_bg,
                         obs_mNND=res.obs_mNND, raw_pvalue=res.p_value))
        end
    end

    # ── Singletons ──────────────────────────────────────────────
    sep_by = build_grouping_columns(:filter_index)
    gdf_filters = groupby(contributions_df, sep_by)
    sorted_keys_s, _, _, _, _ = build_sorted_keys_and_maps(gdf_filters, sep_by)

    for (i, sk) in enumerate(sorted_keys_s)
        intersect_indices = intersect(gdf_filters[sk].data_pt_index, all_indices)
        subset_labels = Vector{Float64}(@view pts.labels[intersect_indices])
        _process_motif!(subset_labels, string(i), "singleton")
    end

    # ── Multi-motifs ────────────────────────────────────────────
    motif_names_local = ["pairs", "triplets", "quadruplets", "quintuplets"]
    for (index, motif_size) in enumerate(motif_sizes)
        df_idx = motif_size - 1
        (df_idx < 1 || df_idx > length(dfs)) && continue
        subdf = dfs[df_idx]
        sep_by_m = build_grouping_columns(:motifs; motif_size=motif_size)
        gdf_by_msyms = groupby(subdf, sep_by_m)
        sorted_keys_m, _, _, _, _ = build_sorted_keys_and_maps(gdf_by_msyms, sep_by_m)
        group_name = motif_size <= 5 ? motif_names_local[min(motif_size - 1, 4)] : "$(motif_size)-motifs"

        for (j, mk) in enumerate(sorted_keys_m)
            intersect_indices = intersect(gdf_by_msyms[mk].data_pt_index, all_indices)
            subset_labels = Vector{Float64}(@view pts.labels[intersect_indices])
            _process_motif!(subset_labels, get_filter_indices_str(mk), group_name)
        end
    end

    # ── BH FDR correction ───────────────────────────────────────
    if isempty(rows)
        @warn "No motifs processed for sensitivity analysis"
        return DataFrame()
    end

    raw_ps = [r.raw_pvalue for r in rows]
    adj_ps = adjust(PValues(raw_ps), BenjaminiHochberg())

    # Build final DataFrame
    df_out = DataFrame(rows)
    df_out.adjusted_pvalue = adj_ps

    csv_path = joinpath(save_path, "nnd_sensitivity.csv")
    CSV.write(csv_path, df_out)
    @info "NND sensitivity analysis saved to $csv_path  ($(nrow(df_out)) rows)"
    
    return df_out
end


"""
    run_nnd_sensitivity_analysis_null(contributions_df, dfs, all_indices, pts, motif_sizes;
                                      ks=SENSITIVITY_KS, save_path="tmp", page_title="n/a", seed=123)

Null-calibration companion to `run_nnd_sensitivity_analysis`.
For each motif, instead of using its actual data-point indices, randomly draws
the **same number** of indices from `all_indices` (guaranteed no clustering by
construction) and runs the batched NND sensitivity test.

This lets you verify that the false-positive rate is ≈ α for your specific
background distribution and sample sizes. Compare the resulting p-value
distribution against the real `nnd_sensitivity.csv` to assess robustness.

# Arguments
Same as `run_nnd_sensitivity_analysis`.

# Additional Keyword Arguments
- `seed::Int=123`: Random seed for reproducible random-subset draws.
  (Distinct from the permutation seed inside `nnd_sensitivity_batch_1d`.)

# Output
Saves `nnd_sensitivity_null.csv` with the same columns as the real analysis:
  dataset, motif_id, motif_type, k, n_subpop, n_background, obs_mNND, raw_pvalue, adjusted_pvalue

# Returns
DataFrame with all null sensitivity analysis results.
"""
function run_nnd_sensitivity_analysis_null(
        contributions_df, dfs, all_indices, pts, motif_sizes;
        ks::Vector{Int} = SENSITIVITY_KS,
        save_path::String = "tmp",
        page_title::String = "n/a",
        seed::Int = 123
    )
    @info "Running NULL NND sensitivity analysis across k ∈ $ks …"
    rng = Random.MersenneTwister(seed)
    background = pts.labels
    N_bg = length(background)
    N_all = length(all_indices)

    RowType = NamedTuple{(:dataset,:motif_id,:motif_type,:k,:n_subpop,:n_background,:obs_mNND,:raw_pvalue),
                          Tuple{String,String,String,Int,Int,Int,Float64,Float64}}
    rows = RowType[]

    # Helper: draw m random indices from all_indices, extract labels, run batched NND
    function _process_null_motif!(m::Int, motif_id::String, motif_type::String)
        m < 2 && return
        m > N_all && return
        # Random draw (without replacement) of m indices from all_indices
        rand_idx = all_indices[Random.randperm(rng, N_all)[1:m]]
        subset_labels = Vector{Float64}(@view pts.labels[rand_idx])

        valid_ks = filter(kv -> kv < m + N_bg, ks)
        isempty(valid_ks) && return

        batch_results = nnd_sensitivity_batch_1d(subset_labels, background; ks=valid_ks)
        for res in batch_results
            push!(rows, (dataset=page_title, motif_id=motif_id, motif_type=motif_type,
                         k=res.k, n_subpop=m, n_background=N_bg,
                         obs_mNND=res.obs_mNND, raw_pvalue=res.p_value))
        end
    end

    # ── Singletons (same grouping, just to get the m per motif) ─
    sep_by = build_grouping_columns(:filter_index)
    gdf_filters = groupby(contributions_df, sep_by)
    sorted_keys_s, _, _, _, _ = build_sorted_keys_and_maps(gdf_filters, sep_by)

    for (i, sk) in enumerate(sorted_keys_s)
        intersect_indices = intersect(gdf_filters[sk].data_pt_index, all_indices)
        m = length(intersect_indices)
        _process_null_motif!(m, string(i), "singleton")
    end

    # ── Multi-motifs ────────────────────────────────────────────
    motif_names_local = ["pairs", "triplets", "quadruplets", "quintuplets"]
    for (index, motif_size) in enumerate(motif_sizes)
        df_idx = motif_size - 1
        (df_idx < 1 || df_idx > length(dfs)) && continue
        subdf = dfs[df_idx]
        sep_by_m = build_grouping_columns(:motifs; motif_size=motif_size)
        gdf_by_msyms = groupby(subdf, sep_by_m)
        sorted_keys_m, _, _, _, _ = build_sorted_keys_and_maps(gdf_by_msyms, sep_by_m)
        group_name = motif_size <= 5 ? motif_names_local[min(motif_size - 1, 4)] : "$(motif_size)-motifs"

        for (j, mk) in enumerate(sorted_keys_m)
            intersect_indices = intersect(gdf_by_msyms[mk].data_pt_index, all_indices)
            m = length(intersect_indices)
            _process_null_motif!(m, get_filter_indices_str(mk), group_name)
        end
    end

    # ── BH FDR correction ───────────────────────────────────────
    if isempty(rows)
        @warn "No motifs processed for null sensitivity analysis"
        return DataFrame()
    end

    raw_ps = [r.raw_pvalue for r in rows]
    adj_ps = adjust(PValues(raw_ps), BenjaminiHochberg())

    df_out = DataFrame(rows)
    df_out.adjusted_pvalue = adj_ps

    csv_path = joinpath(save_path, "nnd_sensitivity_null.csv")
    CSV.write(csv_path, df_out)
    @info "NULL NND sensitivity analysis saved to $csv_path  ($(nrow(df_out)) rows)"

    return df_out
end

export run_nnd_sensitivity_analysis, run_nnd_sensitivity_analysis_null
