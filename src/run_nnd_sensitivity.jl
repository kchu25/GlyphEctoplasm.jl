"""
NND sensitivity analysis: sweep k values across all motifs and apply FDR correction.
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
    rows = Vector{NamedTuple{(:dataset,:motif_id,:motif_type,:k,:n_subpop,:n_background,:obs_mNND,:raw_pvalue),
                              Tuple{String,String,String,Int,Int,Int,Float64,Float64}}}()

    # ── Singletons ──────────────────────────────────────────────
    sep_by = build_grouping_columns(:filter_index)
    gdf_filters = groupby(contributions_df, sep_by)
    sorted_keys_s, _, _, _, _ = build_sorted_keys_and_maps(gdf_filters, sep_by)

    for (i, sk) in enumerate(sorted_keys_s)
        intersect_indices = intersect(gdf_filters[sk].data_pt_index, all_indices)
        subset_labels = Vector{Float64}(@view pts.labels[intersect_indices])
        m = length(subset_labels)
        m < 2 && continue   # need at least 2 points

        motif_id = string(i)
        for kv in ks
            kv >= m + N_bg && continue   # skip if k ≥ pool size
            res = nnd_permutation_test_1d(subset_labels, background; k=kv)
            push!(rows, (dataset=page_title, motif_id=motif_id, motif_type="singleton",
                         k=kv, n_subpop=m, n_background=N_bg,
                         obs_mNND=res.obs_mNND, raw_pvalue=res.p_value))
        end
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
            m = length(subset_labels)
            m < 2 && continue

            motif_id = get_filter_indices_str(mk)
            for kv in ks
                kv >= m + N_bg && continue
                res = nnd_permutation_test_1d(subset_labels, background; k=kv)
                push!(rows, (dataset=page_title, motif_id=motif_id, motif_type=group_name,
                             k=kv, n_subpop=m, n_background=N_bg,
                             obs_mNND=res.obs_mNND, raw_pvalue=res.p_value))
            end
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

export run_nnd_sensitivity_analysis
