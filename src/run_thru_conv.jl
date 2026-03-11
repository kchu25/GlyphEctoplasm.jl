
const motif_names = ["pairs", "triplets", "quadruplets", "quintuplets"]

# Default k values for NND sensitivity analysis
const SENSITIVITY_KS = [1, 3, 5, 7, 10, 15, 20]

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
- `ks::Vector{Int}`: k values to sweep (default: `[1,3,5,7,10,15,20]`)
- `save_path::String`: Directory for saving CSV output
- `page_title::String`: Dataset identifier written into the CSV

# Output
Saves `nnd_sensitivity.csv` with columns:
  dataset, motif_id, motif_type, k, n_subpop, n_background, obs_mNND, raw_pvalue, adjusted_pvalue
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
    for (index, motif_size) in enumerate(motif_sizes)
        df_idx = motif_size - 1
        (df_idx < 1 || df_idx > length(dfs)) && continue
        subdf = dfs[df_idx]
        sep_by_m = build_grouping_columns(:motifs; motif_size=motif_size)
        gdf_by_msyms = groupby(subdf, sep_by_m)
        sorted_keys_m, _, _, _, _ = build_sorted_keys_and_maps(gdf_by_msyms, sep_by_m)
        group_name = motif_size <= 5 ? motif_names[min(motif_size - 1, 4)] : "$(motif_size)-motifs"

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
        return
    end

    raw_ps = [r.raw_pvalue for r in rows]
    adj_ps = adjust(PValues(raw_ps), BenjaminiHochberg())

    # Build final DataFrame
    df_out = DataFrame(rows)
    df_out.adjusted_pvalue = adj_ps

    csv_path = joinpath(save_path, "nnd_sensitivity.csv")
    CSV.write(csv_path, df_out)
    @info "NND sensitivity analysis saved to $csv_path  ($(nrow(df_out)) rows)"
end

# Collect extrema across all DataFrames without allocation
function obtain_xlim(contributions_df_filtered_singletons, dfs)
    xlim = mapreduce(extrema, 
        (a, b) -> (min(a[1], b[1]), max(a[2], b[2])),
        (contributions_df_filtered_singletons.banzhaf, 
            (df.banzhaf for df in dfs)...))

    return xlim
end


function plot_motifs_conv_case(data, m, motif_sizes, 
        contributions_df_filtered_singletons, dfs, pts, pts_test, all_indices;
        interaction_summaries=nothing,
        nav_page_count=4,
        enable_colored_borders = true,
        use_unified=true,
        dpi=65, 
        save_path="tmp", 
        page_title="n/a", 
        rna=false,
        sensitivity_analysis::Bool=false,
        dataset_name::Union{String,Nothing}=nothing
        );

    # motif rendering
    xlim = obtain_xlim(contributions_df_filtered_singletons, dfs)

    config = ConvMotifConfig(data; 
        filter_len=m.hp.pfm_len, dpi=dpi, save_path=save_path, xlim=xlim)
        
    json_motifs = init_json_dict()
    html_dict = init_dict_for_html_render()

    next_idx, sorted_mapping = process_singletons!(
        contributions_df_filtered_singletons, all_indices, pts, config, json_motifs, html_dict; start_idx=1, rna=rna)

    # Remap filter indices in multi-motif DataFrames to use sorted order
    remap_filter_indices!(dfs, sorted_mapping, motif_sizes)
    
    # Remap interaction summaries to use sorted order
    remapped_interaction_summaries = remap_interaction_summaries(interaction_summaries, sorted_mapping, motif_sizes)

    group_ids = [motif_names[min(size-1, 4)] for size in motif_sizes]
    button_texts = ["$(size)-motifs" for size in motif_sizes]

    for (index, (motif_size, group_id, button_text)) in enumerate(zip(motif_sizes, group_ids, button_texts))
        @info "Processing multi-motifs of size: $(motif_size)"
        interaction_summary = remapped_interaction_summaries === nothing ? nothing : remapped_interaction_summaries[index]
        @time next_idx = process_multi_motifs!(dfs, all_indices, pts,
            config, json_motifs, html_dict;             
                interaction_summary=interaction_summary,
                motif_size=motif_size, group_id=group_id, 
                button_text=button_text, start_idx=next_idx, rna=rna                
                )
    end

    # Generate combined panel figure 
    data_pairs = [
        (pts_test.predictions, pts_test.labels, "Predictions", "Labels", "Predictions vs Labels"),
        (pts_test.proc_prod, pts_test.labels, "Learned Predictions", "Labels", "Learned Predictions vs Labels"),
        (pts_test.proc_prod, pts_test.predictions, "Learned Predictions", "Predictions", "Learned Predictions vs Predictions")
    ]
    
    publication_kde_panel(data_pairs, save_path=joinpath(save_path, "generalization.png"))

    # Render generalization page (index2.html)
    render_generalization_page!(save_path; 
        page_title=page_title, 
        nav_page_count=nav_page_count,
        image_filename="generalization.png"
    )

    # ── Optional NND sensitivity analysis ───────────────────────
    if sensitivity_analysis
        ds_name = dataset_name === nothing ? page_title : dataset_name
        run_nnd_sensitivity_analysis(
            contributions_df_filtered_singletons, dfs, all_indices, pts, motif_sizes;
            save_path=save_path, page_title=ds_name
        )
    end

    render_and_save_outputs!(json_motifs, html_dict, 1; 
        html_template=html_template_unified, 
        script_template=script_template,
        css_template=template_css,
        nav_page_count=nav_page_count,
        sequence_paths=[""],
        page_title=page_title,
        save_path=save_path, 
        enable_colored_borders = enable_colored_borders,
        use_unified=use_unified
        )
end