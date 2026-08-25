"""
Multi-motif processing for convolution-based analysis (pairs, triplets, etc.).
"""

"""
    select_distance_configs(sorted_dkeys, counts_map; min_points=NUM_minimum_pts_FOR_BOXPLOT, max_configs=NUM_CONFIG_WITH_DISTANCES)

Filter distance configurations by minimum sample size and select evenly-spaced subset.

# Arguments
- `sorted_dkeys`: Sorted vector of distance configuration keys (NamedTuples)
- `counts_map`: Dict mapping each distance key to its sample count

# Keyword Arguments
- `min_points`: Minimum number of data points required (default: 25)
- `max_configs`: Maximum number of configurations to select (default: 5)

# Returns
- Vector of selected distance keys, evenly spaced across filtered candidates
"""
function select_distance_configs(sorted_dkeys, counts_map; 
        min_points=NUM_minimum_pts_FOR_BOXPLOT, 
        max_configs=NUM_CONFIG_WITH_DISTANCES)
    # If no configurations meet threshold, select top ones by count
    filtered = filter(d_key -> counts_map[d_key] >= min_points, sorted_dkeys)
    
    if isempty(filtered) && !isempty(sorted_dkeys)
        # Fall back to top configurations by count
        sorted_by_count = sort(sorted_dkeys, by=d_key -> counts_map[d_key], rev=true)
        filtered = sorted_by_count[1:min(max_configs, length(sorted_by_count))]
    end
    
    # Select evenly-spaced subset
    n = length(filtered)
    if n <= max_configs
        return filtered
    else
        return [filtered[round(Int, 1 + (i-1)*(n-1)/(max_configs-1))] for i in 1:max_configs]
    end
end

"""
    process_and_register_multi!(json_motifs, html_dict, mode_str, idx, k, d_key, pfm, flat_windows, highlight_region, median_val, count_val, banzhafs, config; kwargs...)

Process one multi-motif variant (pair/triplet/etc): save files with highlighting and register in JSON.
Uses config for rendering parameters (dpi, alpha, use_rna, xlim, filter_len).
"""
function process_and_register_multi!(json_motifs, html_dict, mode_str, idx, k, d_key, 
        pfm, flat_windows, highlight_region, median_val, count_val, banzhafs, config::ConvMotifConfig;
        save_folder_motif, motif_type_subdir, relaxed_median,
        interaction_summary_mode_str = nothing,
        rna=false, nnd_result=nothing, report_location_z::Bool=false)
    
    d_str = get_d_str(d_key)
    paths = build_motif_paths(d_str, save_folder_motif, motif_type_subdir)

    # Save logo and influence plot (using config parameters)
    save_motif_logo(pfm, paths.png.abs, median_val; dpi=config.dpi, alpha=config.alpha, highlighted_regions=highlight_region, rna=rna)
    save_influence_plot(banzhafs, paths.influence.abs; highlighted_regions=highlight_region, xlim=config.xlim)
    
    # Save positional info
    save_positional_info(flat_windows, paths, config.filter_len)
    
    # Save MEME file
    save_as_meme(pfm, paths.meme.abs)
    
    # Build metadata texts (list of strings)
    texts = build_metadata_texts(pfm, paths, median_val, count_val; 
                                interaction_summary_mode_str=interaction_summary_mode_str,
                                use_rna=config.use_rna, 
                                relaxed_median=relaxed_median, nnd_result=nnd_result,
                                report_location_z=report_location_z)
    
    label = get_descriptive_str(k, d_key)
    # Add variant without populating HTML (HTML will be populated once at the end)
    add_motif_variant!(json_motifs, mode_str, paths.png.rel, label, texts)
end

"""
    process_multi_motifs!(df, config, json_motifs, html_dict; kwargs...)

Process multi-motif modes (pairs, triplets, etc.) using a configuration struct.
Builds count matrices and highlight regions, saves motif files, and populates
JSON and HTML dicts sorted by median banzhaf contribution (descending).

# Arguments
- `df`: Vector of DataFrames organized by motif size
- `config::ConvMotifConfig`: Configuration object with all analysis parameters
- `json_motifs`: Dict to populate with saved motif metadata
- `html_dict`: Rendering dict for primary variants

# Keyword Arguments
- `motif_size::Int = 2`: Size of motif (2 for pairs, 3 for triplets, etc.)
- `motif_type::String = "pair_motifs"`: Type identifier for saving paths
- `save_folder = nothing`: Custom save folder (defaults to config.save_path/motif_type)
- `group_id::String = ""`: Namespace for this group (e.g., "pairs_pos")
- `button_text::String = "Multi-Motifs"`: Custom text for the toggle button
- `start_idx::Int = 1`: Starting index for mode numbering

# Returns
- Next available index for mode numbering
"""
function process_multi_motifs!(df, all_indices, pts, config::ConvMotifConfig, json_motifs, html_dict;
        interaction_summary = nothing,
        motif_size::Int = 2,
        motif_type::String = "pair_motifs",
        save_folder = nothing,
        group_id::String = "",
        button_text::String = "Multi-Motifs",
        start_idx::Int = 1,
        rna=false,
        report_location_z::Bool=false,
        show_region_interaction::Bool=false,  # display the interaction line? computed either way
        top_movers_out::Union{Vector{TopMoverEntry}, Nothing}=nothing
    )
    save_folder = save_folder === nothing ? joinpath(config.save_path, motif_type) : save_folder
    df_idx = motif_size - 1
    
    if df_idx < 1 || df_idx > length(df)
        error("motif_size $motif_size not available in input df vector")
    end
    
    subdf = df[df_idx]
    sep_by = build_grouping_columns(:motifs; motif_size=motif_size)
    gdf_by_msyms = groupby(subdf, sep_by)
    sorted_keys, _, _, _, list_of_banzhafs = build_sorted_keys_and_maps(gdf_by_msyms, sep_by)

    @info "Processing multi-motifs of size: $motif_size with $(length(sorted_keys)) motif groups"

    # Build mode prefix with group_id
    mode_prefix = isempty(group_id) ? "mode_" : "mode_$(group_id)_"

    processed_count = 0
    @showprogress for (i, k) in enumerate(sorted_keys)

        # Looked up regardless — it still reaches the CSV through the entry's
        # `epistasis` field below. The flag only gates the popup line.
        interaction_summary_mode_str = if !isnothing(interaction_summary)
            get(interaction_summary, k, nothing)
        else
            nothing
        end
        interaction_display_str = show_region_interaction ? interaction_summary_mode_str : nothing

        idx = start_idx + i - 1
        mode_str = mode_prefix * string(idx)
        k_mode_str = get_k_mode_str(k)
        save_folder_motif = joinpath(save_folder, k_mode_str)
        mkpath(save_folder_motif)
        
        # Calculate relaxed median
        relaxed_median_val = median(list_of_banzhafs[k])
        
        save_box_scatter_distance_relaxed(list_of_banzhafs[k], 
            joinpath(save_folder_motif, "influence_relaxed.png"), WebDisplayMode; xlim=config.xlim)

        # —————————— plot the yy kde indicator plot ———————————————
        intersect_indices = intersect(gdf_by_msyms[k].data_pt_index, all_indices)
        is_in_intersect = all_indices .∈ Ref(Set(intersect_indices))
        fig_intersect = plot_labels_vs_procprod(pts, is_in_intersect; motif_label="Contain motif")
        save(joinpath(save_folder_motif, "yy_kde_intersect.png"), fig_intersect, px_per_unit=1)

        # ————— NND Permutation Test (k=5) ————————
        subpop_pos = findall(is_in_intersect)
        nnd_base = nnd_permutation_test_1d(subpop_pos, pts.labels)
        # Cluster median (median expression of the motif's sequences): augmented
        # onto the NND result so it shows on the card and ranks the top-movers
        # page, matching the mutation case.
        cluster_median = any(is_in_intersect) ? median(@view pts.labels[is_in_intersect]) : NaN
        # Location statistic: closed-form companion to the NND test (see
        # `location_z_1d`). Computed unconditionally — it is free — but only
        # displayed when `report_location_z` is on.
        location_z = location_z_1d(subpop_pos, pts.labels)
        nnd_result = (; nnd_base..., cluster_median, location_z)

        # ——————————————————————— process d_syms variants ————————————————————————————
        gdf_by_dsyms = groupby(gdf_by_msyms[k], build_grouping_columns(:distances; motif_size=motif_size))
        counts_map = Dict(k2 => nrow(gdf_by_dsyms[k2]) for k2 in keys(gdf_by_dsyms))
        count_matrices, highlighted_regions = 
            build_count_matrices_and_highlight(gdf_by_dsyms, config.data, motif_size, config.filter_len, config.float_type)

        # Process all distance variants
        sorted_dkeys = sort(collect(keys(counts_map)), by = distance_key_value)
        selected_dkeys = select_distance_configs(sorted_dkeys, counts_map)

        # Skip this motif if no configurations pass filtering
        if isempty(selected_dkeys)
            @warn "Skipping motif $mode_str: no distance configurations with ≥$NUM_minimum_pts_FOR_BOXPLOT data points"
            continue
        end

        ensure_mode_entry!(json_motifs, mode_str)

        processed_count += 1

        list_of_banzhafs_here = Dict(d_key => gdf_by_dsyms[d_key].banzhaf for d_key in selected_dkeys)

        for d_key in selected_dkeys
            pfm = normalize_countmat(count_matrices[d_key])
            flat_windows = build_motif_windows(gdf_by_dsyms[d_key], motif_size, config.filter_len)
            median_here = median(gdf_by_dsyms[d_key].banzhaf)
            
            process_and_register_multi!(json_motifs, html_dict, mode_str, idx, k, d_key,
                pfm, flat_windows, highlighted_regions[d_key], median_here, counts_map[d_key], 
                list_of_banzhafs_here[d_key], config;
                interaction_summary_mode_str=interaction_display_str,
                save_folder_motif=save_folder_motif, 
                motif_type_subdir=joinpath(motif_type, k_mode_str),
                relaxed_median=relaxed_median_val, rna=rna, nnd_result=nnd_result,
                report_location_z=report_location_z)
        end
        
        # Populate HTML dict with first variant
        filter_indices_str = get_filter_indices_str(k)
        populate_html_dict!(html_dict, idx, json_motifs[mode_str], filter_indices_str, relaxed_median_val, group_id, button_text)

        # Collect a self-contained summary row for the top-movers landing page.
        # Uses the first displayed distance variant's logo/influence as the card,
        # the relaxed median for ranking, and the per-motif interaction summary as
        # the epistasis line. Convolution motifs carry no wild-type track.
        if top_movers_out !== nothing
            first_paths = build_motif_paths(
                get_d_str(selected_dkeys[1]), save_folder_motif, joinpath(motif_type, k_mode_str))
            push!(top_movers_out, TopMoverEntry(
                float(relaxed_median_val),          # Shapley median
                float(cluster_median),              # cluster_median
                float(nnd_result.obs_mNND),         # cluster_nnd
                length(list_of_banzhafs[k]),        # count
                filter_indices_str,                 # display_name
                button_text,                        # group_label
                "",                                 # span (n/a for convolution)
                float(nnd_result.p_value),          # nnd_p
                false,                              # is_singleton
                interaction_summary_mode_str === nothing ? "" : interaction_summary_mode_str,  # epistasis
                WildTypeRegion[],                   # wt_regions (n/a for convolution)
                first_paths.png.rel,                # img (first variant logo)
                "", "",                             # img_reduced/img_region (single-card: no dual view)
                first_paths.influence.rel,          # influence plot
                joinpath(motif_type, k_mode_str, "yy_kde_intersect.png"),  # yy_kde
                collect(String, json_motifs[mode_str][texts_str][1]),       # texts
                # Full distance-variant payload so the popup gets a slider.
                collect(String, json_motifs[mode_str][pwms_str]),           # variant_pwms
                collect(String, json_motifs[mode_str][labels_str]),         # variant_labels
                [collect(String, t) for t in json_motifs[mode_str][texts_str]],  # variant_texts
                float(location_z),                  # location_z (NaN when not computable)
            ))
        end
    end
    
    return start_idx + processed_count
end

"""
    process_multi_motifs!(df, data, json_motifs, html_dict; kwargs...)

Legacy interface for backward compatibility. Creates a temporary config and calls the main function.
Prefer using the config-based interface for new code.
"""
function process_multi_motifs!(df, all_indices, pts, data, json_motifs, html_dict;
        motif_size = 2,
        filter_len = 7,
        motif_type = "pair_motifs",
        SAVE_PATH = "tmp",
        save_folder = nothing,
        float_type = Float32,
        dpi = 65,
        alpha = 1.0,
        use_rna = false,
        start_idx = 1,
        xlim = nothing,
        group_id = "",
        button_text = "Multi-Motifs"
    )
    # Create temporary config
    config = ConvMotifConfig(data; filter_len=filter_len, float_type=float_type,
                            dpi=dpi, alpha=alpha, use_rna=use_rna, xlim=xlim,
                            save_path=save_folder === nothing ? SAVE_PATH : dirname(save_folder))
    
    return process_multi_motifs!(df, all_indices, pts, config, json_motifs, html_dict;
                                motif_size=motif_size, motif_type=motif_type,
                                save_folder=save_folder, group_id=group_id,
                                button_text=button_text, start_idx=start_idx)
end

export process_multi_motifs!, process_and_register_multi!
