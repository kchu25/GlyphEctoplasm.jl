"""
Multi-motif processing for convolution-based analysis (pairs, triplets, etc.).
"""

"""
    process_and_register_multi!(json_motifs, html_dict, mode_str, idx, k, d_key, pfm, flat_windows, highlight_region, median_val, count_val, banzhafs, config; kwargs...)

Process one multi-motif variant (pair/triplet/etc): save files with highlighting and register in JSON.
Uses config for rendering parameters (dpi, alpha, use_rna, xlim, filter_len).
"""
function process_and_register_multi!(json_motifs, html_dict, mode_str, idx, k, d_key, 
        pfm, flat_windows, highlight_region, median_val, count_val, banzhafs, config::ConvMotifConfig;
        save_folder_motif, motif_type_subdir, relaxed_median)
    
    d_str = get_d_str(d_key)
    paths = build_motif_paths(d_str, save_folder_motif, motif_type_subdir)

    # Save logo and influence plot (using config parameters)
    save_motif_logo(pfm, paths.png.abs, median_val; dpi=config.dpi, alpha=config.alpha, highlighted_regions=highlight_region)
    save_influence_plot(banzhafs, paths.influence.abs; highlighted_regions=highlight_region, xlim=config.xlim)
    
    # Save positional info
    save_positional_info(flat_windows, paths, config.filter_len)
    
    # Save MEME file
    save_as_meme(pfm, paths.meme.abs)
    
    # Build metadata texts
    texts = build_metadata_texts(pfm, paths, median_val, count_val; 
                                use_rna=config.use_rna, relaxed_median=relaxed_median)
    
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
function process_multi_motifs!(df, config::ConvMotifConfig, json_motifs, html_dict;
        motif_size::Int = 2,
        motif_type::String = "pair_motifs",
        save_folder = nothing,
        group_id::String = "",
        button_text::String = "Multi-Motifs",
        start_idx::Int = 1
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

    # Build mode prefix with group_id
    mode_prefix = isempty(group_id) ? "mode_" : "mode_$(group_id)_"

    @showprogress for (i, k) in enumerate(sorted_keys)
        idx = start_idx + i - 1
        mode_str = mode_prefix * string(idx)
        k_mode_str = get_k_mode_str(k)
        save_folder_motif = joinpath(save_folder, k_mode_str)
        mkpath(save_folder_motif)
        
        # Calculate relaxed median
        relaxed_median_val = median(list_of_banzhafs[k])
        
        save_box_scatter_distance_relaxed(list_of_banzhafs[k], 
            joinpath(save_folder_motif, "influence_relaxed.png"), WebDisplayMode; xlim=config.xlim)

        gdf_by_dsyms = groupby(gdf_by_msyms[k], build_grouping_columns(:distances; motif_size=motif_size))
        counts_map = Dict(k2 => nrow(gdf_by_dsyms[k2]) for k2 in keys(gdf_by_dsyms))
        count_matrices, highlighted_regions = 
            build_count_matrices_and_highlight(gdf_by_dsyms, config.data, motif_size, config.filter_len, config.float_type)

        ensure_mode_entry!(json_motifs, mode_str)

        # Process all distance variants
        sorted_dkeys = sort(collect(keys(count_matrices)), by = distance_key_value)
        list_of_banzhafs_here = Dict(d_key => gdf_by_dsyms[d_key].banzhaf for d_key in sorted_dkeys)

        for d_key in sorted_dkeys
            pfm = normalize_countmat(count_matrices[d_key])
            flat_windows = build_motif_windows(gdf_by_dsyms[d_key], motif_size, config.filter_len)
            median_here = median(gdf_by_dsyms[d_key].banzhaf)
            
            process_and_register_multi!(json_motifs, html_dict, mode_str, idx, k, d_key,
                pfm, flat_windows, highlighted_regions[d_key], median_here, counts_map[d_key], 
                list_of_banzhafs_here[d_key], config;
                save_folder_motif=save_folder_motif, 
                motif_type_subdir=joinpath(motif_type, k_mode_str),
                relaxed_median=relaxed_median_val)
        end
        
        # Populate HTML dict with first variant
        filter_indices_str = get_filter_indices_str(k)
        populate_html_dict!(html_dict, idx, json_motifs[mode_str], filter_indices_str, relaxed_median_val, group_id, button_text)
    end
    
    return start_idx + length(sorted_keys)
end

"""
    process_multi_motifs!(df, data, json_motifs, html_dict; kwargs...)

Legacy interface for backward compatibility. Creates a temporary config and calls the main function.
Prefer using the config-based interface for new code.
"""
function process_multi_motifs!(df, data, json_motifs, html_dict;
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
    
    return process_multi_motifs!(df, config, json_motifs, html_dict;
                                motif_size=motif_size, motif_type=motif_type,
                                save_folder=save_folder, group_id=group_id,
                                button_text=button_text, start_idx=start_idx)
end

export process_multi_motifs!, process_and_register_multi!
