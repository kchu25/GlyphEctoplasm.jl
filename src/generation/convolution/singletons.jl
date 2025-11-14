"""
Singleton motif processing for convolution-based analysis.
"""

"""
    process_and_register_singleton!(json_motifs, html_dict, idx, k, pfm, gdf_row, config; kwargs...)

Process one singleton motif: save files, build metadata, and register in JSON/HTML dicts.
Uses config for rendering parameters (dpi, alpha, use_rna, xlim, filter_len).
"""
function process_and_register_singleton!(json_motifs, html_dict, idx, k, pfm, gdf_row, config::ConvMotifConfig;
        save_folder, motif_type, median_val, count_val, banzhafs,
        mode_prefix="mode_", group_id="", button_text="Singleton Motifs")
    
    name_base = string(k.filter_index)
    paths = build_motif_paths(name_base, save_folder, motif_type)

    # Save logo and influence plot (using config parameters)
    save_motif_logo(pfm, paths.png.abs, median_val; dpi=config.dpi, alpha=config.alpha, highlighted_regions=nothing)
    save_influence_plot(banzhafs, paths.influence.abs; highlighted_regions=nothing, xlim=config.xlim)
    
    # Save positional info
    save_positional_info(gdf_row, paths, config.filter_len)
    
    # Save MEME file
    save_as_meme(pfm, paths.meme.abs)
    
    # Build metadata texts
    texts = build_metadata_texts(pfm, paths, median_val, count_val; 
                                use_rna=config.use_rna, relaxed_median=nothing)
    
    mode_str = mode_prefix * string(idx)
    label = "pattern $(k.filter_index)"
    filter_indices_str = string(k.filter_index)
    add_motif_entry!(json_motifs, html_dict, mode_str, paths.png.rel, label, texts, idx, filter_indices_str, median_val, group_id, button_text)
end

"""
    process_singletons!(contributions_df, config, json_motifs, html_dict; kwargs...)

Process all singleton motifs using a configuration struct.
Populates JSON and HTML dicts sorted by median banzhaf contribution (descending).

# Arguments
- `contributions_df`: DataFrame with contribution data
- `config::ConvMotifConfig`: Configuration object with all analysis parameters
- `json_motifs`: JSON motif dictionary
- `html_dict`: HTML dictionary

# Keyword Arguments
- `motif_type::String = "singletons"`: Type identifier for saving paths
- `save_folder = nothing`: Custom save folder (defaults to config.save_path/motif_type)
- `group_id::String = ""`: Namespace for this group (e.g., "high_sing")
- `button_text::String = "Singleton Motifs"`: Custom text for the toggle button
- `start_idx::Int = 1`: Starting index for mode numbering
- `pareto_rank = nothing`: Optional Pareto rank filter

# Returns
- Next available index for mode numbering
"""
function process_singletons!(contributions_df, config::ConvMotifConfig, json_motifs, html_dict;
        motif_type::String = "singletons",
        save_folder = nothing,
        group_id::String = "",
        button_text::String = "Singleton Motifs",
        start_idx::Int = 1,
        pareto_rank = nothing
    )
    save_folder = save_folder === nothing ? joinpath(config.save_path, motif_type) : save_folder
    mkpath(save_folder)
    
    sep_by = build_grouping_columns(:filter_index)
    gdf_filters = groupby(contributions_df, sep_by)
    sorted_keys, median_map, _, count_map, list_of_banzhafs =
        build_sorted_keys_and_maps(gdf_filters, sep_by; pareto_rank=pareto_rank)
    count_matrices = build_singleton_count_matrices(gdf_filters, config.data, config.filter_len, config.float_type)

    # Build mode prefix with group_id
    mode_prefix = isempty(group_id) ? "mode_" : "mode_$(group_id)_"

    for (i, k) in enumerate(sorted_keys)
        idx = start_idx + i - 1
        pfm = normalize_countmat(count_matrices[k])
        
        process_and_register_singleton!(json_motifs, html_dict, idx, k, pfm, gdf_filters[k], config;
            save_folder=save_folder, motif_type=motif_type, 
            median_val=median_map[k], count_val=count_map[k], banzhafs=list_of_banzhafs[k],
            mode_prefix=mode_prefix, group_id=group_id, button_text=button_text)
    end
    
    # Return next available index
    return start_idx + length(sorted_keys)
end

"""
    process_singletons!(contributions_df, data, json_motifs, html_dict; kwargs...)

Legacy interface for backward compatibility. Creates a temporary config and calls the main function.
Prefer using the config-based interface for new code.
"""
function process_singletons!(contributions_df, data, json_motifs, html_dict;
        SAVE_PATH = "tmp",
        motif_type = "singletons",
        save_folder = nothing,
        dpi = 65,
        alpha = 1.0,
        use_rna = false,
        filter_len = 7,
        float_type = Float32,
        xlim = nothing,
        group_id = "",
        button_text = "Singleton Motifs",
        start_idx = 1,
        pareto_rank = nothing
    )
    # Create temporary config
    config = ConvMotifConfig(data; filter_len=filter_len, float_type=float_type,
                            dpi=dpi, alpha=alpha, use_rna=use_rna, xlim=xlim,
                            save_path=save_folder === nothing ? SAVE_PATH : dirname(save_folder))
    
    return process_singletons!(contributions_df, config, json_motifs, html_dict;
                              motif_type=motif_type, save_folder=save_folder,
                              group_id=group_id, button_text=button_text,
                              start_idx=start_idx, pareto_rank=pareto_rank)
end

export process_singletons!, process_and_register_singleton!
