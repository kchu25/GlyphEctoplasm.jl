"""
Singleton motif processing for convolution-based analysis.
"""

"""
    process_and_register_singleton!(json_motifs, html_dict, idx, k, pfm, gdf_row, config; kwargs...)

Process one singleton motif: save files, build metadata, and register in JSON/HTML dicts.
Uses config for rendering parameters (dpi, alpha, use_rna, xlim, filter_len).
"""
function process_and_register_singleton!(all_indices, pts, json_motifs, html_dict, idx, k, pfm, gdf_row,      config::ConvMotifConfig;
        save_folder, motif_type, median_val, count_val, banzhafs,
        mode_prefix="mode_", group_id="", button_text="Singleton Motifs", 
        rna=false, is_significant::Bool=true,
        display_index::Union{Int, Nothing}=nothing,
        mad_all_labels::Float64=NaN
        )

    # Use display_index for labels/filenames if provided, otherwise use original k.filter_index
    shown_index = display_index === nothing ? k.filter_index : display_index
    
    name_base = string(shown_index)
    paths = build_motif_paths(name_base, save_folder, motif_type)

    # Save logo and influence plot (using config parameters)
    save_motif_logo(pfm, paths.png.abs, median_val; dpi=config.dpi, alpha=config.alpha, highlighted_regions=nothing, rna=rna)
    save_influence_plot(banzhafs, paths.influence.abs; highlighted_regions=nothing, xlim=config.xlim)
    
    # Save positional info
    save_positional_info(gdf_row, paths, config.filter_len)
    
    # Save MEME file
    save_as_meme(pfm, paths.meme.abs)
    
    # —————————— plot the yy kde indicator plot ————————————————
    intersect_indices = intersect(gdf_row.data_pt_index, all_indices)
    is_in_intersect = all_indices .∈ Ref(Set(intersect_indices))
    fig_intersect = plot_labels_vs_procprod(pts, is_in_intersect; motif_label="Contain motif")
    save(
        joinpath(save_folder, "yy_kde_intersect_$shown_index.png"), 
        fig_intersect, px_per_unit=1)

    # ————— Fligner-Killeen Test for variance comparison ————————
    fk_test_result = FlignerKilleenTest(
            pts.labels, 
            (@view pts.labels[intersect_indices]))
    fk_pvalue = round(pvalue(fk_test_result); digits=4)

    # ————— MAD log ratio (tightness measure) ————————
    subset_labels = @view pts.labels[intersect_indices]
    mad_subset = median(abs.(subset_labels .- median(subset_labels)))
    mad_log_ratio = log((mad_subset + MAD_EPSILON) / (mad_all_labels + MAD_EPSILON))

    # ——————————————————————— end ————————————————————————————

    # Build metadata texts
    texts = build_metadata_texts(pfm, paths, median_val, count_val; 
                                use_rna=config.use_rna, relaxed_median=nothing, 
                                fk_pvalue=fk_pvalue, mad_log_ratio=mad_log_ratio)
    
    mode_str = mode_prefix * string(idx)
    label = "motif $(shown_index)"
    filter_indices_str = string(shown_index)
    add_motif_entry!(json_motifs, html_dict, mode_str, paths.png.rel, label, texts, idx, filter_indices_str, median_val, group_id, button_text; is_significant=is_significant)
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

# Keyword Argumentsgdf_row
- `motif_type::String = "singletons"`: Type identifier for saving paths
- `save_folder = nothing`: Custom save folder (defaults to config.save_path/motif_type)
- `group_id::String = ""`: Namespace for this group (e.g., "high_sing")
- `button_text::String = "Singleton Motifs"`: Custom text for the toggle button
- `start_idx::Int = 1`: Starting index for mode numbering
- `pareto_rank = nothing`: Optional Pareto rank filter

# Returns
- `(next_idx, sorted_mapping)`: Tuple of next available index and Dict mapping original filter_index to sorted order
"""
function process_singletons!(contributions_df, all_indices, pts, config::ConvMotifConfig, json_motifs, html_dict;
        motif_type::String = "singletons",
        save_folder = nothing,
        group_id::String = "",
        button_text::String = "Singleton Motifs",
        start_idx::Int = 1,
        pareto_rank = nothing,
        rna=false,
        get_sorted_mapping_only=false
    )
    save_folder = save_folder === nothing ? joinpath(config.save_path, motif_type) : save_folder
    mkpath(save_folder)
    
    sep_by = build_grouping_columns(:filter_index)
    gdf_filters = groupby(contributions_df, sep_by)
    sorted_keys, median_map, _, count_map, list_of_banzhafs =
        build_sorted_keys_and_maps(gdf_filters, sep_by; pareto_rank=pareto_rank)

    # Build sorted_mapping: original filter_index => sorted order (1-based)
    sorted_mapping = Dict{Int, Int}()
    for (i, k) in enumerate(sorted_keys)
        sorted_mapping[k.filter_index] = i
    end

    if get_sorted_mapping_only
        # do the remapping inplace
        contributions_df.filter_index = [sorted_mapping[idx] for idx in contributions_df.filter_index]
        return nothing, sorted_mapping
    end

    count_matrices = build_singleton_count_matrices(gdf_filters, config.data, config.filter_len, config.float_type)

    # Build mode prefix with group_id
    mode_prefix = isempty(group_id) ? "mode_" : "mode_$(group_id)_"

    # Check if significant column exists
    has_significant_col = hasproperty(contributions_df, :significant)

    # Precompute MAD of all labels (used for log ratio in each motif)
    mad_all_labels = Float64(median(abs.(pts.labels .- median(pts.labels))))

    for (i, k) in enumerate(sorted_keys)
        idx = start_idx + i - 1
        pfm = normalize_countmat(count_matrices[k])
        
        # Extract significance from first row of group (all rows in group have same value)
        is_significant = has_significant_col ? first(gdf_filters[k].significant) : true

        # Use sorted index (i) as display_index for consistent numbering
        process_and_register_singleton!(all_indices, pts, 
            json_motifs, html_dict, idx, k, pfm, gdf_filters[k], config;
            save_folder=save_folder, motif_type=motif_type, 
            median_val=median_map[k], count_val=count_map[k], banzhafs=list_of_banzhafs[k],
            mode_prefix=mode_prefix, group_id=group_id, button_text=button_text, rna=rna,
            is_significant=is_significant, display_index=i, mad_all_labels=mad_all_labels)
    end
    
    # Remap filter_index column in contributions_df to use sorted order (after processing)
    contributions_df.filter_index = [sorted_mapping[idx] for idx in contributions_df.filter_index]

    # Return next available index and the sorted mapping
    return start_idx + length(sorted_keys), sorted_mapping
end

"""
    process_singletons!(contributions_df, data, json_motifs, html_dict; kwargs...)

Legacy interface for backward compatibility. Creates a temporary config and calls the main function.
Prefer using the config-based interface for new code.
"""
function process_singletons!(contributions_df, all_indices, pts, data, json_motifs, html_dict;
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
    
    return process_singletons!(contributions_df, all_indices, pts, config, json_motifs, html_dict;
                              motif_type=motif_type, save_folder=save_folder,
                              group_id=group_id, button_text=button_text,
                              start_idx=start_idx, pareto_rank=pareto_rank)
end

"""
    remap_filter_indices!(dfs, sorted_mapping, motif_sizes)

Remap filter indices in multi-motif DataFrames using the sorted mapping from singletons.
Modifies the m1, m2, m3, ... columns in each DataFrame to use the new sorted indices.

# Arguments
- `dfs`: Vector of DataFrames containing multi-motif data
- `sorted_mapping`: Dict mapping original filter_index to sorted order
- `motif_sizes`: Vector of motif sizes corresponding to each DataFrame
"""
function remap_filter_indices!(dfs::Vector, sorted_mapping::Dict{Int, Int}, motif_sizes::Vector)
    for (df, motif_size) in zip(dfs, motif_sizes)
        m_syms = m_symbols(motif_size)
        for m_sym in m_syms
            if hasproperty(df, m_sym)
                df[!, m_sym] = [get(sorted_mapping, idx, idx) for idx in df[!, m_sym]]
            end
        end
    end
end

"""
    remap_interaction_summaries(interaction_summaries, sorted_mapping, motif_sizes)

Remap filter indices in interaction summaries using the sorted mapping from singletons.
Returns a new vector of dictionaries with remapped keys.

# Arguments
- `interaction_summaries`: Vector of Dicts mapping NamedTuples (m1=x, m2=y, ...) to interaction stats
- `sorted_mapping`: Dict mapping original filter_index to sorted order
- `motif_sizes`: Vector of motif sizes corresponding to each interaction summary

# Returns
- New vector of dictionaries with remapped keys (same types as input)
"""
function remap_interaction_summaries(interaction_summaries, sorted_mapping::Dict{Int, Int}, motif_sizes::Vector)
    interaction_summaries === nothing && return nothing
    
    remapped = similar(interaction_summaries)
    for (i, (summary, motif_size)) in enumerate(zip(interaction_summaries, motif_sizes))
        m_syms = m_symbols(motif_size)
        # Create new dict with same type as original
        new_summary = empty(summary)
        for (key, value) in summary
            # Remap each filter index in the NamedTuple key
            # Use tuple() and typeof(key) to preserve exact key type for Dict lookup
            new_key_values = tuple((get(sorted_mapping, Int(getfield(key, sym)), Int(getfield(key, sym))) for sym in m_syms)...)
            new_key = typeof(key)(new_key_values)
            new_summary[new_key] = value
        end
        remapped[i] = new_summary
    end
    return remapped
end

"""
    compute_sorted_mapping_and_remap!(contributions_df, dfs, motif_sizes;
        interaction_summaries=nothing, pareto_rank=nothing)

Compute the sorted mapping from singleton filter indices to their sorted order
(sorted by median Banzhaf contribution, descending), and remap all DataFrames
and interaction summaries in-place.

This is a lightweight alternative to `process_singletons!` that only performs
the sorting and remapping logic without any rendering, saving, or registration.

# Arguments
- `contributions_df`: DataFrame with `:filter_index` and `:banzhaf` columns (singletons)
- `dfs`: Vector of DataFrames containing multi-motif data (remapped in-place)
- `motif_sizes`: Vector of motif sizes corresponding to each DataFrame in `dfs`
- `interaction_summaries`: Optional vector of interaction summary Dicts (default: `nothing`)
- `pareto_rank`: Optional Pareto rank filter (default: `nothing`)

# Returns
- `(sorted_mapping, remapped_interaction_summaries)` where:
  - `sorted_mapping::Dict{Int,Int}`: original filter_index → 1-based sorted order
  - `remapped_interaction_summaries`: remapped interaction summaries (or `nothing`)

# Side effects
- Mutates `contributions_df.filter_index` to use sorted order
- Mutates columns in `dfs` via `remap_filter_indices!`
"""
function compute_sorted_mapping_and_remap!(contributions_df, dfs, motif_sizes;
        interaction_summaries=nothing, pareto_rank=nothing)
    # Group by filter_index and compute sorted keys
    sep_by = build_grouping_columns(:filter_index)
    gdf_filters = groupby(contributions_df, sep_by)
    sorted_keys, _, _, _, _ =
        build_sorted_keys_and_maps(gdf_filters, sep_by; pareto_rank=pareto_rank)

    # Build sorted_mapping: original filter_index => sorted order (1-based)
    sorted_mapping = Dict{Int, Int}()
    for (i, k) in enumerate(sorted_keys)
        sorted_mapping[k.filter_index] = i
    end

    # Remap contributions_df.filter_index in-place
    contributions_df.filter_index = [sorted_mapping[idx] for idx in contributions_df.filter_index]

    # Remap multi-motif DataFrames in-place
    remap_filter_indices!(dfs, sorted_mapping, motif_sizes)

    # Remap interaction summaries
    remapped_interaction_summaries = remap_interaction_summaries(interaction_summaries, sorted_mapping, motif_sizes)

    return sorted_mapping, remapped_interaction_summaries
end

export process_singletons!, process_and_register_singleton!, remap_filter_indices!, remap_interaction_summaries, compute_sorted_mapping_and_remap!
