"""
Multi-region mutation analysis: collecting, sorting, and registering mutation region motifs.
"""

"""
    prepare_and_collect_mutation_metadata(contributions_df, df_multi_motifs, data, config; kwargs...)

High-level workflow function that handles the complete mutation analysis pipeline:
1. Filters and prepares singleton mutations from contributions_df
2. Filters multi-motif dataframes (pairs, triplets, etc.) by most_common_length_indices
3. Collects metadata for all mutation types

# Arguments
- `contributions_df`: DataFrame with filter-level contributions (for singletons)
- `df_multi_motifs`: Vector of DataFrames indexed by motif_size-1 (e.g., df[1] = pairs, df[2] = triplets)
- `data`: OnehotSEQ2EXP_Dataset
- `config::MutationRegionConfig`: Configuration for analysis

# Keyword Arguments
- `singleton_filter_pareto_rank::Int = 1`: Pareto rank threshold for filtering singleton mutations (only applies to singletons)
- `split_by_sign::Bool = true`: Whether to split contributions by sign when computing Pareto ranks during singleton filtering.
   Note: Final global sorting always splits by sign (positives first, then negatives).
- `motif_sizes::Union{Nothing, Vector{Int}} = nothing`: Which motif sizes to process. 
   If `nothing` (default), automatically infers: [1, 2, ..., length(df_multi_motifs)+1]
- `motif_type_prefix::String = "mutation_regions"`: Prefix for motif type names

# Returns
- Vector of vectors of MotifMetadata, one vector per motif size (in order: singletons, pairs, triplets, ...)

# Example
```julia
config = MutationRegionConfig(data; filter_len=9, xlim=xlim, save_path="tmp")

# Automatic inference (processes all available motif sizes)
all_metadata = prepare_and_collect_mutation_metadata(
    contributions_df, df, data, config
)

# Or specify singleton filtering explicitly
all_metadata = prepare_and_collect_mutation_metadata(
    contributions_df, df, data, config;
    singleton_filter_pareto_rank=2  # Only for singletons
)
```
"""
function prepare_and_collect_mutation_metadata(
    contributions_df, df_multi_motifs, data, config::MutationRegionConfig;
    singleton_filter_pareto_rank::Int = 1,
    split_by_sign::Bool = true,
    motif_sizes::Union{Nothing, Vector{Int}} = nothing,
    motif_type_prefix::String = "mutation_regions"
)
    # Infer motif_sizes if not provided
    # contributions_df → motif_size 1
    # df_multi_motifs[1] → motif_size 2, df_multi_motifs[2] → motif_size 3, etc.
    if motif_sizes === nothing
        motif_sizes = [1; collect(2:(length(df_multi_motifs) + 1))]
    end
    
    # Pre-allocate with proper type
    all_metadata = Vector{Vector{MotifMetadata}}()
    
    for motif_size in motif_sizes
        if motif_size == 1
            # Singletons: filter and rename contributions
            df_filtered = filter_and_rename_for_mutagenesis(
                contributions_df, data;
                pareto_rank=singleton_filter_pareto_rank,
                split_by_sign=split_by_sign
            )
            
            motif_type = "$(motif_type_prefix)_1"
        else
            # Multi-motifs: filter by most_common_length_indices
            df_idx = motif_size - 1
            if df_idx > length(df_multi_motifs)
                @warn "Skipping motif_size=$motif_size: not available in df_multi_motifs"
                continue
            end
            
            if isnothing(data.raw_data.most_common_length_indices)
                df_filtered = df_multi_motifs[df_idx]
            else
                df_filtered = filter(
                    x -> x.data_pt_index ∈ data.raw_data.most_common_length_indices,
                    df_multi_motifs[df_idx]
                )
            end

            motif_type = "$(motif_type_prefix)_$(motif_size)"
        end
        
        # Collect metadata
        metadata = collect_mutation_region_metadata(
            df_filtered, config;
            motif_type=motif_type,
            motif_size=motif_size
        )
        
        push!(all_metadata, metadata)
    end
    
    return all_metadata
end

"""
    build_mutation_aggregates(df_mutated, config, motif_size)

Aggregate mutation data: grouping, sorting, and matrix construction.
Returns a NamedTuple with all computed aggregates.
"""
function build_mutation_aggregates(df_mutated, config::MutationRegionConfig, motif_size::Int)
    sep_by = build_grouping_columns(:mutagenesis; motif_size=motif_size)
    gdf = groupby(df_mutated, sep_by)
    sorted_keys, median_map, mean_map, count_map, list_of_banzhafs = 
        build_sorted_keys_and_maps(gdf, sep_by)
    
    count_matrices_vec, reference_matrices_vec, adjusted_positions_vec = 
        obtain_count_and_reference_matrices(
            sorted_keys, gdf, config.data, config.reference_seq, 
            build_grouping_columns(:motif_positions; motif_size=motif_size), 
            config.filter_len; T=config.float_type, off_region_search=config.off_region_search)
    
    return (
        gdf=gdf, sorted_keys=sorted_keys, median_map=median_map, mean_map=mean_map,
        count_map=count_map, list_of_banzhafs=list_of_banzhafs,
        count_matrices_vec=count_matrices_vec, reference_matrices_vec=reference_matrices_vec,
        adjusted_positions_vec=adjusted_positions_vec
    )
end

"""
    compute_fragment_info(count_mats, ref_pfms, start_positions, reduction_on_ref, motif_size)

Compute fragment count and span for a single motif key.
Returns a NamedTuple with fragment_count, span, group_id, and button_text.

Logic:
- With reduction_on_ref=true: Always use EntroPlots.count_fragments to detect actual fragments
  (can reduce triplets to 2 regions, pairs to 1 region, etc. if some don't match reference)
- Without reduction_on_ref=false: fragment_count = motif_size for multi-motifs, count matrices for singletons
"""
function compute_fragment_info(count_mats, ref_pfms, start_positions, reduction_on_ref::Bool, motif_size::Int)
    if reduction_on_ref
        # Use EntroPlots to compute fragments and span with reference filtering
        # This can reduce the fragment count (e.g., triplets → 2 regions if one is filtered out)
        fragment_count, span_str = EntroPlots.count_fragments(count_mats, ref_pfms, start_positions)
        # Convert dashes to colons for consistency: "36-37, 39-45" → "36:37, 39:45"
        span_str = replace(span_str, "-" => ":")
    else
        # No reference filtering
        if motif_size > 1
            # Multi-motifs: fragment_count = motif_size (pairs=2, triplets=3, etc.)
            fragment_count = motif_size
            # Build span from matrices
            if length(count_mats) == 1
                # Merged into 1 matrix (adjacent merging happened)
                mat_width = size(count_mats[1], 2)
                span_str = "$(start_positions[1]):$(start_positions[1] + mat_width - 1)"
            else
                # Multiple matrices - build span for each
                spans = ["$(pos):$(pos + size(mat, 2) - 1)" for (pos, mat) in zip(start_positions, count_mats)]
                span_str = join(spans, ", ")
            end
        else
            # Singletons without reduction: fragment_count = number of matrices
            fragment_count = length(count_mats)
            spans = ["$(pos):$(pos + size(mat, 2) - 1)" for (pos, mat) in zip(start_positions, count_mats)]
            span_str = join(spans, ", ")
        end
    end
    
    group_id = fragment_count == 1 ? "single_region" : "$(fragment_count)_regions"
    button_text = fragment_count == 1 ? "Single Regions" : "$(fragment_count) Regions"
    
    return (fragment_count=fragment_count, span=span_str, group_id=group_id, button_text=button_text)
end

"""
    collect_mutation_region_metadata(df_mutated, config; kwargs...)

Collect motif metadata for mutation regions using a configuration struct.
Returns a vector of MotifMetadata objects for later sorting and registration.

# Arguments
- `df_mutated`: DataFrame with mutation data
- `config::MutationRegionConfig`: Configuration object with all analysis parameters

# Keyword Arguments
- `motif_type::String = "pair_mutation_regions"`: Type identifier for saving paths
- `motif_size::Int = 2`: Size of motif pairs (1 for singletons, 2 for pairs, etc.)
- `save_folder = nothing`: Custom save folder (defaults to config.save_path/motif_type)

# Returns
- Vector of MotifMetadata objects with all information for rendering
"""
function collect_mutation_region_metadata(df_mutated, config::MutationRegionConfig;
        motif_type::String = "pair_mutation_regions",
        motif_size::Int = 2,
        save_folder = nothing
    )
    save_folder = save_folder === nothing ? joinpath(config.save_path, motif_type) : save_folder
    
    # Step 1: Aggregate data
    agg = build_mutation_aggregates(df_mutated, config, motif_size)
    
    # Step 2: Process each key and build metadata
    metadata = MotifMetadata[]
    for k in agg.sorted_keys
        count_mats = agg.count_matrices_vec[k]
        start_positions = agg.adjusted_positions_vec[k]
        (isempty(count_mats) || any(isempty, count_mats)) && continue
        
        # Build motif data
        motif_data = MotifData(
            k, count_mats, start_positions,
            BitMatrix.(agg.reference_matrices_vec[k]),
            agg.median_map[k], agg.count_map[k],
            agg.list_of_banzhafs[k], agg.gdf[k]
        )
        
        # Compute fragment info
        fragment_info = compute_fragment_info(
            motif_data.count_matrices, 
            motif_data.references, 
            motif_data.positions, 
            config.reduction_on_ref,
            motif_size  # Pass motif_size to distinguish singletons from multi-motifs
        )
        
        # Build complete metadata
        push!(metadata, MotifMetadata(motif_data, config, motif_type, save_folder, motif_size, fragment_info))
    end
    
    return metadata
end

"""
    collect_mutation_region_metadata(df_mutated, data; kwargs...)

Legacy interface for backward compatibility. Creates a temporary config and calls the main function.
Prefer using the config-based interface for new code.

# TODO: Is this used? remove this?
"""
function collect_mutation_region_metadata(df_mutated, data;
        motif_type = "pair_mutation_regions",
        save_folder = nothing,
        filter_len = 9,
        float_type = Float32,
        reference_seq = nothing,
        off_region_search = true,
        motif_size = 2,
        total_length = nothing,
        group_id = "",
        button_text = "Paired Mutation Regions",
        dpi = 65,
        use_rna = false,
        xlim = nothing,
        reduction_on_ref = true
    )
    # Create temporary config
    ref_seq = reference_seq === nothing ? 
        SEQ2EXPdata.consensus_to_bitmatrix_auto(data.raw_data.consensus) : reference_seq
    tot_len = total_length === nothing ? length(data.raw_data.consensus) : total_length
    
    config = MutationRegionConfig(
        data, ref_seq, tot_len, filter_len, float_type,
        off_region_search, reduction_on_ref, dpi, use_rna, xlim, 
        save_folder === nothing ? SAVE_PATH : dirname(save_folder)
    )
    
    return collect_mutation_region_metadata(df_mutated, config; 
                                           motif_type=motif_type, 
                                           motif_size=motif_size, 
                                           save_folder=save_folder)
end

"""
    sort_by_group_and_pareto(metadata_list)

Sort motifs by sign (positive first), then by group_id, then by Pareto rank within each group.
Memory-efficient: computes Pareto ranks in-place per group without excessive allocations.

Returns sorted vector of metadata.
"""
function sort_by_group_and_pareto(metadata_list)
    n = length(metadata_list)
    n == 0 && return metadata_list
    
    # Create augmented entries with sort keys (using NamedTuple wrapper)
    augmented = Vector{NamedTuple}(undef, n)
    for (i, m) in enumerate(metadata_list)
        group_order = m.group_id == "single_region" ? 0 : parse(Int, split(m.group_id, '_')[1])
        sign_order = m.median > 0 ? 0 : 1
        augmented[i] = (
            metadata=m, 
            group_order=group_order, 
            sign_order=sign_order,
            pareto_rank=0  # Will be set later
        )
    end
    
    # Sort by sign and group (in-place, O(n log n))
    sort!(augmented, by = a -> (a.sign_order, a.group_order))
    
    # Compute Pareto ranks within each (sign, group) partition
    i = 1
    while i <= n
        # Find end of current (sign, group) partition
        current_sign = augmented[i].sign_order
        current_group = augmented[i].group_order
        j = i
        while j <= n && 
              augmented[j].sign_order == current_sign && 
              augmented[j].group_order == current_group
            j += 1
        end
        
        # Compute Pareto ranks for this partition [i:j-1]
        partition_ranks = compute_pareto_ranks_subset_wrapped(augmented, i, j-1)
        for (idx, rank) in enumerate(partition_ranks)
            # Update pareto_rank in the NamedTuple
            old = augmented[i+idx-1]
            augmented[i+idx-1] = (
                metadata=old.metadata,
                group_order=old.group_order,
                sign_order=old.sign_order,
                pareto_rank=rank
            )
        end
        
        i = j
    end
    
    # Final sort by (sign, group, pareto_rank, abs(median))
    # For positives: higher abs(median) first (-abs for descending)
    # For negatives: lower abs(median) first (+abs for ascending)
    sort!(augmented, by = a -> (
        a.sign_order, 
        a.group_order, 
        a.pareto_rank, 
        a.sign_order == 0 ? -abs(a.metadata.median) : abs(a.metadata.median)
    ))
    
    # Extract just the metadata (discard sort keys)
    return [a.metadata for a in augmented]
end

"""
    compute_pareto_ranks_subset_wrapped(augmented_list, start_idx, end_idx)

Compute Pareto ranks for a subset of augmented metadata [start_idx:end_idx].
Works with NamedTuple wrappers that have .metadata field.
Returns vector of ranks.
"""
function compute_pareto_ranks_subset_wrapped(augmented_list, start_idx, end_idx)
    n = end_idx - start_idx + 1
    n == 0 && return Int[]
    
    # Extract objectives for this subset (from wrapped metadata)
    objectives = [(abs(augmented_list[i].metadata.median), augmented_list[i].metadata.count) 
                  for i in start_idx:end_idx]
    
    # Compute Pareto ranks
    ranks = zeros(Int, n)
    available = trues(n)
    
    current_rank = 1
    while any(available)
        for i in 1:n
            available[i] || continue
            
            is_dominated = false
            for j in 1:n
                (i == j || !available[j]) && continue
                
                # Check if j dominates i
                dominates = (objectives[j][1] >= objectives[i][1] && objectives[j][2] >= objectives[i][2]) &&
                           (objectives[j][1] > objectives[i][1] || objectives[j][2] > objectives[i][2])
                
                if dominates
                    is_dominated = true
                    break
                end
            end
            
            if !is_dominated
                ranks[i] = current_rank
                available[i] = false
            end
        end
        current_rank += 1
    end
    
    return ranks
end

"""
    register_mutation_region_motifs!(json_motifs, html_dict, motif_metadata_list; 
                                     start_idx=1, sort_globally=true, sort_by_pareto=true)

Save and register collected mutation region motifs. Sorts motifs by sign (positive first),
then by group (single_region first), then by Pareto rank within each group.

Sorting hierarchy:
1. Sign: Positive contributions first, then negative
2. Group: single_region, 2_regions, 3_regions, ... (within each sign)
3. Pareto rank: Based on abs(median) and count (within each group)
4. Ties: Sorted by abs(median) - descending for positives, ascending for negatives

Parameters:
- `json_motifs`: JSON motif dictionary
- `html_dict`: HTML dictionary
- `motif_metadata_list`: Vector (or vector of vectors) of motif metadata from collect_mutation_region_metadata
- `start_idx::Int = 1`: Starting index for mode numbering
- `sort_globally::Bool = true`: If true, sort motifs globally (by sign, group, then Pareto rank)
- `sort_by_pareto::Bool = true`: If true, use Pareto ranking within groups (median, count); else use median only

Returns:
- Next available index for mode numbering
"""
function register_mutation_region_motifs!(json_motifs, html_dict, motif_metadata_list;
        start_idx = 1, sort_globally = true, sort_by_pareto = true)
    
    # Flatten if needed (handles both single vector and vector of vectors)
    # Check if first element is a vector (indicates nested structure)
    all_metadata = if !isempty(motif_metadata_list) && motif_metadata_list[1] isa Vector
        vcat(motif_metadata_list...)
    else
        motif_metadata_list
    end
    
    # Sort globally if requested
    if sort_globally
        if sort_by_pareto
            # Group by sign, then by group_id, compute Pareto ranks within each group, then sort
            all_metadata = sort_by_group_and_pareto(all_metadata)
        else
            # Simple sort: positive first, then by group_id (single_region first), then by abs(median)
            # Positives: high to low abs(median), Negatives: low to high abs(median)
            all_metadata = sort(all_metadata, by=m -> (
                m.median > 0 ? 0 : 1,  # Positive first
                m.group_id == "single_region" ? 0 : parse(Int, split(m.group_id, '_')[1]),
                m.median > 0 ? -abs(m.median) : abs(m.median)  # Desc for pos, asc for neg
            ))
        end
    end
    
    current_idx = start_idx
    
    registered_names = String[]  # Track what we register
    for meta in all_metadata
        mkpath(meta.save_folder)
        
        # Generate file name (with filter prefix for singletons to avoid collisions)
        # and display name (clean span for UI display)
        if !isempty(meta.span) && meta.motif_size == 1 && haskey(meta.key, :m1)
            # Singletons: use filter-prefixed name for file
            # For filename, use simple start:end (first to last position) to avoid commas/spaces
            if meta.group_id != "single_region" && length(meta.positions) > 0
                # Multi-region singleton: use full range for filename
                first_pos = meta.positions[1]
                last_mat = meta.count_matrices[end]
                last_pos = meta.positions[end] + size(last_mat, 2) - 1
                file_name = "$(meta.key.m1)_$(first_pos):$(last_pos)"
            else
                # Single-region singleton: use span as-is
                file_name = "$(meta.key.m1)_$(meta.span)"
            end
            display_name = meta.span  # Always show the actual fragmented span on the card
        else
            # Multi-motifs: include filter indices AND original positions to avoid collisions
            if meta.motif_size > 1
                # Extract filter indices from key (m1, m2, m3, ...)
                filter_indices = [getfield(meta.key, Symbol("m$i")) for i in 1:meta.motif_size]
                filter_str = join(filter_indices, "_")
                
                # Use ORIGINAL positions from the grouping key (m1_position, m2_position, etc.)
                # These are the positions BEFORE any reference filtering/reduction
                original_positions = [getfield(meta.key, Symbol("m$(i)_position")) for i in 1:meta.motif_size]
                position_str = join(original_positions, "_")
                file_name = "$(filter_str)_$(position_str)"
                
                display_name = meta.span  # Show only the reduced span on card, not filters
            else
                # Fallback: use full key string
                file_name = !isempty(meta.span) ? meta.span : "$(get_k_mode_str(meta.key))"
                display_name = file_name
            end
        end
        
        paths = build_motif_paths(file_name, meta.save_folder, meta.motif_type)
        
        try
            EntroPlots.save_logo_with_rect_gaps(
                meta.count_matrices, meta.positions, meta.total_length,
                paths.png.abs; 
                reference_pfms=meta.references, 
                dpi=meta.dpi, 
                rna=meta.use_rna, 
                xrotation=35,
                protein=size(meta.count_matrices[1], 1) == 20,
                uniform_color=true, 
                filter_by_reference=meta.reduction_on_ref
            )
            
            save_influence_plot(meta.banzhafs, paths.influence.abs; xlim=meta.xlim)
            
            flat_windows = build_motif_windows(
                meta.gdf_row, meta.motif_size, meta.filter_len; 
                offset=meta.config.data.prefix_offset
            )
            save_positional_info(flat_windows, paths, meta.filter_len)
            
            texts = build_metadata_texts(
                nothing, paths, meta.median, meta.count; 
                use_rna=meta.use_rna, 
                relaxed_median=nothing, 
                show_meme_and_csv=false
            )
            
            mode_prefix = isempty(meta.group_id) ? "mode_" : "mode_$(meta.group_id)_"
            mode_str = mode_prefix * string(current_idx)
            add_motif_entry!(
                json_motifs, html_dict, mode_str, paths.png.rel, 
                "", texts, current_idx, display_name, meta.median, 
                meta.group_id, meta.button_text
            )
            
            push!(registered_names, display_name)  # Track registration
            current_idx += 1
        catch e
            @warn "Failed to save motif $display_name: $e"
        end
    end
    
    return current_idx
end

"""
    process_mutation_regions!(df_mutated, data, json_motifs, html_dict; kwargs...)

Convenience wrapper that collects and immediately registers mutation region motifs.
For global sorting across multiple dataframes, use collect_mutation_region_metadata 
and register_mutation_region_motifs! separately.

Returns:
- Next available index for mode numbering
"""
function process_mutation_regions!(df_mutated, data, json_motifs, html_dict;
        start_idx = 1, kwargs...)
    
    metadata = collect_mutation_region_metadata(df_mutated, data; kwargs...)
    return register_mutation_region_motifs!(json_motifs, html_dict, metadata; 
                                           start_idx=start_idx, sort_globally=true)
end

export prepare_and_collect_mutation_metadata, build_mutation_aggregates, compute_fragment_info
export collect_mutation_region_metadata, sort_by_group_and_pareto, compute_pareto_ranks_subset_wrapped
export register_mutation_region_motifs!, process_mutation_regions!
