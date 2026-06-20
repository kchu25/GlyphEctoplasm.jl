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
        is_singleton = m.group_id == "single_region"
        group_order = is_singleton ? 0 : parse(Int, split(m.group_id, '_')[1])
        sign_order = m.median > 0 ? 0 : 1
        augmented[i] = (
            metadata=m,
            singleton_order=is_singleton ? 0 : 1,  # Singleton group always first
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
                singleton_order=old.singleton_order,
                group_order=old.group_order,
                sign_order=old.sign_order,
                pareto_rank=rank
            )
        end
        
        i = j
    end
    
    # Final sort by (singleton, sign, group, pareto_rank, abs(median))
    # Singleton group always leads, regardless of sign.
    # For positives: higher abs(median) first (-abs for descending)
    # For negatives: lower abs(median) first (+abs for ascending)
    sort!(augmented, by = a -> (
        a.singleton_order,
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
    motif_filter_indices(meta)

The filter indices `(m1, m2, ...)` of a motif, read from its grouping key.
Returns a tuple of length `meta.motif_size`.
"""
motif_filter_indices(meta) = ntuple(i -> getfield(meta.key, Symbol("m$i")), meta.motif_size)

"""
    motif_filename_and_display(meta) -> (file_name, display_name)

Compute the collision-free output filename and the card display name for a motif:
- single-region singleton: `<m1>_<span>`
- multi-region singleton:  `<m1>_<first>:<last>` (full position range)
- multi-motif:             `<filters>_<positions>` (original, pre-reduction positions)
- fallback:                the span (or the key string)

The card always shows the (reduced) span; only the filename disambiguates.
"""
function motif_filename_and_display(meta)
    if !isempty(meta.span) && meta.motif_size == 1 && haskey(meta.key, :m1)
        # Singletons: filter-prefixed filename, span shown on the card
        if meta.group_id != "single_region" && length(meta.positions) > 0
            # Multi-region singleton: use full range for filename
            first_pos = meta.positions[1]
            last_mat  = meta.count_matrices[end]
            last_pos  = meta.positions[end] + size(last_mat, 2) - 1
            file_name = "$(meta.key.m1)_$(first_pos):$(last_pos)"
        else
            # Single-region singleton: use span as-is
            file_name = "$(meta.key.m1)_$(meta.span)"
        end
        return file_name, meta.span
    elseif meta.motif_size > 1
        # Multi-motifs: filter indices AND original (pre-reduction) positions to avoid collisions
        filter_str   = join(motif_filter_indices(meta), "_")
        position_str = join((getfield(meta.key, Symbol("m$(i)_position")) for i in 1:meta.motif_size), "_")
        return "$(filter_str)_$(position_str)", meta.span
    else
        # Fallback: use full key string
        file_name = !isempty(meta.span) ? meta.span : "$(get_k_mode_str(meta.key))"
        return file_name, file_name
    end
end

"""
    save_indicator_and_nnd(meta, paths, file_name; pts, all_indices, nnd_k) -> nnd_result

Save the per-motif indicator plot (yy-KDE) and run the cluster-tightness (NND)
permutation test, returning its result. The plot is written under the filename
convention the singleton modal derives from the card image
(`<dir>/yy_kde_intersect_<file_name>.png`). May throw; callers wrap this so a
failure (e.g. `pts` not aligned to `all_indices`) only skips this one plot.
"""
function save_indicator_and_nnd(meta, paths, file_name; pts, all_indices, nnd_k)
    is_in_intersect = all_indices .∈ Ref(Set(intersect(meta.gdf_row.data_pt_index, all_indices)))
    kde_fig = plot_labels_vs_procprod(pts, is_in_intersect; motif_label="Contain motif")
    save(joinpath(dirname(paths.png.abs), "yy_kde_intersect_$(file_name).png"), kde_fig, px_per_unit=1)
    nnd = nnd_permutation_test_1d(findall(is_in_intersect), pts.labels; k=nnd_k)
    # Cluster median: where the motif's sequences sit on the expression axis
    # (median of their raw expression labels). Stored alongside the NND result
    # for display and later comparisons.
    cluster_median = any(is_in_intersect) ? median(@view pts.labels[is_in_intersect]) : NaN
    return (; nnd.k, nnd.obs_mNND, nnd.p_value, cluster_median)
end

"""
    motif_cluster_median(meta, pts, all_indices) -> Float64

Median of the expression labels (`pts.labels`) over the sequences that contain
this motif. Returns `NaN` when unavailable (`pts`/`all_indices` missing, or no
overlapping sequences). Cheap (no permutation test) — used to bin motifs for
sorting before the per-motif render loop computes the full NND result.
"""
function motif_cluster_median(meta, pts, all_indices)
    (pts === nothing || all_indices === nothing) && return NaN
    try
        members = Set(intersect(meta.gdf_row.data_pt_index, all_indices))
        mask = all_indices .∈ Ref(members)
        return any(mask) ? median(@view pts.labels[mask]) : NaN
    catch
        return NaN
    end
end

"""
    bin_value(value, lo, hi, bin_count) -> Int

Map `value` to one of `bin_count` equal-width bins spanning `[lo, hi]`, returning
a 0-based bin index in `0:(bin_count-1)`. Degenerate ranges (`hi <= lo`) collapse
to bin `0`; `NaN` values return `-1` (a sentinel that sorts last under
highest-bin-first ordering).
"""
function bin_value(value::Real, lo::Real, hi::Real, bin_count::Int)
    isnan(value) && return -1
    hi <= lo && return 0
    return clamp(floor(Int, (value - lo) / (hi - lo) * bin_count), 0, bin_count - 1)
end

"""
    group_order(meta) -> Int

Region-group rank used as the primary sort key so the group toggles render in a
stable order: `single_region` first (`0`), then `2_regions`, `3_regions`, …
(the integer prefix of `group_id`). Unparseable group ids sort last.
"""
function group_order(meta)
    meta.group_id == "single_region" && return 0
    n = tryparse(Int, first(split(meta.group_id, '_')))
    return n === nothing ? typemax(Int) : n
end

"""
    lookup_interaction(meta, interaction_summaries) -> (is_candidate, interaction_str)

Look up a multi-region motif's interaction summary. `interaction_summaries` is a
vector indexed by motif size (index 1 => size 2, ...), each a Dict keyed by the
`(m1, m2, ...)` filter-index NamedTuple. Returns whether this motif is even a
candidate (multi-region with a dict for its size) and the matched string (or
`nothing`). No remapping needed — the mutation pipeline keeps original indices.
"""
function lookup_interaction(meta, interaction_summaries)
    (interaction_summaries === nothing || meta.motif_size < 2 ||
        meta.motif_size - 1 > length(interaction_summaries)) && return (false, nothing)
    summary_dict = interaction_summaries[meta.motif_size - 1]
    syms = ntuple(i -> Symbol("m$i"), meta.motif_size)
    vals = Int.(motif_filter_indices(meta))
    return (true, get(summary_dict, NamedTuple{syms}(vals), nothing))
end

"""
    register_mutation_region_motifs!(json_motifs, html_dict, motif_metadata_list;
                                     start_idx=1, sort_globally=true, sort_by_pareto=true)

Save and register collected mutation region motifs.

Default sorting (`sort_by_bins=true`) is a binned lexicographic order:
1. group: single_region first, then 2_regions, 3_regions, … (keeps the group toggles ordered)
2. Shapley-median bin (high → low) — `bin_count` equal-width bins over the global range
3. cluster-median bin (high → low) — median expression of the motif's sequences, same binning
4. count (high → low), left un-binned

Binning coarsens near-equal continuous values so the next key can break ties
meaningfully. Set `sort_by_bins=false` to fall back to the older `sort_by_pareto`
(Pareto ranks on |median|/count) or, with both false, the plain lexicographic
sign→group→|median|→count order.

Parameters:
- `json_motifs`: JSON motif dictionary
- `html_dict`: HTML dictionary
- `motif_metadata_list`: Vector (or vector of vectors) of motif metadata from collect_mutation_region_metadata
- `start_idx::Int = 1`: Starting index for mode numbering
- `sort_globally::Bool = true`: If true, sort motifs globally
- `sort_by_bins::Bool = true`: If true, use the binned sort described above (takes precedence)
- `bin_count::Int = 10`: Number of equal-width bins for the Shapley and cluster medians
- `sort_by_pareto::Bool = false`: Fallback Pareto ranking when `sort_by_bins=false`

Returns:
- Next available index for mode numbering
"""
function register_mutation_region_motifs!(json_motifs, html_dict, motif_metadata_list;
        start_idx = 1, sort_globally = true, sort_by_pareto = false,
        sort_by_bins = true, bin_count = 10,
        pts = nothing, all_indices = nothing, interaction_summaries = nothing,
        nnd_k = 15)

    # Flatten if needed (handles both single vector and vector of vectors)
    # Check if first element is a vector (indicates nested structure)
    all_metadata = if !isempty(motif_metadata_list) && motif_metadata_list[1] isa Vector
        vcat(motif_metadata_list...)
    else
        motif_metadata_list
    end

    # Sort globally if requested
    if sort_globally
        if sort_by_bins
            # Binned lexicographic sort. Group order is primary so the group
            # toggles render single_region → 2_regions → 3_regions → … . Within a
            # group, continuous values (Shapley median, cluster median) are
            # coarsened into `bin_count` equal-width bins over their global ranges
            # and compared bin-by-bin — sorting on raw continuous values draws
            # meaningless distinctions between near-equal motifs. Count (un-binned)
            # breaks remaining ties.
            augmented = [(
                meta = m,
                grp = group_order(m),
                shap = float(m.median),
                clus = motif_cluster_median(m, pts, all_indices),
                cnt = m.count,
            ) for m in all_metadata]

            shap_finite = [a.shap for a in augmented if !isnan(a.shap)]
            clus_finite = [a.clus for a in augmented if !isnan(a.clus)]
            shap_lo, shap_hi = isempty(shap_finite) ? (0.0, 0.0) : extrema(shap_finite)
            clus_lo, clus_hi = isempty(clus_finite) ? (0.0, 0.0) : extrema(clus_finite)

            sort!(augmented, by = a -> (
                a.grp,                                                  # single_region, then 2,3,4,… regions
                -bin_value(a.shap, shap_lo, shap_hi, bin_count),        # high Shapley bin first
                -bin_value(a.clus, clus_lo, clus_hi, bin_count),        # high cluster bin first
                -a.cnt,                                                 # higher count first
            ))
            all_metadata = [a.meta for a in augmented]
        elseif sort_by_pareto
            # Group by sign, then by group_id, compute Pareto ranks within each group, then sort
            all_metadata = sort_by_group_and_pareto(all_metadata)
        else
            # Simple lexicographic sort: single_region always first (regardless of
            # sign — a near-zero negative singleton must not be exiled below the
            # multi-region cards), then sign (positive first), then group, then
            # influence, then count.
            # Positives: high to low abs(median); Negatives: low to high abs(median).
            # Count breaks the (rare) influence ties — higher count first.
            all_metadata = sort(all_metadata, by=m -> (
                m.group_id == "single_region" ? 0 : 1,           # Singleton group first
                m.median > 0 ? 0 : 1,  # Positive first
                m.group_id == "single_region" ? 0 : parse(Int, split(m.group_id, '_')[1]),
                m.median > 0 ? -abs(m.median) : abs(m.median),  # Desc for pos, asc for neg
                -m.count                                         # Higher count first
            ))
        end
    end
    
    current_idx = start_idx

    registered_names = String[]  # Track what we register

    # Accumulate failures/diagnostics and report once at the end instead of
    # spamming one log line per motif (keeps stdout quiet on large runs).
    n_failed = 0;            first_error = nothing
    n_indicator_failed = 0;  first_indicator_error = nothing
    n_interaction_candidates = 0; n_interaction_hits = 0

    for meta in all_metadata
        mkpath(meta.save_folder)
        
        file_name, display_name = motif_filename_and_display(meta)

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

            # Per-motif indicator plot + NND test (only when `pts` is supplied).
            # Non-fatal: a failure here (e.g. pts not aligned to all_indices) must
            # not abort the whole card — just skip the indicator plot.
            nnd_result = nothing
            if pts !== nothing && all_indices !== nothing
                try
                    nnd_result = save_indicator_and_nnd(meta, paths, file_name;
                        pts=pts, all_indices=all_indices, nnd_k=nnd_k)
                catch e_ind
                    n_indicator_failed += 1
                    first_indicator_error === nothing && (first_indicator_error = e_ind)
                end
            end

            flat_windows = build_motif_windows(
                meta.gdf_row, meta.motif_size, meta.filter_len; 
                offset=meta.config.data.prefix_offset
            )
            save_positional_info(flat_windows, paths, meta.filter_len)
            
            # Per-motif interaction summary (multi-region motifs only).
            is_interaction_candidate, interaction_str = lookup_interaction(meta, interaction_summaries)
            if is_interaction_candidate
                n_interaction_candidates += 1
                interaction_str !== nothing && (n_interaction_hits += 1)
            end

            texts = build_metadata_texts(
                nothing, paths, meta.median, meta.count;
                use_rna=meta.use_rna,
                relaxed_median=nothing,
                show_meme_and_csv=false,
                interaction_summary_mode_str=interaction_str,
                nnd_result=nnd_result
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
            n_failed += 1
            first_error === nothing && (first_error = e)
        end
    end

    # One summary line each, instead of a per-motif flood.
    n_failed > 0 && @warn "register_mutation_region_motifs!: $n_failed motif(s) failed to render and were skipped" first_error
    n_indicator_failed > 0 && @warn "register_mutation_region_motifs!: indicator plot skipped for $n_indicator_failed motif(s) — check that `pts` is aligned to `all_indices`" first_indicator_error
    n_interaction_candidates > 0 && @info "register_mutation_region_motifs!: interaction summaries matched $n_interaction_hits / $n_interaction_candidates multi-region motifs"

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
