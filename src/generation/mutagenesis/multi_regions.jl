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
        
        # Retain only what downstream rendering actually reads from the group's
        # rows: `:data_pt_index` (population intersects) and the position columns
        # (build_motif_windows). `agg.gdf[k]` is a SubDataFrame *view* and
        # `agg.list_of_banzhafs[k]` is a column *view* — both pin the entire
        # parent DataFrame (all N rows, all columns) alive for the whole render
        # pass. Materialising a slim, owned copy here lets the parent grouped
        # DataFrame be GC'd after collection. Same values ⇒ output unchanged.
        gdf_row_slim = agg.gdf[k][:, [:data_pt_index; m_position_symbols(motif_size)]]
        banzhafs_owned = collect(agg.list_of_banzhafs[k])

        # Build motif data
        motif_data = MotifData(
            k, count_mats, start_positions,
            BitMatrix.(agg.reference_matrices_vec[k]),
            agg.median_map[k], agg.count_map[k],
            banzhafs_owned, gdf_row_slim
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

This is the single point where the mutagenesis path computes a per-motif NND
p-value — for single-region motifs and multi-region ones alike — so it is also
where the companion location z-score is computed. Both come back on the same
NamedTuple: `(k, obs_mNND, p_value, cluster_median, location_z)`.
"""
function save_indicator_and_nnd(meta, paths, file_name; pts, all_indices, nnd_k, bg_max_points=nothing)
    is_in_intersect = all_indices .∈ Ref(Set(intersect(meta.gdf_row.data_pt_index, all_indices)))
    kde_fig = plot_labels_vs_procprod(pts, is_in_intersect; motif_label="Contain motif", bg_max_points=bg_max_points)
    save(joinpath(dirname(paths.png.abs), "yy_kde_intersect_$(file_name).png"), kde_fig, px_per_unit=1)
    carriers = findall(is_in_intersect)
    nnd = nnd_permutation_test_1d(carriers, pts.labels; k=nnd_k)
    # Cluster median: where the motif's sequences sit on the expression axis
    # (median of their raw expression labels). Stored alongside the NND result
    # for display and later comparisons.
    cluster_median = any(is_in_intersect) ? median(@view pts.labels[is_in_intersect]) : NaN
    # Location z: do the carriers behave differently from the population, as
    # opposed to merely coherently (which is what NND asks)? Closed-form, so it
    # adds no measurable cost to this call. NaN when not computable.
    location_z = location_z_1d(carriers, pts.labels)
    return (; nnd.k, nnd.obs_mNND, nnd.p_value, cluster_median, location_z)
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
    motif_cluster_nnd(meta, pts, all_indices; nnd_k) -> Float64

Observed within-group mean k-NN distance (`obs_mNND`) of the sequences that
contain this motif, on the expression axis (`pts.labels`). This is the cheap
half of `nnd_permutation_test_1d` — only the observed statistic, no null
permutation loop — so it can be used to bin motifs for sorting before the
per-motif render loop computes the full (p-valued) NND result. Smaller means a
tighter cluster. Returns `NaN` when unavailable (`pts`/`all_indices` missing, or
fewer than two overlapping sequences).
"""
function motif_cluster_nnd(meta, pts, all_indices; nnd_k)
    (pts === nothing || all_indices === nothing) && return NaN
    try
        members = Set(intersect(meta.gdf_row.data_pt_index, all_indices))
        positions = findall(all_indices .∈ Ref(members))
        length(positions) < 2 && return NaN
        # Gather only the in-group labels directly — avoids allocating a full
        # length-N Float64 copy of `pts.labels` on every motif (this runs once
        # per motif in the pre-sort pass, so the old copy was O(N·M) churn).
        subpop_vals = sort!([Float64(pts.labels[i]) for i in positions])
        k_clamped = min(nnd_k, length(positions) - 1)
        return mean_knn_within_group_1d(subpop_vals, k_clamped)
    catch
        return NaN
    end
end

"""
    parse_span_positions(span) -> Set{Int}

Parse a motif's `span` string into the set of (1-based) residue positions it
marks. Handles the formats `compute_fragment_info` emits — comma-separated
tokens, each a single position (`"28"`) or an inclusive range (`"37:41"`),
optionally wrapped in parentheses, e.g. `"(26, 28, 32:33, 37:41)"` →
`{26,28,32,33,37,38,39,40,41}`. These are the motif's significant (post-reduction)
columns — the residues the WT track marks with arrows.
"""
function parse_span_positions(span::AbstractString)
    positions = Set{Int}()
    s = strip(span, ['(', ')', ' '])
    for tok in split(s, ',')
        t = strip(tok)
        isempty(t) && continue
        if occursin(':', t)
            parts = split(t, ':')
            length(parts) == 2 || continue
            a = tryparse(Int, strip(parts[1])); b = tryparse(Int, strip(parts[2]))
            (a === nothing || b === nothing) && continue
            for p in min(a, b):max(a, b); push!(positions, p); end
        else
            p = tryparse(Int, t)
            p !== nothing && push!(positions, p)
        end
    end
    return positions
end

"""
    motif_alphabet(nrows; use_rna=false) -> Union{Vector{Char},Nothing}

Alphabet matching a count matrix's row order: 20 rows is protein
(`SEQ2EXPdata.AMINO_ACID_LETTERS`, alphabetical — the order the one-hot encoder
uses), 4 rows is nucleotide. `nothing` for anything else.
"""
function motif_alphabet(nrows::Integer; use_rna::Bool=false)
    nrows == 20 && return SEQ2EXPdata.AMINO_ACID_LETTERS
    nrows == 4 && return [(use_rna ? _ind2dna_str_rna : _ind2dna_str_)[i] for i in 1:4]
    return nothing
end

"""
    observed_residue_string(cmat, ref=nothing; use_rna=false) -> String

The residue each column *mutates to*: one character per column, taken as the
mutated letter carrying the most information content — which, in a sequence logo,
is simply the tallest red letter in that column's stack.

Given a reference (backbone) matrix `ref`, the reference's own row is excluded
before taking the argmax. That exclusion is the whole point. The column's overall
argmax is usually the wild-type residue itself — EntroPlots draws that one blue
(`ref_match_color`) and every other residue dark red (`ref_mismatch_color`) — so
reading the plain argmax reports the backbone back as though it were a
substitution. On PIN1 that made half of every generated sentence a no-op
("site 33 to S" where S *is* the wild type).

Ranking by information content is equivalent to ranking by probability within a
column: EntroPlots draws letter `k` at height `ic * col[k]`, where `ic` is a
per-column constant, so the tallest letter is the highest-probability one and no
entropy term has to be recomputed here. Counts are used unnormalised for the same
reason — per-column normalisation is monotone and cannot change the argmax.

Columns with no mutated mass at all (every sequence carries the wild type), and
matrices whose row count matches no known alphabet, yield `_placeholder_char_`,
which [`mutation_description`](@ref) skips rather than inventing a substitution.
With `ref === nothing` this degrades to the plain column argmax.
"""
function observed_residue_string(cmat, ref=nothing; use_rna::Bool=false)
    alphabet = motif_alphabet(size(cmat, 1); use_rna=use_rna)
    alphabet === nothing && return String(fill(_placeholder_char_, size(cmat, 2)))
    nrow = size(cmat, 1)
    chars = Vector{Char}(undef, size(cmat, 2))
    for j in 1:size(cmat, 2)
        col = view(cmat, :, j)
        # Row holding the backbone residue for this column, when we have a
        # reference to compare against; skipped over when picking the argmax.
        skip = (ref !== nothing && j <= size(ref, 2)) ? argmax(view(ref, :, j)) : 0
        best, best_val = 0, zero(eltype(cmat))
        for i in 1:nrow
            i == skip && continue
            if best == 0 || col[i] > best_val
                best, best_val = i, col[i]
            end
        end
        chars[j] = (best > 0 && best_val > 0 && best <= length(alphabet)) ?
                   alphabet[best] : _placeholder_char_
    end
    return String(chars)
end

"""
    motif_wt_regions(meta) -> Vector{WildTypeRegion}

Build the wild-type amino-acid track(s) for a motif's summary card: one
`WildTypeRegion` per count-matrix window, holding the full WT window string
(sliced straight from `data.raw_data.consensus`), a per-column arrow flag, and
the motif's own observed consensus over the same window (what each mutated
column changes *to* — see [`observed_residue_string`](@ref)).

A column is marked only when it is **both** in the motif's `span` — the
significant (post-reduction) positions across *all* regions, so a 4-region motif
gets arrows in every region, not just the first — **and** it actually carries a
mutation, meaning some sequence holds a residue other than the backbone's (the
logo draws such letters dark red, the backbone's own letter blue).

Both halves are needed. Span membership alone marks positions the motif merely
depends on, including ones where every sequence holds the wild type; those are
pure blue and are not substitutions. The mutation test alone would mark
differences anywhere in the window, ignoring the reduction. When the span already
covers the entire window (e.g. `reduction_on_ref=false`, where it is one full
range) the stricter modal-residue test carries it on its own so the whole window
isn't arrowed. With no reference matrix to compare against, span membership is
all there is and is used as-is.

Returns an empty vector if the consensus is unavailable.
"""
function motif_wt_regions(meta)
    regions = WildTypeRegion[]
    consensus = try
        meta.config.data.raw_data.consensus
    catch
        nothing
    end
    consensus === nothing && return regions
    L = length(consensus)
    frag = parse_span_positions(meta.span)
    cmats, refs, pos = meta.count_matrices, meta.references, meta.positions
    n = min(length(cmats), length(pos))
    for r in 1:n
        cmat, startp = cmats[r], pos[r]
        w = size(cmat, 2)
        (startp < 1 || startp + w - 1 > L) && continue
        wt = String(consensus[startp:startp + w - 1])
        in_frag = Bool[(startp + j - 1) in frag for j in 1:w]

        ref = r <= length(refs) ? refs[r] : nothing

        # Two per-column tests against the backbone, both derived from the logo's
        # own colouring (blue = the reference residue, dark red = anything else):
        #
        #   has_mut  any red mass at all — some sequence carries a non-wild-type
        #            residue here. A column failing this is pure backbone and is
        #            not a mutation site however significant the reduction found it.
        #   differs  the column's *tallest* letter is red, i.e. the motif's modal
        #            residue is not the wild type. Stricter than `has_mut`.
        has_mut, differs = if ref === nothing
            nothing, nothing
        else
            hm = falses(w); df = falses(w)
            for j in 1:min(w, size(ref, 2))
                skip = argmax(view(ref, :, j))
                hm[j] = any(i -> i != skip && cmat[i, j] > 0, 1:size(cmat, 1))
                df[j] = argmax(view(cmat, :, j)) != skip
            end
            hm, df
        end

        # A mutation site is a significant position that actually carries a
        # mutation. Span membership alone is not enough: the reduction keeps
        # positions the motif depends on, including ones where every sequence
        # holds the wild type — reporting those as substitutions was wrong (on
        # PIN1 it made half of every sentence a no-op). Note the gate is
        # `has_mut`, not `differs`: a column where most sequences keep the wild
        # type but a mutated subset drives the effect is precisely an interesting
        # site, and `observed_residue_string` names its tallest red letter.
        # The gate is `has_mut` whenever we have a span, whether or not the span
        # covers the whole window. It used to fall back to `differs` when the
        # span covered every column, which is ALWAYS true for a single-region
        # motif — so every single-region motif produced an empty interpretation,
        # and a multi-region motif silently dropped any region whose window the
        # span fully covered. (On the RNA sim, motif `(7:14, 36, 38:41)` reported
        # only sites 36 and 38-41; `7:14` vanished.)
        #
        # `differs` is the wrong test for mutagenesis. It asks whether a column's
        # MODAL residue is non-wild-type, but each variant carries only a handful
        # of substitutions, so the modal residue is the wild type almost
        # everywhere and `differs` is false almost everywhere. `has_mut` asks
        # whether the column carries any mutated mass at all, which is the
        # question the surrounding comment already says we want to ask.
        mutated = if has_mut === nothing
            in_frag                                  # no reference: span is all we have
        elseif !isempty(frag)
            in_frag .& has_mut                       # span known → its columns that carry mutations
        else
            differs                                  # no span at all → modal-residue test
        end
        obs = observed_residue_string(cmat, ref;
            use_rna = hasproperty(meta, :use_rna) ? meta.use_rna : false)
        push!(regions, WildTypeRegion(wt, startp, mutated, obs))
    end
    return regions
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
    bin_value_asc(value, lo, hi, bin_count) -> Int

Like `bin_value`, but for keys sorted ascending (low value first): `NaN` maps to
`bin_count` (one past the last real bin) so missing values still sort last.
"""
function bin_value_asc(value::Real, lo::Real, hi::Real, bin_count::Int)
    b = bin_value(value, lo, hi, bin_count)
    return b < 0 ? bin_count : b
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
    sort_metadata_for_registration(all_metadata; sort_globally, sort_by_bins, sort_by_pareto, bin_count, pts, all_indices, nnd_k)

Return `all_metadata` in render order. Extracted verbatim from
`register_mutation_region_motifs!` — behaviour unchanged:
- `sort_globally=false` → returned as-is.
- `sort_by_bins=true` (default) → binned lexicographic order
  (group → cluster-median bin → cluster-NND bin → Shapley-median bin → count).
- `sort_by_pareto=true` → `sort_by_group_and_pareto`.
- otherwise → plain lexicographic sign → group → |median| → count order.
"""
function sort_metadata_for_registration(all_metadata;
        sort_globally = true, sort_by_bins = true, sort_by_pareto = false,
        bin_count = 10, pts = nothing, all_indices = nothing, nnd_k = 15)
    sort_globally || return all_metadata

    if sort_by_bins
        # Binned lexicographic sort. Group order is primary so the group
        # toggles render single_region → 2_regions → 3_regions → … . Within a
        # group, continuous values (cluster median, cluster NND, Shapley
        # median) are coarsened into `bin_count` equal-width bins over their
        # global ranges and compared bin-by-bin — sorting on raw continuous
        # values draws meaningless distinctions between near-equal motifs.
        # Count (un-binned) breaks remaining ties.
        augmented = [(
            meta = m,
            grp = group_order(m),
            clus = motif_cluster_median(m, pts, all_indices),
            nnd = motif_cluster_nnd(m, pts, all_indices; nnd_k = nnd_k),
            shap = float(m.median),
            cnt = m.count,
        ) for m in all_metadata]

        clus_finite = [a.clus for a in augmented if !isnan(a.clus)]
        nnd_finite = [a.nnd for a in augmented if !isnan(a.nnd)]
        shap_finite = [a.shap for a in augmented if !isnan(a.shap)]
        clus_lo, clus_hi = isempty(clus_finite) ? (0.0, 0.0) : extrema(clus_finite)
        nnd_lo, nnd_hi = isempty(nnd_finite) ? (0.0, 0.0) : extrema(nnd_finite)
        shap_lo, shap_hi = isempty(shap_finite) ? (0.0, 0.0) : extrema(shap_finite)

        sort!(augmented, by = a -> (
            a.grp,                                                  # single_region, then 2,3,4,… regions
            -bin_value(a.clus, clus_lo, clus_hi, bin_count),        # high cluster median bin first
            bin_value_asc(a.nnd, nnd_lo, nnd_hi, bin_count),        # low NND bin first (tighter cluster)
            -bin_value(a.shap, shap_lo, shap_hi, bin_count),        # high Shapley bin first
            -a.cnt,                                                 # higher count first
        ))
        return [a.meta for a in augmented]
    elseif sort_by_pareto
        # Group by sign, then by group_id, compute Pareto ranks within each group, then sort
        return sort_by_group_and_pareto(all_metadata)
    else
        # Simple lexicographic sort: single_region always first (regardless of
        # sign — a near-zero negative singleton must not be exiled below the
        # multi-region cards), then sign (positive first), then group, then
        # influence, then count.
        # Positives: high to low abs(median); Negatives: low to high abs(median).
        # Count breaks the (rare) influence ties — higher count first.
        return sort(all_metadata, by = m -> (
            m.group_id == "single_region" ? 0 : 1,           # Singleton group first
            m.median > 0 ? 0 : 1,  # Positive first
            m.group_id == "single_region" ? 0 : parse(Int, split(m.group_id, '_')[1]),
            m.median > 0 ? -abs(m.median) : abs(m.median),  # Desc for pos, asc for neg
            -m.count                                         # Higher count first
        ))
    end
end

"""
    RegisterDiagnostics()

Mutable, shared per-run counter bag for `register_mutation_region_motifs!`.
`render_one_motif!` bumps these at the exact points the pre-refactor loop body
did, so a counter incremented before a later throw in the same iteration
survives — keeping the end-of-run summary bit-identical to before the extraction.
"""
mutable struct RegisterDiagnostics
    n_failed::Int
    first_error::Any
    n_indicator_failed::Int
    first_indicator_error::Any
    n_interaction_candidates::Int
    n_interaction_hits::Int
end
RegisterDiagnostics() = RegisterDiagnostics(0, nothing, 0, nothing, 0, 0)

"""
    render_one_motif!(json_motifs, html_dict, meta, paths, file_name, display_name, current_idx, diag; kwargs...)

Render one motif's outputs (logo, influence plot, optional indicator/NND plot,
positional CSV) and register its card into `json_motifs`/`html_dict`, pushing a
`TopMoverEntry` into `top_movers_out` when supplied. Extracted verbatim from the
`register_mutation_region_motifs!` loop body — behaviour unchanged.

Non-fatal diagnostics (indicator failures, interaction candidacy/hits) are bumped
in-place on `diag::RegisterDiagnostics` at the same points as the original loop.
A throw propagates to the caller's per-motif `try` (which records `diag.n_failed`).
"""
function render_one_motif!(json_motifs, html_dict, meta, paths, file_name, display_name, current_idx, diag;
        pts = nothing, all_indices = nothing, nnd_k = 15,
        interaction_summaries = nothing, top_movers_out = nothing, bg_max_points = nothing,
        report_location_z::Bool = false,
        show_region_interaction::Bool = false,   # display the region-interaction line? computed either way
        feature_label = nothing)
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

    # The top-movers page shows BOTH logo views side by side: the reduced view
    # (filter_by_reference=true — only the mutated fragments that differ from the
    # backbone, the basis of region interactions) and the full-region view
    # (filter_by_reference=false). Above we rendered the grouped-page default
    # (governed by `reduction_on_ref`); here we render the other view so both are
    # available. `img_reduced`/`img_region` are the relative paths used by the row.
    alt_png_abs = joinpath(meta.save_folder, file_name * "_altview.png")
    alt_png_rel = joinpath(meta.motif_type, file_name * "_altview.png")
    EntroPlots.save_logo_with_rect_gaps(
        meta.count_matrices, meta.positions, meta.total_length,
        alt_png_abs;
        reference_pfms=meta.references,
        dpi=meta.dpi,
        rna=meta.use_rna,
        xrotation=35,
        protein=size(meta.count_matrices[1], 1) == 20,
        uniform_color=true,
        filter_by_reference=!meta.reduction_on_ref
    )
    img_reduced = meta.reduction_on_ref ? paths.png.rel : alt_png_rel
    img_region  = meta.reduction_on_ref ? alt_png_rel : paths.png.rel

    save_influence_plot(meta.banzhafs, paths.influence.abs; xlim=meta.xlim)

    # Per-motif indicator plot + NND test (only when `pts` is supplied).
    # Non-fatal: a failure here (e.g. pts not aligned to all_indices) must
    # not abort the whole card — just skip the indicator plot.
    nnd_result = nothing
    if pts !== nothing && all_indices !== nothing
        try
            nnd_result = save_indicator_and_nnd(meta, paths, file_name;
                pts=pts, all_indices=all_indices, nnd_k=nnd_k, bg_max_points=bg_max_points)
        catch e_ind
            diag.n_indicator_failed += 1
            diag.first_indicator_error === nothing && (diag.first_indicator_error = e_ind)
        end
    end

    flat_windows = build_motif_windows(
        meta.gdf_row, meta.motif_size, meta.filter_len;
        offset=meta.config.data.prefix_offset
    )
    save_positional_info(flat_windows, paths, meta.filter_len)

    # Per-motif interaction summary (multi-region motifs only). Count candidacy
    # here — the exact point the pre-refactor loop did — so a later throw in
    # texts/add_motif_entry! leaves the already-bumped count intact.
    is_interaction_candidate, interaction_str = lookup_interaction(meta, interaction_summaries)
    if is_interaction_candidate
        diag.n_interaction_candidates += 1
        interaction_str !== nothing && (diag.n_interaction_hits += 1)
    end

    # Wild-type context, computed once: it backs both the popup's plain-language
    # interpretation (text slot 7) and the summary row's residue track.
    wt_regions = motif_wt_regions(meta)
    description = mutation_description(wt_regions, float(meta.median), feature_label)

    texts = build_metadata_texts(
        nothing, paths, meta.median, meta.count;
        use_rna=meta.use_rna,
        relaxed_median=nothing,
        show_meme_and_csv=false,
        # Computed above regardless — `show_region_interaction` only decides
        # whether the popup prints it. The CSV keeps the value either way.
        interaction_summary_mode_str=(show_region_interaction ? interaction_str : nothing),
        nnd_result=nnd_result,
        report_location_z=report_location_z,
        description=description
    )

    mode_prefix = isempty(meta.group_id) ? "mode_" : "mode_$(meta.group_id)_"
    mode_str = mode_prefix * string(current_idx)
    add_motif_entry!(
        json_motifs, html_dict, mode_str, paths.png.rel,
        "", texts, current_idx, display_name, meta.median,
        meta.group_id, meta.button_text
    )

    # Collect a self-contained summary row for the top-movers landing page.
    # Uses data already in scope; cluster median/NND/p come from the (cheap
    # half of the) NND result, NaN when `pts` wasn't supplied.
    if top_movers_out !== nothing
        cm = nnd_result === nothing ? NaN : float(nnd_result.cluster_median)
        cn = nnd_result === nothing ? NaN : float(nnd_result.obs_mNND)
        cp = nnd_result === nothing ? NaN : float(nnd_result.p_value)
        lz = nnd_result === nothing ? NaN : float(nnd_result.location_z)
        push!(top_movers_out, TopMoverEntry(
            float(meta.median), cm, cn, meta.count,
            display_name, meta.button_text, meta.span, cp,
            meta.motif_size == 1,
            interaction_str === nothing ? "" : interaction_str,
            wt_regions,
            paths.png.rel, img_reduced, img_region, paths.influence.rel,
            joinpath(meta.motif_type, "yy_kde_intersect_$(file_name).png"),
            collect(String, json_motifs[mode_str][texts_str][1]),
            # No slider here: the mutagenesis page uses a different file
            # layout (per-file yy_kde, no base-folder relaxed plot), so its
            # popup keeps the static singleton modal. (Conv case opts in.)
            String[], String[], Vector{String}[],
            lz,                                     # location_z (NaN when unavailable)
        ))
    end

    return nothing
end

"""
    register_mutation_region_motifs!(json_motifs, html_dict, motif_metadata_list;
                                     start_idx=1, sort_globally=true, sort_by_pareto=true)

Save and register collected mutation region motifs.

Default sorting (`sort_by_bins=true`) is a binned lexicographic order:
1. group: single_region first, then 2_regions, 3_regions, … (keeps the group toggles ordered)
2. cluster-median bin (high → low) — median expression of the motif's sequences, `bin_count` equal-width bins
3. cluster-NND bin (low → high) — observed mean k-NN distance of those sequences; tighter clusters first, same binning
4. Shapley-median bin (high → low) — same binning over the global range
5. count (high → low), left un-binned

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
- `bin_count::Int = 10`: Number of equal-width bins for the cluster median, cluster NND, and Shapley median
- `sort_by_pareto::Bool = false`: Fallback Pareto ranking when `sort_by_bins=false`

Returns:
- Next available index for mode numbering
"""
function register_mutation_region_motifs!(json_motifs, html_dict, motif_metadata_list;
        start_idx = 1, sort_globally = true, sort_by_pareto = false,
        sort_by_bins = true, bin_count = 10,
        pts = nothing, all_indices = nothing, interaction_summaries = nothing,
        nnd_k = 15, top_movers_out = nothing,
        gc_every::Int = 25, bg_max_points = nothing,
        report_location_z::Bool = false,  # also show the location z-score on each card
        show_region_interaction::Bool = false,  # display the region-interaction line? computed either way
        feature_label = nothing)   # assay name (+units) for the interpretation sentence

    # Flatten if needed (handles both single vector and vector of vectors)
    # Check if first element is a vector (indicates nested structure)
    all_metadata = if !isempty(motif_metadata_list) && motif_metadata_list[1] isa Vector
        vcat(motif_metadata_list...)
    else
        motif_metadata_list
    end

    # Sort into render order (binned-lexicographic by default; see helper).
    all_metadata = sort_metadata_for_registration(all_metadata;
        sort_globally = sort_globally, sort_by_bins = sort_by_bins,
        sort_by_pareto = sort_by_pareto, bin_count = bin_count,
        pts = pts, all_indices = all_indices, nnd_k = nnd_k)

    current_idx = start_idx

    registered_names = String[]  # Track what we register

    # Accumulate failures/diagnostics and report once at the end instead of
    # spamming one log line per motif (keeps stdout quiet on large runs). Shared,
    # mutable so `render_one_motif!` can bump a counter at the exact point the
    # old inline loop did and have it survive a later throw in the same iteration.
    diag = RegisterDiagnostics()

    for (iter_no, meta) in enumerate(all_metadata)
        mkpath(meta.save_folder)

        file_name, display_name = motif_filename_and_display(meta)

        paths = build_motif_paths(file_name, meta.save_folder, meta.motif_type)

        try
            render_one_motif!(
                json_motifs, html_dict, meta, paths, file_name, display_name, current_idx, diag;
                pts=pts, all_indices=all_indices, nnd_k=nnd_k,
                interaction_summaries=interaction_summaries, top_movers_out=top_movers_out,
                bg_max_points=bg_max_points, report_location_z=report_location_z,
                show_region_interaction=show_region_interaction,
                feature_label=feature_label
            )
            push!(registered_names, display_name)  # Track registration
            current_idx += 1
        catch e
            diag.n_failed += 1
            diag.first_error === nothing && (diag.first_error = e)
        end

        # Reclaim plotting memory periodically. GR/Plots keeps global figure
        # state and CairoMakie figures back onto C-side Cairo surfaces that
        # Julia's GC does not count toward heap pressure, so a long render loop
        # leaks RSS (~1 GB / few-hundred motifs, measured) unless we tear those
        # down explicitly. `GC.gc(false)` is a cheap *incremental* collection
        # (~1 ms) that still runs the surface finalizers — a full `GC.gc()` here
        # costs ~0.5 s/call and is far too slow to do in a loop. Done every
        # `gc_every` iterations to amortise the cost (default ≈10% slower, ~2×
        # less peak RSS); `gc_every=0` disables it and restores the previous
        # (leakier but marginally faster) behaviour exactly.
        if gc_every > 0 && iter_no % gc_every == 0
            try; EntroPlots.Plots.closeall(); catch; end
            GC.gc(false)
        end
    end

    # One summary line each, instead of a per-motif flood.
    diag.n_failed > 0 && @warn "register_mutation_region_motifs!: $(diag.n_failed) motif(s) failed to render and were skipped" diag.first_error
    diag.n_indicator_failed > 0 && @warn "register_mutation_region_motifs!: indicator plot skipped for $(diag.n_indicator_failed) motif(s) — check that `pts` is aligned to `all_indices`" diag.first_indicator_error
    diag.n_interaction_candidates > 0 && @info "register_mutation_region_motifs!: interaction summaries matched $(diag.n_interaction_hits) / $(diag.n_interaction_candidates) multi-region motifs"

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
