"""
Reference matching utilities for off-region search in mutagenesis analysis.
"""

"""
    find_best_reference_position(mat, reference_seq, s_center, filter_len, search_radius)

Search ±search_radius neighborhood around s_center to find best-matching reference position.
Returns the start position with highest dot product score.
"""
function find_best_reference_position(mat::AbstractMatrix{T}, reference_seq, 
                                     s_center::Int, filter_len::Int, 
                                     search_radius::Int) where T
    ref_cols = size(reference_seq, 2)
    best_score = typemin(T)
    s_best = s_center
    
    for delta in -search_radius:search_radius
        s_candidate = s_center + delta
        e_candidate = s_candidate + filter_len - 1
        
        # Check bounds
        if s_candidate >= 1 && e_candidate <= ref_cols
            ref_window = @view reference_seq[:, s_candidate:e_candidate]
            score = compute_dotproduct_score(mat, ref_window)
            
            if score > best_score
                best_score = score
                s_best = s_candidate
            end
        end
    end
    
    return s_best
end

"""
    accumulate_and_find_reference!(mat, reference_seq, X, g, mp_sym, offset, filter_len, 
                                   idx_col, buf, off_region_search, search_radius)

Accumulate sequence counts for a single motif position and find best reference match.
Returns (buf, s_ref_best) tuple.
"""
function accumulate_and_find_reference!(mat::AbstractMatrix{T}, reference_seq, X, g, 
                                       mp_sym::Symbol, offset::Int, filter_len::Int, 
                                       idx_col, buf, off_region_search::Bool, 
                                       search_radius::Int) where T
    # Extract start position (same for all rows in this group)
    s_raw = Int(g[1, mp_sym])
    
    # Map to X coordinates and accumulate
    s_x = s_raw
    e_x = s_x + filter_len - 1
    buf = accumulate_motif_counts!(mat, X, s_x, e_x, idx_col, buf)
    
    # Find best reference position
    s_ref_best = s_raw + offset
    if off_region_search
        s_ref_best = find_best_reference_position(mat, reference_seq, 
                                                  s_raw + offset, filter_len, 
                                                  search_radius)
    end
    
    return buf, s_ref_best
end

"""
    obtain_count_and_reference_matrices(sorted_keys, gdf, data, reference_seq, mp_syms, filter_len; T=Float32, off_region_search=false)

Build count matrices and reference sequence views for each group in `gdf`.

# Arguments
- `sorted_keys`: Ordered collection of group keys
- `gdf`: GroupedDataFrame with motif position columns
- `data`: Dataset with `.X` field (one-hot encoded sequences)
- `reference_seq`: BitMatrix reference sequence (4 × seq_length)
- `mp_syms`: Vector of motif position column symbols
- `filter_len`: Filter length (window size)
- `T`: Element type for count matrices (default Float32)
- `off_region_search`: If true, search ±3 nt neighborhood for best reference match (default false)

# Returns
- `count_matrices_vec`: Dict mapping keys to Vector{Matrix{T}} (one 4×filter_len matrix per motif position)
- `reference_matrices_vec`: Dict mapping keys to Vector{SubArray} (reference windows)
- `adjusted_positions_vec`: Dict mapping keys to Vector{Int} (best-match reference positions)
"""
function obtain_count_and_reference_matrices(sorted_keys, gdf, data, reference_seq, 
                                             mp_syms, filter_len::Int; T=Float32, 
                                             off_region_search=false)
    # Validate inputs
    filter_len > 0 || throw(ArgumentError("filter_len must be positive"))
    isempty(mp_syms) && throw(ArgumentError("mp_syms cannot be empty"))
    
    # X = data.X
    X = data.onehot_sequences
    offset = data.prefix_offset
    key_type = eltype(sorted_keys)
    n_motifs = length(mp_syms)
    search_radius = 3  # ±3 nucleotides
    
    # Pre-allocate dicts
    count_matrices_vec = Dict{key_type, Vector{Matrix{T}}}()
    reference_matrices_vec = Dict{key_type, Vector{SubArray{Bool, 2}}}()
    adjusted_positions_vec = Dict{key_type, Vector{Int}}()
    
    sizehint!(count_matrices_vec, length(sorted_keys))
    sizehint!(reference_matrices_vec, length(sorted_keys))
    sizehint!(adjusted_positions_vec, length(sorted_keys))

    for k in sorted_keys
        g = gdf[k]
        isempty(g) && continue  # skip empty groups
        
        # Allocate single 3-D array for all motif count matrices (one allocation)
        h_dim = size(X, 1)
        mats_3d = zeros(T, h_dim, filter_len, n_motifs)
        idx_col = g[!, :data_pt_index]
        buf = nothing
        
        # Preallocate outputs
        reference_here = Vector{SubArray{Bool, 2}}(undef, n_motifs)
        adjusted_positions = Vector{Int}(undef, n_motifs)

        for (j, mp_sym) in enumerate(mp_syms)
            mat = view(mats_3d, :, :, j)
            
            # Accumulate counts and find best reference position
            buf, s_ref_best = accumulate_and_find_reference!(mat, reference_seq, X, g, 
                                                             mp_sym, offset, filter_len, 
                                                             idx_col, buf, off_region_search, 
                                                             search_radius)
            
            # Store results
            adjusted_positions[j] = s_ref_best
            e_ref_best = s_ref_best + filter_len - 1
            reference_here[j] = @view reference_seq[:, s_ref_best:e_ref_best]
        end

        # Convert 3-D array to vector of matrices for consistent interface
        count_matrices_vec[k] = [mats_3d[:, :, j] for j in 1:n_motifs]
        reference_matrices_vec[k] = reference_here
        adjusted_positions_vec[k] = adjusted_positions
        
        # Post-processing: merge overlapping matrices if off_region_search is enabled
        if off_region_search && n_motifs > 1
            merged_mats, merged_positions = 
                merge_overlapping_matrices(count_matrices_vec[k], adjusted_positions, filter_len)
            
            # Update with merged results
            count_matrices_vec[k] = merged_mats
            adjusted_positions_vec[k] = merged_positions
            
            # Update reference views to match merged positions
            reference_matrices_vec[k] = [
                @view reference_seq[:, pos:pos+size(mat,2)-1] 
                for (pos, mat) in zip(merged_positions, merged_mats)
            ]
        end
    end
    
    return count_matrices_vec, reference_matrices_vec, adjusted_positions_vec
end

export find_best_reference_position, accumulate_and_find_reference!, obtain_count_and_reference_matrices
