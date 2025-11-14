"""
Matrix operations for mutagenesis analysis: merging overlapping regions and accumulating counts.
"""

"""
    merge_overlapping_matrices(mats::Vector{Matrix{T}}, positions::Vector{Int}, 
                              filter_len::Int) where T

Merge count matrices that have overlapping reference regions by concatenating 
non-overlapping parts. Returns merged matrices and updated positions.

Example: positions=[41, 49], filter_len=9
  - Matrix 1 covers 41:49
  - Matrix 2 covers 49:57 (overlap at position 49)
  - Merged: mat1[:, 1:8] concatenated with mat2[:, 2:9] → covers 41:57
"""
function merge_overlapping_matrices(mats::Vector{Matrix{T}}, positions::Vector{Int}, 
                                   filter_len::Int) where T
    n = length(mats)
    n == 0 && return mats, positions
    n == 1 && return mats, positions
    
    # Sort by position
    perm = sortperm(positions)
    sorted_positions = positions[perm]
    sorted_mats = mats[perm]
    
    merged_mats = Matrix{T}[]
    merged_positions = Int[]
    
    i = 1
    while i <= n
        # Start a new merged matrix
        s_start = sorted_positions[i]
        s_end = s_start + filter_len - 1
        current_mat = sorted_mats[i]
        
        # Check for consecutive overlapping matrices
        j = i + 1
        while j <= n
            s_next = sorted_positions[j]
            
            # Check if next matrix overlaps OR is adjacent (touching) with current range
            # Adjacent means s_next == s_end + 1 (e.g., 35:43 and 44:52)
            if s_next <= s_end + 1
                # Calculate overlap
                overlap_start = s_next
                overlap_cols = s_end - overlap_start + 1  # how many columns overlap
                
                # Take non-overlapping part from next matrix
                mat_next = sorted_mats[j]
                if overlap_cols < filter_len
                    non_overlap_part = mat_next[:, overlap_cols+1:end]
                    current_mat = hcat(current_mat, non_overlap_part)
                    s_end = s_next + filter_len - 1  # extend range
                else
                    # Fully overlapped, skip this matrix
                end
                j += 1
            else
                break
            end
        end
        
        # Store merged result
        push!(merged_mats, current_mat)
        push!(merged_positions, s_start)
        
        i = j
    end
    
    return merged_mats, merged_positions
end

"""
    accumulate_motif_counts!(mat, X, s, e, idx_col, buf)

Accumulate sequence windows into count matrix `mat` for all data points in `idx_col`.
"""
function accumulate_motif_counts!(mat::AbstractMatrix{T}, X, s::Int, e::Int, 
                                  idx_col, buf) where T
    @inbounds for i in eachindex(idx_col)
        buf = accumulate_window!(mat, X, s, e, idx_col[i], buf)
    end
    return buf
end

"""
    compute_dotproduct_score(mat, ref_window)

Compute dot product similarity score between count matrix and reference window.
"""
function compute_dotproduct_score(mat::AbstractMatrix{T}, ref_window) where T
    score = zero(T)
    @inbounds for col in axes(mat, 2)
        for row in axes(mat, 1)
            score += mat[row, col] * T(ref_window[row, col])
        end
    end
    return score
end

export merge_overlapping_matrices, accumulate_motif_counts!, compute_dotproduct_score
