"""
    trim_dash(arr::Vector{Char}) -> Vector{Char}

Remove leading and trailing placeholder characters from array.

# Arguments
- `arr`: Character array potentially containing placeholder dashes

# Returns
Trimmed array, or ['−','−','−'] if all elements are placeholders

# Example
```julia
trim_dash(['−', 'A', 'C', 'G', '−'])  # ['A', 'C', 'G']
trim_dash(['−', '−', '−'])            # ['-', '-', '-']
```
"""
function trim_dash(arr::Vector{Char})
    start_idx = findfirst(x -> x != _placeholder_char_, arr)
    end_idx = findlast(x -> x != _placeholder_char_, arr)
    
    if start_idx === nothing || end_idx === nothing
        return ['-', '-', '-']  # All elements are placeholders
    else
        return arr[start_idx:end_idx]
    end
end

"""
    get_relaxed_consensus_str(pfm; dash = _placeholder_char_, prob_thresh = 0.5, rna = false) -> String

Generate consensus sequence string from position frequency matrix.

# Arguments
- `pfm`: Position frequency matrix (nucleotides × positions)
- `dash`: Placeholder character for low-confidence positions
- `prob_thresh`: Probability threshold for keeping base calls (default: 0.5)
- `rna`: If true, use RNA alphabet (U) instead of DNA (T)

# Returns
Consensus sequence string with low-confidence positions replaced by dashes

# Algorithm
1. Find the nucleotide with highest probability at each position
2. Replace positions with probability < threshold with placeholder
3. Trim leading and trailing placeholders
4. Return as string

# Example
```julia
pfm = [0.8 0.1 0.2;
       0.1 0.7 0.3;
       0.05 0.1 0.4;
       0.05 0.1 0.1]
get_relaxed_consensus_str(pfm)  # "AC−" (if third position < threshold)
```
"""
function get_relaxed_consensus_str(pfm; dash = _placeholder_char_, prob_thresh = 0.5, rna = false)
    argmax_inds = reshape(argmax(pfm, dims = 1), (size(pfm, 2),))
    
    char_array = if rna
        [_ind2dna_str_rna[i[1]] for i in argmax_inds]
    else
        [_ind2dna_str_[i[1]] for i in argmax_inds]
    end
    
    # Replace low-confidence positions with placeholder
    char_array[findall((@view pfm[argmax_inds]) .< prob_thresh)] .= dash
    
    return join(trim_dash(char_array))
end

"""
    countmat2consensus(countmat::AbstractMatrix; pseudo_count = float_type(0.1), rna = false) -> String

Convert count matrix to consensus sequence.

# Arguments
- `countmat`: Count matrix (nucleotides × positions)
- `pseudo_count`: Pseudocount to add for smoothing (default: 0.1)
- `rna`: If true, use RNA alphabet

# Returns
Consensus sequence string

# Algorithm
1. Add pseudocount to all positions for smoothing
2. Normalize to create position frequency matrix
3. Generate consensus string using threshold filtering

# Example
```julia
counts = [10 2 5;
          2 15 3;
          1 1 8;
          1 1 2]
countmat2consensus(counts)  # "ACG"
```
"""
function pfm2consensus(pfm; rna = false)
    consensus_str = get_relaxed_consensus_str(pfm; rna = rna)
    return consensus_str
end