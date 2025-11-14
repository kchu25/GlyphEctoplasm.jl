"""
Symbol generation utilities for motif column naming.

These functions generate column name symbols used throughout the analysis pipeline
for motif identities, positions, and distances.
"""

"""
    m_symbols(n::Integer) -> Vector{Symbol}

Generate motif identity column symbols: [:m1, :m2, ..., :mn]

Used to identify which filters/motifs are involved in a multi-motif pattern.

# Arguments
- `n`: Number of motifs in the pattern (must be positive)

# Returns
Vector of symbols representing motif identity columns

# Examples
```julia
m_symbols(1)  # Returns [:m1]
m_symbols(2)  # Returns [:m1, :m2]
m_symbols(3)  # Returns [:m1, :m2, :m3]
```
"""
function m_symbols(n::Integer)
    @assert n > 0 "n must be a positive integer"
    return [Symbol("m$i") for i in 1:n]
end

"""
    m_position_symbols(n::Integer) -> Vector{Symbol}

Generate motif position column symbols: [:m1_position, :m2_position, ..., :mn_position]

Used to track the sequence positions where each motif occurs.

# Arguments
- `n`: Number of motifs in the pattern (must be positive)

# Returns
Vector of symbols representing motif position columns

# Examples
```julia
m_position_symbols(1)  # Returns [:m1_position]
m_position_symbols(2)  # Returns [:m1_position, :m2_position]
m_position_symbols(3)  # Returns [:m1_position, :m2_position, :m3_position]
```
"""
function m_position_symbols(n::Integer)
    @assert n > 0 "n must be a positive integer"
    return [Symbol("m$(i)_position") for i in 1:n]
end

"""
    d_symbols(n::Integer) -> Vector{Symbol}

Generate distance column symbols between consecutive motifs.

For n motifs, there are (n-1) distance measurements between consecutive pairs.

# Arguments
- `n`: Number of motifs in the pattern (must be ≥ 1)

# Returns
Vector of symbols representing inter-motif distance columns

# Examples
```julia
d_symbols(1)  # Returns [] (no distances for single motif)
d_symbols(2)  # Returns [:d12] (distance between m1 and m2)
d_symbols(3)  # Returns [:d12, :d23] (distances between m1-m2 and m2-m3)
d_symbols(4)  # Returns [:d12, :d23, :d34]
```

# Note
Distance symbols follow the pattern :d{i}{i+1} where:
- d12 = distance from motif 1 to motif 2
- d23 = distance from motif 2 to motif 3
- etc.
"""
function d_symbols(n::Integer)
    @assert n ≥ 1 "n must be >= 1"
    n == 1 && return Symbol[]  # No distances for single motif
    return [Symbol("d$i$(i+1)") for i in 1:(n-1)]
end

export m_symbols, m_position_symbols, d_symbols
