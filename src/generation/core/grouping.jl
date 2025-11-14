"""
Grouping column specification builders for motif analysis.
"""

"""
    build_grouping_columns(criterion::Symbol; motif_size::Union{Nothing,Integer}=nothing)

Build grouping column specification based on analysis criterion.

# Criteria
- `:filter_index` - Group by convolutional filter index (singletons)
- `:filter_and_position` - Group by filter and sequence position (mutagenesis)
- `:motif_positions` - Group by motif occurrence positions only (multi-motifs, ignores which filters)
- `:motifs` - Group by motif filter identities (which filters fired)
- `:distances` - Group by inter-motif distances
- `:motif_identity_and_distance` - Group by filter identities and inter-motif distances (multi-motifs)
- `:mutagenesis` - Complete specification for mutagenesis: filters, distances, AND positions

# Arguments
- `criterion`: Symbol specifying grouping strategy
- `motif_size`: Required for multi-motif criteria (2 for pairs, 3 for triplets, etc.)

# Returns
Vector of column symbols to use with `groupby(df, columns)`

# Examples
```julia
# Singleton analysis: group by which filter activated
build_grouping_columns(:filter_index)  # => [:filter_index]

# Mutagenesis: group by filter AND position
build_grouping_columns(:filter_and_position)  # => [:filter_index, :position]

# Paired motifs: group by which two filters, ignoring position
build_grouping_columns(:motifs, motif_size=2)  # => [:m1, :m2]

# Paired motifs: group by filter pair and distance
build_grouping_columns(:motif_identity_and_distance, motif_size=2)  # => [:m1, :m2, :d12]

# Paired mutations: complete specification
build_grouping_columns(:mutagenesis, motif_size=2)  
# => [:m1, :m2, :d12, :m1_position, :m2_position]
```
"""
function build_grouping_columns(criterion::Symbol; motif_size::Union{Nothing,Integer}=nothing)
    if criterion == :filter_index
        return [:filter_index]
        
    elseif criterion == :filter_and_position
        return [:filter_index, :position]
        
    elseif criterion == :motifs
        isnothing(motif_size) && error("motif_size required for :motifs criterion")
        return m_symbols(motif_size)

    elseif criterion == :distances
        isnothing(motif_size) && error("motif_size required for :distances criterion")
        return d_symbols(motif_size)
        
    elseif criterion == :motif_identity_and_distance
        isnothing(motif_size) && error("motif_size required for :motif_identity_and_distance criterion")
        m_syms = m_symbols(motif_size)
        d_syms = d_symbols(motif_size)
        return vcat(m_syms, d_syms)
        
    elseif criterion == :mutagenesis
        isnothing(motif_size) && error("motif_size required for :mutagenesis criterion")
        if motif_size == 1
            return vcat(m_symbols(1), m_position_symbols(1))
        else
            m_syms = m_symbols(motif_size)
            d_syms = d_symbols(motif_size)
            mp_syms = m_position_symbols(motif_size)
            return vcat(m_syms, d_syms, mp_syms)
        end
        
    elseif criterion == :motif_positions
        isnothing(motif_size) && error("motif_size required for :motif_positions criterion")
        return m_position_symbols(motif_size)

    else
        error("Unknown grouping criterion: $criterion. Valid options: " *
              ":filter_index, :filter_and_position, :motifs, :distances, " *
              ":motif_identity_and_distance, :mutagenesis, :motif_positions")
    end
end

export build_grouping_columns
