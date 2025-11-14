"""
Pareto ranking and sorting utilities for motif prioritization.
"""

"""
    filter_pareto_rank(gdf::GroupedDataFrame; 
                       rank::Int=1,
                       split_by_sign::Bool=true,
                       maximize_objectives::Vector{Symbol}=[:abs_median_banzhaf, :count])

Filter grouped DataFrame to keep only groups up to a specified Pareto rank based on multiple objectives.

# Arguments
- `gdf`: GroupedDataFrame to filter
- `rank`: Maximum Pareto rank to keep (1 = best, 2 = second tier, etc.)
- `split_by_sign`: If true, apply Pareto filtering separately to positive and negative median groups
- `maximize_objectives`: Vector of objectives to maximize. Default is [:abs_median_banzhaf, :count]

# Returns
New GroupedDataFrame containing only groups with Pareto rank ≤ specified rank

# Pareto Rank Definition
- Rank 1: Groups that are not dominated by any other group
- Rank 2: Groups dominated only by rank-1 groups
- Rank k: Groups dominated only by groups of rank < k

Group A dominates B if: A is better than or equal to B on all objectives, 
and strictly better on at least one.

# Examples
```julia
# Keep only best Pareto-optimal groups
gdf_filtered = filter_pareto_rank(gdf; rank=1)

# Keep top 2 tiers of Pareto-optimal groups
gdf_filtered = filter_pareto_rank(gdf; rank=2)

# Apply separately to positive and negative groups
gdf_filtered = filter_pareto_rank(gdf; rank=1, split_by_sign=true)

# Don't split by sign, treat all together
gdf_filtered = filter_pareto_rank(gdf; rank=1, split_by_sign=false)
```
"""
function filter_pareto_rank(gdf::GroupedDataFrame; 
                            rank::Int=1,
                            split_by_sign::Bool=true,
                            maximize_objectives::Vector{Symbol}=[:abs_median_banzhaf, :count])
    
    @assert rank >= 1 "Pareto rank must be >= 1"
    
    # Compute summary statistics for all groups
    df_summary = combine(gdf,
        :banzhaf => median => :median_banzhaf,
        :banzhaf => (x -> abs(median(x))) => :abs_median_banzhaf,
        nrow => :count
    )
    
    # Add group keys
    keycols = keys(gdf)[1] |> keys |> collect
    for col in keycols
        df_summary[!, col] = [k[col] for k in keys(gdf)]
    end
    
    if split_by_sign
        # Split into positive and negative groups
        df_pos = filter(row -> row.median_banzhaf > 0, df_summary)
        df_neg = filter(row -> row.median_banzhaf < 0, df_summary)
        
        # Find Pareto ranks in each
        pareto_pos_indices = find_pareto_ranks(df_pos, maximize_objectives, rank)
        pareto_neg_indices = find_pareto_ranks(df_neg, maximize_objectives, rank)
        
        # Combine the keys
        pareto_keys_pos = [NamedTuple(row[keycols]) for row in eachrow(df_pos[pareto_pos_indices, :])]
        pareto_keys_neg = [NamedTuple(row[keycols]) for row in eachrow(df_neg[pareto_neg_indices, :])]
        pareto_keys = vcat(pareto_keys_pos, pareto_keys_neg)
    else
        # Find Pareto ranks across all groups
        pareto_indices = find_pareto_ranks(df_summary, maximize_objectives, rank)
        pareto_keys = [NamedTuple(row[keycols]) for row in eachrow(df_summary[pareto_indices, :])]
    end
    
    # Filter original grouped dataframe to keep only Pareto-optimal groups
    # Create a new dataframe with only the Pareto-optimal rows
    pareto_rows = DataFrame()
    for k in pareto_keys
        append!(pareto_rows, DataFrame(gdf[k]))
    end
    
    # Return as grouped dataframe with same grouping
    return groupby(pareto_rows, keycols)
end

# Backward compatibility alias
filter_pareto_rank1(gdf::GroupedDataFrame; kwargs...) = filter_pareto_rank(gdf; rank=1, kwargs...)

"""
    find_pareto_ranks(df::DataFrame, objectives::Vector{Symbol}, max_rank::Int) -> Vector{Int}

Find indices of rows with Pareto rank ≤ max_rank in a DataFrame based on given objectives.
All objectives are assumed to be maximized.

# Arguments
- `df`: DataFrame with objective columns
- `objectives`: Vector of column symbols representing objectives to maximize
- `max_rank`: Maximum rank to include (1 = non-dominated only, 2 = first two tiers, etc.)

# Returns
Vector of row indices that have Pareto rank ≤ max_rank
"""
function find_pareto_ranks(df::DataFrame, objectives::Vector{Symbol}, max_rank::Int)
    n = nrow(df)
    if n == 0
        return Int[]
    end
    
    ranks = zeros(Int, n)
    available = trues(n)  # Track which rows haven't been assigned a rank yet
    
    for current_rank in 1:max_rank
        # Find non-dominated rows among the available ones
        rank_indices = Int[]
        
        for i in 1:n
            !available[i] && continue
            
            is_dominated = false
            
            # Check if i is dominated by any other available row
            for j in 1:n
                (i == j || !available[j]) && continue
                
                # Check if j dominates i
                all_geq = true
                any_greater = false
                
                for obj in objectives
                    val_i = df[i, obj]
                    val_j = df[j, obj]
                    
                    if val_j < val_i
                        all_geq = false
                        break
                    elseif val_j > val_i
                        any_greater = true
                    end
                end
                
                if all_geq && any_greater
                    is_dominated = true
                    break
                end
            end
            
            if !is_dominated
                push!(rank_indices, i)
                ranks[i] = current_rank
            end
        end
        
        # Mark assigned rows as unavailable for next iteration
        for idx in rank_indices
            available[idx] = false
        end
        
        # If no rows were assigned this rank and we still have available rows,
        # we need to continue to higher ranks, but we've reached max_rank
        if isempty(rank_indices) && any(available)
            break
        end
    end
    
    return findall(x -> 1 <= x <= max_rank, ranks)
end

"""
    build_sorted_keys_and_maps(gdf, keycols; sort_rev=true, pareto_rank=nothing)

Build sorted keys with median banzhaf values and occurrence counts from grouped DataFrame.
Returns `(sorted_keys, median_map, mean_map, count_map, list_of_banzhafs)` sorted by median banzhaf (descending by default).
"""
function build_sorted_keys_and_maps(gdf, keycols; sort_rev=true, pareto_rank=nothing)

    if !isnothing(pareto_rank)
        gdf = filter_pareto_rank(gdf; rank=pareto_rank, split_by_sign=true,
                                 maximize_objectives=[:abs_median_banzhaf, :count])
    end

    # create a summary df with median banzhaf and count
    df_summary = combine(gdf, 
        :banzhaf => median => :banzhaf_median, 
        :banzhaf => mean => :banzhaf_mean,
        nrow => :count)
    sort!(df_summary, :banzhaf_median, rev=sort_rev)
    # build a small DataFrame with only the key columns (works for Symbol or Vector{Symbol})
    keycols_vec = isa(keycols, Symbol) ? [keycols] : keycols
    keydf = df_summary[:, keycols_vec]
    sorted_keys = [NamedTuple(r) for r in eachrow(keydf)]
    median_map = Dict(zip(sorted_keys, df_summary.banzhaf_median))
    mean_map = Dict(zip(sorted_keys, df_summary.banzhaf_mean))
    count_map = Dict(zip(sorted_keys, df_summary.count))
    list_of_banzhafs = Dict(k=>gdf[k].banzhaf for k in sorted_keys);

    return sorted_keys, median_map, mean_map, count_map, list_of_banzhafs
end

export filter_pareto_rank, filter_pareto_rank1, find_pareto_ranks, build_sorted_keys_and_maps
