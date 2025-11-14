"""
Single mutation region analysis and filtering.
"""

"""
    filter_and_rename_for_mutagenesis(contributions_df, data; pareto_rank=3, split_by_sign=true)

Simple wrapper to filter by Pareto rank and rename columns for mutagenesis analysis.
Just automates these 4 lines:
1. Filter to common length
2. Group and apply Pareto filtering  
3. Convert back to DataFrame
4. Rename columns (:filter_index → :m1, :position → :m1_position)

# Example
```julia
df_ready = filter_and_rename_for_mutagenesis(contributions_df, data; pareto_rank=3)
```
"""
function filter_and_rename_for_mutagenesis(contributions_df::DataFrame, data; 
                                           pareto_rank::Int=3, 
                                           split_by_sign::Bool=true)
    # Filter to most common length
    if isnothing(data.raw_data.most_common_length_indices)
        df_filtered = contributions_df
    else
        df_filtered = filter(x -> x.data_pt_index ∈ data.raw_data.most_common_length_indices, 
                        contributions_df)
    end
    
    # Group and apply Pareto filtering
    gdf = groupby(df_filtered, [:filter_index, :position])
    gdf_pareto = filter_pareto_rank(gdf; rank=pareto_rank, split_by_sign=split_by_sign)
    
    # Convert back to DataFrame
    df_result = DataFrame(gdf_pareto)
    
    # Rename columns
    rename!(df_result, Dict(:filter_index => :m1, :position => :m1_position))
    
    return df_result
end

export filter_and_rename_for_mutagenesis
