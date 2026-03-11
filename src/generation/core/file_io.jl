""" 
File I/O operations for motif data: saving logos, influence plots, positional info, and MEME files.
"""

"""
    precompute_knn_dists_1d(sorted_pool, k_max)

For every position in the sorted array, precompute the distances to its
1st, 2nd, …, k_max-th nearest neighbors by scanning left/right.
Returns a matrix `D` of size `(N, k_max)` where `D[i, j]` is the
j-th nearest-neighbor distance for position `i`.

Since the pool is sorted, neighbors are always adjacent — this is O(N * k_max).
"""
function precompute_knn_dists_1d(sorted_pool::Vector{Float64}, k_max::Int)
    N = length(sorted_pool)
    D = zeros(Float64, N, k_max)
    for i in 1:N
        li, ri = i - 1, i + 1
        for j in 1:k_max
            dl = li >= 1 ? sorted_pool[i] - sorted_pool[li] : Inf  # sorted, so no abs needed
            dr = ri <= N ? sorted_pool[ri] - sorted_pool[i] : Inf
            if dl <= dr
                D[i, j] = dl
                li -= 1
            else
                D[i, j] = dr
                ri += 1
            end
        end
    end
    return D
end

"""
    nnd_permutation_test_1d(subpop, background; k=5, B=10_000, seed=42)

Test whether `subpop` points are more tightly clustered than expected under
random labeling, using k-NN distances and a permutation test.
Optimized for 1D: sort once, precompute all k-NN distances, then each
permutation is O(m) table lookups.

# Arguments
- `subpop`: Subpopulation labels (any real-valued vector)
- `background`: Full background labels (any real-valued vector)

# Keyword Arguments
- `k::Int=5`: Number of nearest neighbors.
- `B::Int=10_000`: Number of permutations
- `seed::Int=42`: Random seed for reproducibility

# Returns
NamedTuple `(k, obs_mNND, p_value)` where:
- `k`: The k value used
- `obs_mNND`: Observed mean k-NN distance for the subpopulation
- `p_value`: Fraction of permutations with mean k-NN distance ≤ observed
"""
function nnd_permutation_test_1d(
    subpop::AbstractVector{<:Real},
    background::AbstractVector{<:Real};
    k::Int = 5,
    B::Int = 10_000,
    seed::Int = 42
)
    rng = Random.MersenneTwister(seed)
    pool = Vector{Float64}(vcat(subpop, background))
    sp = sortperm(pool)
    sorted_pool = pool[sp]
    N = length(pool)
    m = length(subpop)

    k_use = min(k, N - 1)

    # Precompute k-NN distance table for every position: O(N * k_use)
    D = precompute_knn_dists_1d(sorted_pool, k_use)

    # Map subpopulation to sorted-array positions
    rank_of = zeros(Int, N)
    for (rank, orig) in enumerate(sp)
        rank_of[orig] = rank
    end
    sub_ranks = sort([rank_of[i] for i in 1:m])

    # Helper: mean of k-th NN distance for given indices
    function compute_mNND(indices::Vector{Int})
        s = 0.0
        for idx in indices
            s += D[idx, k_use]
        end
        return s / m
    end

    # Observed statistic
    obs_mNND = compute_mNND(sub_ranks)

    # Permutation null
    count_leq = 0
    full_perm = Vector{Int}(undef, N)

    for _ in 1:B
        Random.randperm!(rng, full_perm)
        perm_indices = sort!(collect(@view(full_perm[1:m])))
        null_mNND = compute_mNND(perm_indices)
        if null_mNND <= obs_mNND
            count_leq += 1
        end
    end

    p_value = count_leq / B
    return (; k=k_use, obs_mNND, p_value)
end

"""
    build_motif_paths(name_base::AbstractString, save_folder::AbstractString, motif_type::AbstractString)

Build absolute and relative file paths for motif output files (PNG logo, CSV positions, MEME format).

# Arguments
- `name_base`: Base filename (without extension)
- `save_folder`: Absolute path to save folder
- `motif_type`: Motif type subdirectory name (e.g., "singletons", "pairs_positive")

# Returns
Named tuple with absolute (.abs) and relative (.rel) paths for:
- `png`: Logo image file
- `influence`: Influence plot image file  
- `csv`: Positional information CSV
- `meme`: MEME format motif file

# Example
```julia
paths = build_motif_paths("filter_42", "/path/to/tmp3/singletons", "singletons")
# paths.png.abs => "/path/to/tmp3/singletons/filter_42.png"
# paths.png.rel => "singletons/filter_42.png"
```
"""
function build_motif_paths(name_base::AbstractString, save_folder::AbstractString, motif_type::AbstractString)
    png_fn = name_base * ".png"
    influence_fn = name_base * "_influence.png"
    csv_fn = name_base * ".csv"
    meme_fn = name_base * ".meme"
    
    return (
        png = (abs = joinpath(save_folder, png_fn), rel = joinpath(motif_type, png_fn)),
        csv = (abs = joinpath(save_folder, csv_fn), rel = joinpath(motif_type, csv_fn)),
        meme = (abs = joinpath(save_folder, meme_fn), rel = joinpath(motif_type, meme_fn)),
        influence = (abs = joinpath(save_folder, influence_fn), rel = joinpath(motif_type, influence_fn))
    )
end

"""
    save_motif_logo(pfm, png_path, median_val; dpi=65, alpha=1.0, highlighted_regions=nothing)

Save motif logo plot as PNG.
"""
function save_motif_logo(pfm, png_path, median_val; dpi=65, alpha=1.0, highlighted_regions=nothing, 
    rna=false)
    if highlighted_regions === nothing
        save_logoplot(pfm, png_path; 
            dpi=dpi, alpha=alpha, uniform_color=false, pos=median_val > 0, 
            _margin_=0(EntroPlots.Plots.mm), tight=true, yaxis=false, xaxis=false, 
            rna=rna
            )
    else
        save_logoplot(pfm, png_path; dpi=dpi, alpha=alpha, uniform_color=false, 
            pos=median_val > 0, highlighted_regions=highlighted_regions,
            _margin_=0(EntroPlots.Plots.mm), tight=true, yaxis=false, xaxis=false, rna=rna
            )
    end
end

"""
    save_influence_plot(banzhafs, influence_path; highlighted_regions=nothing, xlim=nothing)

Save influence box/scatter plot.
"""
function save_influence_plot(banzhafs, influence_path; highlighted_regions=nothing, xlim=nothing)
    if highlighted_regions === nothing
        BanzhafPlots.save_box_scatter_distance_default(banzhafs, influence_path, BanzhafPlots.WebDisplayMode; xlim=xlim)
    else
        BanzhafPlots.save_box_scatter_distance_fixed(banzhafs, influence_path, BanzhafPlots.WebDisplayMode; xlim=xlim)
    end
end

"""
    save_positional_info(pos_data, paths, filter_len)

Save positional information CSV. Handles both DataFrame (singletons) and 
Vector of tuples (multi-motifs) formats.
"""
function save_positional_info(pos_data, paths, filter_len)
    if pos_data === nothing
        return
    end
    
    if isa(pos_data, Vector)  # flat_windows for multi-motifs
        save_pos_info_as_csv(pos_data, paths.csv.abs)
    else  # DataFrame/SubDataFrame for singletons
        save_pos_info_as_csv(pos_data, filter_len, paths.csv.abs)
    end
end

"""
    build_metadata_texts(pfm, paths, median_val, count_val; 
                        use_rna=false, relaxed_median=nothing, show_meme_and_csv=true)

Build text entries for JSON metadata display.
Returns array of formatted strings for influence median, construction info, 
consensus, and file links.
"""
function build_metadata_texts(pfm, paths, median_val, count_val; 
                             interaction_summary_mode_str=nothing,
                             use_rna=false, relaxed_median=nothing,
                             show_meme_and_csv=true,
                             nnd_result=nothing)

    @assert !isnothing(count_val) "number of counts used to construct the logo must be provided"
    if !isnothing(pfm)                             
        pfm_len = size(pfm, 2)
        construct_str = string("The PWM was constructed from ", count_val, " sequences of length ", pfm_len, ".")
    else
        construct_str = string("The PWM was constructed from ", count_val, " sequences")
    end

    if !isnothing(pfm)
        consensus_str = "Consensus: "*pfm2consensus(pfm; rna=use_rna)
    else
        consensus_str = ""
    end

    if show_meme_and_csv
        meme_link = string("<a href=\"", paths.meme.rel, "\">.meme file</a>")
        csv_str = fill_csv_link(paths.csv.rel)
        meme_csv_combined = string(meme_link, " | ", csv_str)
    else
        meme_csv_combined = ""
    end

    # Build influence median string(s)
    if relaxed_median !== nothing
        # Multi-motif case: show both relaxed and fixed distance medians
        influence_median = string(
            "Influence Median (Relaxed): <strong>", round(relaxed_median, digits=2),
            "</strong> | (Fixed): <strong>", round(median_val, digits=2), "</strong>"
        )
    else
        # Singleton case: just show the median
        influence_median = string("Influence Median: <strong>", round(median_val, digits=2), "</strong>")
    end

    interact_part = begin 
        if isnothing(interaction_summary_mode_str)
            ""
        else 
            interaction_summary_mode_str
        end
    end
    
    # Build NND test display row
    variance_row = if isnothing(nnd_result)
        ""
    else
        k_used = nnd_result.k
        # Build NND observed mean distance string
        nnd_str = string("mean ", k_used, "-NN dist: <strong>", @sprintf("%.3f", nnd_result.obs_mNND), "</strong>")

        # Build NND p-value string with significance coloring
        pv = nnd_result.p_value
        pval_str = if pv > NND_PVALUE_CUTOFF
            string("<span style=\"color: gray;\">p-value (NND): ", @sprintf("%.2g", pv), "</span>")
        else
            string("p-value (NND): <strong>", @sprintf("%.2g", pv), "</strong>")
        end

        string(nnd_str, " &nbsp;|&nbsp; ", pval_str)
    end

    return [influence_median, construct_str, consensus_str, meme_csv_combined, interact_part, variance_row]
end

export build_motif_paths, save_motif_logo, save_influence_plot, save_positional_info, build_metadata_texts, nnd_permutation_test_1d
