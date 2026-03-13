""" 
File I/O operations for motif data: saving logos, influence plots, positional info, and MEME files.
"""

"""
    mean_knn_within_group_1d(sorted_vals, k)

Compute the mean k-NN distance **within** a sorted group of values.
For each point, scan left/right among its group-mates to find the k
nearest neighbors, then average the k-th NN distance across all points.

Since the group is sorted, neighbors are always adjacent — this is O(m * k).
"""
function mean_knn_within_group_1d(sorted_vals::AbstractVector{Float64}, k::Int)
    m = length(sorted_vals)
    k_clamped = min(k, m - 1)
    sum_knn_dists = 0.0
    for i in 1:m
        left_idx, right_idx = i - 1, i + 1
        knn_dist = 0.0
        for _ in 1:k_clamped
            dist_left = left_idx >= 1 ? sorted_vals[i] - sorted_vals[left_idx] : Inf
            dist_right = right_idx <= m ? sorted_vals[right_idx] - sorted_vals[i] : Inf
            if dist_left <= dist_right
                knn_dist = dist_left
                left_idx -= 1
            else
                knn_dist = dist_right
                right_idx += 1
            end
        end
        sum_knn_dists += knn_dist
    end
    return sum_knn_dists / m
end

"""
    nnd_permutation_test_1d(subpop_positions, background; k=10, B=1_000, seed=42)

Test whether the points at `subpop_positions` within `background` are more
tightly clustered **among themselves** than a random draw of the same size,
using within-group k-NN distances and a permutation test.

**Observed statistic**: sort the m subpop values → compute mean k-NN
distance within that sorted group.

**Null distribution**: draw m values at random from background → sort →
compute the same within-group mean k-NN distance.  Repeat B times.

p-value = fraction of null draws whose within-group mNND ≤ observed.

# Arguments
- `subpop_positions`: Integer vector of **positional indices** (1-based)
  into `background` identifying the subpopulation members.
- `background`: Full background labels (real-valued vector, length N).

# Keyword Arguments
- `k::Int=10`: Number of nearest neighbors.
- `B::Int=1_000`: Number of permutations.
- `seed::Int=42`: Random seed for reproducibility.

# Returns
NamedTuple `(k, obs_mNND, p_value)`.
"""
function nnd_permutation_test_1d(
    subpop_positions::AbstractVector{<:Integer},
    background::AbstractVector{<:Real};
    k::Int = 10,
    B::Int = 1_000,
    seed::Int = 42
)
    rng = Random.MersenneTwister(seed)
    bg_vals = Vector{Float64}(background)
    n_bg = length(bg_vals)
    n_subpop = length(subpop_positions)
    k_clamped = min(k, n_subpop - 1)   # can only have m-1 neighbors within a group of m

    # Observed: sort subpop values, compute within-group mNND
    subpop_vals = sort!([bg_vals[i] for i in subpop_positions])
    obs_mNND = mean_knn_within_group_1d(subpop_vals, k_clamped)

    # Null: draw m random values from background, compute within-group mNND
    count_significant = 0
    perm_indices = Vector{Int}(undef, n_bg)
    null_vals = Vector{Float64}(undef, n_subpop)

    for _ in 1:B
        Random.randperm!(rng, perm_indices)
        @inbounds for i in 1:n_subpop
            null_vals[i] = bg_vals[perm_indices[i]]
        end
        sort!(null_vals)
        null_mNND = mean_knn_within_group_1d(null_vals, k_clamped)
        if null_mNND <= obs_mNND
            count_significant += 1
        end
    end

    p_value = count_significant / B
    return (; k=k_clamped, obs_mNND, p_value)
end

"""
    nnd_sensitivity_batch_1d(subpop_positions, background; ks, B=1_000, seed=42)

Batched NND permutation test across multiple k values for a single motif.
Each k is tested using **within-group** k-NN distances — the subpopulation's
mean k-NN distance is computed among its own m members, then compared to
null draws of m random background values.

Shares the same random draws across all k values for efficiency.

# Arguments
- `subpop_positions`: Integer vector of **positional indices** (1-based)
  into `background` identifying the subpopulation members.
- `background`: Full background labels (real-valued vector, length N).

# Keyword Arguments
- `ks::Vector{Int}`: k values to evaluate.
- `B::Int=1_000`: Number of permutations.
- `seed::Int=42`: Random seed for reproducibility.

# Returns
Vector of NamedTuples `(k, obs_mNND, p_value)`, one per entry in `ks`.
"""
function nnd_sensitivity_batch_1d(
    subpop_positions::AbstractVector{<:Integer},
    background::AbstractVector{<:Real};
    ks::Vector{Int},
    B::Int = 1_000,
    seed::Int = 42
)
    rng = Random.MersenneTwister(seed)
    bg_vals = Vector{Float64}(background)
    n_bg = length(bg_vals)
    n_subpop = length(subpop_positions)

    # Clamp k values to m-1 (max neighbors within a group of m)
    ks_clamped = [min(kv, n_subpop - 1) for kv in ks]
    n_ks = length(ks)

    # Observed: sort subpop values, compute within-group mNND for each k
    subpop_vals = sort!([bg_vals[i] for i in subpop_positions])
    obs_mNNDs = [mean_knn_within_group_1d(subpop_vals, ku) for ku in ks_clamped]

    # Null: draw m random values, sort, compute within-group mNND for all k
    counts_significant = zeros(Int, n_ks)
    perm_indices = Vector{Int}(undef, n_bg)
    null_vals = Vector{Float64}(undef, n_subpop)

    for _ in 1:B
        Random.randperm!(rng, perm_indices)
        @inbounds for i in 1:n_subpop
            null_vals[i] = bg_vals[perm_indices[i]]
        end
        sort!(null_vals)

        # Evaluate all k values for this draw
        for (j, ku) in enumerate(ks_clamped)
            null_mNND = mean_knn_within_group_1d(null_vals, ku)
            if null_mNND <= obs_mNNDs[j]
                counts_significant[j] += 1
            end
        end
    end

    inv_B = 1.0 / B
    return [(; k=ks_clamped[j], obs_mNND=obs_mNNDs[j], p_value=counts_significant[j] * inv_B) for j in 1:n_ks]
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

export build_motif_paths, save_motif_logo, save_influence_plot, save_positional_info, build_metadata_texts, nnd_permutation_test_1d, nnd_sensitivity_batch_1d
