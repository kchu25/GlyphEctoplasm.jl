"""
File I/O operations for motif data: saving logos, influence plots, positional info, and MEME files.
"""

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
function save_motif_logo(pfm, png_path, median_val; dpi=65, alpha=1.0, highlighted_regions=nothing)
    if highlighted_regions === nothing
        save_logoplot(pfm, png_path; dpi=dpi, alpha=alpha, uniform_color=false, pos=median_val > 0)
    else
        save_logoplot(pfm, png_path; dpi=dpi, alpha=alpha, uniform_color=false, 
                     pos=median_val > 0, highlighted_regions=highlighted_regions)
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
                             use_rna=false, relaxed_median=nothing,
                             show_meme_and_csv=true)

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

    return [influence_median, construct_str, consensus_str, meme_csv_combined, "", ""]
end

export build_motif_paths, save_motif_logo, save_influence_plot, save_positional_info, build_metadata_texts
