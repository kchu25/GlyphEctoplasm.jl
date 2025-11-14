
"""
    box_scatter_template(banzhaf_values; 
                         mode=DefaultMode, 
                         title_str="", 
                         box_color=:lightgray,
                         xlim=nothing)

Create a box-scatter figure with standard styling. Returns a Figure object.

# Arguments
- `banzhaf_values`: Vector of values to visualize

# Keyword Arguments
- `mode::DisplayMode`: Display mode preset (DefaultMode, WebDisplayMode, PublicationMode)
- `title_str::String`: Y-axis label/title (default: "")
- `box_color`: Color for the box plot (default: :lightgray)
- `xlim::Union{Nothing,Tuple{Real,Real}}`: X-axis limits (default: nothing for auto)

# Examples
```julia
# Default mode with title
f = box_scatter_template(banzhaf_values; title_str="Filter 1")

# Web display mode (no legend, no labels)
f = box_scatter_template(banzhaf_values; mode=WebDisplayMode)

# Web mode with uniform x-axis limits
f = box_scatter_template(banzhaf_values; mode=WebDisplayMode, xlim=(-0.1, 1.0))
```
"""
function box_scatter_template(banzhaf_values; 
                             mode::DisplayMode=DefaultMode,
                             title_str::String="", 
                             box_color=:lightgray,
                             xlim::Union{Nothing,Tuple{<:Real,<:Real}}=nothing)
    # Small wrapper that returns a Figure for a given dataset.
    f, _ = create_simple_boxplot_figure(banzhaf_values;
        mode=mode,
        figsize=(725, 115),
        ylabel=title_str,
        ylabelsize=12,
        figure_padding=(5, 5, 2, 3),
        xlim=xlim,
        box_config=BoxPlotConfig(
            box_height=0.1,
            box_color=box_color,
            whisker_cap_height=0.1
        ))
    return f
end

"""
    save_box_scatter(banzhaf_values, fp; 
                     mode=DefaultMode,
                     title="", 
                     box_color=:lightgray,
                     xlim=nothing,
                     px_per_unit=1)

Render and save a box-scatter figure to the specified file path.

# Arguments
- `banzhaf_values`: Vector of values to visualize
- `fp::AbstractString`: File path to save the figure

# Keyword Arguments
- `mode::DisplayMode`: Display mode preset (DefaultMode, WebDisplayMode, PublicationMode)
- `title::String`: Y-axis label/title (default: "")
- `box_color`: Color for the box plot (default: :lightgray)
- `xlim::Union{Nothing,Tuple{Real,Real}}`: X-axis limits (default: nothing for auto)
- `px_per_unit::Real`: Pixels per unit for output resolution (default: 1)
"""
function save_box_scatter(
    banzhaf_values, fp::AbstractString; 
    mode::DisplayMode=DefaultMode,
    title::String="", 
    box_color=:lightgray,
    xlim::Union{Nothing,Tuple{<:Real,<:Real}}=nothing,
    px_per_unit::Real=1)
    # Generalized saver: render a box-scatter figure and save to `fp`.
    f = box_scatter_template(banzhaf_values; 
                            mode=mode, 
                            title_str=title, 
                            box_color=box_color,
                            xlim=xlim)
    save(fp, f, px_per_unit=px_per_unit)
end

# Convenience wrappers kept for backwards compatibility

"""
    save_box_scatter_distance_default(banzhaf_values, fp, mode=DefaultMode; xlim=nothing)

Save a box-scatter with no title and an orange box color.

This lightweight wrapper accepts an optional third (positional) argument
`mode::DisplayMode` which will be forwarded to `save_box_scatter`. Use
`WebDisplayMode` to produce the minimal web-friendly output (no legend,
no axis labels).

Examples
    save_box_scatter_distance_default(vals, "out.png")
    save_box_scatter_distance_default(vals, "out_web.png", WebDisplayMode)
    save_box_scatter_distance_default(vals, "out_web.png", WebDisplayMode; xlim=(-0.1, 1.0))
"""
save_box_scatter_distance_default(banzhaf_values, fp::AbstractString, _mode_; xlim=nothing) = 
    save_box_scatter(banzhaf_values, fp; mode=_mode_, box_color=(:orange, 0.98), xlim=xlim)

"""
    save_box_scatter_distance_relaxed(banzhaf_values, fp, mode=DefaultMode; xlim=nothing)

Save a box-scatter labeled "Distance Relaxed" using a light-gray box color.

Accepts an optional `mode::DisplayMode` (third positional argument) which is
forwarded to `save_box_scatter`. This lets you request the web-friendly
variant when generating assets for the UI.

Examples
    save_box_scatter_distance_relaxed(vals, "relaxed.png")
    save_box_scatter_distance_relaxed(vals, "relaxed_web.png", WebDisplayMode)
    save_box_scatter_distance_relaxed(vals, "relaxed_web.png", WebDisplayMode; xlim=(-0.1, 1.0))
"""
save_box_scatter_distance_relaxed(banzhaf_values, fp::AbstractString, _mode_; xlim=nothing) = 
    save_box_scatter(banzhaf_values, fp; mode=_mode_, title="Distance Relaxed", xlim=xlim)

"""
    save_box_scatter_distance_fixed(banzhaf_values, fp, mode=DefaultMode; xlim=nothing)

Save a box-scatter labeled "Distance Fixed" using an orange box color.

Accepts an optional `mode::DisplayMode` (passed as the third positional
argument) which is forwarded to `save_box_scatter` for producing web-friendly
or publication variants.

Examples
    save_box_scatter_distance_fixed(vals, "fixed.png")
    save_box_scatter_distance_fixed(vals, "fixed_web.png", WebDisplayMode)
    save_box_scatter_distance_fixed(vals, "fixed_web.png", WebDisplayMode; xlim=(-0.1, 1.0))
"""
save_box_scatter_distance_fixed(banzhaf_values, fp::AbstractString, _mode_; xlim=nothing) = 
    save_box_scatter(banzhaf_values, fp; mode=_mode_, title="Distance Fixed", box_color=(:orange, 0.98), xlim=xlim)
