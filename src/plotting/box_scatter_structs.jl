# Display mode enumeration for different use cases
@enum DisplayMode begin
    DefaultMode      # Standard display with all labels and legend
    WebDisplayMode   # Minimal display for web (no legend, no labels)
    PublicationMode  # Reserved for future publication-ready styling
end

# Configuration struct for box plot styling
Base.@kwdef struct BoxPlotConfig
    box_height::Float64 = 0.165
    box_offset_y::Float64 = -0.075
    box_color = (:orange, 0.98)
    line_width::Float64 = 1.5
    whisker_cap_height::Union{Float64,Nothing} = nothing  # If nothing, uses box_height * 1.2
end

# Configuration struct for scatter plot styling
Base.@kwdef struct ScatterConfig
    jitter_scale::Float64 = 0.05  # Reduced from 0.25 to prevent dots touching borders
    color = :maroon
    alpha::Float64 = 0.6
    markersize::Float64 = 2.5
    marker = :circle
end

# Configuration struct for legend styling
Base.@kwdef struct LegendConfig
    margin::Tuple{Int,Int,Int,Int} = (10, 10, 35, 10)  # (top, right, bottom, left)
    backgroundcolor = (:white, 0.9)
    framecolor = (:black, 0.0)
    framewidth::Float64 = 0.6
    labelsize::Int = 13
    patchsize::Tuple{Int,Int} = (0, 0)
    padding::Tuple{Int,Int,Int,Int} = (6, 6, 4, 4)
    halign::Symbol = :center  # Horizontal alignment (:left, :center, :right)
    valign::Symbol = :top     # Vertical alignment (:top, :center, :bottom)
    tellheight::Bool = false  # Whether legend controls its row height
    tellwidth::Bool = false   # Whether legend controls its column width
    enabled::Bool = true      # Whether to show the legend at all
end

# Configuration struct for layout
Base.@kwdef struct LayoutConfig
    row_legend::Float64 = 0.2
    row_box::Float64 = 0.7
    row_scatter::Float64 = 0.5
    row_spacer::Float64 = 0.0001
    row_gap::Union{Float64,Nothing} = nothing  # Gap between rows (uses rowgap! if set)
    xlabel::String = "Influence (Banzhaf indices)"
    xlabelsize::Int = 17
    show_xlabel::Bool = true  # Whether to show x-axis label
end

"""
    get_mode_configs(mode::DisplayMode; xlim=nothing)

Get preset configurations for a specific display mode.

# Arguments
- `mode::DisplayMode`: The display mode (DefaultMode, WebDisplayMode, or PublicationMode)
- `xlim::Union{Nothing,Tuple{Real,Real}}`: Optional x-axis limits

# Returns
- `NamedTuple` with fields: `legend_config`, `layout_config`, `show_ylabel`, `xlim`

# Examples
```julia
configs = get_mode_configs(DefaultMode)
configs = get_mode_configs(WebDisplayMode; xlim=(-0.1, 1.1))
```
"""
function get_mode_configs(mode::DisplayMode; xlim::Union{Nothing,Tuple{<:Real,<:Real}}=nothing)
    if mode == DefaultMode
        return (
            legend_config = LegendConfig(
                labelsize=13,
                halign=:center,
                valign=:top,
                margin=(10, 10, 10, 0),
                tellheight=false,
                enabled=true
            ),
            layout_config = LayoutConfig(
                row_legend=0.22,
                row_gap=2.0,
                row_box=0.48,
                row_scatter=0.30,
                xlabel="Influence (Banzhaf indices)",
                xlabelsize=14,
                show_xlabel=true
            ),
            show_ylabel = true,
            xlim = xlim
        )
    elseif mode == WebDisplayMode
        return (
            legend_config = LegendConfig(
                enabled=false  # No legend in web display mode
            ),
            layout_config = LayoutConfig(
                row_legend=0.0001,  # Minimal space since no legend
                row_gap=2.0,
                row_box=0.70,
                row_scatter=0.30,
                xlabel="",
                xlabelsize=14,
                show_xlabel=false  # No x-axis label
            ),
            show_ylabel = false,  # No y-axis label
            xlim = xlim
        )
    elseif mode == PublicationMode
        # Placeholder for future publication styling
        return (
            legend_config = LegendConfig(
                labelsize=14,
                halign=:center,
                valign=:top,
                margin=(10, 10, 10, 0),
                tellheight=false,
                enabled=true
            ),
            layout_config = LayoutConfig(
                row_legend=0.22,
                row_gap=2.0,
                row_box=0.48,
                row_scatter=0.30,
                xlabel="Influence (Banzhaf indices)",
                xlabelsize=16,
                show_xlabel=true
            ),
            show_ylabel = true,
            xlim = xlim
        )
    else
        error("Unknown DisplayMode: $mode")
    end
end
