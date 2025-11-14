
# Calculate box plot statistics
function calculate_boxplot_stats(data)
    q1, q3 = quantile(data, [0.25, 0.75])
    iqr = q3 - q1
    med = median(data)
    whisker_lo = max(minimum(data), q1 - 1.5 * iqr)
    whisker_hi = min(maximum(data), q3 + 1.5 * iqr)
    
    return (q1=q1, q3=q3, median=med, whisker_lo=whisker_lo, whisker_hi=whisker_hi)
end

# Draw custom box plot
function draw_boxplot!(ax, data; 
                       box_height=0.125, 
                       box_offset_y=0.0,
                       box_color=(:orange, 0.98),
                       line_width=1.5,
                       whisker_cap_height=nothing,  # If nothing, uses box_height * 1.2
                       label="")
    
    stats = calculate_boxplot_stats(data)
    half_height = box_height / 2
    
    # Calculate whisker cap height
    if isnothing(whisker_cap_height)
        cap_half = half_height * 1.2
    else
        cap_half = whisker_cap_height / 2
    end
    
    # Draw box (IQR) - return this for legend
    box_plot = poly!(ax, 
          [stats.q1, stats.q3, stats.q3, stats.q1], 
          [box_offset_y - half_height, box_offset_y - half_height, 
           box_offset_y + half_height, box_offset_y + half_height];
          color=box_color, strokecolor=:black, strokewidth=line_width,
          label=label)
    
    # Draw whiskers
    lines!(ax, [stats.whisker_lo, stats.q1], [box_offset_y, box_offset_y]; 
           color=:black, linewidth=line_width)
    lines!(ax, [stats.q3, stats.whisker_hi], [box_offset_y, box_offset_y]; 
           color=:black, linewidth=line_width)
    
    # Draw whisker caps (only if cap_half > 0)
    if cap_half > 0
        lines!(ax, [stats.whisker_lo, stats.whisker_lo], 
               [box_offset_y - cap_half, box_offset_y + cap_half]; 
               color=:black, linewidth=line_width)
        lines!(ax, [stats.whisker_hi, stats.whisker_hi], 
               [box_offset_y - cap_half, box_offset_y + cap_half]; 
               color=:black, linewidth=line_width)
    end
    
    # Draw median line
    lines!(ax, [stats.median, stats.median], 
           [box_offset_y - half_height, box_offset_y + half_height]; 
           color=:black, linewidth=2.0)
    
    return stats, box_plot
end


"""
    setup_box_axis!(ax, stats, box_offset_y, data_range)

Configure the box plot axis with appropriate limits and hide decorations.
Y-limits are set to ±0.25 to provide buffer for scatter markers (increased from ±0.22).
"""
function setup_box_axis!(ax, stats, box_offset_y, data_range)
    xmin = minimum([stats.whisker_lo]) - 0.02 * data_range
    xmax = maximum([stats.whisker_hi]) + 0.5 * data_range
    xlims!(ax, xmin, xmax)
    ylims!(ax, box_offset_y - 0.25, box_offset_y + 0.25)  # Increased from 0.22 to 0.25 for more buffer
    
    hideydecorations!(ax)
    hidexdecorations!(ax)
    hidespines!(ax)
end

"""
    create_median_legend!(gl, box_plot, median_value, format_median, config::LegendConfig)

Create a legend displaying the median value with customizable styling.
Uses proper Unicode minus sign (−) for negative values for better typography.
"""
function create_median_legend!(gl, box_plot, median_value, format_median, config::LegendConfig)
    # Skip legend creation if disabled
    if !config.enabled
        return nothing
    end
    
    formatted_value = format_median(median_value)
    # Replace hyphen-minus with proper Unicode minus sign for negative numbers
    formatted_value = replace(string(formatted_value), r"^-" => "−")
    med_str = "median: $(formatted_value)"
    
    # Place legend in row 1 (the dedicated legend row)
    Legend(gl[1, 1], 
           [box_plot], 
           [med_str],
           tellwidth=config.tellwidth,
           tellheight=config.tellheight,
           halign=config.halign,
           valign=config.valign,
           margin=config.margin,
           backgroundcolor=config.backgroundcolor,
           framecolor=config.framecolor,
           framewidth=config.framewidth,
           labelsize=config.labelsize,
           patchsize=config.patchsize,
           padding=config.padding)
end

"""
    get_adaptive_markersize(n_points, base_size)

Compute adaptive marker size based on number of points.
Increases size when fewer than 100 points for better visibility.
Scaling is conservative to avoid dots touching plot borders.
"""
function get_adaptive_markersize(n_points::Int, base_size::Real)
    if n_points < 100
        # Conservative scaling: 100 points → 1.0x, 50 points → 1.3x, 10 points → 1.6x
        scale_factor = 1.0 + (100 - n_points) / 90.0  # Even more conservative
        return base_size * scale_factor
    else
        return base_size
    end
end

"""
    draw_scatter_with_jitter!(ax, values, config::ScatterConfig)

Draw a scatter plot with vertical jitter for better visualization of point density.
Jitter is clamped to ±0.12 to ensure dots never touch the plot borders (y-limits are ±0.25).
"""
function draw_scatter_with_jitter!(ax, values, config::ScatterConfig)
    # Apply jitter using the configured scale
    jittered_y = config.jitter_scale .* randn(length(values))
    
    adaptive_size = get_adaptive_markersize(length(values), config.markersize)
    
    # Clamp jitter more aggressively to account for large adaptive marker sizes
    # Y-axis limits are ±0.25, clamp to ±0.12 to leave 0.13 buffer for marker radius
    max_jitter = 0.12
    jittered_y = clamp.(jittered_y, -max_jitter, max_jitter)
    
    # Debug: print the actual markersize being used (using println to ensure it prints)
    # println("DEBUG: Scatter plot with $(length(values)) points, base size=$(config.markersize), adaptive size=$(adaptive_size), jitter_scale=$(config.jitter_scale), max_jitter=$(maximum(abs.(jittered_y)))")
    
    # Only add stroke outline when there are few points (< 100) for better visibility
    if length(values) < 100
        scatter!(ax, values, jittered_y;
                 color=config.color, 
                 alpha=config.alpha, 
                 markersize=adaptive_size, 
                 marker=config.marker,
                 strokecolor=:black,
                 strokewidth=0.8)
    else
        scatter!(ax, values, jittered_y;
                 color=config.color, 
                 alpha=config.alpha, 
                 markersize=adaptive_size, 
                 marker=config.marker)
    end
    
    # CRITICAL: Set y-limits explicitly to prevent auto-scaling from defeating our jitter clamp
    # Must match the buffer we designed for (±0.25 with jitter clamped to ±0.12)
    ylims!(ax, -0.25, 0.25)
    
    hidexdecorations!(ax)
    hideydecorations!(ax)
    hidespines!(ax, :t, :r, :l, :b)
end

"""
    setup_xlabel_axis!(ax, config::LayoutConfig)

Configure the spacer axis that displays the x-axis label.
"""
function setup_xlabel_axis!(ax, config::LayoutConfig)
    hideydecorations!(ax)
    hidespines!(ax, :t, :r, :l)
    
    # Only show x-axis label if enabled
    if config.show_xlabel
        ax.xlabel = config.xlabel
        ax.xlabelsize = config.xlabelsize
    else
        ax.xlabel = ""
    end
end

"""
    create_boxplot_scatter(gl::GridLayout, banzhaf_values, format_median;
                          box_config=BoxPlotConfig(),
                          scatter_config=ScatterConfig(),
                          legend_config=LegendConfig(),
                          layout_config=LayoutConfig())

Create a box plot with scatter visualization in the provided GridLayout.

# Arguments
- `gl::GridLayout`: The grid layout to populate with the visualization
- `banzhaf_values`: Vector of values to visualize
- `format_median`: Function to format the median value for display

# Keyword Arguments
- `box_config::BoxPlotConfig`: Configuration for box plot styling
  - `box_height`: Height of the box (IQR)
  - `whisker_cap_height`: Height of the whisker caps (nothing = auto, 0 = no caps)
  - `box_offset_y`: Vertical offset of the box
  - `box_color`: Color of the box
  - `line_width`: Width of box and whisker lines
- `scatter_config::ScatterConfig`: Configuration for scatter plot styling
- `legend_config::LegendConfig`: Configuration for legend styling
  - `labelsize`: Font size of the legend text
  - `halign`: Horizontal alignment (:left, :center, :right)
  - `valign`: Vertical alignment (:top, :center, :bottom)
  - `margin`: (top, right, bottom, left) margins around legend
  - `padding`: Internal padding within legend box
- `layout_config::LayoutConfig`: Configuration for layout and labels
  - `row_legend`, `row_box`, `row_scatter`, `row_spacer`: Relative row heights
  - `row_gap`: Gap between rows (uses rowgap! if set)

# Returns
- Named tuple containing the axes and statistics

# Examples
```julia
# Adjust box height and remove whisker caps
result = create_boxplot_scatter(gl, data, fmt;
    box_config=BoxPlotConfig(box_height=0.2, whisker_cap_height=0.0))

# Larger legend font positioned at left
result = create_boxplot_scatter(gl, data, fmt;
    legend_config=LegendConfig(labelsize=18, halign=:left, margin=(5, 10, 15, 10)))

# Adjust spacing between boxplot and scatter
result = create_boxplot_scatter(gl, data, fmt;
    layout_config=LayoutConfig(row_gap=10.0))
```
"""
function create_boxplot_scatter(gl::GridLayout, banzhaf_values, format_median;
                                box_config=BoxPlotConfig(),
                                scatter_config=ScatterConfig(),
                                legend_config=LegendConfig(),
                                layout_config=LayoutConfig())
    # Create axes in the grid layout
    ax_box = Axis(gl[2, 1])
    ax_scatter = Axis(gl[3, 1])
    ax_spacer = Axis(gl[4, 1])
    
    # Draw box plot with configurable whisker cap height
    stats, box_plot = draw_boxplot!(ax_box, banzhaf_values; 
                                    box_height=box_config.box_height,
                                    box_offset_y=box_config.box_offset_y,
                                    box_color=box_config.box_color,
                                    line_width=box_config.line_width,
                                    whisker_cap_height=box_config.whisker_cap_height)
    
    # Calculate data range for axis limits
    # Use a minimum range of 1.0 to handle zero-variance data
    data_range = max(stats.whisker_hi - stats.whisker_lo, 1.0)
    
    # Setup box plot axis
    setup_box_axis!(ax_box, stats, box_config.box_offset_y, data_range)
    
    # Create legend with median value
    create_median_legend!(gl, box_plot, stats.median, format_median, legend_config)
    
    # Draw scatter plot with jitter
    draw_scatter_with_jitter!(ax_scatter, banzhaf_values, scatter_config)
    
    # Setup x-axis label
    setup_xlabel_axis!(ax_spacer, layout_config)
    
    # Link x-axes for synchronized zooming/panning
    linkxaxes!(ax_box, ax_scatter, ax_spacer)
    
    # Adjust row sizes
    rowsize!(gl, 1, Relative(layout_config.row_legend))
    rowsize!(gl, 2, Relative(layout_config.row_box))
    rowsize!(gl, 3, Relative(layout_config.row_scatter))
    rowsize!(gl, 4, Relative(layout_config.row_spacer))
    
    # Set row gap if specified
    if !isnothing(layout_config.row_gap)
        rowgap!(gl, layout_config.row_gap)
    end
    
    return (axes=(box=ax_box, scatter=ax_scatter, spacer=ax_spacer), 
            stats=stats, 
            box_plot=box_plot)
end