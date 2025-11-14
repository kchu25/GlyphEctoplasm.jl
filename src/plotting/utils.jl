
"""
    create_simple_boxplot_figure(banzhaf_values; 
                                  mode=DefaultMode,
                                  ylabel="",
                                  ylabelsize=14,
                                  xlabel="Influence (Banzhaf indices)",
                                  xlabelsize=17,
                                  figsize=(725, 115),
                                  figure_padding=(5, 5, 5, 2),
                                  xlim=nothing,
                                  box_config=BoxPlotConfig(),
                                  legend_config=LegendConfig(),
                                  layout_config=LayoutConfig())

Create a simple boxplot-scatter figure with customizable options or use preset display modes.

# Arguments
- `banzhaf_values`: Vector of values to visualize

# Keyword Arguments
- `mode::DisplayMode`: Preset display mode (DefaultMode, WebDisplayMode, or PublicationMode). 
                       When specified, overrides legend_config, layout_config, ylabel display, 
                       and xlim if provided to mode.
- `ylabel::String`: Y-axis label (default: "")
- `ylabelsize::Int`: Y-axis label font size (default: 14)
- `xlabel::String`: X-axis label (default: "Influence (Banzhaf indices)") - ignored if using mode
- `xlabelsize::Int`: X-axis label font size (default: 17) - ignored if using mode
- `figsize::Tuple{Int,Int}`: Figure size in pixels (default: (725, 115))
- `figure_padding::Tuple{Int,Int,Int,Int}`: Figure padding (left, right, bottom, top) (default: (5, 5, 5, 2))
- `xlim::Union{Nothing,Tuple{Real,Real}}`: X-axis limits (e.g., (-0.1, 1.1)). Default: nothing (auto)
- `box_config::BoxPlotConfig`: Box plot configuration (used in all modes)
- `legend_config::LegendConfig`: Legend configuration (overridden by mode)
- `layout_config::LayoutConfig`: Layout configuration (overridden by mode)

# Returns
- `(figure, result)`: Tuple containing the Figure object and the result from create_boxplot_scatter

# Example
```julia
# Use preset web display mode (no legend, no labels)
f, result = create_simple_boxplot_figure(banzhaf_values; mode=WebDisplayMode)

# Use default mode with custom ylabel
f, result = create_simple_boxplot_figure(banzhaf_values; 
    mode=DefaultMode,
    ylabel="Filter 1",
    ylabelsize=16)

# Manual configuration (original behavior)
f, result = create_simple_boxplot_figure(banzhaf_values; 
    ylabel="Filter 1",
    ylabelsize=16,
    xlabel="Custom X Label",
    xlabelsize=18,
    xlim=(-0.1, 1.0),
    box_config=BoxPlotConfig(box_color=:lightgray))
```
"""
function create_simple_boxplot_figure(banzhaf_values; 
                                      mode::Union{Nothing,DisplayMode}=nothing,
                                      ylabel="",
                                      ylabelsize=14,
                                      xlabel="Influence (Banzhaf indices)",
                                      xlabelsize=17,
                                      figsize=(725, 115),
                                      figure_padding=(5, 5, 1, 2), # left, right, bottom, top
                                      xlim::Union{Nothing,Tuple{<:Real,<:Real}}=nothing,
                                      box_config=BoxPlotConfig(
                                          box_height=0.1,
                                          box_color=:lightgray,
                                          whisker_cap_height=0.1
                                      ),
                                      scatter_config=ScatterConfig(),  # Add scatter_config parameter
                                      legend_config=LegendConfig(
                                          labelsize=13,
                                          halign=:center,
                                          valign=:top,
                                          margin=(10, 10, 10, 0),
                                          tellheight=false
                                      ),
                                      layout_config=LayoutConfig(
                                          row_legend=0.22,
                                          row_gap=2.0,
                                          row_box=0.48,
                                          row_scatter=0.30
                                      ))
    
    # If mode is specified, use preset configurations
    if !isnothing(mode)
        mode_configs = get_mode_configs(mode; xlim=xlim)
        legend_config = mode_configs.legend_config
        layout_config = mode_configs.layout_config
        show_ylabel = mode_configs.show_ylabel
        xlim = mode_configs.xlim
    else
        # Manual mode: use provided configs and show ylabel if it's not empty
        show_ylabel = !isempty(ylabel)
    end
    
    # Create a custom layout config with the provided xlabel and xlabelsize
    # (unless mode already set these via layout_config)
    if isnothing(mode)
        custom_layout_config = LayoutConfig(
            row_legend=layout_config.row_legend,
            row_box=layout_config.row_box,
            row_scatter=layout_config.row_scatter,
            row_spacer=layout_config.row_spacer,
            row_gap=layout_config.row_gap,
            xlabel=xlabel,
            xlabelsize=xlabelsize,
            show_xlabel=layout_config.show_xlabel
        )
    else
        custom_layout_config = layout_config
    end
    
    # Create figure and grid layout
    f = Figure(size=figsize, figure_padding=figure_padding)
    gl = f[1, 1] = GridLayout()
    
    # Create the visualization
    result = create_boxplot_scatter(gl, banzhaf_values, format_median;
        box_config=box_config,
        scatter_config=scatter_config,  # Pass scatter_config through
        legend_config=legend_config,
        layout_config=custom_layout_config)
    
    # Apply xlim if specified (to both box and scatter axes)
    if !isnothing(xlim)
        xlims!(result.axes.box, xlim...)
        xlims!(result.axes.scatter, xlim...)
        xlims!(result.axes.spacer, xlim...)
    end
    
    # Add ylabel if provided and enabled by mode
    if show_ylabel && !isempty(ylabel)
        Label(f[1, 1, Left()], ylabel, 
              rotation=π/2, 
              tellheight=false,
              fontsize=ylabelsize,
              valign=:top,  # In Left() protrusion, :top moves label down toward top of figure
              padding=(0, 15, 0, 0))  # Add 15px spacing between ylabel and plot
    end
    
    return f, result
end

"""
    create_logo_boxplot_figure(logo_path, banzhaf_values;
                                figsize=(1020, 230),
                                figure_padding=15,
                                logo_width=343,
                                colgap=0,
                                logo_config=LogoConfig(),
                                box_config=BoxPlotConfig(),
                                layout_config=LayoutConfig())

Create a figure with a logo on the left and boxplot-scatter on the right.

# Arguments
- `logo_path::String`: Path to the logo image file (PNG)
- `banzhaf_values`: Vector of values to visualize in the boxplot-scatter

# Keyword Arguments
- `figsize::Tuple{Int,Int}`: Figure size in pixels (default: (1020, 230))
- `figure_padding::Int`: Figure padding on all sides (default: 15)
- `logo_width::Int`: Fixed width for the logo column in pixels (default: 343)
- `colgap::Int`: Gap between logo and boxplot columns in pixels (default: 0)
- `logo_config::LogoConfig`: Configuration for logo display
- `box_config::BoxPlotConfig`: Configuration for boxplot styling
- `layout_config::LayoutConfig`: Configuration for boxplot-scatter layout

# Returns
- `(figure, result)`: Tuple containing the Figure object and the result from create_boxplot_scatter

# Example
```julia
f, result = create_logo_boxplot_figure("tmp/singletons/3.png", banzhaf_values;
    logo_width=300,
    box_config=BoxPlotConfig(box_height=0.065, whisker_cap_height=0.1))
```
"""
function create_logo_boxplot_figure(logo_path, banzhaf_values;
                                    figsize=(1020, 230),
                                    figure_padding=15,
                                    logo_width=343,
                                    colgap=0,
                                    logo_config=LogoConfig(
                                        rotate=true,
                                        preserve_aspect=true,
                                        row_span=1:5,
                                        x_autolimit_margin=(0.0, 0.0),
                                        y_autolimit_margin=(0.0, 0.0)
                                    ),
                                    box_config=BoxPlotConfig(
                                        box_height=0.065,
                                        whisker_cap_height=0.1
                                    ),
                                    scatter_config=ScatterConfig(),  # Add scatter_config parameter
                                    layout_config=LayoutConfig(
                                        row_gap=2.0,
                                        row_box=0.50,
                                        row_scatter=0.125
                                    ))
    
    # Create figure
    f = Figure(size=figsize, figure_padding=figure_padding)
    
    # Create two grid layouts: one for logo (left), one for boxplot-scatter (right)
    logo_gl = f[1, 1] = GridLayout()
    plot_gl = f[1, 2] = GridLayout()
    
    # Set alignmode for logo grid to flush it to the left
    logo_gl.alignmode = Mixed(left=0, right=65)
    
    # Plot the logo on the left - spans all rows with preserved aspect
    plot_logo!(logo_gl, logo_path; config=logo_config)
    
    # Set row sizes AFTER the content exists
    rowsize!(logo_gl, 1, Relative(0.25))   # Smaller top
    rowsize!(logo_gl, 2, Relative(0.85))   # Much larger middle (logo area)
    rowsize!(logo_gl, 3, Relative(0.6))    # Medium bottom
    rowsize!(logo_gl, 4, Relative(0.001))  # Tiny spacer
    
    # Plot the boxplot-scatter on the right
    result = create_boxplot_scatter(plot_gl, banzhaf_values, format_median;
        box_config=box_config,
        scatter_config=scatter_config,  # Pass scatter_config through
        layout_config=layout_config)
    
    # Adjust column widths - logo has fixed width, plot gets remaining space
    colsize!(f.layout, 1, Fixed(logo_width))
    
    # Set the gap between the two columns
    colgap!(f.layout, colgap)
    
    return f, result
end


"""
    save_figure(fig, path; overwrite=true)

Save a `Figure` to disk. Creates parent directories if necessary. If the file
already exists and `overwrite=false` an error is thrown.

# Keyword Arguments
- `overwrite::Bool` : whether to overwrite an existing file (default: true)

# Returns
- the `path` string on success
"""
function save_figure(fig, path::AbstractString; overwrite::Bool=true)
    # Ensure parent directory exists
    parent = dirname(path)
    if !isempty(parent) && !isdir(parent)
        mkpath(parent)
    end

    if isfile(path) && !overwrite
        error("File exists and overwrite=false: $path")
    end

    # Try the standard Makie/CairoMakie save call. `save(path, fig)` is the
    # usual form when `using CairoMakie` is active in the calling script.
    try
        save(path, fig)
    catch err
        # Some environments accept save(fig, path) instead — try that as fallback.
        try
            save(fig, path)
        catch err2
            rethrow(err2)
        end
    end

    return path
end
