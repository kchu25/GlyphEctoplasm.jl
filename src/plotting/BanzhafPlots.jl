module BanzhafPlots

using CairoMakie
using StatsBase
using FileIO

# Include all component files
include("helper.jl")
include("logo_plotting_structs.jl")
include("box_scatter_structs.jl")
include("logo_plotting.jl")
include("box_scatter.jl")
include("utils.jl")
include("wrapper.jl")

# Export helper functions
export format_median

# Export structs
export BoxPlotConfig, ScatterConfig, LegendConfig, LayoutConfig
export LogoConfig
export DisplayMode, DefaultMode, WebDisplayMode, PublicationMode
export get_mode_configs

# Export main plotting functions
export create_boxplot_scatter
export plot_logo!

# Export high-level convenience functions
export create_simple_boxplot_figure
export create_logo_boxplot_figure
export save_figure

# Export wrapper functions
export box_scatter_template
export save_box_scatter
export save_box_scatter_distance_default
export save_box_scatter_distance_relaxed
export save_box_scatter_distance_fixed

end # module
