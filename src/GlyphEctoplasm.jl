"""
# GlyphEctoplasm Module

A comprehensive module for motif analysis visualization and rendering.

## Usage

```julia
using GlyphEctoplasm

# All functionality is now available:
# - BanzhafPlots submodule for plotting
# - Generation functions for processing motifs
# - Template constants for HTML/CSS/JS
```

## Exported Functions

### Data Structures
- `init_json_dict()` - Initialize JSON dictionary
- `init_dict_for_html_render()` - Initialize HTML rendering dictionary

### Convolution Processing
- `process_singletons!()` - Process singleton motifs
- `process_multi_motifs!()` - Process paired/multi motifs

### Rendering
- `render_and_save_outputs!()` - Render and save HTML outputs

### Configuration
- `ConvMotifConfig` - Configuration for convolution-based motif analysis

### Templates
- `html_template_unified` - Unified HTML template
- `script_template` - JavaScript template
- `template_css` - CSS template
"""

# TODO: make the template files constants

module GlyphEctoplasm

# External dependencies
using Statistics: median, mean
using SEQ2EXPdata
using ProgressMeter
using DataFrames
using Mustache
using EntroPlots
using CSV
using JSON3
using CairoMakie
using StatsBase
using FileIO
using Printf
using Colors
using FixedPointNumbers: N0f8
using PNGFiles
using IndirectArrays
using KernelDensity: kde
using Random
using MultipleTesting: adjust, PValues, BenjaminiHochberg

# =============================================================================
# BanzhafPlots Submodule
# =============================================================================
include(joinpath("plotting", "BanzhafPlots.jl"))
using .BanzhafPlots

# Export BanzhafPlots submodule
export BanzhafPlots

# =============================================================================
# Core Rendering Utilities
# =============================================================================
include(joinpath("core", "logo_saving.jl"))
include(joinpath("core", "points_dump.jl"))
include(joinpath("core", "json_html_dict.jl"))
include(joinpath("core", "png_optimize.jl"))
include(joinpath("core", "constants.jl"))
include(joinpath("core", "consensus.jl"))
include(joinpath("core", "html_generation.jl"))
include(joinpath("core", "templates.jl"))
include(joinpath("core", "path_utils.jl"))

# =============================================================================
# Generation System
# =============================================================================
include(joinpath("generation", "generation.jl"))

# =============================================================================
# Rendering 
# =============================================================================
include("generalization.jl")
include("run_nnd_sensitivity.jl")
include("run_thru_conv.jl")
include("run_thru_mut.jl")

# =============================================================================
# Constants
# =============================================================================
const SENSITIVITY_KS = [1, 3, 5, 7, 10, 15, 20, 25, 30, 35, 40]

# =============================================================================
# Exports
# =============================================================================

# Data structure initialization
export init_json_dict, init_dict_for_html_render

# Configuration types
export ConvMotifConfig, MutMotifConfig

# Processing functions
export process_singletons!, process_multi_motifs!
export process_single_mut_region!, process_multi_mut_regions!

# Rendering functions
export render_and_save_outputs!, render_generalization_page!,
       render_statistics_page!, render_readme_page!

# PNG size optimization
export optimize_pngs!, quantize_png!

# Indicator-plot point cloud (rebuilds the yy-KDE figures from files alone)
export save_indicator_points

# Sensitivity analysis
export run_nnd_sensitivity_analysis, run_nnd_sensitivity_analysis_null, SENSITIVITY_KS

# Template constants
export html_template_unified, html_template_generalization,
       html_template_statistics, html_template_readme,
       script_template, template_css

# Rendering
export plot_motifs_conv_case, plot_motifs_mut_case

println("✓ Render module loaded successfully")

end
