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
include(joinpath("core", "json_html_dict.jl"))
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
include("run_thru_conv.jl")
include("run_thru_mut.jl")

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
export render_and_save_outputs!

# Template constants
export html_template_unified, script_template, template_css

# Rendering
export plot_motifs_conv_case, plot_motifs_mut_case

println("✓ Render module loaded successfully")

end
