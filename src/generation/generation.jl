"""
Main entry point for the generation module.

This file includes all generation submodules in the correct order for both
convolution-based and mutagenesis-based motif analysis.

# Usage

```julia
include("generation/generation.jl")

# Now all generation functions are available:
# - Core: types, sorting, file I/O, rendering
# - Convolution: singletons, multi-motifs
# - Mutagenesis: single regions, multi-regions
# - Utils: distance keys, consensus
```

# Module Organization

```
generation/
├── core/              # Shared functionality
│   ├── types.jl              # Config structs & data types
│   ├── data_structures.jl    # JSON/HTML dict management
│   ├── grouping.jl           # Grouping column builders
│   ├── sorting.jl            # Pareto ranking & sorting
│   ├── file_io.jl            # File I/O operations
│   └── rendering.jl          # HTML/JS/CSS rendering
├── convolution/       # Convolution-based analysis
│   ├── singletons.jl         # Singleton motif processing
│   └── multi_motifs.jl       # Pairs, triplets processing
├── mutagenesis/       # Mutation-based analysis
│   ├── singletons.jl         # Single mutation regions
│   ├── multi_regions.jl      # Paired/triple mutations
│   ├── matrix_operations.jl  # Matrix merging & accumulation
│   └── reference_matching.jl # Off-region search utilities
└── utils/             # General utilities
    └── distance_keys.jl      # Distance key extraction
```
"""

using Statistics: median, mean
using EntroPlots
using CSV
using JSON3

# Get the directory where this file is located
const GENERATION_DIR = @__DIR__

# =============================================================================
# Core Module (required by all other modules)
# =============================================================================

include(joinpath(GENERATION_DIR, "core", "types.jl"))
include(joinpath(GENERATION_DIR, "core", "data_structures.jl"))
include(joinpath(GENERATION_DIR, "core", "grouping.jl"))
include(joinpath(GENERATION_DIR, "core", "sorting.jl"))
include(joinpath(GENERATION_DIR, "core", "file_io.jl"))
include(joinpath(GENERATION_DIR, "core", "rendering.jl"))

# =============================================================================
# Utilities Module
# =============================================================================

include(joinpath(GENERATION_DIR, "utils", "symbols.jl"))
include(joinpath(GENERATION_DIR, "utils", "distance_keys.jl"))

# =============================================================================
# Convolution Module
# =============================================================================

include(joinpath(GENERATION_DIR, "convolution", "singletons.jl"))
include(joinpath(GENERATION_DIR, "convolution", "multi_motifs.jl"))

# =============================================================================
# Mutagenesis Module
# =============================================================================

include(joinpath(GENERATION_DIR, "mutagenesis", "matrix_operations.jl"))
include(joinpath(GENERATION_DIR, "mutagenesis", "reference_matching.jl"))
include(joinpath(GENERATION_DIR, "mutagenesis", "singletons.jl"))
include(joinpath(GENERATION_DIR, "mutagenesis", "multi_regions.jl"))

println("✓ Generation module loaded successfully")
