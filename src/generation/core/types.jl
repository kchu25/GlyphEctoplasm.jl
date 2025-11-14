"""
Configuration types for motif analysis and rendering.

This module defines configuration structs and data types used across the motif
generation pipeline, separating concerns between data processing, rendering,
and analysis configuration.
"""

# =============================================================================
# Configuration Structs
# =============================================================================

"""
    MutationRegionConfig

Configuration for mutation region analysis and rendering.
Groups all parameters needed for mutation region metadata collection and visualization.

# Fields
- `data::OnehotSEQ2EXP_Dataset`: Dataset containing sequences and expression data
- `reference_seq::BitMatrix`: Reference sequence as a bit matrix
- `total_length::Int`: Total length of sequences
- `filter_len::Int`: Length of convolutional filters
- `float_type::Type`: Floating point precision (Float32 or Float64)
- `off_region_search::Bool`: Whether to search for off-region mutations
- `reduction_on_ref::Bool`: Whether to filter by reference sequence
- `dpi::Int`: DPI for saved images
- `use_rna::Bool`: Whether to use RNA (U) instead of DNA (T)
- `xlim::Union{Nothing, Tuple{Float64, Float64}}`: X-axis limits for plots
- `save_path::String`: Base path for saving outputs
"""
struct MutationRegionConfig
    # Data parameters
    data::OnehotSEQ2EXP_Dataset
    reference_seq::BitMatrix
    total_length::Int
    filter_len::Int
    float_type::Type
    
    # Analysis parameters
    off_region_search::Bool
    reduction_on_ref::Bool
    
    # Rendering parameters
    dpi::Int
    use_rna::Bool
    xlim::Union{Nothing, Tuple{Float64, Float64}}
    
    # Path configuration
    save_path::String
end

"""
    MutationRegionConfig(data; kwargs...)

Convenience constructor with sensible defaults.

# Arguments
- `data::OnehotSEQ2EXP_Dataset`: Required dataset

# Keyword Arguments
- `reference_seq`: Reference sequence (defaults to consensus from data)
- `total_length`: Total sequence length (defaults to consensus length)
- `filter_len::Int = 9`: Filter length
- `float_type::Type = Float32`: Floating point precision
- `off_region_search::Bool = true`: Enable off-region search
- `reduction_on_ref::Bool = true`: Enable reference filtering
- `dpi::Int = 65`: Image DPI
- `use_rna::Bool = false`: Use RNA notation
- `xlim = nothing`: Plot x-axis limits
- `save_path::String = "tmp"`: Output directory
"""
function MutationRegionConfig(data;
    reference_seq = SEQ2EXPdata.consensus_to_bitmatrix_auto(data.raw_data.consensus),
    total_length = length(data.raw_data.consensus),
    filter_len::Int = 9,
    float_type::Type = Float32,
    off_region_search::Bool = true,
    reduction_on_ref::Bool = true,
    dpi::Int = 65,
    use_rna::Bool = false,
    xlim = nothing,
    save_path::String = "tmp"
)
    MutationRegionConfig(
        data, reference_seq, total_length, filter_len, float_type,
        off_region_search, reduction_on_ref, dpi, use_rna, xlim, save_path
    )
end

"""
    ConvMotifConfig

Configuration for convolution-based motif analysis (singletons, pairs, triplets).

# Fields
- `data::OnehotSEQ2EXP_Dataset`: Dataset containing sequences and expression data
- `filter_len::Int`: Length of convolutional filters
- `float_type::Type`: Floating point precision (Float32 or Float64)
- `dpi::Int`: DPI for saved images
- `alpha::Float64`: Alpha transparency for visualizations
- `use_rna::Bool`: Whether to use RNA (U) instead of DNA (T)
- `xlim::Union{Nothing, Tuple{Float64, Float64}}`: X-axis limits for plots
- `save_path::String`: Base path for saving outputs
"""
struct ConvMotifConfig
    # Data parameters
    data::OnehotSEQ2EXP_Dataset
    filter_len::Int
    float_type::Type
    
    # Rendering parameters
    dpi::Int
    alpha::Float64
    use_rna::Bool
    xlim::Union{Nothing, Tuple{Float64, Float64}}
    
    # Path configuration
    save_path::String
end

"""
    ConvMotifConfig(data; kwargs...)

Convenience constructor with sensible defaults.

# Arguments
- `data::OnehotSEQ2EXP_Dataset`: Required dataset

# Keyword Arguments
- `filter_len::Int = 7`: Filter length
- `float_type::Type = Float32`: Floating point precision
- `dpi::Int = 65`: Image DPI
- `alpha::Float64 = 1.0`: Alpha transparency
- `use_rna::Bool = false`: Use RNA notation
- `xlim = nothing`: Plot x-axis limits
- `save_path::String = "tmp"`: Output directory
"""
function ConvMotifConfig(data;
    filter_len::Int = 7,
    float_type::Type = Float32,
    dpi::Int = 65,
    alpha::Float64 = 1.0,
    use_rna::Bool = false,
    xlim = nothing,
    save_path::String = "tmp"
)
    ConvMotifConfig(data, filter_len, float_type, dpi, alpha, use_rna, xlim, save_path)
end

# =============================================================================
# Data Structs
# =============================================================================

"""
    MotifData

Core motif data extracted from mutation or convolution analysis.
Contains the essential numerical and sequence data for a motif.

# Fields
- `key`: Grouping key (filter indices, mutation positions, etc.)
- `count_matrices::Vector`: Count matrices for the motif
- `positions::Vector`: Positions of motif occurrences
- `references::Vector{BitMatrix}`: Reference sequence matrices
- `median::Float64`: Median contribution value
- `count::Int`: Number of occurrences
- `banzhafs::Vector`: Banzhaf contribution values
- `gdf_row`: Corresponding row from grouped dataframe
"""
struct MotifData
    key::Any
    count_matrices::Vector
    positions::Vector
    references::Vector{BitMatrix}
    median::Float64
    count::Int
    banzhafs::Vector
    gdf_row::Any
end

"""
    MotifMetadata

Complete metadata for a motif, combining data, configuration, and computed grouping info.
Used for registration and rendering in the HTML output.

# Fields
- `data::MotifData`: Core motif data
- `config::MutationRegionConfig`: Configuration used for analysis
- `motif_type::String`: Type identifier (e.g., "pair_mutation_regions")
- `save_folder::String`: Folder where outputs are saved
- `motif_size::Int`: Size of the motif (1 for singletons, 2 for pairs, etc.)
- `group_id::String`: Group identifier (e.g., "single_region", "2_regions")
- `button_text::String`: Display text for UI buttons
- `span::String`: Span string for position ranges
"""
struct MotifMetadata
    data::MotifData
    config::MutationRegionConfig
    
    # Motif-specific settings
    motif_type::String
    save_folder::String
    motif_size::Int
    
    # Computed grouping info
    group_id::String
    button_text::String
    span::String
end

"""
    MotifMetadata(data, config, motif_type, save_folder, motif_size, fragment_info)

Convenience constructor using fragment_info tuple.
"""
function MotifMetadata(data::MotifData, config::MutationRegionConfig, 
                       motif_type::String, save_folder::String, motif_size::Int, 
                       fragment_info::NamedTuple)
    MotifMetadata(
        data, config, motif_type, save_folder, motif_size,
        fragment_info.group_id, fragment_info.button_text, fragment_info.span
    )
end

# =============================================================================
# Utility Functions for Accessing Struct Fields
# =============================================================================

"""
    Base.getproperty(meta::MotifMetadata, sym::Symbol)

Custom property accessor for MotifMetadata that provides convenient access
to commonly used fields from nested structs (for backward compatibility).
"""
function Base.getproperty(meta::MotifMetadata, sym::Symbol)
    # Handle direct fields
    if sym in (:data, :config, :motif_type, :save_folder, :motif_size, :group_id, :button_text, :span)
        return getfield(meta, sym)
    end
    
    # Convenience accessors for config fields
    if sym in (:dpi, :use_rna, :xlim, :filter_len, :total_length, :reduction_on_ref, :off_region_search, :float_type, :reference_seq)
        return getfield(getfield(meta, :config), sym)
    end
    
    # Convenience accessors for data fields
    if sym in (:key, :count_matrices, :positions, :references, :median, :count, :banzhafs, :gdf_row)
        return getfield(getfield(meta, :data), sym)
    end
    
    error("type MotifMetadata has no field $sym")
end

"""
    Base.propertynames(meta::MotifMetadata)

List all accessible properties including nested ones.
"""
function Base.propertynames(meta::MotifMetadata)
    return (:data, :config, :motif_type, :save_folder, :motif_size, :group_id, :button_text, :span,
            :dpi, :use_rna, :xlim, :filter_len, :total_length, :reduction_on_ref, :off_region_search,
            :float_type, :reference_seq, :key, :count_matrices, :positions, :references, 
            :median, :count, :banzhafs, :gdf_row)
end

# Export types and constructors
export MutationRegionConfig, ConvMotifConfig, MotifData, MotifMetadata
