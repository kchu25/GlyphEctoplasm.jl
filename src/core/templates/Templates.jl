# ============================================================================
# TEMPLATES MODULE
# ============================================================================

"""
Main module for motif visualization templates.

This module provides modular templates for generating interactive HTML pages
for visualizing biological sequence motifs with various features including:
- Responsive grid layouts
- Interactive modals
- Sequence highlighting
- Dynamic styling based on Banzhaf values

# Structure
The templates are organized into separate files:
- `css_template.jl`: Stylesheet for all components
- `html_templates.jl`: HTML page structures
- `js_core.jl`: Core JavaScript (navigation, data loading, sliders)
- `js_sequence.jl`: Sequence highlighting functionality
- `js_modals.jl`: Modal window management
- `js_styling.jl`: Dynamic card styling

# Usage
```julia
include("templates/Templates.jl")

# Access individual templates
css = template_css
html = html_template
js = build_script_template(j=1, upto=6, mode_counts=10, sequence_paths=["path/to/seq.fasta"])
```

# Exported Templates
- `template_css`: Complete CSS stylesheet
- `html_template`: HTML template for multi-motif pages (pairs, triplets, etc.)
- `html_template_singleton`: HTML template for singleton motif pages
- `html_hover_default`: Hover window template for metadata
- `html_end`: Closing HTML tags
- `build_script_template()`: Function to build complete JavaScript from components
"""

# Include all template components
include("css_template.jl")
include("html_templates.jl")
include("js_core.jl")
include("js_sequence.jl")
include("js_modals.jl")
include("js_styling.jl")

# Export CSS templates
export template_css

# Export HTML templates
export html_template
export html_template_singleton
export html_hover_default
export html_end

# Export JavaScript component builders
export script_core
export script_sequence
export script_modals
export script_styling

"""
    build_script_template(; j, upto, mode_counts, sequence_file_paths)

Build complete JavaScript template by combining all components.

# Arguments
- `j::Int`: Current page number
- `upto::Int`: Total number of pages
- `mode_counts::Int`: Number of motif modes to display
- `sequence_file_paths::Vector{String}`: Array of FASTA file paths

# Returns
- `String`: Complete JavaScript code ready for Mustache rendering

# Example
```julia
script = build_script_template(
    j=2,
    upto=6,
    mode_counts=15,
    sequence_file_paths=["data/sequences.fasta", "data/sequences2.fasta"]
)
```
"""
function build_script_template(; j::Int, upto::Int, mode_counts::Int, sequence_file_paths::Vector{String})
    # Note: The Mustache variables are embedded in each component template
    # They will be filled when the template is rendered, not here
    # This function just combines the component templates
    
    # Build complete script by concatenating components
    # The order matters for functionality
    return script_core * "\n\n" * 
           script_sequence * "\n\n" * 
           script_modals * "\n\n" * 
           script_styling
end

export build_script_template

# For backward compatibility, provide the original combined template
"""
    script_template

Complete JavaScript template (backward compatibility).

This is the original monolithic template. For new code, prefer using
`build_script_template()` which combines the modular components.
"""
const script_template = mt"""
""" * 
script_core * 
script_sequence * 
script_modals * 
script_styling * 
"""
"""

export script_template
