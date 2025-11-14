# ============================================================================
# TEMPLATES - Refactored Modular Version
# ============================================================================

"""
This file provides a unified interface to all template components.
The actual templates are now organized in the templates/ subdirectory.

For the modular structure, see:
- templates/css_template.jl
- templates/html_templates.jl
- templates/js_core.jl
- templates/js_sequence.jl
- templates/js_modals.jl
- templates/js_styling.jl

This file maintains backward compatibility by re-exporting all templates.
"""

# Include the modular templates directly (not as a module)
include("templates/css_template.jl")
include("templates/html_templates.jl")
include("templates/js_core.jl")
include("templates/js_sequence.jl")
include("templates/js_modals.jl")
include("templates/js_styling.jl")

# Combine all JavaScript components into one template
# Build the combined script by concatenating the raw strings
const _combined_script = script_core_str * "\n\n" * 
                         script_sequence_str * "\n\n" * 
                         script_modals_str * "\n\n" * 
                         script_styling_str

# Convert to Mustache template
const script_template = Mustache.parse(_combined_script)

