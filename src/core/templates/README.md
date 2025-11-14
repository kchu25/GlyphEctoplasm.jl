# Templates Module - Refactored Structure

This directory contains the refactored, modular template system for motif visualization pages.

## File Structure

```
templates/
├── README.md              # This file
├── Templates.jl           # Main module file (optional, for module-based usage)
├── css_template.jl        # Complete CSS stylesheet
├── html_templates.jl      # HTML page templates (multi-motif and singleton)
├── js_core.jl            # Core JavaScript (navigation, data loading, sliders)
├── js_sequence.jl        # Sequence highlighting functionality
├── js_modals.jl          # Modal window management
└── js_styling.jl         # Dynamic card styling based on Banzhaf values
```

## Purpose

The original `templates.jl` file was a large monolithic file (~2000 lines) containing:
- CSS styles
- HTML templates
- JavaScript code for multiple features

This has been refactored into smaller, focused modules for better:
- **Maintainability**: Each file has a clear, single responsibility
- **Readability**: Easier to find and understand specific functionality
- **Extensibility**: Can add new features without modifying unrelated code
- **Testing**: Individual components can be tested in isolation

## Components

### CSS Template (`css_template.jl`)
Contains all styling for:
- Grid layouts (multi-motif and singleton)
- Modal windows
- Navigation bars
- Interactive elements (sliders, buttons, hover effects)
- Animations and transitions

**Exports**: `template_css`

### HTML Templates (`html_templates.jl`)
Contains page structure templates:
- `html_template`: Multi-motif pages (pairs, triplets, quadruplets)
- `html_template_singleton`: Singleton motif pages
- `html_hover_default`: Hover window for metadata
- `html_end`: Closing HTML tags

**Exports**: `html_template`, `html_template_singleton`, `html_hover_default`, `html_end`

### JavaScript - Core (`js_core.jl`)
Core JavaScript functionality:
- Page navigation generation
- JSON data loading
- Slider initialization
- Image/text updates

**Exports**: `script_core`

### JavaScript - Sequence Highlighting (`js_sequence.jl`)
Sequence highlighting features:
- FASTA file parsing
- CSV position data parsing
- Sequence highlighting with color coding
- Modal display for highlighted sequences

**Exports**: `script_sequence`

### JavaScript - Modals (`js_modals.jl`)
Modal window management for:
- Image display modal
- Text display modal (with copy functionality)
- Cluster view modal (image + text)
- Singleton motif modal
- Multi-motif modal
- Keyboard/click event handlers

**Exports**: `script_modals`

### JavaScript - Styling (`js_styling.jl`)
Dynamic styling based on Banzhaf values:
- Card border styling (thickness and color)
- Gradient effects
- Color coding (positive=red, negative=blue)

**Exports**: `script_styling`

## Usage

### Simple Usage (Backward Compatible)
The main `templates.jl` file in the parent directory now includes all modular components:

```julia
include("core/templates.jl")

# Use templates as before
css = template_css
html = html_template
js = script_template
```

### Modular Usage
You can include individual components as needed:

```julia
include("core/templates/css_template.jl")
include("core/templates/html_templates.jl")
include("core/templates/js_core.jl")
# ... include only what you need

# Or use the builder function
include("core/templates.jl")
js = build_script_template(
    j=2,
    upto=6,
    mode_counts=15,
    sequence_file_paths=["data/seq.fasta"]
)
```

### Custom Combinations
Build custom JavaScript by combining only needed components:

```julia
include("core/templates/js_core.jl")
include("core/templates/js_modals.jl")

# Custom script without sequence highlighting or styling
custom_js = script_core * "\n\n" * script_modals
```

## Template Variables

All templates use Mustache.jl syntax (`mt"..."`). Common variables:

### CSS Template
- No variables (static stylesheet)

### HTML Templates
- `protein_name`: Page title
- `j`: Current page number
- `upto`: Total number of pages
- `DF`: DataFrame with motif data

### JavaScript Templates
- `j`: Current page number
- `upto`: Total number of pages
- `mode_counts`: Number of motif modes
- `sequence_file_paths`: Array of FASTA file paths

## Migration Notes

No code changes required for existing scripts. The refactored structure:
- Maintains all original variable names
- Preserves all functionality
- Provides same exports
- Fully backward compatible

The original monolithic file is preserved as `templates_original.jl` for reference.

## Future Improvements

Potential enhancements enabled by modular structure:
- Add new modal types without affecting existing code
- Swap CSS themes by replacing `css_template.jl`
- Add alternative visualizations (e.g., D3.js) as new JS components
- Create variant HTML templates for different layouts
- Unit test individual components
