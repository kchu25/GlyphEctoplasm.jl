# Generation Module

Modular architecture for motif generation and visualization across both convolution-based and mutagenesis-based analyses.

## Directory Structure

```
generation/
├── generation.jl              # Main entry point - load this file
├── core/                      # Shared functionality (required by all)
│   ├── types.jl              # Configuration structs & data types
│   ├── data_structures.jl    # JSON/HTML dictionary management
│   ├── grouping.jl           # Grouping column specifications
│   ├── sorting.jl            # Pareto ranking & sorting utilities
│   ├── file_io.jl            # File I/O operations (save logos, CSVs, MEME)
│   └── rendering.jl          # HTML/JavaScript/CSS rendering
├── convolution/              # Convolution-based motif analysis
│   ├── singletons.jl         # Single filter motif processing
│   └── multi_motifs.jl       # Paired/triplet filter motif processing
├── mutagenesis/              # Mutation-based motif analysis
│   ├── singletons.jl         # Single mutation region filtering
│   ├── multi_regions.jl      # Multi-region mutation analysis
│   ├── matrix_operations.jl  # Matrix merging & accumulation
│   └── reference_matching.jl # Off-region search & reference matching
└── utils/                    # General utilities
    └── distance_keys.jl      # Distance key value extraction
```

## Usage

### Basic Setup

```julia
# Load all generation functionality
include("generation/generation.jl")

# Now you have access to all generation functions
```

### Convolution-Based Analysis

```julia
using DataFrames

# Create configuration
config = ConvMotifConfig(data;
    filter_len = 7,
    dpi = 65,
    alpha = 1.0,
    use_rna = false,
    xlim = (x_min, x_max),
    save_path = "output_dir"
)

# Initialize output dictionaries
json_motifs = init_json_dict()
html_dict = init_dict_for_html_render()

# Process singleton motifs
next_idx = process_singletons!(
    contributions_df, config, json_motifs, html_dict;
    motif_type = "singletons_high",
    group_id = "high_sing",
    button_text = "High Influence Singletons",
    start_idx = 1
)

# Process multi-motifs (pairs, triplets, etc.)
next_idx = process_multi_motifs!(
    df, config, json_motifs, html_dict;
    motif_size = 2,
    motif_type = "pairs_positive",
    group_id = "pairs_pos",
    button_text = "Positive Paired Motifs",
    start_idx = next_idx
)

# Render HTML output
render_and_save_outputs!(json_motifs, html_dict, 1;
    html_template = html_template_unified,
    script_template = script_template,
    css_template = template_css,
    save_path = "output_dir",
    nav_page_count = 3,
    sequence_paths = sequence_paths,
    page_title = "Motif Analysis",
    use_unified = true,
    enable_colored_borders = true
)
```

### Mutagenesis-Based Analysis

```julia
# Create configuration for mutagenesis
config = MutationRegionConfig(data;
    filter_len = 9,
    off_region_search = true,
    reduction_on_ref = true,
    dpi = 65,
    use_rna = false,
    xlim = xlim,
    save_path = "output_dir"
)

# Collect metadata for all motif sizes (singletons, pairs, triplets)
all_metadata = prepare_and_collect_mutation_metadata(
    contributions_df, df_multi_motifs, data, config;
    singleton_filter_pareto_rank = 1,
    split_by_sign = true
)

# Register all motifs with global sorting
next_idx = register_mutation_region_motifs!(
    json_motifs, html_dict, all_metadata;
    start_idx = 1,
    sort_globally = true,
    sort_by_pareto = true
)
```

## Key Features

### 1. Configuration Objects

Two main configuration types encapsulate all analysis parameters:

- **`ConvMotifConfig`**: For convolution-based analysis
- **`MutationRegionConfig`**: For mutagenesis analysis

These eliminate the need to pass dozens of individual parameters.

### 2. Modular Architecture

The codebase is organized by **concern**:

- **Core**: Functionality shared across all analysis types
- **Convolution**: Specific to filter-based convolution analysis
- **Mutagenesis**: Specific to mutation region analysis
- **Utils**: General-purpose utilities

### 3. Sorting & Prioritization

Multiple sorting strategies:

- **Median Banzhaf**: Sort by contribution magnitude
- **Pareto Ranking**: Multi-objective optimization (median + occurrence count)
- **Group-based**: Sort within groups (e.g., single regions, 2-regions, 3-regions)
- **Sign-based**: Positive contributions first, then negative

### 4. Flexible Grouping

The `build_grouping_columns()` function supports various grouping criteria:

- `:filter_index` - Group by which filter activated
- `:motifs` - Group by filter identities (multi-motifs)
- `:distances` - Group by inter-motif distances
- `:mutagenesis` - Complete specification for mutations
- `:motif_positions` - Group by occurrence positions

### 5. File I/O

Standardized file operations:

- **Logo plots**: PWM visualizations with optional highlighting
- **Influence plots**: Box/scatter plots of Banzhaf contributions
- **Positional CSVs**: Occurrence positions for each motif
- **MEME format**: Standard motif format for external tools

### 6. HTML Rendering

Dynamic HTML generation with:

- Multiple motif groups with toggle buttons
- Colored borders based on contribution sign/magnitude
- JSON data embedding for client-side interactivity
- Responsive layout with navigation

## Migration from Old Code

### Before (Old Structure)

```julia
include("generation_types.jl")
include("generation_helpers.jl")
include("generation_adders.jl")
include("generation_save_and_render.jl")
include("generation_singletons.jl")
include("generation_multi.jl")
include("generation_mutations.jl")
include("generation_mutations_subroutines.jl")
include("generation_mutagenesis_singletons.jl")
```

### After (New Structure)

```julia
include("generation/generation.jl")
```

That's it! All functionality is loaded in the correct order.

### API Compatibility

All public functions maintain backward compatibility. The refactoring only:

1. **Reorganizes** code into logical modules
2. **Improves** naming and documentation
3. **Adds** new features (e.g., Pareto ranking)
4. **Maintains** all existing functionality

## Design Principles

1. **Separation of Concerns**: Each file has a single, clear responsibility
2. **DRY (Don't Repeat Yourself)**: Shared code is in `core/`
3. **Configuration over Parameters**: Use config objects instead of many kwargs
4. **Explicit Dependencies**: Clear module loading order
5. **Backward Compatibility**: Existing code continues to work

## Testing

After refactoring, test that your existing scripts work:

```julia
# Should work exactly as before
include("generation/generation.jl")

# Your existing analysis code here...
```

## Future Enhancements

Possible improvements:

1. **Module System**: Convert to proper Julia modules with explicit exports
2. **Unit Tests**: Add comprehensive test suite
3. **Performance**: Profile and optimize hot paths
4. **Documentation**: Add more examples and tutorials
5. **Validation**: Input validation and error handling

## Questions?

See the original files for more details:

- Original: `generation_*.jl` files in project root
- Refactored: `generation/` directory with organized structure

The refactoring maintains all functionality while improving maintainability and extensibility.
