# Generation Module - Quick Reference

## 🚀 Quick Start

```julia
include("generation/generation.jl")  # That's it!
```

## 📁 Module Structure

```
generation/
├── core/              # Shared by all analysis types
├── convolution/       # Conv-motif specific (singletons, pairs, triplets)
├── mutagenesis/       # Mutation-region specific
└── utils/             # General utilities
```

## ✅ What You Get

- **Types**: `ConvMotifConfig`, `MutationRegionConfig`, `MotifData`, `MotifMetadata`
- **Convolution**: `process_singletons!`, `process_multi_motifs!`
- **Mutagenesis**: `collect_mutation_region_metadata`, `register_mutation_region_motifs!`
- **Utilities**: `build_grouping_columns`, `filter_pareto_rank`, `build_sorted_keys_and_maps`
- **I/O**: `save_motif_logo`, `save_influence_plot`, `build_metadata_texts`
- **Rendering**: `render_and_save_outputs!`

## 🔧 Common Tasks

### Convolution Analysis (Singletons + Multi-motifs)

```julia
# Create config
config = ConvMotifConfig(data; filter_len=7, dpi=65, save_path="tmp")

# Initialize
json_motifs = init_json_dict()
html_dict = init_dict_for_html_render()

# Process motifs
next_idx = process_singletons!(contributions_df, config, json_motifs, html_dict; start_idx=1)
next_idx = process_multi_motifs!(df, config, json_motifs, html_dict; motif_size=2, start_idx=next_idx)

# Render
render_and_save_outputs!(json_motifs, html_dict, 1; 
    html_template=html_template_unified, 
    script_template=script_template,
    save_path="tmp", 
    use_unified=true)
```

### Mutagenesis Analysis (Mutation Regions)

```julia
# Create config
config = MutationRegionConfig(data; filter_len=9, off_region_search=true, save_path="tmp")

# Collect metadata from all sizes
all_metadata = prepare_and_collect_mutation_metadata(contributions_df, df, data, config)

# Register with Pareto sorting
next_idx = register_mutation_region_motifs!(json_motifs, html_dict, all_metadata; 
    sort_by_pareto=true)

# Render
render_and_save_outputs!(...)
```

## 📊 Sorting Behavior

| Analysis Type | Default Sort | Pareto Option |
|---------------|-------------|---------------|
| Conv Singletons | Median Banzhaf ↓ | Yes (optional) |
| Conv Multi-motifs | Median Banzhaf ↓ | No |
| Mutation Regions | Pareto (auto) | Yes (default) |

## 🔍 Finding Functions

```
core/types.jl           → Config structs
core/sorting.jl         → Pareto ranking, build_sorted_keys_and_maps
core/grouping.jl        → build_grouping_columns
core/file_io.jl         → save_motif_logo, build_motif_paths
core/rendering.jl       → render_and_save_outputs!
convolution/singletons  → process_singletons!
convolution/multi       → process_multi_motifs!
mutagenesis/multi       → collect/register mutation regions
utils/symbols.jl        → m_symbols, d_symbols, m_position_symbols
```

## 📚 Documentation

- **README.md** - Full usage guide with examples
- **MIGRATION_GUIDE.md** - Step-by-step migration from old code
- **REFACTORING_SUMMARY.md** - What changed and why
- **MODULE_INDEPENDENCE.md** - Dependency management

## ✨ Key Features

- ✅ **Self-contained** - No external dependencies on old_refactored/
- ✅ **Modular** - Clear separation of concerns
- ✅ **Backward compatible** - Existing code works unchanged
- ✅ **Well-documented** - Comprehensive inline and external docs
- ✅ **Flexible** - Config objects OR individual parameters
- ✅ **Extensible** - Easy to add new analysis types

## 🎯 Updated Scripts

- ✅ `try_home.jl` - Convolution analysis
- ✅ `try_home4_updated.jl` - Mutagenesis analysis

## 🐛 Troubleshooting

**"Function not found"**
```julia
# Make sure you include the main file:
include("generation/generation.jl")
```

**"Still depends on helpers.jl"**
```julia
# Remove this line (no longer needed):
# include("old_refactored/banzhaf_conv_assign/helpers.jl")
```

**"Missing core functions"**
```julia
# These are still needed (project-wide utilities):
include("core/logo_saving.jl")
include("core/json_html_dict.jl")
# ... etc
```

## 📞 Need Help?

Check the comprehensive documentation:
1. **README.md** - Start here
2. **MIGRATION_GUIDE.md** - Updating existing scripts
3. **Inline docs** - Function-level documentation

---

**Happy analyzing! 🎉**
