# Migration Guide: Using the Refactored Generation Module

This guide helps you update existing scripts to use the new modular generation system.

## Quick Start

### Before (Old Way)
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

### After (New Way)
```julia
include("generation/generation.jl")
```

That's the only change needed! All functions remain available and work identically.

## File-by-File Migration

### Scripts Using Convolution Analysis

**Files to update**: Scripts with conv-motif analysis (singletons, pairs, triplets)

**Example**: `try_home.jl`, `try_tune_*.jl`

#### Changes Required

1. **Remove old includes**:
```julia
# DELETE these lines:
include("generation_types.jl")
include("generation_helpers.jl")
include("generation_adders.jl")
include("generation_save_and_render.jl")
include("generation_singletons.jl")
include("generation_multi.jl")
```

2. **Add single new include**:
```julia
# ADD this line:
include("generation/generation.jl")
```

3. **Everything else stays the same**:
```julia
# This code works unchanged:
config = ConvMotifConfig(data; filter_len=7, dpi=65, ...)
json_motifs = init_json_dict()
html_dict = init_dict_for_html_render()

next_idx = process_singletons!(contributions_df, config, json_motifs, html_dict; ...)
next_idx = process_multi_motifs!(df, config, json_motifs, html_dict; ...)

render_and_save_outputs!(json_motifs, html_dict, 1; ...)
```

### Scripts Using Mutagenesis Analysis

**Files to update**: Scripts with mutation region analysis

#### Changes Required

1. **Remove old includes**:
```julia
# DELETE these lines:
include("generation_types.jl")
include("generation_helpers.jl")
include("generation_mutations.jl")
include("generation_mutations_subroutines.jl")
include("generation_mutagenesis_singletons.jl")
```

2. **Add single new include**:
```julia
# ADD this line:
include("generation/generation.jl")
```

3. **Everything else stays the same**:
```julia
# This code works unchanged:
config = MutationRegionConfig(data; filter_len=9, ...)

all_metadata = prepare_and_collect_mutation_metadata(
    contributions_df, df_multi_motifs, data, config; ...
)

next_idx = register_mutation_region_motifs!(
    json_motifs, html_dict, all_metadata; ...
)
```

## Common Patterns

### Pattern 1: Mixed Singletons and Multi-Motifs

```julia
include("generation/generation.jl")

# Initialize
json_motifs = init_json_dict()
html_dict = init_dict_for_html_render()

# Process different groups
next_idx = process_singletons!(...; group_id="high_sing", start_idx=1)
next_idx = process_multi_motifs!(...; group_id="pairs_pos", start_idx=next_idx)
next_idx = process_multi_motifs!(...; group_id="triplets_pos", start_idx=next_idx)

# Render unified page
render_and_save_outputs!(...; use_unified=true)
```

### Pattern 2: Multiple Mutation Regions

```julia
include("generation/generation.jl")

# Collect metadata for all sizes
all_metadata = prepare_and_collect_mutation_metadata(
    contributions_df, df_multi_motifs, data, config;
    singleton_filter_pareto_rank=1
)

# Register with global sorting
register_mutation_region_motifs!(
    json_motifs, html_dict, all_metadata;
    sort_globally=true,
    sort_by_pareto=true
)
```

### Pattern 3: Legacy Parameters (Still Supported)

```julia
# You can still use the old interface with individual parameters:
process_singletons!(
    contributions_df, data, json_motifs, html_dict;
    SAVE_PATH="tmp",
    filter_len=7,
    dpi=65,
    alpha=1.0,
    use_rna=false,
    xlim=nothing
)

# But the new config-based interface is cleaner:
config = ConvMotifConfig(data; filter_len=7, dpi=65, alpha=1.0, use_rna=false, xlim=nothing, save_path="tmp")
process_singletons!(contributions_df, config, json_motifs, html_dict)
```

## Checklist for Each Script

- [ ] Replace all `generation_*.jl` includes with `include("generation/generation.jl")`
- [ ] Keep all other includes (core, plotting, helpers) unchanged
- [ ] Test that the script runs without errors
- [ ] Verify output files are identical to before
- [ ] (Optional) Refactor to use config objects instead of individual parameters

## Scripts to Update

Search your project for files including generation modules:

```bash
# Find all scripts that need updating
grep -r "include.*generation_" *.jl
```

Common files that may need updates:
- `try_home.jl` ✅ (already updated)
- `try_home3.jl`
- `try_home4.jl`
- `try_home4_updated.jl`
- `try_tune_lac.jl`
- `try_tune_rna_compete.jl`
- `develop.jl`
- Any custom analysis scripts

## Troubleshooting

### "Function not found" errors

**Problem**: Getting `UndefVarError` for generation functions

**Solution**: Make sure you're including the new file:
```julia
include("generation/generation.jl")  # Not generation_*.jl
```

### "File not found" errors

**Problem**: `cannot open file generation/generation.jl`

**Solution**: Make sure you're in the project root directory:
```julia
cd("/home/shane/Desktop/academia/code/experimental/bio_seq_cnn_tmp5")
include("generation/generation.jl")
```

### Dependencies still missing

**Problem**: Errors about missing functions like `init_json_dict` or `m_symbols`

**Solution**: Keep these includes (they're external dependencies):
```julia
include("core/logo_saving.jl")
include("core/json_html_dict.jl")
include("core/consensus.jl")
include("core/html_generation.jl")
include("core/templates.jl")
include("core/path_utils.jl")
include("old_refactored/banzhaf_conv_assign/helpers.jl")
```

## Testing Your Migration

After updating a script, test it:

```bash
julia your_script.jl
```

Expected results:
1. Script runs without errors
2. Output files are generated in expected locations
3. HTML visualizations render correctly
4. Output is identical to before refactoring

## Benefits of Migration

Once migrated, you get:

✅ **Simpler includes**: One line instead of 9  
✅ **Better organization**: Know where to find specific functionality  
✅ **Improved documentation**: Comprehensive README and inline docs  
✅ **Future-proof**: Easier to extend and maintain  
✅ **Same functionality**: Everything works exactly as before  

## Questions?

See:
- `generation/README.md` - Full documentation
- `generation/REFACTORING_SUMMARY.md` - What changed and why
- Original files - Still available in project root for reference

---

**Happy migrating! 🚀**

The refactored system maintains full backward compatibility while providing better organization and maintainability.
