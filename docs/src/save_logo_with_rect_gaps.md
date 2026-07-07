# `save_logo_with_rect_gaps`

Renders a multi-region sequence logo — several count matrices laid out along a
shared coordinate axis with visual gaps (rectangles) between them — and writes it
to an image file.

This function lives in the external **EntroPlots** package (v0.2.1,
`src/old/plot_logo_w_arr_gaps.jl:210`), not in GlyphEctoplasm. It's the call that
produces every mutagenesis / multi-region logo PNG. In this repo it's invoked from
`render_one_motif!` in
[`src/generation/mutagenesis/multi_regions.jl`](../../src/generation/mutagenesis/multi_regions.jl#L761),
using fields off the `meta::MotifMetadata` object.

## Signature

```julia
save_logo_with_rect_gaps(
    count_matrices,
    starting_indices,
    total_length,
    save_name::String;
    arrow_shape_scale_ratio::Real = 1.0,
    height_top::Real            = 2.0,
    dpi                         = 65,
    rna                         = false,
    protein                     = false,
    uniform_color               = true,
    basic_fcn                   = get_rectangle_basic,
    xrotation                   = 0,
    reference_pfms::Union{Nothing, Vector{BitMatrix}} = nothing,
    filter_by_reference::Bool   = true,
    filter_tolerance::Real      = 1e-9,
)
```

Returns the underlying `Plots.jl` plot object `p` and, as a side effect, saves it
to `save_name` via `savefig`.

## Positional arguments

| Argument | Type | Meaning |
|---|---|---|
| `count_matrices` | `Vector` of matrices | One count matrix per region. Rows = alphabet (4 for DNA/RNA, 20 for protein), columns = positions. Each column must sum to `> 0` and hold no negatives (asserted). Normalized internally to a PFM via `counts ./ sum(counts, dims=1)`. |
| `starting_indices` | `Vector{Int}` | The start coordinate of each matrix on the shared axis — one integer per matrix, **must be sorted ascending** and non-overlapping. Region `i` occupies columns `starting_indices[i] : starting_indices[i] + size(count_matrices[i], 2) - 1`. Drives the x-tick labels and the gap sizes. |
| `total_length` | `Int` | Total width (in columns) of the coordinate space the regions live in. Used to compute the gap increments between regions so the layout spans the full span. |
| `save_name` | `String` | Output file path. Extension picks the format (`.png` here). |

In GlyphEctoplasm these map to `meta.count_matrices`, `meta.positions`,
`meta.total_length`, and `paths.png.abs` respectively (see
[multi_regions.jl:761-771](../../src/generation/mutagenesis/multi_regions.jl#L761-L771)).
Note `meta.positions` is a `Vector{Int}` of start positions — not ranges.

## Keyword arguments

| Keyword | Default | Meaning |
|---|---|---|
| `arrow_shape_scale_ratio` | `1.0` | Scales the gap/rectangle shapes drawn between regions. |
| `height_top` | `2.0` | Top of the y-axis / logo height in bits. |
| `dpi` | `65` | Output resolution. GlyphEctoplasm passes `meta.dpi`. |
| `rna` | `false` | Use RNA alphabet (U instead of T). GlyphEctoplasm passes `meta.use_rna`. |
| `protein` | `false` | Protein (20-row) logo. GlyphEctoplasm passes `size(count_matrices[1], 1) == 20`. |
| `uniform_color` | `true` | Color all glyphs uniformly rather than per-base coloring. |
| `basic_fcn` | `get_rectangle_basic` | Shape function for the inter-region gap marks. |
| `xrotation` | `0` | Rotation (degrees) of x-tick labels. GlyphEctoplasm passes `35`. |
| `reference_pfms` | `nothing` | Optional `Vector{BitMatrix}`, one per region. When given, each region's logo is drawn relative to its reference, and its length must match `count_matrices`. GlyphEctoplasm passes `meta.references`. |
| `filter_by_reference` | `true` | When a reference is supplied, drop columns that match the reference before plotting (keeps only the mutated/informative columns). GlyphEctoplasm passes `meta.reduction_on_ref`. Errors if *every* column matches ("nothing to plot after filtering"). |
| `filter_tolerance` | `1e-9` | Tolerance for the reference-match test used by `filter_by_reference`. |

## What it does internally

1. **Validate** each count matrix (no negatives, no zero-sum columns) and, if given,
   check `length(reference_pfms) == length(count_matrices)`.
2. Delegate to `logoplot_with_rect_gaps`, which:
   - Optionally **filters** columns against the reference (`apply_count_filter`).
   - **Normalizes** count matrices to PFMs.
   - Computes per-region offsets and the adjusted total width via
     `get_offset_from_start(starting_indices, pfms, total_length)` — the gaps
     between regions come from the difference between consecutive
     `starting_indices` and the matrix widths.
   - Draws each region's logo with `logoplot!` at its computed x-offset, and the
     rectangle gap shapes with `arrowplot!`.
   - Builds x-tick labels as the real coordinates
     `starting_indices[i] + 0 : starting_indices[i] + width - 1`.
3. **Save** with `savefig(p, save_name)` and return `p`.

## Invariants / gotchas

- `starting_indices` must be **sorted** and the region ranges must **not overlap**
  (`check_overlap` errors otherwise).
- With `filter_by_reference = true`, a motif whose columns all equal the reference
  raises an error rather than emitting an empty logo — callers in this repo wrap
  the render in a per-motif `try` so one such failure doesn't abort the batch.
- `protein` is inferred from the matrix row count (20) at the call site; RNA vs DNA
  is an explicit flag.

## Call site in this repo

```julia
EntroPlots.save_logo_with_rect_gaps(
    meta.count_matrices, meta.positions, meta.total_length,
    paths.png.abs;
    reference_pfms = meta.references,
    dpi            = meta.dpi,
    rna            = meta.use_rna,
    xrotation      = 35,
    protein        = size(meta.count_matrices[1], 1) == 20,
    uniform_color  = true,
    filter_by_reference = meta.reduction_on_ref,
)
```
