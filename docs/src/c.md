# Tutorial: `obtain_count_and_reference_matrices`

This is the function that turns grouped mutagenesis rows into the three arrays the
renderer needs for each motif: the **count matrices**, the **reference windows**,
and the **positions**. Those three become `meta.count_matrices`,
`meta.references`, and `meta.positions`, which are what
`save_logo_with_rect_gaps` draws (see [`save_logo_with_rect_gaps.md`](save_logo_with_rect_gaps.md)).

Source: [`src/generation/mutagenesis/reference_matching.jl:87`](../../src/generation/mutagenesis/reference_matching.jl#L87).

## Where it sits in the pipeline

```
df_mutated
  └─ build_mutation_aggregates(df_mutated, config, motif_size)   multi_regions.jl:112
        ├─ groupby(:mutagenesis columns)                    → gdf, sorted_keys
        └─ obtain_count_and_reference_matrices(...)         ← THIS FUNCTION
              → count_matrices_vec, reference_matrices_vec, adjusted_positions_vec
                    │
                    ▼  (per key k, packed into MotifData)
              collect_mutation_region_metadata → MotifMetadata (`meta`)
                    │
                    ▼
              render_one_motif! → save_logo_with_rect_gaps(meta.count_matrices,
                                                            meta.positions, ...)
```

The call site ([multi_regions.jl:118-122](../../src/generation/mutagenesis/multi_regions.jl#L118-L122)):

```julia
count_matrices_vec, reference_matrices_vec, adjusted_positions_vec =
    obtain_count_and_reference_matrices(
        sorted_keys, gdf, config.data, config.reference_seq,
        build_grouping_columns(:motif_positions; motif_size=motif_size),
        config.filter_len; T=config.float_type,
        off_region_search=config.off_region_search)
```

## Input arguments

| Argument | Type | What it is |
|---|---|---|
| `sorted_keys` | ordered collection | The group keys to iterate, in render order. Each key identifies one motif (a `(m1,…)`/position tuple from the `:mutagenesis` grouping). Output dicts are keyed by these. |
| `gdf` | `GroupedDataFrame` | The mutated rows grouped by motif. For a key `k`, `gdf[k]` is the sub-DataFrame of all data points that hit that motif. Must contain `:data_pt_index` and the motif-position columns. |
| `data` | dataset object (`config.data`) | Holds `data.onehot_sequences` — the one-hot encoded sequences `X`, a 4×L×1×N array (4 bases, L positions, N sequences) — and `data.prefix_offset`, an integer added to raw positions to map sequence coordinates into reference coordinates. |
| `reference_seq` | `BitMatrix` (4 × seq_length) | The reference (consensus/wild-type) one-hot sequence. Windows of it are sliced out as each motif's reference, and (when `off_region_search`) it's the target the counts are aligned against. |
| `mp_syms` | `Vector{Symbol}` | The motif-position column names, one per region — e.g. `[:m1_pos, :m2_pos]` for a pair. `n_motifs = length(mp_syms)` is how many regions this motif has. Built by `build_grouping_columns(:motif_positions; motif_size)`. |
| `filter_len` | `Int` | Window width (columns) of each count matrix / reference window. Must be `> 0`. |
| `T` | type (kw, default `Float32`) | Element type of the count matrices (`config.float_type`). |
| `off_region_search` | `Bool` (kw, default `false`) | If `true`, don't trust the recorded position — search ±3 nt around it for the best-matching reference window, then merge regions that end up overlapping. |

## What it returns

Three `Dict`s, each keyed by the same `sorted_keys`:

| Output | Type | Contents |
|---|---|---|
| `count_matrices_vec` | `Dict{key, Vector{Matrix{T}}}` | Per key, one `4 × filter_len` count matrix **per region** — how often each base appeared at each window column across all data points in that group. |
| `reference_matrices_vec` | `Dict{key, Vector{SubArray{Bool,2}}}` | Per key, one reference window (`@view` into `reference_seq`) per region, aligned to the chosen start. |
| `adjusted_positions_vec` | `Dict{key, Vector{Int}}` | Per key, the reference **start position** of each region (a plain `Int`, not a range). This is what becomes `meta.positions`. |

## How it builds them — step by step

For each key `k` in `sorted_keys` (skipping empty groups):

**1. Allocate one 3-D scratch array for all regions.**
```julia
h_dim   = size(X, 1)                     # 4 for DNA/RNA
mats_3d = zeros(T, h_dim, filter_len, n_motifs)   # one slice per region
idx_col = g[!, :data_pt_index]           # which sequences are in this group
```
Doing all regions in a single allocation (`mats_3d`) instead of `n_motifs` separate
matrices is a deliberate memory optimization.

**2. For each region `j` (`mp_sym` in `mp_syms`):** call
`accumulate_and_find_reference!` on the `j`-th slice `view(mats_3d, :, :, j)`. That
helper does two things ([reference_matching.jl:44](../../src/generation/mutagenesis/reference_matching.jl#L44)):

- **Accumulate counts.** It reads the region's raw start `s_raw = g[1, mp_sym]`
  (identical for every row in the group), computes the window `s_x : s_x+filter_len-1`,
  and sums the one-hot slices of every data point in the group into the count
  matrix via `accumulate_motif_counts!` → `accumulate_window!`
  ([logo_saving.jl:7](../../src/core/logo_saving.jl#L7)):
  ```julia
  dest[row, col+1] += X[row, s+col, 1, didx]   # for each data point didx
  ```
  So column `c` of the count matrix is a histogram of bases seen at reference
  position `s+c` across the group.

- **Pick the reference start `s_ref_best`.** Default is just `s_raw + offset`
  (`offset = data.prefix_offset`). If `off_region_search`, it instead scans
  `delta ∈ -3:+3` around `s_raw + offset` and keeps the start whose reference
  window maximizes the dot-product score against the accumulated counts
  (`find_best_reference_position` → `compute_dotproduct_score`). This corrects for
  small positional drift between the recorded position and the reference.

**3. Store per-region results.**
```julia
adjusted_positions[j] = s_ref_best
reference_here[j]     = @view reference_seq[:, s_ref_best : s_ref_best+filter_len-1]
```

**4. Materialize the count matrices** by splitting the 3-D array into a vector:
```julia
count_matrices_vec[k]     = [mats_3d[:, :, j] for j in 1:n_motifs]
reference_matrices_vec[k] = reference_here
adjusted_positions_vec[k] = adjusted_positions
```

**5. (only if `off_region_search` and `n_motifs > 1`) Merge overlapping regions.**
Because the ±3 search can nudge two adjacent regions into overlapping (or touching)
windows, `merge_overlapping_matrices` fuses them into a single wider matrix and
recomputes the positions ([matrix_operations.jl](../../src/generation/mutagenesis/matrix_operations.jl)).
The reference views are then rebuilt to match the merged widths:
```julia
reference_matrices_vec[k] = [
    @view reference_seq[:, pos : pos+size(mat,2)-1]
    for (pos, mat) in zip(merged_positions, merged_mats)
]
```
This is why, after this stage, a region's matrix width may exceed `filter_len` and
why its span is always reconstructed as `pos : pos + size(mat,2) - 1` rather than
`pos : pos + filter_len - 1`.

## Worked mini-example

A 2-region motif (`n_motifs = 2`, `filter_len = 9`), one key `k` whose group has
40 data points, `off_region_search = false`, `prefix_offset = 5`:

- Region 1 raw start `s_raw = 36` → counts window `36:44` accumulated over the 40
  sequences → a `4×9` matrix; reference start `36 + 5 = 41`.
- Region 2 raw start `s_raw = 55` → counts window `55:63` → another `4×9` matrix;
  reference start `60`.

Results for key `k`:
```julia
count_matrices_vec[k]     == [mat1 (4×9), mat2 (4×9)]
adjusted_positions_vec[k] == [41, 60]
reference_matrices_vec[k] == [ref[:,41:49], ref[:,60:68]]
```

Downstream, `save_logo_with_rect_gaps(count_matrices, [41,60], total_length, …)`
places region 1 at column 41, region 2 at column 60, and draws a rectangle gap for
the columns 50–59 in between.

## Gotchas

- **Counts, not frequencies.** Matrices hold raw base counts; normalization to a
  PFM (`counts ./ sum(counts, dims=1)`) happens later inside EntroPlots.
- **`adjusted_positions` are scalar starts**, reference-coordinate (raw +
  `prefix_offset`), and must stay sorted/non-overlapping for the logo layout.
- **Empty groups are skipped** — a key with no rows simply gets no dict entry, so
  don't assume every `sorted_key` appears in the outputs.
- **`off_region_search` can change matrix widths** via merging; anything reading the
  outputs must derive spans from `size(mat, 2)`, not from `filter_len`.
