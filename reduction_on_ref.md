# The `reduction_on_ref` option

This is about a single toggle in the **mutation-region** pipeline that decides
whether the *reference* (consensus) sequence is used to **filter and reduce** what
each motif shows. It touches two things: (1) how regions are *counted and labelled*,
and (2) how the *logo* is drawn. Both live in
[`src/generation/mutagenesis/multi_regions.jl`](src/generation/mutagenesis/multi_regions.jl).

## The idea in one sentence

A mutation-region motif is built from raw count matrices over some positions. With
`reduction_on_ref=true`, those matrices are compared against the reference and the
columns/fragments that *don't match the reference* are dropped — so the count, the
span, and the logo all reflect the **reference-matching** part only. With
`reduction_on_ref=false`, nothing is dropped: you see the raw geometry.

## Where you set it

It's a keyword on the entry point, and it flows into the config object:

```julia
plot_motifs_mut_case(data, m, contributions_df_filtered, dfs;
    reduction_on_ref = false,   # <- this option
    ...)
```

```julia
# run_thru_mut.jl — passed straight into the config
m_config = MutationRegionConfig(data;
    ...,
    reduction_on_ref = reduction_on_ref
)
```

Heads up on defaults — they disagree depending on the door you come through:

- `plot_motifs_mut_case` defaults to **`false`** ([run_thru_mut.jl:8](src/run_thru_mut.jl#L8)).
- The `MutationRegionConfig` struct itself defaults to **`true`** ([types.jl:79](src/generation/core/types.jl#L79)).

Since `plot_motifs_mut_case` always passes its own value through, the effective
default *for the normal entry point* is `false`. Just don't assume "the config
default" is what you get.

## Effect 1 — fragment count and span

This is `compute_fragment_info` ([multi_regions.jl:143](src/generation/mutagenesis/multi_regions.jl#L143)).
The toggle picks one of two completely different branches:

```julia
function compute_fragment_info(count_mats, ref_pfms, start_positions, reduction_on_ref::Bool, motif_size::Int)
    if reduction_on_ref
        # Ask EntroPlots to detect the *actual* fragments against the reference.
        # This can REDUCE the count (a pair where one region doesn't match the
        # reference collapses to a single region, etc.)
        fragment_count, span_str = EntroPlots.count_fragments(count_mats, ref_pfms, start_positions)
        span_str = replace(span_str, "-" => ":")
    else
        if motif_size > 1
            fragment_count = motif_size           # <- hardcoded: pairs=2, triplets=3, ...
            if length(count_mats) == 1
                mat_width = size(count_mats[1], 2)
                span_str = "$(start_positions[1]):$(start_positions[1] + mat_width - 1)"
            else
                span_str = join(["$(pos):$(pos + size(mat,2) - 1)" for (pos, mat) in zip(start_positions, count_mats)], ", ")
            end
        else
            ...
        end
    end

    group_id    = fragment_count == 1 ? "single_region" : "$(fragment_count)_regions"
    button_text = fragment_count == 1 ? "Single Regions" : "$(fragment_count) Regions"
    return (fragment_count=fragment_count, span=span_str, group_id=group_id, button_text=button_text)
end
```

The crucial difference:

- **`true`** → `count_fragments` derives **both** `fragment_count` *and* `span` from the
  same reference comparison, so the group label ("2 regions") and the span string
  always agree.
- **`false`** → `fragment_count` is just `motif_size` (hardcoded), while the span is
  read off the raw matrices. These two can **disagree** — that's how you get a card
  filed under "2 regions" that displays a single contiguous span like `1:44` (the two
  regions merged into one matrix, but the label stayed at `motif_size = 2`).

## Effect 2 — the logo

The same flag is handed to the logo renderer as `filter_by_reference`
([multi_regions.jl:513](src/generation/mutagenesis/multi_regions.jl#L513)):

```julia
EntroPlots.save_logo_with_rect_gaps(
    meta.count_matrices, meta.positions, meta.total_length, paths.png.abs;
    reference_pfms = meta.references,
    ...,
    filter_by_reference = meta.reduction_on_ref   # <- same toggle
)
```

So with `true`, the drawn logo only keeps the reference-matching columns (gaps where
it doesn't match); with `false`, every column of the raw count matrix is drawn.

## Putting it together

| | `reduction_on_ref = false` (default) | `reduction_on_ref = true` |
|---|---|---|
| Region count | `motif_size`, fixed | recomputed from the reference (can shrink) |
| Span string | raw matrix extents | reference-matching fragments |
| Label vs span | can disagree (`2 regions` ↔ `1:44`) | always consistent |
| Logo | full raw matrix | reference-filtered columns |

## Which one do I want?

- Use **`false`** when you want the **raw** picture — every position the motif covers,
  with the group buckets just reflecting how many filters went into the motif
  (pair/triplet/…). Simple and literal, at the cost of the label sometimes not
  matching the span.
- Use **`true`** when you want the report to reflect what **actually matches the
  reference** — counts, spans, and logos all reduced consistently. A "pair" that only
  half-matches will honestly move into the single-region bucket.

See also [[pareto]] for how, once these groups are decided, the motifs are ordered
within each region-count bucket.
