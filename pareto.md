# How the mutation-region motifs are sorted (Pareto)

This is about the order motifs appear in the **mutation-region** report — i.e. what
`register_mutation_region_motifs!` does before it renders cards. The interesting
bit is `sort_by_group_and_pareto` in
[`src/generation/mutagenesis/multi_regions.jl`](src/generation/mutagenesis/multi_regions.jl).

## The idea

Each motif has two numbers we care about:

- **`abs(median)`** — how *strong* its effect is (magnitude of the median Banzhaf
  contribution, sign ignored).
- **`count`** — how *often* it shows up (number of occurrences).

We'd like both to be high, but they trade off: a rare motif can have a huge effect,
a common one a modest effect. Rather than blend them into one score (which forces an
arbitrary weighting), we rank by **Pareto dominance**.

## What "Pareto rank" means here

Motif A *dominates* motif B if A is **≥ B on both** numbers and **strictly greater on
at least one**. In other words, A is just plain better — no trade-off to argue about.

- **Rank 1** = the motifs nobody dominates. These are the best available trade-offs
  between "strong" and "frequent" (the Pareto front).
- Remove rank 1, look at what's left, and the new front becomes **rank 2**. Repeat
  until everything has a rank.

So lower rank = closer to the best-of-both-worlds frontier. (`compute_pareto_ranks_subset_wrapped`
does exactly this peel-the-front loop.)

## The full ordering

Pareto rank isn't the only thing — it's one tier in a hierarchy. Motifs are sorted by,
in order:

1. **Sign** — positive-contribution motifs first, then negative.
2. **Group** — `single_region` first, then `2_regions`, `3_regions`, … (by region count).
3. **Pareto rank** — within each sign+group bucket, rank 1 first, then 2, …
4. **Tie-break** — within the same rank, by `abs(median)`:
   - positives: **strongest first** (descending),
   - negatives: **weakest first** (ascending).

Worth noting: Pareto ranks are computed **per (sign, group) bucket**, not globally — so
a 3-region motif is only ever compared against other 3-region motifs of the same sign.

## The toggle

`register_mutation_region_motifs!` takes `sort_by_pareto`:

- `true` (default) → the Pareto scheme above.
- `false` → a simpler sort: same sign-then-group hierarchy, but ordered directly by
  `abs(median)` (count is ignored entirely).

There's also `sort_globally`; set it to `false` to skip sorting and keep whatever order
the motifs were collected in.
