# Multi-run: how several models become one displayed result

Read this before changing anything about runs, seeds or the run folder layout.

## The layout

```
<save_path>/
  models/          model_<seed>.jld2, processor_<seed>_pp_<i>.jld2   shared, seed-keyed
  json/            trial_seed_<seed>.json                            shared, from tuning
  results_*.csv    the tuning score table                            shared
  run_1/           renderings_*/  top_motifs.csv  motifs_cache_*  cache/  run_info.json
  run_2/           the same, from the next-best seed
  top_motifs.csv   a copy of the chosen run's, written by combine_runs
  runs.html        the run index
  runs.json        the same facts, machine-readable
```

**With `n_runs = 1` the `run_*` level does not exist at all.** Outputs go
straight into `<save_path>` exactly as they always did. That is the whole
backward-compatibility story, and it is worth preserving.

## Why only four things move

Model and tuning artifacts are already keyed by seed, so several runs coexist in
one `models/` and one `json/`. Tuning happens **once**: `tune_hyperparameters`
returns a row per trial, and every trial's config is on disk, so training the
runner-up seed reloads its config instead of re-tuning.

Only these are run-specific, and only these are addressed through
`MotifInference.output_path(trc)`:

| artifact | why it cannot be shared |
|---|---|
| `motifs_cache_output_<i>.jld2` | derived from this model |
| `cache/motifs_size_*.arrow` | derived from `m` and `processor` |
| `renderings_<feature>/` | this model's figures |
| `top_motifs.csv` | this model's ranking |

## What a combiner may and may not assume

**May assume**, because it was measured across two independent runs of 33 RNA
datasets:

* carrier ids are comparable across runs — `seq_index` in
  `mutation_regions_*/*.csv` and `data_index` in `indicator_points.csv` index
  the dataset, not the model, and the index universe was identical in 12/12
  datasets checked;
* reported positions are in **original untrimmed coordinates** in every run;
* window widths match across runs for a given dataset;
* three of the four top-mover ranking keys — `cluster_median`, `cluster_nnd`,
  `count` — are computed from labels and carriers, so they are on the assay's
  own scale and comparable between models. Only the Shapley median is on a
  model-internal scale.

**May not assume:**

* `filter_index` means anything across runs. The 50 filters have no canonical
  order. Match on carriers or positions, never on filter id.
* runs are of equal quality. Held-out R² differed by more than 0.1 between two
  runs on 8 of 33 datasets, once by 0.48. A combiner that pools runs blindly
  lets a bad model outvote a good one.
* agreement implies correctness. A dataset whose model learned nothing (R² 0.02)
  had two runs agree perfectly, because they agreed on the same nothing.

## Adding a policy

Subtype `RunCombiner` and add one method:

```julia
struct SupportFiltered <: RunCombiner
    tau::Float64
end
choose_runs(c::SupportFiltered, runs::Vector{RunRecord}) = ...
```

`combine_runs` handles discovery, the manifest and the top-level view.

A known cost, before anyone reaches for a stability filter: cross-run support is
easier to achieve for motifs with many carriers, so filtering on it is biased
against rare combinations — which may be exactly the findings of interest.
Measure the attrition on your own result before adopting a threshold.

## Convolution

Nothing in `multirun.jl` is mutagenesis-specific, and `plot_motifs_conv_case`
produces the same folder shape. To extend:

1. give `run_thru_conv.jl` the same `output_path(trc)` treatment the mut path has
   — it already goes through the shared `render_html`, so this may need nothing;
2. note the conv generalization page is rendered unconditionally, whereas the
   mut one appears only when held-out points exist. A combiner that assumes
   `index2.html` exists will be right for conv and wrong for mut.

---

# Consensus top movers

One page above the runs, showing the findings the runs agree on. Built by
`consensus_top_movers(save_path; combiner=ConsensusTopMovers())`.

## What it writes into `<save_path>`

| file | contents |
|---|---|
| `index.html` | the consensus top-movers page: 5 increasing and 5 decreasing findings, each a real motif from a real run |
| `index4.html` | readme placeholder |
| `consensus_top_motifs.csv` | one row per displayed finding, with full provenance |
| `consensus.json` | the selection and its parameters |

Nothing inside the run folders is modified, so re-running with different
parameters is free and reversible.

## The two steps

Scoring uses cross-run agreement, so it must happen **before** duplicates are
collapsed — otherwise the evidence it needs is gone.

1. **Score.** `motif_support` gives each motif the mean, over the other runs, of
   its best match there. Matching (`motif_omega`) requires the same direction,
   overlapping windows, and shared carriers. A run with no candidate contributes
   `0`, not a skipped term.
2. **Collapse.** `dedup_motifs` greedily emits the highest-support motif and drops
   everything of the same direction whose window overlaps it by at least `rho`.
   Greedy rather than connected components, which chain one window to the next
   and fuse a whole region into a single finding.

Full specification, with the reasoning for each condition and the ground-truth
validation: `multirun_findings/12_method.md`.

## Tracing a displayed row

`consensus_top_motifs.csv` carries `from_run`, `logo_path` and `carriers_csv`,
each relative to the folder holding the CSV. So a row resolves directly to the
model that produced it, the image on the page, and the variants behind it:

```
from_run      run_3
logo_path     run_3/renderings_1/mutation_regions_1/12_15:22.png
carriers_csv  run_3/renderings_1/mutation_regions_1/12_15:22.csv
```

`support`, `n_runs` and `n_motifs` say how much corroboration the finding has.

## Requirements

Each run folder needs `top_motifs.json`, written beside `top_motifs.csv` on every
mutagenesis render. Runs produced before that dump existed must be re-rendered;
with the seed supplied this reloads the cached model rather than retraining.

## Convolution

Nothing in the combiner is mutagenesis-specific. `plot_motifs_conv_case` produces
the same folder shape, so the conv path needs only `write_top_movers_json` wiring
in beside its `write_top_movers_csv` call to use all of this unchanged.
