"""
    save_indicator_points(pts, all_indices, save_path; file_name="indicator_points.csv")

Dump the per-datapoint coordinates behind the indicator (yy-KDE) plots.

Every `yy_kde_intersect_*.png` in a rendering folder draws the *same* point
cloud — `pts.proc_prod` on x, `pts.labels` on y — and differs only in which
points are highlighted. That membership is already on disk (each motif's
`<name>.csv` lists its `seq_index`es), but the coordinates themselves only ever
existed inside the running process: `pts` is recomputed from the cached model
and processor on every run and was never serialised. Writing it once here makes
every indicator plot in the folder reproducible from files alone — no model, no
GPU, no package stack.

One row per population point, in `all_indices` order — the ordering the
plotting code already assumes — so row `i` is the point for sequence
`all_indices[i]`. To rebuild a motif's highlight mask:

    pop  = CSV.read("<rendering>/indicator_points.csv", DataFrame)
    mem  = Set(CSV.read("<rendering>/mutation_regions_1/<name>.csv", DataFrame).seq_index)
    mask = [i in mem for i in pop.data_index]   # matches `is_in_intersect`

Returns the path written, or `nothing` when there is nothing to write (`pts` or
`all_indices` missing, or the two are not aligned). A misalignment warns rather
than throws, matching the best-effort handling the indicator plots themselves
get — a failed dump must never take down a render.
"""
function save_indicator_points(pts, all_indices, save_path;
        file_name::AbstractString="indicator_points.csv")
    (pts === nothing || all_indices === nothing) && return nothing

    idx       = vec(collect(all_indices))
    labels    = vec(pts.labels)
    proc_prod = vec(pts.proc_prod)

    n = length(idx)
    if !(length(labels) == length(proc_prod) == n)
        @warn "indicator points not saved: `pts` is not aligned to `all_indices`" n_all_indices = n n_labels = length(labels) n_proc_prod = length(proc_prod)
        return nothing
    end

    df = DataFrame(data_index = idx, proc_prod = proc_prod, label = labels)

    # `proc_prod` is what the indicator plots put on the x-axis; `predictions`
    # is the model's own output (what the generalization page plots against the
    # labels). Keep both when available — one extra column, and it saves a
    # second forward pass if the generalization scatter is ever wanted as data.
    if hasproperty(pts, :predictions)
        preds = vec(pts.predictions)
        length(preds) == n && (df.prediction = preds)
    end

    mkpath(save_path)
    out = joinpath(save_path, file_name)
    CSV.write(out, df)
    @info "Saved indicator-plot points" path = out n = n
    return out
end
