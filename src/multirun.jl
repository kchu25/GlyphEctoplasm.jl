# ─────────────────────────────────────────────────────────────────────────────
# Combining several runs into one displayed result.
#
# A multi-run pipeline writes each seed's outputs into its own folder:
#
#   <save_path>/            models/   json/   results_*.csv     shared
#   <save_path>/run_1/      renderings_*/  top_motifs.csv  motifs_cache_*  run_info.json
#   <save_path>/run_2/      the same, from the next-best seed
#
# Every run is a complete, self-consistent rendering. Nothing here merges or
# synthesises a motif: a combiner only DECIDES which run's material is shown,
# and records why. That keeps every displayed logo a real logo from a real
# model, with real carriers behind it.
#
# EXTENDING THIS. Add a struct subtyping `RunCombiner` and one method:
#
#     choose_runs(c::MyCombiner, runs::Vector{RunRecord}) -> Vector{RunRecord}
#
# `combine_runs` handles discovery, the manifest and the top-level view. Two
# facts make richer policies possible without retraining anything:
#
#   * carrier ids (`seq_index` in mutation_regions_*/*.csv, `data_index` in
#     indicator_points.csv) index the DATASET, not the model, so they are
#     directly comparable between runs;
#   * reported positions are in original untrimmed coordinates in every run.
#
# `filter_index` is NOT comparable across runs -- the filters have no canonical
# order -- so match on carriers or positions, never on filter id.
#
# CONVOLUTION. Nothing here is mutagenesis-specific; the layout is produced by
# `plot_motifs_conv_case` in the same shape. When conv gains multi-run support
# the only extra work is that its generalization page is unconditional, whereas
# the mutagenesis page is rendered only when held-out points exist.
# ─────────────────────────────────────────────────────────────────────────────

"""One run's folder, as found on disk."""
struct RunRecord
    dir::String          # absolute path to the run folder
    name::String         # "run_1"
    rank::Int            # 1 = best-scoring seed; 0 when unknown
    seed::Int            # -1 when unknown
    has_summary::Bool    # a top_motifs.csv exists
end

"""How to pick what gets displayed. See the note at the top of this file."""
abstract type RunCombiner end

"""
    UseSingleRun(index = 1)

Show one run and nothing else. `index` is the RANK, so the default shows the
best-scoring seed. This is the default policy: it changes no numbers and keeps
every page internally consistent, which a mixed page would not.
"""
struct UseSingleRun <: RunCombiner
    index::Int
end
UseSingleRun() = UseSingleRun(1)

"""
    discover_runs(save_path) -> Vector{RunRecord}

Every `run_*` folder under `save_path`, ordered by rank. Reads the `run_info.json`
the pipeline stamps into each; falls back to the trailing digits of the folder
name when that file is absent.
"""
function discover_runs(save_path::AbstractString)
    isdir(save_path) || return RunRecord[]
    out = RunRecord[]
    for f in sort(readdir(save_path))
        dir = joinpath(save_path, f)
        (isdir(dir) && occursin(r"^run_\d+$", f)) || continue
        rank, seed = 0, -1
        info = joinpath(dir, "run_info.json")
        if isfile(info)
            txt = read(info, String)
            m = match(r"\"run\"\s*:\s*(\d+)", txt);  m === nothing || (rank = parse(Int, m[1]))
            m = match(r"\"seed\"\s*:\s*(-?\d+)", txt); m === nothing || (seed = parse(Int, m[1]))
        end
        rank == 0 && (rank = parse(Int, match(r"(\d+)$", f)[1]))
        push!(out, RunRecord(dir, f, rank, seed, isfile(joinpath(dir, "top_motifs.csv"))))
    end
    sort!(out, by = r -> r.rank)
    return out
end

"""
    choose_runs(combiner, runs) -> Vector{RunRecord}

Which runs the display is built from. The only method a new combiner must add.
"""
function choose_runs(c::UseSingleRun, runs::Vector{RunRecord})
    isempty(runs) && return RunRecord[]
    i = findfirst(r -> r.rank == c.index, runs)
    i === nothing && throw(ArgumentError(
        "UseSingleRun($(c.index)): no run with rank $(c.index); found ranks $(getfield.(runs, :rank))"))
    return [runs[i]]
end

"""
    combine_runs(save_path; combiner=UseSingleRun()) -> NamedTuple

Build the top-level view of a multi-run result folder. Returns
`(; runs, chosen, summary, index)`.

What it writes into `save_path`:

  * `top_motifs.csv` — a copy of the chosen run's, so any existing
    `run_N/*/top_motifs.csv` glob keeps working unchanged;
  * `runs.html` — a small index naming every run, its seed, and which one is
    displayed, linking into each run's own pages;
  * `runs.json` — the same facts, for programs.

Nothing inside the run folders is modified, so re-running with a different
combiner is free and reversible. A folder with no `run_*` subfolders is a
single-run result and is left completely alone.
"""
function combine_runs(save_path::AbstractString; combiner::RunCombiner=UseSingleRun())
    runs = discover_runs(save_path)
    if isempty(runs)
        @debug "combine_runs: no run_* folders under $save_path; single-run layout, nothing to do"
        return (; runs, chosen=RunRecord[], summary=nothing, index=nothing)
    end
    chosen = choose_runs(combiner, runs)
    policy = string(nameof(typeof(combiner)))

    summary = nothing
    src = joinpath(first(chosen).dir, "top_motifs.csv")
    if isfile(src)
        summary = joinpath(save_path, "top_motifs.csv")
        cp(src, summary; force=true)
    else
        @warn "combine_runs: chosen run has no top_motifs.csv" run=first(chosen).name
    end

    open(joinpath(save_path, "runs.json"), "w") do io
        print(io, "{\"policy\":\"", policy, "\",\"chosen\":[",
              join(("\"$(r.name)\"" for r in chosen), ","), "],\"runs\":[")
        print(io, join(("{\"name\":\"$(r.name)\",\"rank\":$(r.rank),\"seed\":$(r.seed)," *
                        "\"has_summary\":$(r.has_summary)}" for r in runs), ","))
        print(io, "]}")
    end

    index = joinpath(save_path, "runs.html")
    open(index, "w") do io
        print(io, "<!DOCTYPE html><html><head><meta charset=\"UTF-8\">",
              "<title>Runs</title><link rel=\"stylesheet\" href=\"styles.css\">",
              "<style>body{font-family:Helvetica,Arial,sans-serif;margin:2rem auto;max-width:52rem}",
              "table{border-collapse:collapse;width:100%}td,th{padding:.4rem .7rem;",
              "border-bottom:1px solid #ddd;text-align:left}.on{font-weight:700}",
              "</style></head><body><h1>Runs</h1><p>Policy: <code>", policy,
              "</code>. The displayed result comes from the marked run; the others are ",
              "kept for comparison and are complete renderings in their own right.</p>",
              "<table><tr><th>run</th><th>rank</th><th>seed</th><th>pages</th><th></th></tr>")
        for r in runs
            on = any(c -> c.name == r.name, chosen)
            print(io, "<tr", on ? " class=\"on\"" : "", "><td>", r.name, "</td><td>", r.rank,
                  "</td><td>", r.seed < 0 ? "?" : string(r.seed), "</td><td>")
            pages = filter(x -> startswith(x, "renderings"), readdir(r.dir))
            print(io, isempty(pages) ? "&mdash;" :
                  join(("<a href=\"$(r.name)/$(p)/index.html\">$(p)</a>" for p in pages), " "))
            print(io, "</td><td>", on ? "displayed" : "", "</td></tr>")
        end
        print(io, "</table></body></html>")
    end
    @info "combine_runs" policy chosen=[r.name for r in chosen] of=length(runs) index
    return (; runs, chosen, summary, index)
end

export RunCombiner, UseSingleRun, RunRecord, discover_runs, choose_runs, combine_runs
