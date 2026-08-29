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



# ─────────────────────────────────────────────────────────────────────────────
# Consensus top movers.
#
# One page above the runs, showing the findings the runs AGREE on. Each row is a
# real motif from a real model -- nothing is merged or synthesised; the runs
# only decide which motif is worth showing and how confident to be about it.
#
# Two steps, in this order. Scoring uses cross-run agreement, so collapsing
# duplicates first would destroy the evidence it needs.
#
#   1. score    supp(g) = mean over the other runs of the best matching motif,
#               where matching needs the same direction, overlapping windows and
#               shared carriers.
#   2. collapse greedily emit the highest-support motif and drop everything whose
#               window overlaps it by >= rho.
#
# Full specification and the reasoning for each condition: docs/multirun.md and
# multirun_findings/12_method.md.
#
# CONVOLUTION: nothing here is mutagenesis-specific except where the run folders
# come from. `plot_motifs_conv_case` writes the same top_motifs.json, so the conv
# path needs only that dump wiring in to use this unchanged.
# ─────────────────────────────────────────────────────────────────────────────

"The consensus page has two destinations and no numbered pages of its own."
const NAV_CONSENSUS = "[[\"index.html\",\"Top movers\"],[\"index4.html\",\"Readme\"]]"

"One candidate motif, lifted out of a finished run folder."
struct RunMotif
    entry::TopMoverEntry
    run::String            # "run_3"
    dir::Int               # +1 / -1
    win::UnitRange{Int}    # filter window, original coordinates
    carriers::Set{Int}
    carriers_path::String  # run-relative path of the carrier CSV, for tracing
end

"""
    ConsensusTopMovers(; n=5, rho=0.25, min_runs=1)

Select the top `n` increasing and `n` decreasing findings across all runs.

`rho` is the window-overlap threshold at which two motifs are considered the
same finding: `rho = 0.25` keeps features 3 nt apart distinct while collapsing
near-identical windows. `min_runs` drops findings supported by fewer than that
many runs; the default keeps everything, since support is meant to be reported
rather than used as a gate.
"""
struct ConsensusTopMovers <: RunCombiner
    n::Int
    rho::Float64
    min_runs::Int
end
ConsensusTopMovers(; n::Int=5, rho::Real=0.25, min_runs::Int=1) =
    ConsensusTopMovers(n, Float64(rho), min_runs)

choose_runs(::ConsensusTopMovers, runs::Vector{RunRecord}) = runs   # all of them

"Window of a carrier CSV: every row shares one, so the first row is enough."
function _carrier_window(path::AbstractString)
    isfile(path) || return nothing
    ls = readlines(path); length(ls) < 2 && return nothing
    h = split(ls[1], ','); si = findfirst(==("start_position"), h); ei = findfirst(==("end_position"), h)
    (si === nothing || ei === nothing) && return nothing
    lo, hi = typemax(Int), typemin(Int); ids = Set{Int}()
    xi = findfirst(==("seq_index"), h)
    for l in ls[2:end]
        p = split(l, ','); length(p) < max(si, ei, something(xi, 0)) && continue
        a = tryparse(Int, p[si]); b = tryparse(Int, p[ei])
        (a === nothing || b === nothing) && continue
        lo = min(lo, a); hi = max(hi, b)
        xi === nothing || (v = tryparse(Int, p[xi]); v === nothing || push!(ids, v))
    end
    lo > hi ? nothing : (lo:hi, ids)
end

"Load every run's top movers as `RunMotif`s, with asset paths rewritten."
function collect_run_motifs(save_path::AbstractString, runs::Vector{RunRecord})
    out = RunMotif[]
    for rec in runs
        js = joinpath(rec.dir, "top_motifs.json")
        isfile(js) || (@warn "no top_motifs.json; run skipped" run=rec.name; continue)
        # assets are stored relative to the RENDERING folder, so find it
        rends = filter(x -> startswith(x, "renderings") && isdir(joinpath(rec.dir, x)),
                       readdir(rec.dir))
        isempty(rends) && continue
        prefix = joinpath(rec.name, first(rends))
        pos, neg, car = read_top_movers_json(js; prefix=prefix)
        for (i, entry) in enumerate(vcat(pos, neg))
            dir  = i <= length(pos) ? 1 : -1
            info = _carrier_window(joinpath(save_path, car[i]))
            info === nothing && continue
            push!(out, RunMotif(entry, rec.name, dir, info[1], info[2], car[i]))
        end
    end
    out
end

# ---- the two steps ----------------------------------------------------------

_ov(a::Set{Int}, b::Set{Int}) =
    (isempty(a) || isempty(b)) ? 0.0 : length(intersect(a, b)) / min(length(a), length(b))

"Pairwise agreement: same direction, overlapping windows, shared carriers."
function motif_omega(g::RunMotif, h::RunMotif)
    g.dir == h.dir || return 0.0
    isempty(intersect(g.win, h.win)) && return 0.0
    _ov(g.carriers, h.carriers)
end

"Fraction of the window two motifs share; 1 = identical, 0 = disjoint."
function window_rho(g::RunMotif, h::RunMotif)
    w = min(length(g.win), length(h.win)); w == 0 && return 0.0
    length(intersect(g.win, h.win)) / w
end

"""
    motif_support(ms) -> Vector{Float64}

Step 1. For each motif, the mean over the OTHER runs of its best match there. A
run with no matching candidate contributes `0`, not a skipped term — skipping
would reward a motif for being unmatched.
"""
function motif_support(ms::Vector{RunMotif})
    runs = unique(m.run for m in ms)
    length(runs) < 2 && return fill(NaN, length(ms))
    out = Vector{Float64}(undef, length(ms))
    for (i, g) in enumerate(ms)
        tot = 0.0; cnt = 0
        for r in runs
            r == g.run && continue
            best = 0.0
            for h in ms
                h.run == r || continue
                o = motif_omega(g, h)
                o > best && (best = o)
            end
            tot += best; cnt += 1
        end
        out[i] = cnt == 0 ? NaN : tot / cnt
    end
    return out
end

"""
    dedup_motifs(ms, supp; rho=0.25, n=5) -> (positives, negatives)

Step 2. Greedy: emit the highest-support motif, drop everything of the same
direction whose window overlaps it by at least `rho`, repeat. Returns up to `n`
findings per direction as `(motif, support, n_members, n_runs)` tuples.

Greedy rather than connected components: single linkage chains one window to the
next and fuses a whole region into a single finding.
"""
function dedup_motifs(ms::Vector{RunMotif}, supp::Vector{Float64};
                      rho::Real=0.25, n::Int=5, min_runs::Int=1)
    out = Dict(1 => NamedTuple[], -1 => NamedTuple[])
    for dir in (1, -1)
        idx = [i for i in eachindex(ms) if ms[i].dir == dir && !isnan(supp[i])]
        left = Set(idx)
        while !isempty(left) && length(out[dir]) < n
            i = argmax(j -> supp[j], collect(left))
            members = [j for j in left if window_rho(ms[i], ms[j]) >= rho]
            nruns = length(unique(ms[j].run for j in members))
            nruns >= min_runs &&
                push!(out[dir], (motif=ms[i], support=supp[i],
                                 members=length(members), runs=nruns,
                                 region=minimum(minimum(ms[j].win) for j in members):
                                        maximum(maximum(ms[j].win) for j in members)))
            setdiff!(left, members)
        end
    end
    out[1], out[-1]
end

# ---- rendering --------------------------------------------------------------

"""
    _consensus_note(f, n_runs)

The one-line provenance note under each consensus row.

Reports how many runs FOUND the region, as a fraction of all runs, rather than
how many motifs landed in the cluster. The motif count is dominated by
within-run duplication -- one run routinely describes the same region several
times -- so it measures filter redundancy, not reproducibility, and reads as
though it were the latter.
"""
function _consensus_note(f, n_runs::Int)
    string("<span class=\"consensus-note\">support ", round(f.support, digits=2),
           " &middot; found in ", f.runs, "/", n_runs, " runs",
           " &middot; region ", first(f.region), ":", last(f.region),
           " &middot; shown from <code>", f.motif.run, "</code></span>")
end

function _runs_nav(runs::Vector{RunRecord})
    isempty(runs) && return ""
    items = map(runs) do r
        rends = filter(x -> startswith(x, "renderings") && isdir(joinpath(r.dir, x)), readdir(r.dir))
        tgt = isempty(rends) ? "" : joinpath(r.name, first(rends), "index.html")
        isempty(tgt) ? "" : string("<a href=\"", tgt, "\">", r.name,
                                   r.seed < 0 ? "" : " (seed $(r.seed))", "</a>")
    end
    string("<div class=\"runs-nav\"><span class=\"runs-nav-label\">Individual runs</span>",
           join(filter(!isempty, items), " "), "</div>")
end

"""
    consensus_top_movers(save_path; combiner=ConsensusTopMovers(), page_title="n/a", ...)

Build the top-level consensus page from the `run_*` folders under `save_path`.

Writes `index.html` (the consensus top movers) and `index1.html` (a readme
placeholder) at `save_path`, plus `consensus.json` recording what was selected
and why. Nothing inside the run folders is touched, so re-running with different
settings is free.

Returns `(; positives, negatives, runs)` where each finding carries its
representative motif, support, member count and run count.
"""
function consensus_top_movers(save_path::AbstractString;
        combiner::ConsensusTopMovers = ConsensusTopMovers(),
        page_title::AbstractString = "n/a",
        protein_name = nothing,
        protein_length = nothing,
        wild_type = nothing,
        feature_label = nothing)

    runs = discover_runs(save_path)
    length(runs) < 2 && throw(ArgumentError(
        "consensus_top_movers needs at least 2 runs under $save_path; found $(length(runs))"))

    ms = collect_run_motifs(save_path, runs)
    isempty(ms) && throw(ArgumentError("no motifs could be read from the run folders"))
    sup = motif_support(ms)
    pos, neg = dedup_motifs(ms, sup; rho=combiner.rho, n=combiner.n, min_runs=combiner.min_runs)

    # The rows are real entries from real runs; only the provenance line is new.
    note(f) = _consensus_note(f, length(runs))
    # Annotated, not inferred. A dataset can easily have findings in one
    # direction only -- Pitt_2010_ribozyme yields negatives and no positives --
    # and an empty comprehension gives `Vector{Any}`, which the typed
    # `positives::AbstractVector{TopMoverEntry}` keyword below rejects with a
    # TypeError. Naming the element type makes the empty case a valid empty
    # vector rather than a crash.
    pe = TopMoverEntry[f.motif.entry for f in pos]
    ne = TopMoverEntry[f.motif.entry for f in neg]
    for (fs, es) in ((pos, pe), (neg, ne))
        for (f, e) in zip(fs, es)
            isempty(e.texts) || push!(e.texts, note(f))
        end
    end

    # The page links styles.css relative to itself, and only the RENDERING folders
    # get one. Without this the consensus page loads with no stylesheet at all and
    # is unreadable. Same build artifact as every other render, written the same way.
    open(joinpath(save_path, "styles.css"), "w") do io
        print(io, Mustache.render(template_css))
    end

    render_top_movers_page!(save_path;
        positives = pe, negatives = ne,
        page_title = page_title, nav_page_count = 1,
        protein_name = protein_name, protein_length = protein_length,
        wild_type = wild_type, feature_label = feature_label,
        show_epistasis = false, modal_scroll_fix = true,
        consensus_rows = [note(f) for f in vcat(pos, neg)],
        runs_nav = _runs_nav(runs),
        nav_override = NAV_CONSENSUS)

    render_readme_page!(save_path; page_title=page_title, has_summary=true,
                        nav_override=NAV_CONSENSUS)

    write_consensus_csv(joinpath(save_path, "consensus_top_motifs.csv"), pos, neg;
                        protein_name=protein_name, label=feature_label)

    open(joinpath(save_path, "consensus.json"), "w") do io
        JSON3.write(io, Dict("rho"=>combiner.rho, "n"=>combiner.n,
            "min_runs"=>combiner.min_runs, "runs"=>[r.name for r in runs],
            "findings"=>[Dict("direction"=>d, "support"=>f.support, "members"=>f.members,
                              "runs"=>f.runs, "region"=>[first(f.region), last(f.region)],
                              "from"=>f.motif.run, "span"=>f.motif.entry.span)
                         for (d, fs) in (("positive", pos), ("negative", neg)) for f in fs]))
    end
    @info "consensus_top_movers" runs=length(runs) motifs=length(ms) positives=length(pos) negatives=length(neg)
    return (; positives=pos, negatives=neg, runs)
end

export RunCombiner, UseSingleRun, ConsensusTopMovers, RunRecord, RunMotif,
       discover_runs, choose_runs, combine_runs, consensus_top_movers,
       collect_run_motifs, motif_support, dedup_motifs, motif_omega, window_rho,
       write_consensus_csv

"""
    write_consensus_csv(path, positives, negatives; protein_name=nothing, label=nothing)

The consensus selection as one flat table, so every displayed row can be traced
back to the run, the logo file and the carrier list it came from.

Columns, in order:

    protein_name label direction rank support n_runs n_motifs region_start
    region_end from_run display_name group_label span is_singleton count
    shapley_median cluster_median cluster_nnd nnd_pvalue location_z wt_seq
    mut_positions description logo_path carriers_csv

`logo_path` and `carriers_csv` are relative to the folder holding this file, so
they resolve directly: `run_3/renderings_1/mutation_regions_1/12_15:22.png`.
"""
function write_consensus_csv(path::AbstractString, positives, negatives;
                             protein_name=nothing, label=nothing)
    mkpath(dirname(path))
    esc(x) = (t = string(x); occursin(r"[\",\n]", t) ? "\"" * replace(t, "\"" => "\"\"") * "\"" : t)
    open(path, "w") do io
        println(io, join(["protein_name","label","direction","rank","support","n_runs",
            "n_motifs","region_start","region_end","from_run","display_name","group_label",
            "span","is_singleton","count","shapley_median","cluster_median","cluster_nnd",
            "nnd_pvalue","location_z","wt_seq","mut_positions","description",
            "logo_path","carriers_csv"], ","))
        for (dir, fs) in (("positive", positives), ("negative", negatives))
            for (rank, f) in enumerate(fs)
                e = f.motif.entry
                mut = Int[]
                for reg in e.wt_regions, (j, flag) in enumerate(reg.mutated)
                    flag && push!(mut, reg.start + j - 1)
                end
                println(io, join(esc.([
                    something(protein_name, ""), something(label, ""), dir, rank,
                    round(f.support, digits=4), f.runs, f.members,
                    first(f.region), last(f.region), f.motif.run,
                    e.display_name, e.group_label, e.span, e.is_singleton, e.count,
                    e.median, e.cluster_median, e.cluster_nnd, e.nnd_p, e.location_z,
                    join((r.wt for r in e.wt_regions), ";"), join(mut, ";"),
                    mutation_description(e.wt_regions, e.median, label),
                    e.img, f.motif.carriers_path]), ","))
            end
        end
    end
    return path
end
