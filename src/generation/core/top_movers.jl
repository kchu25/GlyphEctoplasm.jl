"""
Top-movers summary page.

A provider-agnostic landing page that surfaces the strongest positive and
negative motifs as a row-by-row list: a clickable card on the left (default card
styling) and a metadata panel on the right (Shapley median, cluster median,
cluster NND, mutation position, count, group). Clicking a card opens a
self-contained detail modal embedded in the page — it does not depend on the
grouped page's index conventions, so the same `TopMoverEntry` + render API works
for both the mutagenesis and (later) the convolution case.
"""

"""
    WildTypeRegion(wt, start, mutated, obs)

One region's wild-type amino-acid window for the summary card: the full WT string
`wt` (un-abbreviated), the 1-based `start` position of its first residue, a
per-column `mutated` flag (true where the motif's observed consensus differs from
the wild type), and `obs`, the motif's own observed consensus over the same
window. The panel draws one up-arrow under each mutated residue.

`obs` is what the motif mutates *to*: column `j` sits at sequence position
`start + j - 1`, carries wild-type residue `wt[j]`, and the motif substitutes
`obs[j]` there. Together with `mutated` that is everything
[`mutation_description`](@ref) needs to state the substitution in words.
"""
struct WildTypeRegion
    wt::String
    start::Int
    mutated::Vector{Bool}
    obs::String
end

"""
    region_window_span(regions) -> String

The span of sequence the motif's FILTER windows cover, e.g. `"(17:24, 36:43)"`,
derived from each region's start and window length.

This is not the same thing as `TopMoverEntry.span`, which lists the individual
positions that carry mutations (`"(17:18, 20:23, 36:40, 42)"`). The region-view
logo draws whole windows, so labelling it with scattered mutation positions
mislabels what the picture shows. One region renders bare (`"17:24"`); two or
more are parenthesised, matching how `span` is written.

Returns `""` for an empty vector, which callers use to fall back to `span`.
"""
function region_window_span(regions::AbstractVector{WildTypeRegion})
    isempty(regions) && return ""
    wins = [string(r.start, ":", r.start + length(r.wt) - 1) for r in regions]
    return length(wins) == 1 ? wins[1] : string("(", join(wins, ", "), ")")
end

"""
    mutation_description(regions, median, label) -> String

Plain-language summary of what a mutation motif does, e.g.

    Mutations at sites 66, 69, 72 to K, L, P respectively increase ddG_ML_float (kcal/mol).

Sites are the arrowed (mutated) columns of `regions`, in ascending position
order across every region — so a multi-region motif contributes all of its
sites, not just the first region's. Each site's residue is the motif's observed
consensus there (`WildTypeRegion.obs`), i.e. what the position mutates *to*.

Direction comes from the sign of `median` (the motif's Shapley median): positive
reads "increase", negative "decrease". `label` is the feature/assay name shown
verbatim, units included; when it is `nothing` or blank the sentence falls back
to "the measured value".

Returns `""` when there is nothing to describe — no regions, no arrowed columns,
or a zero/NaN median — so callers can treat empty as "omit the section".

Grammar is singular for a lone site ("Mutation at site 66 to K increases …") and
drops the "respectively", which only makes sense for a list.

!!! note "Where the residues come from"
    The observed residue is `argmax` over the motif's count-matrix column, mapped
    through the alphabet the sequences were encoded with — `AMINO_ACID_LETTERS`
    (20 rows, alphabetical) for protein, A/C/G/T(U) for nucleotide. It is the
    motif's *modal* residue, not a per-sequence call, so a motif whose column is
    split across several residues is summarised by its most common one.
"""
function mutation_description(regions::AbstractVector{WildTypeRegion},
        median::Real, label::Union{AbstractString,Nothing})
    isempty(regions) && return ""
    (isnan(median) || median == 0) && return ""

    sites = Int[]
    residues = Char[]
    for reg in regions
        for j in 1:min(length(reg.mutated), length(reg.obs))
            reg.mutated[j] || continue
            ch = reg.obs[j]
            ch == _placeholder_char_ && continue   # no confident call at this column
            push!(sites, reg.start + j - 1)
            push!(residues, ch)
        end
    end
    isempty(sites) && return ""

    # Ascending by position, so a multi-region motif reads left to right along
    # the sequence rather than in region-construction order.
    ord = sortperm(sites)
    sites, residues = sites[ord], residues[ord]

    what = (label === nothing || isempty(strip(String(label)))) ?
           "the measured value" : strip(String(label))
    verb_up = median > 0

    # The magnitude is the motif's median Banzhaf value. That number is already
    # in the assay's own units: BanzhafInference composes the scale-back function
    # into the computation, so a run trained on log or z-scored labels still
    # reports influences on the original scale.
    mag = _format_magnitude(abs(float(median)))

    if length(sites) == 1
        return string("Mutation at site ", sites[1], " to ", residues[1], " ",
                      verb_up ? "increases" : "decreases", " ", what,
                      " by about ", mag, ".")
    end
    # Positions collapse to ranges and the residues collapse the SAME way, so
    # the two lists stay token-for-token aligned:
    #   "7-14, 36, 38-41 to GGUGUUGC, U, CGAC respectively"
    runs = integer_runs(sites)
    pos_toks = _run_position_tokens(sites, runs)
    res_toks = _run_residue_tokens(residues, runs)
    # "respectively" only earns its place when there are two or more tokens to
    # pair up. One run collapses to a single token on each side ("sites 7-14 to
    # CUCUCGUC"), where the word is noise.
    resp = length(pos_toks) > 1 ? " respectively " : " "
    return string("Mutations at sites ", join(pos_toks, ", "),
                  " to ", join(res_toks, ", "),
                  resp, verb_up ? "increase" : "decrease", " ", what,
                  " by about ", mag, ".")
end

"""
    _format_magnitude(x) -> String

`x` rendered with 3 significant digits, and without an exponent for values a
reader can read directly. Assay units span many orders of magnitude -- ddG is
around 1, a ribozyme rate can be 1e-3 -- so a fixed number of decimal places
would print `0.00` for a real effect.
"""
function _format_magnitude(x::Real)
    (isnan(x) || !isfinite(x)) && return "?"
    x == 0 && return "0"
    ax = abs(x)
    if ax >= 0.01 && ax < 1e5
        digits = ax >= 10 ? 1 : (ax >= 1 ? 2 : 4)
        return string(round(x; digits=digits))
    end
    return string(round(x; sigdigits=3))
end

"""
    integer_runs(xs) -> Vector{UnitRange{Int}}

Index ranges into `xs` covering each maximal run of consecutive integers.
`[11,12,13,21,23,24]` gives `[1:3, 4:4, 5:6]`.

Positions and their residues are collapsed against the same run list, which is
what keeps the two halves of the interpretation sentence aligned. Input is
assumed sorted ascending.
"""
function integer_runs(xs::AbstractVector{<:Integer})
    isempty(xs) && return UnitRange{Int}[]
    runs = UnitRange{Int}[]
    i, n = 1, length(xs)
    while i <= n
        j = i
        while j < n && xs[j + 1] == xs[j] + 1
            j += 1
        end
        push!(runs, i:j)
        i = j + 1
    end
    return runs
end

# A run is worth collapsing only at length 3+. "4-5" is no shorter than "4, 5"
# and reads worse, and the residue side would turn a two-letter run into a
# two-letter token that looks like a single substitution.
_collapsible(r) = length(r) >= 3

"""
    _run_position_tokens(sites, runs) -> Vector{String}

One token per run: `"11-17"` for a collapsible run, otherwise one token per
position.
"""
function _run_position_tokens(sites, runs)
    toks = String[]
    for r in runs
        if _collapsible(r)
            push!(toks, string(sites[first(r)], "-", sites[last(r)]))
        else
            for k in r
                push!(toks, string(sites[k]))
            end
        end
    end
    return toks
end

"""
    _run_residue_tokens(residues, runs) -> Vector{String}

The residue counterpart of [`_run_position_tokens`](@ref): a collapsible run
becomes one concatenated string (`"GGUGUUGC"`), so it lines up with the `"7-14"`
on the position side. Non-collapsible runs contribute one token per residue.
"""
function _run_residue_tokens(residues, runs)
    toks = String[]
    for r in runs
        if _collapsible(r)
            push!(toks, String(residues[r]))
        else
            for k in r
                push!(toks, string(residues[k]))
            end
        end
    end
    return toks
end

"""
    collapse_runs(xs) -> String

A sorted integer vector with consecutive runs written as ranges:
`[11,12,13,14,15,16,17,21,23,24,25,26,27]` becomes `"11-17, 21, 23-27"`.
Runs shorter than 3 stay listed individually.
"""
collapse_runs(xs::AbstractVector{<:Integer}) =
    join(_run_position_tokens(xs, integer_runs(xs)), ", ")

"""
    TopMoverEntry

One summary row. Carries the ranking keys, the metadata shown in the right-hand
panel, and a self-contained modal payload (image paths + text fields) so the
detail modal can be populated without the main page's `jsonData`.

Path fields (`img`, `influence`, `yy_kde`) are relative to the page's save
folder, matching the relative paths used by the grouped page.

!!! note "`location_z` is the last field, on purpose"
    It was added after the other 20 and is appended at the end so that every
    pre-existing 20-argument positional construction still works — a
    20-argument call is forwarded to the 21-argument one with `location_z =
    NaN`. Put it anywhere else and silent argument shifting becomes possible.
"""
struct TopMoverEntry
    # Ranking keys (group-less binned lexicographic order)
    median::Float64          # Shapley median (sign decides positive/negative column)
    cluster_median::Float64  # median expression of the motif's sequences
    cluster_nnd::Float64     # observed mean k-NN distance (tighter = smaller)
    count::Int
    # Right-panel metadata
    display_name::String
    group_label::String
    span::String             # mutation position(s), e.g. "36:45" or "36:45, 50:60"
    nnd_p::Float64           # NND permutation p-value (NaN when unavailable)
    is_singleton::Bool
    epistasis::String        # interaction-coefficient summary (caller HTML; "" if none)
    # Wild-type region(s): full WT string + per-residue mutation flags (drawn as
    # an amino-acid track with up-arrows under each mutated residue).
    wt_regions::Vector{WildTypeRegion}
    # Self-contained modal payload (paths relative to save folder)
    img::String              # motif logo (grouped-page default view)
    # Mutagenesis dual-view logos for the top-movers page: `img_reduced` is the
    # reduced view (only the mutated fragments that differ from the backbone —
    # the basis of region interactions), `img_region` the full-region view. Both
    # `""` ⇒ single-card layout (convolution case).
    img_reduced::String
    img_region::String
    influence::String        # Shapley/influence plot
    yy_kde::String           # indicator (yy-KDE) plot
    texts::Vector{String}    # metadata text fields (rendered as HTML in the modal)
    # Multi-motif distance variants (empty for singletons). When non-empty the
    # detail popup shows the multi-motif modal with an inter-motif-distance slider
    # instead of the static singleton modal, mirroring the grouped Motifs page.
    variant_pwms::Vector{String}            # per-variant logo paths (one per distance)
    variant_labels::Vector{String}          # per-variant distance labels
    variant_texts::Vector{Vector{String}}   # per-variant metadata text fields
    # Location statistic: (mean carrier label − mean population label) /
    # (sd(population) / sqrt(n_carriers)). NaN when unavailable — see
    # `location_z_1d`. Trailing field so 20-argument calls keep working.
    location_z::Float64
end

# Legacy 20-argument constructor: everything that built a `TopMoverEntry` before
# `location_z` existed keeps compiling and keeps producing the same entry, with
# the new statistic marked unavailable (NaN).
TopMoverEntry(median, cluster_median, cluster_nnd, count, display_name, group_label,
              span, nnd_p, is_singleton, epistasis, wt_regions, img, img_reduced,
              img_region, influence, yy_kde, texts, variant_pwms, variant_labels,
              variant_texts) =
    TopMoverEntry(median, cluster_median, cluster_nnd, count, display_name, group_label,
                  span, nnd_p, is_singleton, epistasis, wt_regions, img, img_reduced,
                  img_region, influence, yy_kde, texts, variant_pwms, variant_labels,
                  variant_texts, NaN)

"""
    select_top_movers(entries; n=5, bin_count=10, min_count=1) -> (positives, negatives)

Split `entries` by the sign of `median` and return the top `n` of each, ordered
by the group-less binned lexicographic key — the same key the grouped page sorts
on, minus the region-group component. The two columns are **mirror images** on
the sign-bearing keys so each leads with its strongest signal:

  positives: cluster-median bin (high → low) → cluster-NND bin (low → high) →
             Shapley-median bin (high → low) → count (high → low)
  negatives: cluster-median bin (low → high) → cluster-NND bin (low → high) →
             Shapley-median bin (low → high) → count (high → low)

So the most positive cluster median sits at the top of the positive column and
the most negative cluster median at the top of the negative column. The
NND-tightness (tighter first) and count (higher first) tiebreakers stay the same
in both. Continuous keys are coarsened into `bin_count` equal-width bins over the
global range of `entries` (mirrors the display sort), so near-equal motifs aren't
split on meaningless differences. Motifs with `median == 0` are dropped.

`min_count` gates the summary on cluster size: motifs whose cluster holds fewer
than `min_count` points are dropped before the sign split. Tiny clusters make the
lead ranking keys degenerate — `cluster_median` is a one- or two-point "median"
and `cluster_nnd` is `NaN`/a single pairwise distance below ~3 points, with no
permutation-test power — so a single lucky sequence could otherwise top a column.
The default `1` keeps every motif (no gating); the runners pass a higher floor.
Note this is cluster size (`count`), not the `is_singleton` flag (which marks a
single-region motif vs. a multi-motif combination — a different axis).
"""
function select_top_movers(entries::AbstractVector{TopMoverEntry}; n::Int=5, bin_count::Int=10, min_count::Int=1)
    isempty(entries) && return (TopMoverEntry[], TopMoverEntry[])

    if min_count > 1
        entries = [e for e in entries if e.count >= min_count]
        isempty(entries) && return (TopMoverEntry[], TopMoverEntry[])
    end

    clus = [e.cluster_median for e in entries if !isnan(e.cluster_median)]
    nnd  = [e.cluster_nnd    for e in entries if !isnan(e.cluster_nnd)]
    shap = [e.median         for e in entries if !isnan(e.median)]
    clus_lo, clus_hi = isempty(clus) ? (0.0, 0.0) : extrema(clus)
    nnd_lo,  nnd_hi  = isempty(nnd)  ? (0.0, 0.0) : extrema(nnd)
    shap_lo, shap_hi = isempty(shap) ? (0.0, 0.0) : extrema(shap)

    # Positives: strongest (highest) cluster median / Shapley first.
    pos_key(e) = (
        -bin_value(e.cluster_median, clus_lo, clus_hi, bin_count),   # high cluster median first
        bin_value_asc(e.cluster_nnd, nnd_lo, nnd_hi, bin_count),     # low NND (tighter) first
        -bin_value(e.median, shap_lo, shap_hi, bin_count),           # high Shapley first
        -e.count,                                                    # higher count first
    )
    # Negatives: mirror — most negative (lowest) cluster median / Shapley first.
    # `bin_value_asc` puts low bins first and keeps NaN last.
    neg_key(e) = (
        bin_value_asc(e.cluster_median, clus_lo, clus_hi, bin_count), # low cluster median first
        bin_value_asc(e.cluster_nnd, nnd_lo, nnd_hi, bin_count),      # low NND (tighter) first
        bin_value_asc(e.median, shap_lo, shap_hi, bin_count),         # low Shapley first
        -e.count,                                                     # higher count first
    )
    positives = sort([e for e in entries if e.median > 0], by = pos_key)
    negatives = sort([e for e in entries if e.median < 0], by = neg_key)
    return (first(positives, n), first(negatives, n))
end

# Minimal HTML escaping for text injected server-side into the summary rows.
_tm_esc(s) = replace(string(s), "&" => "&amp;", "<" => "&lt;", ">" => "&gt;")

"""
    parse_epistasis(s) -> (coef, pval)

Pull the interaction coefficient and p-value out of a caller-built epistasis
string such as `"β_interaction: +0.64, se: 0.0172, p-value: 1.60e-268"`. HTML
tags are stripped first. Returns `(coef::String|nothing, pval::String|nothing)`;
the coefficient is returned with an explicit leading sign. When the coefficient
can't be found the caller should fall back to showing the original string.
"""
function parse_epistasis(s::AbstractString)
    txt = replace(string(s), r"<[^>]*>" => "")
    num = raw"([+-]?\d*\.?\d+(?:[eE][+-]?\d+)?)"
    cm = match(Regex("interaction\\s*:?\\s*" * num), txt)
    pm = match(Regex("p[-\\s]?value\\s*:?\\s*" * num), txt)
    coef = cm === nothing ? nothing : String(cm.captures[1])
    pval = pm === nothing ? nothing : String(pm.captures[1])
    if coef !== nothing && !startswith(coef, "+") && !startswith(coef, "-")
        coef = "+" * coef
    end
    return (coef, pval)
end

"""
    wt_track_html(regions) -> String

Render the wild-type amino-acid track(s): one row per region — a small span label
plus a monospace track of residues, each mutated residue marked by an up-arrow and
its position number underneath. Columns are stacked (residue / arrow / position) so
they stay aligned without relying on glyph widths. Returns `""` when empty.
"""
function wt_track_html(regions::AbstractVector{WildTypeRegion})
    isempty(regions) && return ""
    region_html = String[]
    for reg in regions
        cols = String[]
        for (j, ch) in enumerate(reg.wt)
            mut = j <= length(reg.mutated) && reg.mutated[j]
            cls    = mut ? "wt-col is-mut" : "wt-col"
            arrow  = mut ? "▲" : ""
            posnum = mut ? string(reg.start + j - 1) : ""
            push!(cols, string(
                "<div class=\"", cls, "\">",
                "<span class=\"wt-aa\">", _tm_esc(string(ch)), "</span>",
                "<span class=\"wt-arrow\">", arrow, "</span>",
                "<span class=\"wt-pos\">", posnum, "</span></div>"))
        end
        span_lbl = "$(reg.start)–$(reg.start + length(reg.wt) - 1)"
        push!(region_html, string(
            "<div class=\"wt-region\"><div class=\"wt-span\">", span_lbl,
            "</div><div class=\"wt-track\">", join(cols), "</div></div>"))
    end
    return "<div class=\"wt-block\">" * join(region_html, "\n") * "</div>"
end

"""
    top_mover_row_html(entry, rowid) -> String

Render one summary row: a clickable card (left) wired to `openTopMover(rowid)`
and a metadata table (right). `rowid` is the 0-based index into the page's
`topMoverData` JS array.
"""
function top_mover_row_html(e::TopMoverEntry, rowid::Int; show_epistasis::Bool=true)
    name = _tm_esc(e.display_name)
    # Mutagenesis rows: the name is a mutation-position span — prefix it with a
    # muted "position(s)" label so the bare coordinates read clearly. (Conv rows
    # have an empty span and keep the name as-is.)
    name_html = if isempty(strip(e.span))
        "<div class=\"top-mover-name\">$name</div>"
    else
        word = occursin(',', e.span) ? "positions" : "position"
        "<div class=\"top-mover-name\"><span class=\"top-mover-name-label\">$word</span>$name</div>"
    end
    epi_block = if !show_epistasis
        ""                                          # caller opted out (e.g. convolution case)
    else
        epi = if isempty(strip(e.epistasis))
            "<div class=\"epi-line epi-none\">region interaction: n/a</div>"
        else
            coef, pval = parse_epistasis(e.epistasis)
            if coef === nothing
                "<div class=\"epi-line\">$(e.epistasis)</div>"   # unknown format: show as-is
            else
                pv = pval === nothing ? "" : "<div class=\"epi-line epi-sub\">p-value: $(pval)</div>"
                "<div class=\"epi-line\">region interaction: <strong>$(coef)</strong></div>$pv"
            end
        end
        "\n            <div class=\"top-mover-epistasis\">$epi</div>"
    end
    # Dual-view (mutagenesis): two logo cards — reduced view then region view.
    # Both open the same motif modal, each showing its own view. When the view
    # paths are empty (convolution) fall back to the original single card.
    dual = !isempty(e.img_reduced) && !isempty(e.img_region)
    # The region-view logo draws whole filter windows, so it is labelled with
    # those windows rather than with the scattered mutation positions in `span`.
    # The reduced view, the metadata panel and the CSV keep `span` unchanged.
    region_name = let rs = region_window_span(e.wt_regions)
        isempty(rs) ? name : _tm_esc(rs)
    end
    cards_html, row_class = if dual
        (string(
            "<div class=\"top-mover-card\" data-median=\"$(e.median)\" onclick=\"openTopMover($rowid, 'reduced')\">",
            "<img src=\"$(e.img_reduced)\" alt=\"$name\" class=\"singleton-img\">",
            "<span class=\"singleton-filter-overlay\">$name</span></div>",
            "<div class=\"top-mover-card\" data-median=\"$(e.median)\" onclick=\"openTopMover($rowid, 'region')\">",
            "<img src=\"$(e.img_region)\" alt=\"$region_name\" class=\"singleton-img\">",
            "<span class=\"singleton-filter-overlay\">$region_name</span></div>"),
         "top-mover-row dual")
    else
        (string(
            "<div class=\"top-mover-card\" data-median=\"$(e.median)\" onclick=\"openTopMover($rowid)\">",
            "<img src=\"$(e.img)\" alt=\"$name\" class=\"singleton-img\">",
            "<span class=\"singleton-filter-overlay\">$name</span></div>"),
         "top-mover-row")
    end
    """
    <div class="$row_class">
        $cards_html
        <div class="top-mover-meta">
            $name_html$epi_block
        </div>
        <div class="top-mover-wt">$(wt_track_html(e.wt_regions))</div>
    </div>"""
end

"""
    top_mover_header_row_html() -> String

A grid-aligned header row placed above a dual-view list so the two logo columns
are self-labelled ("Reduced view" / "Region view").
"""
top_mover_header_row_html() = """
    <div class="top-mover-row dual top-mover-header">
        <div class="top-mover-colhead">Reduced view</div>
        <div class="top-mover-colhead">Region view</div>
        <div></div>
        <div></div>
    </div>"""

"""
    render_top_movers_page!(save_path; positives, negatives, page_title, nav_page_count, file)

Write the top-movers landing page to `save_path/file` (default `index.html`).
Reuses `styles.css` (already emitted by the grouped-page render). `positives`
and `negatives` are the vectors returned by [`select_top_movers`](@ref).

`feature_label` is the assay/feature name (units included, e.g.
`"ddG_ML_float (kcal/mol)"`). When given, each row's popup gains a plain-language
interpretation line built by [`mutation_description`](@ref) — *"Mutations at
sites 66, 69, 72 to K, L, P respectively increase ddG_ML_float (kcal/mol)."* The
sentence is omitted for motifs with no wild-type context, so the convolution case
simply never shows it.
"""
function render_top_movers_page!(save_path::AbstractString;
        positives::AbstractVector{TopMoverEntry},
        negatives::AbstractVector{TopMoverEntry},
        page_title::AbstractString="n/a",
        nav_page_count::Integer=4,
        feature_label::Union{AbstractString,Nothing}=nothing,
        protein_name::Union{AbstractString,Nothing}=nothing,
        protein_length::Union{Integer,Nothing}=nothing,
        wild_type::Union{AbstractString,Nothing}=nothing,  # full WT sequence → adds a copy button by the title
        show_epistasis::Bool=true,      # set false to hide the epistasis line (convolution case)
        modal_scroll_fix::Bool=false,   # set true to give the detail popup its own scrollbar
        generalization_warning::AbstractString="",  # "" => nothing is shown (default)
        transform_note::AbstractString="",  # rendered at the BOTTOM of the page; only used
                                            # when there is no generalization page to carry it
        runs_nav::AbstractString="",        # consensus page: links to the individual runs
        nav_override::AbstractString="null",# JS array of [href,label]; "null" = numbered nav
        consensus_rows::Union{Nothing,AbstractVector}=nothing,
                                            # one provenance line per row, positives then
                                            # negatives, appended under each card
        file::AbstractString="index.html")

    mkpath(save_path)

    # Optional <head> injection. The detail popup reuses the singleton modal,
    # whose overlay scrolls at the page level; `modal_scroll_fix` adds a popup-local
    # scrollbar so tall content is reachable without scrolling the dimmed backdrop.
    extra_head = modal_scroll_fix ?
        "<style>#singletonModal .singleton-modal-content,#multiMotifModal .multi-modal-content{max-height:90vh;overflow-y:auto;}</style>" : ""

    # Slick protein masthead: name (when supplied) plus an optional length note.
    has_name = protein_name !== nothing && !isempty(strip(protein_name))
    has_len  = protein_length !== nothing
    protein_header = if !has_name && !has_len
        ""
    else
        name_html = has_name ?
            string("<span class=\"protein-title-eyebrow\">Protein</span>",
                   "<span class=\"protein-title-name\">", _tm_esc(protein_name), "</span>") : ""
        sep_html = (has_name && has_len) ? "<span class=\"protein-title-sep\">·</span>" : ""
        len_html = has_len ?
            string("<span class=\"protein-title-meta\">length: ", protein_length, " amino acids</span>") : ""
        # One-click copy of the full wild-type sequence, placed beside the title.
        wt_html = (wild_type !== nothing && !isempty(strip(wild_type))) ?
            string("<button type=\"button\" class=\"wt-copy-btn\" data-wt=\"",
                   _tm_esc(wild_type), "\" onclick=\"copyWildType(this)\">Copy wild-type</button>") : ""
        string("<div class=\"protein-title\">", name_html, sep_html, len_html, wt_html, "</div>")
    end

    # rowid space: positives first, then negatives — matches the JS data array.
    all_rows = vcat(collect(positives), collect(negatives))
    payload = [(img=e.img, img_reduced=e.img_reduced, img_region=e.img_region,
                influence=e.influence, yy_kde=e.yy_kde,
                title=e.display_name, texts=e.texts,
                description=mutation_description(e.wt_regions, e.median, feature_label),
                variant_pwms=e.variant_pwms, variant_labels=e.variant_labels,
                variant_texts=e.variant_texts) for e in all_rows]
    data_js = "const topMoverData = " * JSON3.write(payload) * ";"

    # Prepend the "Reduced view / Region view" column header once per non-empty
    # list when any row uses the dual-view (mutagenesis) layout.
    is_dual(e) = !isempty(e.img_reduced) && !isempty(e.img_region)
    header = (any(is_dual, positives) || any(is_dual, negatives)) ?
        top_mover_header_row_html() * "\n" : ""

    # Consensus pages append a provenance line under each card; `nothing` (every
    # ordinary render) leaves the rows byte-identical to before.
    note(k) = consensus_rows === nothing || k > length(consensus_rows) ? "" : String(consensus_rows[k])
    pos_html = isempty(positives) ? "<p class=\"top-mover-empty\">No motifs increase the measured value.</p>" :
        header * join((top_mover_row_html(e, i - 1; show_epistasis=show_epistasis) * note(i) for (i, e) in enumerate(positives)), "\n")
    neg_html = isempty(negatives) ? "<p class=\"top-mover-empty\">No motifs decrease the measured value.</p>" :
        header * join((top_mover_row_html(e, length(positives) + i - 1; show_epistasis=show_epistasis) * note(length(positives) + i) for (i, e) in enumerate(negatives)), "\n")

    html_rendered = Mustache.render(html_template_top_movers;
        protein_name=page_title,
        protein_header=protein_header,
        generalization_warning=generalization_warning,
        transform_note=transform_note,
        runs_nav=runs_nav,
        nav_override=nav_override,
        extra_head=extra_head,
        upto=nav_page_count,
        top_mover_data=data_js,
        positive_rows=pos_html,
        negative_rows=neg_html)

    open(joinpath(save_path, file), "w") do io
        print(io, html_rendered)
    end
end

# Numeric coercion for the CSV columns: `parse_epistasis` hands back strings (it
# feeds the HTML panel), and absent fields must stay empty rather than become 0.
_tm_num(::Nothing) = missing
_tm_num(s::AbstractString) = something(tryparse(Float64, s), missing)

"""
    top_movers_dataframe(positives, negatives; protein_name=nothing, label=nothing) -> DataFrame

Flatten the two top-mover columns into one tidy table, one row per motif, at full
precision — the HTML panel rounds for display, this does not.

Rows carry `direction` (`"positive"`/`"negative"`) and a 1-based `rank` within
that direction, so the two columns stay distinguishable after datasets are
concatenated. `protein_name` and `label` are stamped onto every row so that
`vcat`-ing the per-dataset files yields a self-describing run-level table.

The epistasis string is split back into numeric `epistasis_coef` /
`epistasis_pvalue` (`missing` for singletons, which have no interaction term).
Wild-type context is flattened to `wt_seq` (region strings, `;`-joined) and
`mut_positions` (absolute 1-based positions whose consensus differs from wild
type, `;`-joined) — enough to line up recurring positions across datasets.

`description` carries the same plain-language sentence the card's popup shows
(see [`mutation_description`](@ref)) — e.g. *"Mutations at sites 66, 69, 72 to
K, L, P respectively increase ddG_ML_float (kcal/mol)."* — so the interpretation
travels with the table instead of living only in the rendered HTML. It is `""`
for motifs with no wild-type context (the convolution case) and for any motif
whose columns yield no confident residue call.

Note the table can hold fewer than `n` rows per direction, or none: motifs whose
cluster falls under `select_top_movers`'s `min_count` are dropped before ranking.

`report_location_z=true` adds one further column, `location_z`, immediately after
`nnd_pvalue`: the motif's location z-score (see [`location_z_1d`](@ref)), i.e.
how far its carriers' mean label sits from the population mean in standard
errors. `NaN` where it could not be computed. The column is **absent** by
default, so an existing consumer of this table sees the same 18 columns in the
same order it always did.
"""
function top_movers_dataframe(positives::AbstractVector{TopMoverEntry},
                              negatives::AbstractVector{TopMoverEntry};
                              protein_name::Union{AbstractString,Nothing}=nothing,
                              label::Union{AbstractString,Nothing}=nothing,
                              report_location_z::Bool=false)
    # Columns are declared up front (rather than inferred from pushed rows) so a
    # dataset with no surviving motifs still writes a full header — otherwise the
    # empty file would derail a glob-and-vcat over the run.
    df = DataFrame(
        protein_name     = String[],
        label            = String[],
        direction        = String[],
        rank             = Int[],
        display_name     = String[],
        group_label      = String[],
        span             = String[],
        is_singleton     = Bool[],
        count            = Int[],
        shapley_median   = Float64[],
        cluster_median   = Float64[],
        cluster_nnd      = Float64[],
        nnd_pvalue       = Float64[],
        epistasis_coef   = Union{Missing,Float64}[],
        epistasis_pvalue = Union{Missing,Float64}[],
        wt_seq           = String[],
        mut_positions    = String[],
        description      = String[],
    )
    pname = protein_name === nothing ? "" : String(protein_name)
    lbl   = label === nothing ? "" : String(label)
    locz  = Float64[]   # parallel to the rows; only spliced in when opted into
    for (direction, entries) in (("positive", positives), ("negative", negatives))
        for (rank, e) in enumerate(entries)
            push!(locz, e.location_z)
            coef, pval = parse_epistasis(e.epistasis)
            mut_pos = Int[]
            for reg in e.wt_regions, (j, flag) in enumerate(reg.mutated)
                flag && push!(mut_pos, reg.start + j - 1)
            end
            push!(df, (pname, lbl, direction, rank, e.display_name, e.group_label,
                       e.span, e.is_singleton, e.count, e.median, e.cluster_median,
                       e.cluster_nnd, e.nnd_p, _tm_num(coef), _tm_num(pval),
                       join((reg.wt for reg in e.wt_regions), ";"), join(mut_pos, ";"),
                       mutation_description(e.wt_regions, e.median, label)))
        end
    end
    # Spliced in right after `nnd_pvalue` so the tightness p-value and the
    # location z read side by side. Left out entirely when not opted into, which
    # is what keeps the default file byte-identical to the pre-change one.
    if report_location_z
        insertcols!(df, columnindex(df, :nnd_pvalue) + 1, :location_z => locz)
    end
    return df
end

"""
    write_top_movers_csv(path, positives, negatives; protein_name=nothing, label=nothing, append=false)

Write [`top_movers_dataframe`](@ref) to `path`, creating parent directories as
needed. `append=true` adds rows without repeating the header, so a multi-output
run can accumulate every label into one per-dataset file. Returns the DataFrame.

`report_location_z=true` adds the `location_z` column (see
[`top_movers_dataframe`](@ref)). Keep it consistent across the writes that share
one appended file, or the header and the later rows will disagree.
"""
function write_top_movers_csv(path::AbstractString,
                              positives::AbstractVector{TopMoverEntry},
                              negatives::AbstractVector{TopMoverEntry};
                              protein_name::Union{AbstractString,Nothing}=nothing,
                              label::Union{AbstractString,Nothing}=nothing,
                              report_location_z::Bool=false,
                              append::Bool=false)
    df = top_movers_dataframe(positives, negatives; protein_name=protein_name, label=label,
                              report_location_z=report_location_z)
    mkpath(dirname(path))
    CSV.write(path, df; append=append, writeheader=!append)
    return df
end

export TopMoverEntry, select_top_movers, render_top_movers_page!,
       top_movers_dataframe, write_top_movers_csv

# ─────────────────────────────────────────────────────────────────────────────
# Machine-readable dump of the top movers.
#
# `write_top_movers_csv` loses everything a page needs to be REDRAWN: the image
# paths, the metadata texts, the wild-type regions. A consensus page built from
# several runs has to rebuild real `TopMoverEntry` values from finished run
# folders, so this dumps enough to do that, plus the path of each motif's carrier
# CSV — which is what `multirun.jl` matches motifs on across runs.
#
# Paths are stored exactly as the entry holds them: relative to the rendering
# folder. A reader that mounts the run elsewhere prefixes them; see
# `read_top_movers_json`.
# ─────────────────────────────────────────────────────────────────────────────

_wt_dict(r::WildTypeRegion) = Dict("wt"=>r.wt, "start"=>r.start,
                                   "mutated"=>collect(r.mutated), "obs"=>r.obs)

"""
    top_movers_payload(positives, negatives; protein_name=nothing, label=nothing) -> Dict

The data behind a rendered top-movers page, as plain `Dict`s. See
[`write_top_movers_json`](@ref).
"""
function top_movers_payload(positives::AbstractVector{TopMoverEntry},
                            negatives::AbstractVector{TopMoverEntry};
                            protein_name=nothing, label=nothing)
    entries = Dict{String,Any}[]
    for (direction, es) in (("positive", positives), ("negative", negatives))
        for (rank, e) in enumerate(es)
            push!(entries, Dict{String,Any}(
                "direction"=>direction, "rank"=>rank,
                "median"=>e.median, "cluster_median"=>e.cluster_median,
                "cluster_nnd"=>e.cluster_nnd, "count"=>e.count,
                "display_name"=>e.display_name, "group_label"=>e.group_label,
                "span"=>e.span, "nnd_p"=>e.nnd_p, "is_singleton"=>e.is_singleton,
                "epistasis"=>e.epistasis, "location_z"=>e.location_z,
                "wt_regions"=>[_wt_dict(r) for r in e.wt_regions],
                "img"=>e.img, "img_reduced"=>e.img_reduced, "img_region"=>e.img_region,
                "influence"=>e.influence, "yy_kde"=>e.yy_kde,
                "texts"=>collect(e.texts),
                "variant_pwms"=>collect(e.variant_pwms),
                "variant_labels"=>collect(e.variant_labels),
                "variant_texts"=>[collect(v) for v in e.variant_texts],
                # the motif's carrier list lives beside its logo, same stem
                "carriers_csv"=>replace(e.img, r"\.png$"=>".csv"),
            ))
        end
    end
    Dict{String,Any}("protein_name"=>something(protein_name, ""),
                     "label"=>something(label, ""), "entries"=>entries)
end

"""
    write_top_movers_json(path, positives, negatives; protein_name=nothing, label=nothing)

Write [`top_movers_payload`](@ref) to `path`. Written beside `top_motifs.csv` on
every render; the CSV stays the human/analysis artifact and this is the one a
combiner reads.
"""
function write_top_movers_json(path::AbstractString,
                               positives::AbstractVector{TopMoverEntry},
                               negatives::AbstractVector{TopMoverEntry};
                               protein_name=nothing, label=nothing)
    mkpath(dirname(path))
    open(path, "w") do io
        JSON3.write(io, top_movers_payload(positives, negatives;
                                           protein_name=protein_name, label=label))
    end
    return path
end

"""
    read_top_movers_json(path; prefix="") -> (positives, negatives, carriers)

Rebuild `TopMoverEntry` values from a dump. `prefix` is prepended to every asset
path, so a page written one level above the run folder can point at
`run_3/renderings_x/...` while the run's own page keeps working unchanged.

`carriers` is a parallel vector of carrier-CSV paths (also prefixed), aligned to
`vcat(positives, negatives)`.
"""
function read_top_movers_json(path::AbstractString; prefix::AbstractString="")
    d = JSON3.read(read(path, String))
    pos = TopMoverEntry[]; neg = TopMoverEntry[]; car = String[]
    pre(p) = isempty(prefix) || isempty(p) ? String(p) : joinpath(prefix, String(p))
    for e in d.entries
        wt = [WildTypeRegion(String(r.wt), Int(r.start), Bool.(collect(r.mutated)), String(r.obs))
              for r in e.wt_regions]
        entry = TopMoverEntry(
            Float64(e.median), Float64(e.cluster_median), Float64(e.cluster_nnd), Int(e.count),
            String(e.display_name), String(e.group_label), String(e.span), Float64(e.nnd_p),
            Bool(e.is_singleton), String(e.epistasis), wt,
            pre(e.img), pre(e.img_reduced), pre(e.img_region), pre(e.influence), pre(e.yy_kde),
            String.(collect(e.texts)), String.(collect(e.variant_pwms)),
            String.(collect(e.variant_labels)),
            [String.(collect(v)) for v in e.variant_texts], Float64(e.location_z))
        push!(String(e.direction) == "positive" ? pos : neg, entry)
        push!(car, pre(e.carriers_csv))
    end
    return pos, neg, car
end

export top_movers_payload, write_top_movers_json, read_top_movers_json
