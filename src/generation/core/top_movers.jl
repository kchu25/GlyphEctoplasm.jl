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
    WildTypeRegion(wt, start, mutated)

One region's wild-type amino-acid window for the summary card: the full WT string
`wt` (un-abbreviated), the 1-based `start` position of its first residue, and a
per-column `mutated` flag (true where the motif's observed consensus differs from
the wild type). The panel draws one up-arrow under each mutated residue.
"""
struct WildTypeRegion
    wt::String
    start::Int
    mutated::Vector{Bool}
end

"""
    TopMoverEntry

One summary row. Carries the ranking keys, the metadata shown in the right-hand
panel, and a self-contained modal payload (image paths + text fields) so the
detail modal can be populated without the main page's `jsonData`.

Path fields (`img`, `influence`, `yy_kde`) are relative to the page's save
folder, matching the relative paths used by the grouped page.
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
    img::String              # motif logo
    influence::String        # Shapley/influence plot
    yy_kde::String           # indicator (yy-KDE) plot
    texts::Vector{String}    # metadata text fields (rendered as HTML in the modal)
    # Multi-motif distance variants (empty for singletons). When non-empty the
    # detail popup shows the multi-motif modal with an inter-motif-distance slider
    # instead of the static singleton modal, mirroring the grouped Motifs page.
    variant_pwms::Vector{String}            # per-variant logo paths (one per distance)
    variant_labels::Vector{String}          # per-variant distance labels
    variant_texts::Vector{Vector{String}}   # per-variant metadata text fields
end

"""
    select_top_movers(entries; n=5, bin_count=10) -> (positives, negatives)

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
"""
function select_top_movers(entries::AbstractVector{TopMoverEntry}; n::Int=5, bin_count::Int=10)
    isempty(entries) && return (TopMoverEntry[], TopMoverEntry[])

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

# Format a numeric stat for the summary chips: "—" for NaN, optional leading
# "+" so positive Shapley values read as gains.
function _tm_num(x::Real; signed::Bool=false, digits::Int=2)
    isnan(x) && return "—"
    s = string(round(x; digits=digits))
    (signed && x > 0) ? "+" * s : s
end

# Compact p-value: blank when unavailable, "<0.001" below the floor.
_tm_pval(p::Real) = isnan(p) ? "" : (p < 0.001 ? "&lt;0.001" : string(round(p; digits=3)))

# One labeled stat chip (small caps label over a bold value).
_tm_chip(label::AbstractString, value::AbstractString) = string(
    "<div class=\"tm-chip\"><span class=\"tm-chip-label\">", label,
    "</span><span class=\"tm-chip-val\">", value, "</span></div>")

"""
    stats_chips_html(entry) -> String

Render the right-hand stat strip from the entry's metadata: Shapley median
(sign-colored), expression median, mean NND (with NND p-value when available),
mutation position (when present), and sequence count.
"""
function stats_chips_html(e::TopMoverEntry)
    mcls = e.median > 0 ? "tm-pos" : (e.median < 0 ? "tm-neg" : "")
    mval = "<span class=\"$mcls\">$(_tm_num(e.median; signed=true))</span>"

    nnd_val = _tm_num(e.cluster_nnd)
    pstr = _tm_pval(e.nnd_p)
    if nnd_val != "—" && !isempty(pstr)
        nnd_val = string(nnd_val, " <span class=\"tm-chip-sub\">p ", pstr, "</span>")
    end

    chips = String[
        _tm_chip("Shapley median", mval),
        _tm_chip("Expression median", _tm_num(e.cluster_median)),
        _tm_chip("Mean NND", nnd_val),
    ]
    isempty(strip(e.span)) || push!(chips, _tm_chip("Position", _tm_esc(e.span)))
    push!(chips, _tm_chip("Count", string(e.count)))
    return "<div class=\"top-mover-stats\">" * join(chips) * "</div>"
end

"""
    top_mover_row_html(entry, rowid) -> String

Render one summary row: a clickable card (left) wired to `openTopMover(rowid)`
and an info column (right) holding a group badge + name, a labeled stat strip,
and — when present — the epistasis line and wild-type track. `rowid` is the
0-based index into the page's `topMoverData` JS array.
"""
function top_mover_row_html(e::TopMoverEntry, rowid::Int; show_epistasis::Bool=true)
    name = _tm_esc(e.display_name)
    group = _tm_esc(e.group_label)
    epi_block = if !show_epistasis
        ""                                          # caller opted out (e.g. convolution case)
    else
        epi = if isempty(strip(e.epistasis))
            "<div class=\"epi-line epi-none\">epistasis: n/a</div>"
        else
            coef, pval = parse_epistasis(e.epistasis)
            if coef === nothing
                "<div class=\"epi-line\">$(e.epistasis)</div>"   # unknown format: show as-is
            else
                pv = pval === nothing ? "" : "<div class=\"epi-line epi-sub\">p-value: $(pval)</div>"
                "<div class=\"epi-line\">epistasis: <strong>$(coef)</strong></div>$pv"
            end
        end
        "\n            <div class=\"top-mover-epistasis\">$epi</div>"
    end
    wt = wt_track_html(e.wt_regions)
    wt_block = isempty(wt) ? "" : "\n            <div class=\"top-mover-wt\">$wt</div>"
    group_badge = isempty(strip(group)) ? "" : "<span class=\"top-mover-group-badge\">$group</span>"
    """
    <div class="top-mover-row">
        <div class="top-mover-card" data-median="$(e.median)" onclick="openTopMover($rowid)">
            <img src="$(e.img)" alt="$name" class="singleton-img">
            <span class="singleton-filter-overlay">$name</span>
        </div>
        <div class="top-mover-info">
            <div class="top-mover-head">$group_badge<span class="top-mover-name">$name</span></div>
            $(stats_chips_html(e))$epi_block$wt_block
        </div>
    </div>"""
end

"""
    render_top_movers_page!(save_path; positives, negatives, page_title, nav_page_count, file)

Write the top-movers landing page to `save_path/file` (default `index.html`).
Reuses `styles.css` (already emitted by the grouped-page render). `positives`
and `negatives` are the vectors returned by [`select_top_movers`](@ref).
"""
function render_top_movers_page!(save_path::AbstractString;
        positives::AbstractVector{TopMoverEntry},
        negatives::AbstractVector{TopMoverEntry},
        page_title::AbstractString="n/a",
        nav_page_count::Integer=4,
        protein_name::Union{AbstractString,Nothing}=nothing,
        protein_length::Union{Integer,Nothing}=nothing,
        show_epistasis::Bool=true,      # set false to hide the epistasis line (convolution case)
        modal_scroll_fix::Bool=false,   # set true to give the detail popup its own scrollbar
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
        string("<div class=\"protein-title\">", name_html, sep_html, len_html, "</div>")
    end

    # rowid space: positives first, then negatives — matches the JS data array.
    all_rows = vcat(collect(positives), collect(negatives))
    payload = [(img=e.img, influence=e.influence, yy_kde=e.yy_kde,
                title=e.display_name, texts=e.texts,
                variant_pwms=e.variant_pwms, variant_labels=e.variant_labels,
                variant_texts=e.variant_texts) for e in all_rows]
    data_js = "const topMoverData = " * JSON3.write(payload) * ";"

    pos_html = isempty(positives) ? "<p class=\"top-mover-empty\">No positive motifs.</p>" :
        join((top_mover_row_html(e, i - 1; show_epistasis=show_epistasis) for (i, e) in enumerate(positives)), "\n")
    neg_html = isempty(negatives) ? "<p class=\"top-mover-empty\">No negative motifs.</p>" :
        join((top_mover_row_html(e, length(positives) + i - 1; show_epistasis=show_epistasis) for (i, e) in enumerate(negatives)), "\n")

    html_rendered = Mustache.render(html_template_top_movers;
        protein_name=page_title,
        protein_header=protein_header,
        extra_head=extra_head,
        upto=nav_page_count,
        top_mover_data=data_js,
        positive_rows=pos_html,
        negative_rows=neg_html)

    open(joinpath(save_path, file), "w") do io
        print(io, html_rendered)
    end
end

export TopMoverEntry, select_top_movers, render_top_movers_page!
