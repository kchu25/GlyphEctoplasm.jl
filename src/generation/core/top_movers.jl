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
    # Self-contained modal payload (paths relative to save folder)
    img::String              # motif logo
    influence::String        # Shapley/influence plot
    yy_kde::String           # indicator (yy-KDE) plot
    texts::Vector{String}    # metadata text fields (rendered as HTML in the modal)
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
    top_mover_row_html(entry, rowid) -> String

Render one summary row: a clickable card (left) wired to `openTopMover(rowid)`
and a metadata table (right). `rowid` is the 0-based index into the page's
`topMoverData` JS array.
"""
function top_mover_row_html(e::TopMoverEntry, rowid::Int)
    fmt(x) = isnan(x) ? "—" : @sprintf("%.3f", x)
    pval_row = isnan(e.nnd_p) ? "" :
        "<tr><td>cluster NND p-value</td><td>$(@sprintf("%.3f", e.nnd_p))</td></tr>"
    name = _tm_esc(e.display_name)
    """
    <div class="top-mover-row">
        <div class="top-mover-card" data-median="$(e.median)" onclick="openTopMover($rowid)">
            <img src="$(e.img)" alt="$name" class="singleton-img">
            <span class="singleton-filter-overlay">$name</span>
        </div>
        <div class="top-mover-meta">
            <div class="top-mover-name">$name</div>
            <table class="top-mover-table">
                <tr><td>Shapley median</td><td>$(fmt(e.median))</td></tr>
                <tr><td>cluster median</td><td>$(fmt(e.cluster_median))</td></tr>
                <tr><td>cluster NND</td><td>$(fmt(e.cluster_nnd))</td></tr>
                $pval_row
                <tr><td>position</td><td>$(_tm_esc(e.span))</td></tr>
                <tr><td>count</td><td>$(e.count)</td></tr>
                <tr><td>group</td><td>$(_tm_esc(e.group_label))</td></tr>
            </table>
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
        file::AbstractString="index.html")

    mkpath(save_path)

    # rowid space: positives first, then negatives — matches the JS data array.
    all_rows = vcat(collect(positives), collect(negatives))
    payload = [(img=e.img, influence=e.influence, yy_kde=e.yy_kde,
                title=e.display_name, texts=e.texts) for e in all_rows]
    data_js = "const topMoverData = " * JSON3.write(payload) * ";"

    pos_html = isempty(positives) ? "<p class=\"top-mover-empty\">No positive motifs.</p>" :
        join((top_mover_row_html(e, i - 1) for (i, e) in enumerate(positives)), "\n")
    neg_html = isempty(negatives) ? "<p class=\"top-mover-empty\">No negative motifs.</p>" :
        join((top_mover_row_html(e, length(positives) + i - 1) for (i, e) in enumerate(negatives)), "\n")

    html_rendered = Mustache.render(html_template_top_movers;
        protein_name=page_title,
        upto=nav_page_count,
        top_mover_data=data_js,
        positive_rows=pos_html,
        negative_rows=neg_html)

    open(joinpath(save_path, file), "w") do io
        print(io, html_rendered)
    end
end

export TopMoverEntry, select_top_movers, render_top_movers_page!
