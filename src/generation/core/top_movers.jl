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

    if length(sites) == 1
        return string("Mutation at site ", sites[1], " to ", residues[1], " ",
                      verb_up ? "increases" : "decreases", " ", what, ".")
    end
    return string("Mutations at sites ", join(sites, ", "), " to ",
                  join(residues, ", "), " respectively ",
                  verb_up ? "increase" : "decrease", " ", what, ".")
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
end

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
    cards_html, row_class = if dual
        (string(
            "<div class=\"top-mover-card\" data-median=\"$(e.median)\" onclick=\"openTopMover($rowid, 'reduced')\">",
            "<img src=\"$(e.img_reduced)\" alt=\"$name\" class=\"singleton-img\">",
            "<span class=\"singleton-filter-overlay\">$name</span></div>",
            "<div class=\"top-mover-card\" data-median=\"$(e.median)\" onclick=\"openTopMover($rowid, 'region')\">",
            "<img src=\"$(e.img_region)\" alt=\"$name\" class=\"singleton-img\">",
            "<span class=\"singleton-filter-overlay\">$name</span></div>"),
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

    pos_html = isempty(positives) ? "<p class=\"top-mover-empty\">No positive motifs.</p>" :
        header * join((top_mover_row_html(e, i - 1; show_epistasis=show_epistasis) for (i, e) in enumerate(positives)), "\n")
    neg_html = isempty(negatives) ? "<p class=\"top-mover-empty\">No negative motifs.</p>" :
        header * join((top_mover_row_html(e, length(positives) + i - 1; show_epistasis=show_epistasis) for (i, e) in enumerate(negatives)), "\n")

    html_rendered = Mustache.render(html_template_top_movers;
        protein_name=page_title,
        protein_header=protein_header,
        generalization_warning=generalization_warning,
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
"""
function top_movers_dataframe(positives::AbstractVector{TopMoverEntry},
                              negatives::AbstractVector{TopMoverEntry};
                              protein_name::Union{AbstractString,Nothing}=nothing,
                              label::Union{AbstractString,Nothing}=nothing)
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
    for (direction, entries) in (("positive", positives), ("negative", negatives))
        for (rank, e) in enumerate(entries)
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
    return df
end

"""
    write_top_movers_csv(path, positives, negatives; protein_name=nothing, label=nothing, append=false)

Write [`top_movers_dataframe`](@ref) to `path`, creating parent directories as
needed. `append=true` adds rows without repeating the header, so a multi-output
run can accumulate every label into one per-dataset file. Returns the DataFrame.
"""
function write_top_movers_csv(path::AbstractString,
                              positives::AbstractVector{TopMoverEntry},
                              negatives::AbstractVector{TopMoverEntry};
                              protein_name::Union{AbstractString,Nothing}=nothing,
                              label::Union{AbstractString,Nothing}=nothing,
                              append::Bool=false)
    df = top_movers_dataframe(positives, negatives; protein_name=protein_name, label=label)
    mkpath(dirname(path))
    CSV.write(path, df; append=append, writeheader=!append)
    return df
end

export TopMoverEntry, select_top_movers, render_top_movers_page!,
       top_movers_dataframe, write_top_movers_csv
