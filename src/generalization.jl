
"""Compute R² coefficient"""
function _compute_r2(y_true::AbstractVector{T}, y_pred::AbstractVector{T}) where T<:AbstractFloat
    ss_res = sum((y_true .- y_pred).^2)
    ss_tot = sum((y_true .- mean(y_true)).^2)
    return T(1) - ss_res / ss_tot
end

# Shared axis spine/grid style
function _style_axis!(ax)
    ax.xgridvisible      = false
    ax.ygridvisible      = false
    ax.topspinevisible   = false
    ax.rightspinevisible = false
    ax.spinewidth        = 0.8
    ax.xtickwidth        = 0.8
    ax.ytickwidth        = 0.8
end

"""
    publication_kde_panel(data_pairs;
        figsize=(1200, 400), nlevels=12, colormap=:BuPu,
        hist_bins=40, bandwidth=nothing, dpi=300, save_path=nothing)

Create a publication-quality multi-panel 2D KDE density plot with
marginal histograms. Uses smooth kernel density estimation (via
KernelDensity.jl + contourf!) instead of hexbins, giving a continuous,
artifact-free density cloud — the standard in genomics expression papers.

`data_pairs`         vector of (x, y, xlabel, ylabel, title) tuples.
`nlevels`            number of filled contour levels; more = smoother gradient.
`colormap`           colormap; default :BuPu (white→purple, perceptually uniform).
`bandwidth`          2D KDE bandwidth (nothing = automatic Silverman's rule).
`marginal_bandwidth` 1D KDE bandwidth for marginals (nothing = automatic).
`dpi`                dots-per-inch for raster output; 300 for print-ready.
"""
function publication_kde_panel(data_pairs;
        figsize=(1200, 400), nlevels=12,
        colormap=:BuPu,
        bandwidth=nothing, marginal_bandwidth=nothing, dpi=300, save_path=nothing)

    n = length(data_pairs)

    fig = Figure(
        size           = figsize,
        fontsize       = 11,
        figure_padding = (12, 20, 12, 12),
    )

    main_axes = Axis[]

    for (i, (x, y, xlabel, ylabel, title)) in enumerate(data_pairs)
        xv = vec(Float64.(x))
        yv = vec(Float64.(y))
        r2 = _compute_r2(yv, xv)  # y_true=y (labels), y_pred=x (predictions)

        col = fig[1, i] = GridLayout()

        # ── top marginal ────────────────────────────────────────────────────
        ax_top = Axis(col[1, 1];
            xticksvisible      = false,
            xticklabelsvisible = false,
            yticksvisible      = false,
            yticklabelsvisible = false,
            leftspinevisible   = false,
            rightspinevisible  = false,
            topspinevisible    = false,
            bottomspinevisible = false,
            xgridvisible       = false,
            ygridvisible       = false,
            height             = 50,
        )
        density!(ax_top, xv;
            color       = (:gray55, 0.45),
            strokecolor = (:gray40, 0.8),
            strokewidth = 1.0,
            (isnothing(marginal_bandwidth) ? () : (bandwidth = marginal_bandwidth,))...
        )

        # ── main KDE axis ───────────────────────────────────────────────────
        ax = Axis(col[2, 1];
            xlabel         = xlabel,
            ylabel         = ylabel,
            title          = title,
            titlesize      = 12,
            xlabelsize     = 11,
            ylabelsize     = 11,
            xticklabelsize = 9,
            yticklabelsize = 9,
            aspect         = DataAspect(),
        )
        _style_axis!(ax)
        push!(main_axes, ax)

        # 2D KDE — evaluate on a regular grid then filled contour
        kd = isnothing(bandwidth) ?
            kde((xv, yv)) :
            kde((xv, yv); bandwidth=(bandwidth, bandwidth))

        # Power-transform density to expand low-density tail regions for better visibility
        # sqrt and higher powers compress the high-density core while stretching tails
        density_transformed = kd.density.^0.3

        # contourf gives the smooth filled gradient
        cf = contourf!(ax,
            kd.x, kd.y, density_transformed;
            levels   = nlevels,
            colormap = colormap,
            extendhigh = :auto,
        )

        # ── right marginal ──────────────────────────────────────────────────
        ax_right = Axis(col[2, 2];
            xticksvisible      = false,
            xticklabelsvisible = false,
            yticksvisible      = false,
            yticklabelsvisible = false,
            leftspinevisible   = false,
            rightspinevisible  = false,
            topspinevisible    = false,
            bottomspinevisible = false,
            xgridvisible       = false,
            ygridvisible       = false,
            width              = 50,
        )
        density!(ax_right, yv;
            color       = (:gray55, 0.45),
            strokecolor = (:gray40, 0.8),
            strokewidth = 1.0,
            direction   = :y,
            (isnothing(marginal_bandwidth) ? () : (bandwidth = marginal_bandwidth,))...
        )

        # ── colorbar ────────────────────────────────────────────────────────
        Colorbar(col[2, 3], cf;
            label         = "density",
            labelsize     = 9,
            ticklabelsize = 8,
            tickwidth     = 0.6,
            width         = 8,
            spinewidth    = 0.6,
        )

        linkyaxes!(ax, ax_right)
        linkxaxes!(ax, ax_top)

        # identity line — solid black
        lims = (min(minimum(xv), minimum(yv)), max(maximum(xv), maximum(yv)))
        lines!(ax, [lims[1], lims[2]], [lims[1], lims[2]];
            color     = :black,
            linewidth = 1.5,
        )

        # R² annotation
        text!(ax, 0.05, 0.93;
            text     = @sprintf("R² = %.3f", r2),
            space    = :relative,
            align    = (:left, :top),
            fontsize = 10,
            font     = :bold,
            color    = :black,
        )

        rowgap!(col, 1, 2)
        colgap!(col, 1, 2)
        colgap!(col, 2, 4)
        rowsize!(col, 1, Relative(0.15))
        colsize!(col, 2, Relative(0.15))
    end

    length(main_axes) > 1 && linkyaxes!(main_axes...)
    colgap!(fig.layout, 28)

    if !isnothing(save_path)
        save(save_path, fig; px_per_unit = dpi / 96)
    end

    return fig
end

# ────────────────────────────────────────────────────────────────────────────
# Low-generalization warning banner
#
# Why this exists: the motif pages ALWAYS report groups. The final selection is
# `top_and_bot_counts` = top 8 + bottom 8 by median Banzhaf, which is a RANKING,
# not a significance test. Measured on null data (labels replaced with pure
# noise, and separately with the real labels shuffled), the pipeline still
# returned exactly 16 singleton groups in 9 out of 9 runs, with positions, signs
# and magnitudes that look identical to real findings. In one of six null runs
# the magnitudes were inside the range seen on genuine signal.
#
# What DOES separate the two cleanly is held-out R²: ~0.87 on real signal,
# ~0.00 on both nulls. So R² is the honest gate, and this banner surfaces it on
# the page rather than leaving it in a log nobody reads.
# ────────────────────────────────────────────────────────────────────────────

"""Default held-out R² below which reported motifs are flagged as unreliable."""
const DEFAULT_GENERALIZATION_WARN_THRESHOLD = 0.15

"""
    generalization_warning_html(test_r2; threshold=DEFAULT_GENERALIZATION_WARN_THRESHOLD)

HTML banner warning that the motifs on this page may be meaningless, or `""`
when the model generalizes well enough (or when `test_r2` is unavailable).

Returns `""` for `nothing` and `NaN` so callers can pass a possibly-missing
value without branching.
"""
function generalization_warning_html(test_r2; threshold=DEFAULT_GENERALIZATION_WARN_THRESHOLD)
    test_r2 === nothing && return ""
    (isa(test_r2, Real) && isfinite(test_r2)) || return ""
    test_r2 >= threshold && return ""
    r2s = string(round(test_r2, digits=3))
    ths = string(round(threshold, digits=3))
    return string(
        "<div class=\"generalization-warning\" role=\"alert\">",
        "<span class=\"generalization-warning-badge\">Low generalization</span>",
        "<div class=\"generalization-warning-body\">",
        "<strong>These motifs may be meaningless.</strong> ",
        "Held-out R&sup2; is <code>", r2s, "</code>, below the ", ths, " threshold, ",
        "so the model does not predict unseen variants.",
        "<br>The motif list is a ranking, not a significance test: the top and bottom ",
        "groups are always reported, even on data with no signal at all. ",
        "On null data (pure noise, and shuffled labels) this pipeline still returned ",
        "a full set of groups with plausible-looking positions and magnitudes. ",
        "Treat everything below as unverified until R&sup2; improves.",
        "</div></div>")
end

# ────────────────────────────────────────────────────────────────────────────
# Label-transform note
#
# The pipeline can now choose the label normalization from the dataset itself
# (MotifInference `:auto`). That choice changes the meaning of every number on
# the page: under a log transform the predictions, the Shapley medians and the
# indicator axis are all on the log scale. A reader who does not know which
# scale they are looking at cannot read the page correctly.
#
# So whenever the choice was made automatically, it is stated here.
# ────────────────────────────────────────────────────────────────────────────

"""
    transform_note_html(note) -> String

An informational banner naming the label transform, or `""` when `note` is
empty or `nothing`.

`note` is a plain sentence produced by the caller, for example
*"log, chosen automatically: every label is positive (min 0.001) and the labels
are right skewed (skew 5.99). Predictions and motif values are on the log
scale."*

This is not a warning. It renders in the neutral informational style, because a
transform having been applied is normal and expected.

`bottom=true` adds the `transform-note-bottom` class, which is the placement the
motif and summary pages use: the note is a footnote about the axis units, not
something a reader needs before they can look at the figures. The low-generalization
banner keeps the top of the page, because that one IS a warning and changes whether
the page should be read at all.

Where it actually appears is decided by the caller — see `plot_motifs_mut_case`.
Every transform the pipeline chooses maps its numbers back to the assay's own
units, so the note is a record of what happened during training, not a caveat on
what is displayed. It therefore sits at the bottom of the generalization page
alone, and is left off the motif and summary pages entirely.
"""
function transform_note_html(note; bottom::Bool=false)
    note === nothing && return ""
    txt = strip(string(note))
    isempty(txt) && return ""
    esc = replace(txt, "&" => "&amp;", "<" => "&lt;", ">" => "&gt;")
    cls = bottom ? "transform-note transform-note-bottom" : "transform-note"
    return string(
        "<div class=\"", cls, "\" role=\"note\">",
        "<span class=\"transform-note-badge\">Label transform</span>",
        "<div class=\"transform-note-body\">", esc, "</div></div>")
end

