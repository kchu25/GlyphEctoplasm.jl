
"""
    estimate_point_density(x, y)

Estimate per-point local density using 2D KDE, returning a density value for each (x,y) pair.
"""
function estimate_point_density(x, y)
    k = kde((x, y))
    # Interpolate KDE grid back to each point's location
    # Find nearest grid indices for each point
    densities = Vector{Float64}(undef, length(x))
    for i in eachindex(x)
        ix = searchsortedlast(k.x, x[i])
        iy = searchsortedlast(k.y, y[i])
        ix = clamp(ix, 1, length(k.x))
        iy = clamp(iy, 1, length(k.y))
        densities[i] = k.density[ix, iy]
    end
    return densities
end


"""
    _sync_xy_ticks!(ax; n_ticks=5)

Force an Axis to have identical x and y limits and identical tick positions,
which is appropriate for any prediction-vs-observed plot where both axes
represent the same quantity.  Computes a shared range from the current data
limits, generates `n_ticks` nicely rounded ticks, and applies them to both axes.
"""
function _sync_xy_ticks!(ax; n_ticks=5)
    # Read current auto limits
    xl = ax.xaxis.attributes.limits[]
    yl = ax.yaxis.attributes.limits[]
    lo = min(xl[1], yl[1])
    hi = max(xl[2], yl[2])

    # Compute nice round tick positions using Makie's own algorithm
    ticks = Makie.get_tickvalues(Makie.LinearTicks(n_ticks), lo, hi)

    # Extend limits slightly beyond outermost ticks
    margin = (ticks[end] - ticks[1]) * 0.04
    shared_lims = (ticks[1] - margin, ticks[end] + margin)

    ax.xticks = ticks
    ax.yticks = ticks
    xlims!(ax, shared_lims)
    ylims!(ax, shared_lims)
    return ax
end


"""
    plot_labels_vs_procprod(pts, is_in_intersect; show_density=true, show_r2=false, motif_label="Contain motif")

Plot predicted values vs. observed labels, highlighting which sequences contain a given motif.

The goal is to see where motif-containing sequences sit relative to the full population. With 
thousands of points, a flat scatter becomes illegible, so this function uses **four visual layers** 
to make density patterns clear:

## Visual Layers (back to front)

**Layer 1 — Background population scatter**  
All non-motif sequences plotted as medium gray-blue circles (markersize=8, alpha=0.85). These form 
the "baseline cloud" that motif-containing points are compared against. Stroke width is 0.7 for 
subtle definition.

**Layer 2 — Filled density contours (motif points only)**  
A 2D kernel density estimate (KDE) is computed over motif-containing points using a Gaussian kernel 
with bandwidth chosen by Scott's rule (h ∝ n^{-1/6}). The density surface is rendered as filled 
contours (`contourf!`) with 8 smoothly graded levels using a warm colormap:
- Low density → nearly transparent pale yellow (α ≈ 0.02)
- Medium density → translucent orange (α ≈ 0.15–0.35)
- High density → semi-opaque rust-red (α ≈ 0.55)

This creates a visual "heatmap wash" where dense motif clusters glow warmer and more intense, while 
individual markers remain visible through the transparency.

**Layer 3 — Foreground scatter (motif points)**  
Individual motif-containing points drawn as uniform orange circles (markersize=7, alpha=0.90). All 
points have the same color and transparency — no per-point density scaling — so the density structure 
comes entirely from the contours beneath them. Stroke width is 0.7 with a dark outline for crisp 
edges.

**Layer 4 — Contour lines (drawn last, on top)**  
Bold contour lines (linewidth=2.5) marking the 3 highest density isosurfaces, drawn in dark brown-red 
(RGBA 0.30, 0.05, 0.00, 0.95). By placing these *after* the scatter, they remain visible even over 
dense point clusters, clearly delineating the core concentration zones. Only 3 levels are shown to 
avoid clutter — these mark the densest regions where motif-containing points concentrate most strongly.

## Technical Details

- **KDE bandwidth**: Scott's rule, which uses σ̂ · n^{-1/(d+4)} where σ̂ is the standard deviation 
  and d=2 for bivariate data. This gives adaptive smoothing that works across varying sample sizes.
  
- **Contour placement**: The `contourf!` call creates the filled regions in Layer 2, then `contour!` 
  draws just the lines in Layer 4 (same KDE, different rendering). This ensures lines aren't hidden 
  by scatter.

- **Legend markers**: Both background and foreground get separate invisible scatter points 
  (markersize=12) for the legend, so legend entries are large and clear regardless of actual plot 
  marker size.

- **Minimum point threshold**: Density features (contours + filled regions) only appear when 
  `n_motif > 10`. Below that, motif points render as flat scatter to avoid spurious density artifacts.

## Arguments
- `pts::NamedTuple`: Must contain `proc_prod` (predicted values, x-axis) and `labels` (observed 
  values, y-axis) as numeric vectors of equal length.
  
- `is_in_intersect::BitVector`: Boolean mask with `true` for sequences containing the motif. Length 
  must match `pts.proc_prod` and `pts.labels`.
  
- `show_density::Bool = true`: Enable/disable density visualization (contours + filled regions). 
  When `false`, motif points render as plain scatter with no KDE overlay.
  
- `show_r2::Bool = false`: Annotate the plot with R² (coefficient of determination) between 
  predictions and labels, displayed in the bottom-right corner.
  
- `motif_label::String = "Contain motif"`: Legend label for the motif-containing group.

## Returns
- `fig::Figure`: A Makie Figure (800×800 px) with axes configured for publication: no top/right 
  spines, large tick labels (28pt), large axis labels (32pt), clean white background.

## Example
```julia
# pts has fields: proc_prod (predictions), labels (observed)
# is_motif is a BitVector marking which sequences contain the motif
fig = plot_labels_vs_procprod(pts, is_motif; 
    show_density=true, show_r2=true, motif_label="Contains TATA box")
save("motif_enrichment.png", fig)
```
"""
function plot_labels_vs_procprod(pts, is_in_intersect; show_density=true, show_r2=false, motif_label="Contain motif")
    fig = Figure(size=(800, 800))
    ax = Axis(fig[1, 1], 
        xlabel="Predicted values", 
        ylabel="Labels", 
        xlabelsize=32,
        ylabelsize=32,
        xticklabelsize=28,
        yticklabelsize=28,
        topspinevisible=false,
        rightspinevisible=false,
        xticksize=8,
        yticksize=8,
        spinewidth=2,
        xgridvisible=false,
        ygridvisible=false)

    # ── Indices ──
    bg_mask = .!is_in_intersect
    fg_mask = is_in_intersect

    bg_x = pts.proc_prod[bg_mask]
    bg_y = pts.labels[bg_mask]
    fg_x = pts.proc_prod[fg_mask]
    fg_y = pts.labels[fg_mask]

    n_fg = sum(fg_mask)

    # ═══════════════════════════════════════════════════════════════════════
    # Layer 1 (back): Scatter for non-motif population
    # ═══════════════════════════════════════════════════════════════════════
    if length(bg_x) > 0
        scatter!(ax, bg_x, bg_y,
            color=RGBA(0.70, 0.70, 0.80, 0.85), 
            markersize=1, 
            marker=:circle,
            strokewidth=0.7, 
            strokecolor=RGBA(0.5, 0.5, 0.6, 0.5))
        
        # Legend entry with matching size as motif legend
        scatter!(ax, [NaN], [NaN],
            color=RGBA(0.70, 0.70, 0.80, 0.85), markersize=12, marker=:circle,
            strokewidth=0.5, strokecolor=RGBA(0.5, 0.5, 0.6, 0.5),
            label="Other sequences")
    end

    # ═══════════════════════════════════════════════════════════════════════
    # Layer 2 (mid): Filled density contours for motif-containing points
    # ═══════════════════════════════════════════════════════════════════════
    if show_density && n_fg > 10
        kde_result = kde((fg_x, fg_y))

        # Warm filled contours: light wash → strong red-orange
        warm_cmap = cgrad([RGBA(1.0, 0.88, 0.55, 0.02),
                           RGBA(0.99, 0.55, 0.15, 0.15),
                           RGBA(0.92, 0.25, 0.08, 0.35),
                           RGBA(0.80, 0.08, 0.02, 0.55)])

        contourf!(ax, kde_result.x, kde_result.y, kde_result.density,
            levels=8, colormap=warm_cmap)
    end

    # ═══════════════════════════════════════════════════════════════════════
    # Layer 3 (front): Motif-containing scatter with uniform transparency
    # ═══════════════════════════════════════════════════════════════════════
    if n_fg > 0
        scatter!(ax, fg_x, fg_y,
            color=RGBA(0.98, 0.40, 0.10, 0.90), markersize=6, marker=:circle, 
            strokewidth=0.7, strokecolor=RGBA(0.1, 0.1, 0.1, 0.6))
        
        # Legend entry with matching orange color
        scatter!(ax, [NaN], [NaN],
            color=RGBA(0.98, 0.40, 0.10, 0.90), markersize=12, marker=:circle,
            strokewidth=0.8, strokecolor=RGBA(0.1, 0.1, 0.1, 0.9),
            label=motif_label)
    end

    # ═══════════════════════════════════════════════════════════════════════
    # Layer 4 (top): Contour lines on top so they're always visible
    # ═══════════════════════════════════════════════════════════════════════
    if show_density && n_fg > 10
        # Bold contour lines for clearest structure — just the densest regions
        contour!(ax, kde_result.x, kde_result.y, kde_result.density,
            levels=3, color=RGBA(0.30, 0.05, 0.00, 0.95), linewidth=2.5)
    end

    # ═══════════════════════════════════════════════════════════════════════
    # R² annotation
    # ═══════════════════════════════════════════════════════════════════════
    if show_r2
        y_true = pts.labels
        y_pred = pts.proc_prod
        ss_res = sum((y_true .- y_pred).^2)
        ss_tot = sum((y_true .- mean(y_true)).^2)
        r2 = 1 - ss_res / ss_tot
        
        text!(ax, 0.98, 0.02, text=@sprintf("R² = %.3f", r2), 
            align=(:right, :bottom), fontsize=22, color=:black,
            space=:relative)
    end

    axislegend(ax, position=:lt, labelsize=28, markersize=18, 
               framevisible=false, backgroundcolor=:white)

    _sync_xy_ticks!(ax)
    
    return fig
end