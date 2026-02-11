
"""
    plot_labels_vs_procprod(pts, is_in_intersect; show_density=true, show_r2=false, motif_label="Contain motif")

Generate a scatter plot comparing predicted values (x-axis) against observed labels (y-axis).
Points are distinguished by whether they contain the specified motif pattern.

The density contours are estimated using two-dimensional kernel density estimation (KDE) with a 
Gaussian kernel. The bandwidth is automatically selected using Scott's rule (bandwidth ∝ n^(-1/6)
where n is the sample size). Each contour line represents an isoline of the estimated probability 
density function f̂(x,y), with 6 equally-spaced levels shown spanning from the minimum to maximum 
estimated density values (Note: density values represent the estimated probability per unit area 
in the (prediction, label) space, with higher values indicating greater concentration of points). 
The contours are computed using a bivariate Gaussian kernel: K_h(x,y) = (1/(2πh²))exp(-(x²+y²)/(2h²)), 
where h is the bandwidth parameter. Regions enclosed by inner contours indicate higher local density 
of motif-containing points, revealing clustered patterns in the prediction-label space.

# Arguments
- `pts`: NamedTuple containing `proc_prod` (predictions) and `labels` fields
- `is_in_intersect`: BitVector indicating which points contain the motif
- `show_density`: Boolean flag to display density contours (default: true)
- `show_r2`: Boolean flag to display R² value (default: false)
- `motif_label`: String label for motif-containing points (default: "Contain motif")

# Returns
- `fig`: Makie Figure object
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

    # Use a modern, ergonomic color scheme
    nonmotif_color = RGBA(0.88, 0.88, 0.95, 0.85)  # very light gray-blue
    motif_color = RGBA(0.98, 0.45, 0.15, 0.95)    # vivid orange

    # Add density contours for motif-containing points only
    if show_density && sum(is_in_intersect) > 2
        motif_x = pts.proc_prod[is_in_intersect]
        motif_y = pts.labels[is_in_intersect]
        kde_result = kde((motif_x, motif_y))
        contour!(ax, kde_result.x, kde_result.y, kde_result.density,
            levels=6, color=:gray30, linewidth=2.5, linestyle=:dash, alpha=0.7)
    end

    # Plot non-motif points first (background)
    scatter!(ax, pts.proc_prod[.!is_in_intersect], pts.labels[.!is_in_intersect],
        color=nonmotif_color, markersize=13, marker=:circle, label="Other sequences",
        strokewidth=0.8, strokecolor=RGBA(0.2,0.2,0.2,0.25))

    # Plot motif-containing points last (foreground)
    scatter!(ax, pts.proc_prod[is_in_intersect], pts.labels[is_in_intersect],
        color=motif_color, markersize=15, marker=:circle, label=motif_label,
        strokewidth=1.5, strokecolor=RGBA(0.1,0.1,0.1,0.9))

    # Calculate and display R² if requested
    if show_r2
        y_true = pts.labels
        y_pred = pts.proc_prod
        ss_res = sum((y_true .- y_pred).^2)
        ss_tot = sum((y_true .- mean(y_true)).^2)
        r2 = 1 - ss_res / ss_tot
        
        # Add R² text to plot
        text!(ax, 0.98, 0.02, text=@sprintf("R² = %.3f", r2), 
            align=(:right, :bottom), fontsize=22, color=:black,
            space=:relative)
    end

    axislegend(ax, position=:lt, labelsize=28, framevisible=false, backgroundcolor=:white)
    return fig
end
