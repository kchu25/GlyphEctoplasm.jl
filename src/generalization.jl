
"""Compute R² coefficient"""
function _compute_r2(y_true::AbstractVector{T}, y_pred::AbstractVector{T}) where T<:AbstractFloat
    ss_res = sum((y_true .- y_pred).^2)
    ss_tot = sum((y_true .- mean(y_true)).^2)
    return T(1) - ss_res / ss_tot
end

"""
    publication_scatter_panel(data_pairs; 
        figsize=(1200, 400), markersize=6, alpha=0.5, save_path=nothing)

Create a publication-quality multi-panel scatter plot with aligned axes.
`data_pairs` is a vector of (x, y, xlabel, ylabel, title) tuples.
"""
function publication_scatter_panel(data_pairs; 
        figsize=(1200, 400), markersize=6, alpha=0.5, save_path=nothing)
    
    n_panels = length(data_pairs)
    
    fig = Figure(size=figsize, fontsize=12)
    
    axes = Axis[]
    
    for (i, (x, y, xlabel, ylabel, title)) in enumerate(data_pairs)
        r2 = _compute_r2(y, x)
        
        ax = Axis(fig[1, i],
            xlabel=xlabel,
            ylabel=ylabel,
            title=title,
            titlesize=14,
            xlabelsize=12,
            ylabelsize=12,
            aspect=DataAspect()
        )
        push!(axes, ax)
        
        scatter!(ax, vec(x), vec(y), 
            markersize=markersize, 
            color=(:steelblue, alpha),
            strokewidth=0.5,
            strokecolor=:black
        )
        
        # Add identity line
        lims = (min(minimum(x), minimum(y)), max(maximum(x), maximum(y)))
        lines!(ax, [lims[1], lims[2]], [lims[1], lims[2]], 
            color=:red, linewidth=1.5, linestyle=:dash)
        
        # Add R² annotation
        text!(ax, 0.05, 0.95, 
            text=@sprintf("R² = %.3f", r2),
            space=:relative,
            fontsize=11,
            font=:bold
        )
    end
    
    # Link y-axes for alignment if they share the same ylabel
    linkyaxes!(axes...)
    
    # Adjust column gaps
    colgap!(fig.layout, 20)
    
    if !isnothing(save_path)
        save(save_path, fig, px_per_unit=3)
    end
    
    return fig
end
