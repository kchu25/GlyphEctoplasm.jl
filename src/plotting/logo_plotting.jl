

###########################################################
# LOGO PLOTTING
###########################################################

"""
    plot_logo!(gl::GridLayout, logo_path::String; config::LogoConfig=LogoConfig())

Load and plot a logo image (PNG) in the provided GridLayout.

# Arguments
- `gl::GridLayout`: The grid layout to populate with the logo
- `logo_path::String`: Path to the logo image file

# Keyword Arguments
- `config::LogoConfig`: Configuration for logo display styling

# Returns
- The Axis containing the logo image
"""
function plot_logo!(gl::GridLayout, logo_path::String; config::LogoConfig=LogoConfig())
    # Create axis with optional aspect ratio preservation
    if config.preserve_aspect
        ax = Axis(gl[config.row_span, 1], aspect=DataAspect())
    else
        ax = Axis(gl[config.row_span, 1])
    end
    
    # Set axis margins to reduce whitespace
    ax.xautolimitmargin = config.x_autolimit_margin
    ax.yautolimitmargin = config.y_autolimit_margin
    
    # Load the logo image
    logo = load(logo_path)
    
    # Rotate if requested (similar to rain_cloud_plt.jl which uses rotr90)
    if config.rotate
        logo = rotr90(logo)
    end
    
    # Display the image
    image!(ax, logo)
    
    # Hide decorations
    if config.hide_spines
        hidespines!(ax)
    end
    if config.hide_decorations
        hidedecorations!(ax)
    end
    
    return ax
end
