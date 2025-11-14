# Configuration struct for logo plotting
Base.@kwdef struct LogoConfig
    rotate::Bool = true  # Whether to rotate the logo 90° counterclockwise
    hide_decorations::Bool = true
    hide_spines::Bool = true
    preserve_aspect::Bool = true  # Whether to preserve the logo's aspect ratio
    row_span = 1:1  # Which rows to span in the GridLayout
    x_autolimit_margin::Tuple{Float64,Float64} = (0.0, 0.0)  # (left, right) margin for x-axis
    y_autolimit_margin::Tuple{Float64,Float64} = (0.0, 0.0)  # (top, bottom) margin for y-axis
end