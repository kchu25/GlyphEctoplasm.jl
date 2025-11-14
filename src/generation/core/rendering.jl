"""
HTML, JavaScript, and CSS rendering for motif visualization.
"""

"""
    render_and_save_outputs!(json_motifs, html_dict, j; kwargs...)

Render and write all visualization files: HTML, JavaScript, CSS, and JSON data.
Creates server-free static pages with embedded data.

# Parameters
- `json_motifs::Dict`: Dictionary containing motif data
- `html_dict`: HTML rendering dictionary
- `j::Integer`: Current page number
- `html_template`: Mustache template for HTML
- `script_template`: Mustache template for JavaScript
- `css_template`: Mustache template for CSS
- `save_path::AbstractString`: Directory to save output files
- `nav_page_count::Integer`: Number of navigation links to display at top of page (e.g., 4 for Pattern influence, Generalization, Readme, Statistics)
- `sequence_paths::Vector{String}`: Paths to sequence files
- `page_title::AbstractString = "n/a"`: Custom title for the page (displayed at top and in browser tab)
- `use_unified::Bool = false`: Use unified template that can mix singletons and multi-motifs
- `enable_colored_borders::Bool = true`: Apply dynamic colored borders to cards based on median values
"""
function render_and_save_outputs!(json_motifs::Dict, html_dict, j::Integer;
        html_template,
        script_template,
        css_template,
        save_path::AbstractString,
        nav_page_count::Integer,
        sequence_paths::Vector{String},
        page_title::AbstractString="n/a",
        use_unified::Bool=false,
        enable_colored_borders::Bool=true)
    
    # Ensure save_path directory exists
    mkpath(save_path)
    
    # Save JSON as JS file with const jsonData declaration
    open(joinpath(save_path, "data$j.js"), "w") do io
        print(io, "const jsonData = ")
        JSON3.pretty(io, json_motifs)
        println(io, ";")
    end
    
    # Copy CSS (once)
    css_path = joinpath(save_path, "styles.css")
    if !isfile(css_path)
        out = Mustache.render(css_template)
        open(css_path, "w") do io
            print(io, out)
        end
    end
    
    # Choose template: if use_unified is true, use html_template_unified
    # Otherwise use the provided html_template
    template_to_use = use_unified ? html_template_unified : html_template
    
    # Render HTML
    df = html_dict |> DataFrame
    html_rendered = Mustache.render(template_to_use, DF=df, protein_name=page_title, j=j)
    html_rendered = html_rendered * html_end
    
    # Render JavaScript with border coloring flag
    script_rendered = Mustache.render(script_template,
        mode_counts=size(df, 1), j=j, upto=nav_page_count, sequence_file_paths=sequence_paths)
    
    # Add global variable to control border coloring
    border_control_script = "\n// Control colored borders\nwindow.ENABLE_COLORED_BORDERS = $(enable_colored_borders ? "true" : "false");\n"
    script_rendered = border_control_script * script_rendered
    
    # Write HTML and JS files
    open(joinpath(save_path, "index$j.html"), "w") do io
        print(io, html_rendered)
    end
    open(joinpath(save_path, "scripts$j.js"), "w") do io
        print(io, script_rendered)
    end
end

export render_and_save_outputs!
