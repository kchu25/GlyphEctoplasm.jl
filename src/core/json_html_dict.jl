# Nicely organized initializer for the JSON/HTML motif structure

"""
init_json_dict(; keys=String[]) -> Dict{String,Dict{String,JSONValue}}

Create an empty JSON-like dictionary used by the renderer. The expected
shape is

  Dict{String, Dict{String, JSONValue}}

where JSONValue is one of:
  - String
  - Vector{String}
  - Vector{Vector{String}}

Arguments
- keys: optional vector of top-level keys (Strings) to pre-create.

Returns
- A Dict of the described shape, optionally pre-populated with the
  provided top-level keys.
"""

const JSONValue = Union{String, Vector{String}, Vector{Vector{String}}}
const JSONInner = Dict{String, JSONValue}

function init_json_dict(; keys::AbstractVector{<:AbstractString}=String[])
    json_x_motifs = Dict{String, JSONInner}()
    # pre-create top-level keys with empty inner dicts if requested
    for k in keys
        json_x_motifs[string(k)] = JSONInner()
    end
    return json_x_motifs
end

const HTMLValue = Vector{String}
const HTMLDict = Dict{String, HTMLValue}

"""
init_dict_for_html_render(; extra_keys=String[]) -> HTMLDict

Create and return a Dict mapping HTML tag identifiers (strings) to an empty
vector of strings. The returned dictionary has type `Dict{String, Vector{String}}`
so downstream code can push! strings into each entry.

Keyword arguments
- extra_keys: optional collection of additional top-level keys (strings) to
  include alongside the standard set defined by the module constants.

Examples
- init_dict_for_html_render()
- init_dict_for_html_render(extra_keys=["custom1","custom2"])
"""
function init_dict_for_html_render()
    default_keys = (
        tag_div_img_id,
        tag_i,
        tag_img_src,
        tag_img_alt,
        tag_div_text_id,
        tag_p_id1_default,
        tag_p_id2_default,
        tag_p_id3_default,
        tag_p_id4_default,
        tag_p_id5_default,
        tag_p_id6_default,
        tag_div_slide_id,
        tag_max_comb,
        tag_filter_indices,
        tag_median_val,
        tag_group_id,
        tag_button_text,
    )

    html_dict = HTMLDict()
    for k in default_keys
        html_dict[string(k)] = String[]
    end

    return html_dict
end

"""
    populate_html_dict!(html_dict, ind, json_mode, filter_indices, median_val, group_id, button_text)

Populate HTML dictionary with data from JSON mode.

Transfers the first entry from each JSON array to the corresponding
HTML tag field. This provides default values for initial page load.

# Arguments
- `html_dict`: HTML dictionary to populate
- `ind`: Mode index
- `json_mode`: JSON mode data
- `filter_indices`: Filter indices to display (e.g., "42" or "(8, 112)")
- `median_val`: Median Banzhaf value for dynamic styling
- `group_id`: Group identifier for grouping motifs (e.g., "g1", "g2")
- `button_text`: Custom text for the toggle button

# Populated Fields
- Image container ID and source
- Text container ID
- Six paragraph fields (contribution, count, links, etc.)
- Slide container ID
- Maximum combination index
- Filter indices overlay
- Median value for styling
- Group ID and button text

# Example
```julia
populate_html_dict!(html_dict, 1, json_mode, "42", 0.523, "g1", "High Confidence Motifs")
# html_dict["img_src"] now contains first PWM path
# html_dict["filter_indices"] contains "42"
# html_dict["median_val"] contains "0.523"
# html_dict["group_id"] contains "g1"
# html_dict["button_text"] contains "High Confidence Motifs"
```
"""
function populate_html_dict!(html_dict, idx, json_mode, filter_indices::AbstractString, median_val::Real, 
                            group_id::AbstractString="", button_text::AbstractString="Singleton Motifs")
    push!(html_dict[tag_i], "$idx")
    push!(html_dict[tag_div_img_id], "imageContainer$idx")
    push!(html_dict[tag_img_src], json_mode[pwms_str][1])
    push!(html_dict[tag_img_alt], json_mode[labels_str][1])
    push!(html_dict[tag_div_text_id], "textContainer$idx")
    push!(html_dict[tag_p_id1_default], json_mode[texts_str][1][1])
    push!(html_dict[tag_p_id2_default], json_mode[texts_str][1][2])
    push!(html_dict[tag_p_id3_default], json_mode[texts_str][1][3])
    push!(html_dict[tag_p_id4_default], json_mode[texts_str][1][4])
    push!(html_dict[tag_p_id5_default], json_mode[texts_str][1][5])
    push!(html_dict[tag_p_id6_default], json_mode[texts_str][1][6])
    push!(html_dict[tag_div_slide_id], "slideContainer$idx")
    push!(html_dict[tag_max_comb], "$(length(json_mode[pwms_str])-1)")
    push!(html_dict[tag_filter_indices], filter_indices)
    push!(html_dict[tag_median_val], string(median_val))
    push!(html_dict[tag_group_id], group_id)
    push!(html_dict[tag_button_text], button_text)
end