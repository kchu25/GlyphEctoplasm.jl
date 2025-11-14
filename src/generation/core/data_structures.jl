"""
JSON and HTML data structure management for motif rendering.
"""

"""
    ensure_mode_entry!(json_motifs, mode_str)

Initialize mode entry in JSON dict if it doesn't exist. Returns the mode dict.
"""
function ensure_mode_entry!(json_motifs::Dict, mode_str::AbstractString)
    if !haskey(json_motifs, mode_str)
        json_motifs[mode_str] = Dict{String, Any}()
        json_motifs[mode_str][pwms_str] = String[]
        json_motifs[mode_str][labels_str] = String[]
        json_motifs[mode_str][texts_str] = Vector{Vector{String}}()
    end
    return json_motifs[mode_str]
end

"""
    add_motif_entry!(json_x_motifs, html_dict, mode_str, png_rel, label, texts, idx, filter_indices, median_val, group_id, button_text)

Add motif data to JSON dict and update HTML rendering dict (for primary motif display).
"""
function add_motif_entry!(json_x_motifs::Dict, html_dict, mode_str::AbstractString,
        png_rel::AbstractString, label::AbstractString, texts::Vector{String}, idx::Integer, 
        filter_indices::AbstractString, median_val::Real, group_id::AbstractString="", button_text::AbstractString="Singleton Motifs")
    ensure_mode_entry!(json_x_motifs, mode_str)
    push!(json_x_motifs[mode_str][pwms_str], png_rel)
    push!(json_x_motifs[mode_str][labels_str], label)
    push!(json_x_motifs[mode_str][texts_str], texts)
    populate_html_dict!(html_dict, idx, json_x_motifs[mode_str], filter_indices, median_val, group_id, button_text)
end

"""
    add_motif_variant!(json_motifs, mode_str, png_rel, label, texts)

Add motif variant to JSON dict without updating HTML (for multi-motif distance variants).
"""
function add_motif_variant!(json_motifs::Dict, mode_str::AbstractString,
        png_rel::AbstractString, label::AbstractString, texts::Vector{String})
    ensure_mode_entry!(json_motifs, mode_str)
    push!(json_motifs[mode_str][pwms_str], png_rel)
    push!(json_motifs[mode_str][labels_str], label)
    push!(json_motifs[mode_str][texts_str], texts)
end

export ensure_mode_entry!, add_motif_entry!, add_motif_variant!
