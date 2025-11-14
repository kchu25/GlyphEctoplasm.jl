
"""
    process_singleton_motif!(json_mode, k, x_motifs, x_avg_contribs, rna)

Process singleton motif data and populate JSON structure.

# Arguments
- `json_mode`: JSON dictionary for this mode
- `k`: Motif key
- `x_motifs`: Motif data
- `x_avg_contribs`: Average contributions
- `rna`: Use RNA alphabet

# Generated Links
- PWM image
- Box plot with contribution
- MEME file download
- CSV for sequence highlighting
- Consensus sequence

# Example Output Structure
```julia
json_mode["pwms"] = ["path/to/pwm.png"]
json_mode["labels"] = ["pattern 1"]
json_mode["texts"] = [[
    "<a>contribution: 0.45</a>",
    "# occurrences: 1234",
    "<a>.meme file</a>",
    "<a>string highlight</a>",
    "ACGTACGT",
    ""
]]
```
"""
function process_singleton_motif!(json_mode, k, x_motifs, x_avg_contribs, rna)
    push!(json_mode[pwms_str], get_pwm_singleton_save_path(k))
    push!(json_mode[labels_str], "pattern $(k)")
    
    counts_here = (@view x_motifs[k][:, 1]) |> sum |> Int
    meme_link = "<a href=\"$(x_dir[1])/$k.meme\">.meme file</a>"
    consensus_str = Consensus.countmat2consensus(x_motifs[k]; rna = rna)
    csv_str = HtmlGeneration.fill_csv_link(get_csv_singleton_save_path(k))
    
    push!(
        json_mode[texts_str],
        [
            HtmlGeneration.fill_box_plot_kd_link(get_box_singleton_save_path(k), x_avg_contribs[k]),
            "# occurrences: $counts_here",
            meme_link,
            csv_str,
            consensus_str,
            ""
        ],
    )
end