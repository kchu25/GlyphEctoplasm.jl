const motif_names = ["pairs", "triplets", "quadruplets", "quintuplets"]

function plot_motifs_conv_case(data, m, motif_sizes, 
        contributions_df_filtered, dfs; 
        dpi=65, save_path="tmp", xlim=(-2,2), 
        page_title="n/a"
        );

    config = ConvMotifConfig(data; 
        filter_len=m.hp.pfm_len, dpi=dpi, save_path=save_path, xlim=xlim)
        
    json_motifs = init_json_dict()
    html_dict = init_dict_for_html_render()

    next_idx = process_singletons!(
        contributions_df_filtered, config, json_motifs, html_dict; start_idx=1)

    group_ids = [motif_names[min(size-1, 4)] for size in motif_sizes]
    button_texts = ["$(size)-motifs" for size in motif_sizes]

    for (motif_size, group_id, button_text) in zip(motif_sizes, group_ids, button_texts)
        @info "Processing multi-motifs of size: $(motif_size)"
        next_idx = process_multi_motifs!(dfs, 
            config, json_motifs, html_dict; 
                motif_size=motif_size, group_id=group_id, 
                button_text=button_text, start_idx=next_idx
                )
    end
    
    render_and_save_outputs!(json_motifs, html_dict, 1; 
        html_template=html_template_unified, 
        script_template=script_template,
        css_template=template_css,
        nav_page_count=4,
        sequence_paths=[""],
        page_title=page_title,
        save_path=save_path, 
        enable_colored_borders = true,
        use_unified=true)
end
