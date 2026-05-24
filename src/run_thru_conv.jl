
const motif_names = ["pairs", "triplets", "quadruplets", "quintuplets"]

# Collect extrema across all DataFrames without allocation
function obtain_xlim(contributions_df_filtered_singletons, dfs)
    xlim = mapreduce(extrema, 
        (a, b) -> (min(a[1], b[1]), max(a[2], b[2])),
        (contributions_df_filtered_singletons.banzhaf, 
            (df.banzhaf for df in dfs)...))

    return xlim
end


function plot_motifs_conv_case(data, m, motif_sizes, 
        contributions_df_filtered_singletons, dfs, pts, pts_test, all_indices;
        interaction_summaries=nothing,
        nav_page_count=4,
        enable_colored_borders = true,
        use_unified=true,
        dpi=65, 
        save_path="tmp", 
        page_title="n/a", 
        rna=false,
        sensitivity_analysis::Bool=false,
        dataset_name::Union{String,Nothing}=nothing
        );

    # ── Sensitivity-analysis-only mode (skip all rendering) ────
    if sensitivity_analysis
        ds_name = dataset_name === nothing ? page_title : dataset_name
        run_nnd_sensitivity_analysis(
            contributions_df_filtered_singletons, dfs, all_indices, pts, motif_sizes;
            save_path=save_path, page_title=ds_name
        )
        run_nnd_sensitivity_analysis_null(
            contributions_df_filtered_singletons, dfs, all_indices, pts, motif_sizes;
            save_path=save_path, page_title=ds_name
        )
        return nothing
    end

    # ── Full motif rendering pipeline ───────────────────────────
    xlim = obtain_xlim(contributions_df_filtered_singletons, dfs)

    config = ConvMotifConfig(data; 
        filter_len=m.hp.pfm_len, dpi=dpi, save_path=save_path, xlim=xlim)
        
    json_motifs = init_json_dict()
    html_dict = init_dict_for_html_render()

    next_idx, sorted_mapping = process_singletons!(
        contributions_df_filtered_singletons, all_indices, pts, config, json_motifs, html_dict; start_idx=1, rna=rna)

    # Remap filter indices in multi-motif DataFrames to use sorted order
    remap_filter_indices!(dfs, sorted_mapping, motif_sizes)
    
    # Remap interaction summaries to use sorted order
    remapped_interaction_summaries = remap_interaction_summaries(interaction_summaries, sorted_mapping, motif_sizes)

    group_ids = [motif_names[min(size-1, 4)] for size in motif_sizes]
    button_texts = ["$(size)-motifs" for size in motif_sizes]

    for (index, (motif_size, group_id, button_text)) in enumerate(zip(motif_sizes, group_ids, button_texts))
        @info "Processing multi-motifs of size: $(motif_size)"
        interaction_summary = remapped_interaction_summaries === nothing ? nothing : remapped_interaction_summaries[index]
        @time next_idx = process_multi_motifs!(dfs, all_indices, pts,
            config, json_motifs, html_dict;             
                interaction_summary=interaction_summary,
                motif_size=motif_size, group_id=group_id, 
                button_text=button_text, start_idx=next_idx, rna=rna                
                )
    end

    # Generate combined panel figure 
    data_pairs = [
        (pts_test.predictions, pts_test.labels, "Predictions", "Labels", "Predictions vs Labels"),
        (pts_test.proc_prod, pts_test.labels, "Learned Predictions", "Labels", "Learned Predictions vs Labels"),
        (pts_test.proc_prod, pts_test.predictions, "Learned Predictions", "Predictions", "Learned Predictions vs Predictions")
    ]
    
    publication_kde_panel(data_pairs, save_path=joinpath(save_path, "generalization.png"))

    # Render generalization page (index2.html)
    render_generalization_page!(save_path;
        page_title=page_title,
        nav_page_count=nav_page_count,
        image_filename="generalization.png"
    )

    # Render statistics (index3.html) and readme (index4.html) docs pages
    render_statistics_page!(save_path; page_title=page_title, nav_page_count=nav_page_count)
    render_readme_page!(save_path;     page_title=page_title, nav_page_count=nav_page_count)

    render_and_save_outputs!(json_motifs, html_dict, 1; 
        html_template=html_template_unified, 
        script_template=script_template,
        css_template=template_css,
        nav_page_count=nav_page_count,
        sequence_paths=[""],
        page_title=page_title,
        save_path=save_path, 
        enable_colored_borders = enable_colored_borders,
        use_unified=use_unified
        )
end