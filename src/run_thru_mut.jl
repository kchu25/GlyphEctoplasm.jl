
function plot_motifs_mut_case(data, m,
        contributions_df_filtered, dfs;
        pts=nothing, pts_test=nothing, all_indices=nothing,
        interaction_summaries=nothing,
        nnd_k=15,
        dpi=65, save_path="tmp2", xlim=(-2,2),
        page_title="Mutation Regions Analysis",
        protein_name::Union{String,Nothing}=nothing,  # shown atop the summary page when given
        reduction_on_ref=false,
        float_type=Float32,
        use_rna=false,
        off_region_search=true,
        split_by_sign=true,
        sort_globally=true,
        sort_by_bins=true,      # default: single_region → Shapley-bin → cluster-bin → count
        bin_count=10,           # equal-width bins for Shapley median and cluster median
        sort_by_pareto=false,   # fallback when sort_by_bins=false (Pareto on |median|/count)
        nav_page_count=4,
        use_unified=true,
        enable_colored_borders=true,
        optimize_pngs::Bool=true,
        png_colors::Int=64,
        gc_every::Int=25       # incremental GC + Plots.closeall() every N motifs to cap RSS; 0 disables
    )
    # Population index basis for the per-motif indicator plots. Must share the
    # ordering of `pts` (same contract as the convolution case). Defaults to the
    # subset the mutation pipeline already filters on.
    all_idx = all_indices === nothing ? data.raw_data.most_common_length_indices : all_indices
    # Single configuration object for all mutation region analysis
    m_config = MutationRegionConfig(data;
        filter_len = m.receptive_field,
        float_type = float_type,
        use_rna = use_rna,
        off_region_search = off_region_search,
        xlim = xlim,
        save_path = save_path,
        dpi = dpi,
        reduction_on_ref = reduction_on_ref
    )

    all_metadata = prepare_and_collect_mutation_metadata(
        contributions_df_filtered, dfs, data, m_config;
        singleton_filter_pareto_rank = 1,
        split_by_sign = split_by_sign  # Splits by sign when computing Pareto ranks for singletons
    );

    json_motifs = init_json_dict()
    html_dict = init_dict_for_html_render()

    # Accumulates one self-contained row per motif for the top-movers landing page.
    top_movers = TopMoverEntry[]

    register_mutation_region_motifs!(
        json_motifs, html_dict, all_metadata;
        start_idx = 1,
        sort_globally = sort_globally,    # Enable hierarchical sorting
        sort_by_bins = sort_by_bins,       # Binned: single_region → Shapley-bin → cluster-bin → count
        bin_count = bin_count,             # Equal-width bins for the two medians
        sort_by_pareto = sort_by_pareto,   # Fallback Pareto ranking when sort_by_bins=false
        pts = pts,                         # Per-datapoint preds/labels for indicator plots (nothing => skip)
        all_indices = all_idx,             # Population basis aligned to `pts`
        interaction_summaries = interaction_summaries,  # Per-motif interaction text (string-valued, conv format)
        nnd_k = nnd_k,                     # k for the cluster-tightness (NND) permutation test
        top_movers_out = top_movers,       # Collect summary rows for the landing page
        gc_every = gc_every                # Periodic incremental GC to bound peak RSS on large runs
    )

    render_and_save_outputs!(json_motifs, html_dict, 1;
        html_template = html_template_unified,
        script_template = script_template,
        css_template = template_css,
        save_path = save_path,
        nav_page_count = nav_page_count,  # Show navigation for 4 pages: Motif influence, Generalization, Readme, Statistics
        sequence_paths = [""],
        page_title = page_title,
        use_unified = use_unified,
        enable_colored_borders = enable_colored_borders,
        has_summary = true                # Show the Summary (top-movers) nav link
    )

    # Generalization page (index2.html) — only when held-out predictions are
    # supplied. Mirrors the convolution case; closes the nav_page_count=4 gap.
    if pts_test !== nothing
        data_pairs = [
            (pts_test.predictions, pts_test.labels, "Predictions", "Labels", "Predictions vs Labels"),
            (pts_test.proc_prod, pts_test.labels, "Learned Predictions", "Labels", "Learned Predictions vs Labels"),
            (pts_test.proc_prod, pts_test.predictions, "Learned Predictions", "Predictions", "Learned Predictions vs Predictions")
        ]
        publication_kde_panel(data_pairs, save_path=joinpath(save_path, "generalization.png"))
        render_generalization_page!(save_path;
            page_title=page_title,
            nav_page_count=nav_page_count,
            image_filename="generalization.png",
            has_summary=true
        )
    end

    # Render statistics (index3.html) and readme (index4.html) docs pages
    render_statistics_page!(save_path; page_title=page_title, nav_page_count=nav_page_count, has_summary=true)
    render_readme_page!(save_path;     page_title=page_title, nav_page_count=nav_page_count, has_summary=true)

    # Top-movers summary becomes the landing page (index.html). Ranking is the
    # group-less binned lexicographic order; positives/negatives split by sign.
    positives, negatives = select_top_movers(top_movers; n=5, bin_count=bin_count)
    protein_length = try
        length(data.raw_data.consensus)
    catch
        nothing
    end
    wild_type = try
        String(data.raw_data.consensus)
    catch
        nothing
    end
    render_top_movers_page!(save_path;
        positives=positives, negatives=negatives,
        page_title=page_title, nav_page_count=nav_page_count,
        protein_name=protein_name, protein_length=protein_length,
        wild_type=wild_type
    )

    # Shrink emitted PNGs in place (filenames unchanged; HTML/JS references intact)
    optimize_pngs && optimize_pngs!(save_path; ncolors=png_colors)
end
