PDF_PLOTS = [
    "results_table",
    # "results_table_scaled",
    # "cartesian",
    # "mean_all_metrics_hbarplot",
    # "map_scores_barplot",
    # "all_metrics_barplot",
    # "full_panel",
]

plots_pattern = f"outputs/{scenario}/plots/{{plot}}.{{ext}}"
umap_baseline_pattern = f"outputs/{scenario}/projection/{{workflow}}_umap.parquet"
umap_pattern = f"outputs/{scenario}/projection/{{workflow}}_{{method}}_umap.parquet"



rule pivot_scores:
    input:
        path="outputs/{scenario}/metrics/" + config["preproc"] + "_all_metrics.parquet",
    output:
        path="outputs/{scenario}/metrics/" + config["preproc"] + "_all_metrics_pivoted.parquet",
        # dirty but then I don't have to change the plotting stuff
        # mock_path_a="outputs/{scenario}/plots/data/pivot_scores.parquet", 
        # mock_path_b="outputs/{scenario}/plots/data/tidy_scores.parquet", 
    run:
        plot.data.pivot_scores(input.path, output.path)


# rule prepare_embeddings:
#     input:
#         embd_files=expand(umap_pattern, workflow=WORKFLOWS, method=METHODS)
#         + expand(umap_baseline_pattern, workflow=WORKFLOWS),
#     output:
#         "outputs/{scenario}/plots/data/embeddings.parquet",
#     run:
#         plot.data.prepare_embeddings(input.embd_files, *output)


rule results_table:
    input:
        "outputs/{scenario}/metrics/" + config["preproc"] + "_all_metrics_pivoted.parquet",
    output:
        "outputs/{scenario}/plots/results_table.{ext}",
    run:
        plot.panel.results_table(*input, *output)


# rule results_table_scaled:
#     input:
#         "outputs/{scenario}/plots/data/pivot_scores.parquet",
#     output:
#         "outputs/{scenario}/plots/results_table_scaled.{ext}",
#     run:
#         plot.figures.results_table(*input, *output, min_max_scale=True)


# rule cartesian_plane:
#     input:
#         "outputs/{scenario}/plots/data/tidy_scores.parquet",
#     output:
#         "outputs/{scenario}/plots/cartesian.{ext}",
#     params:
#         min_cvar=0.01,
#     run:
#         plot.legacy.cartesian_plane(*input, *params, *output)


# rule barplot_all_metrics:
#     input:
#         "outputs/{scenario}/plots/data/tidy_scores.parquet",
#     output:
#         "outputs/{scenario}/plots/all_metrics_barplot.{ext}",
#     run:
#         plot.bar.all_metrics(*input, *output)


# rule barplot_map_scores:
#     input:
#         "outputs/{scenario}/plots/data/tidy_scores.parquet",
#     output:
#         "outputs/{scenario}/plots/map_scores_barplot.{ext}",
#     run:
#         plot.bar.map_scores(*input, *output)


# rule hbarplot_all_metrics:
#     input:
#         "outputs/{scenario}/plots/data/pivot_scores.parquet",
#     output:
#         "outputs/{scenario}/plots/mean_all_metrics_hbarplot.{ext}",
#     run:
#         plot.bar.all_metrics_h(*input, *output)


# rule umap_batch:
#     input:
#         "outputs/{scenario}/plots/data/embeddings.parquet",
#         "outputs/{scenario}/plots/data/pivot_scores.parquet",
#     output:
#         "outputs/{scenario}/plots/umap_batch.{ext}",
#     run:
#           plot.figures.umap_batch(*input, *output)


# rule umap_source:
#     input:
#         "outputs/{scenario}/plots/data/embeddings.parquet",
#         "outputs/{scenario}/plots/data/pivot_scores.parquet",
#     output:
#         "outputs/{scenario}/plots/umap_source.{ext}",
#     run:
#         plot.figures.umap_source(*input, *output)


# rule umap_compound:
#     input:
#         "outputs/{scenario}/plots/data/embeddings.parquet",
#         "outputs/{scenario}/plots/data/pivot_scores.parquet",
#     output:
#         "outputs/{scenario}/plots/umap_compound.{ext}",
#     run:
#         plot.figures.umap_compound(*input, *output)


# rule umap_microscope:
#     input:
#         "outputs/{scenario}/plots/data/embeddings.parquet",
#         "outputs/{scenario}/plots/data/pivot_scores.parquet",
#     output:
#         "outputs/{scenario}/plots/umap_microscope.{ext}",
#     run:
#         plot.figures.umap_microscope(*input, *output)


# rule full_panel:
#     input:
#         "outputs/{scenario}/plots/data/embeddings.parquet",
#         "outputs/{scenario}/plots/data/pivot_scores.parquet",
#     output:
#         "outputs/{scenario}/plots/full_panel.{ext}",
#     run:
#         plot.panel.full_panel(*input, *output, scenario)
