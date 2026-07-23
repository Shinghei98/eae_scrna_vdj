################################################################################
# Figure 1
#
# Figure 1a: external experimental-design schematic
# Figure 1b: cell-type marker heatmap
# Figure 1c: all-cell UMAP
# Figures 1d1-1d4: M1, M2, L1, and L2 all-cell UMAPs
# Figures 1e1-1e3: T-cell, macrophage/microglia, and DC UMAPs
# Figure 1f: broad cell-type composition by sample
# Figure 1g: scCODA broad cell-type composition
# Figure 1h: T-cell cluster composition heatmap and Fisher tests
################################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(ggplot2)
  library(cowplot)
  library(dplyr)
  library(tidyr)
  library(Matrix)
  library(grid)
  library(ComplexHeatmap)
  library(circlize)
})

################################################################################
# Shared Figure 1 paths and annotated object
################################################################################

base_dir <- "/Users/shingheimok/Desktop/Phd_Study/research/Wu_lab/b_cells"
figure1_dir <- file.path(base_dir, "figures", "figure_1")
dir.create(figure1_dir, recursive = TRUE, showWarnings = FALSE)

# Figure 1a is the external experimental-design schematic and is not generated
# by this R script.

rds_file <- file.path(base_dir, "srt_fullannot_with_bcell_dbscan_eps017_min5_min31.rds")
qs_file <- sub("\\.rds$", ".qs", rds_file)
output_dir <- file.path(base_dir, "outs", "output", "all_cells_umap_celltype_major_template_fit")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

obj <- if (file.exists(qs_file) && requireNamespace("qs", quietly = TRUE)) {
  qs::qread(qs_file)
} else {
  readRDS(rds_file)
}

################################################################################
# Figure 1h: T-cell cluster composition across M1, M2, L1, and L2
################################################################################

set.seed(1234)

figure1h_object_path <- file.path(
  base_dir,
  "tcell_annotation",
  "tcell_validated_obj_pc10_k30_res04_labeled.rds"
)
figure1h_output_dir <- file.path(
  base_dir,
  "outs",
  "output",
  "tcell_cluster_composition_pc10_k30_res04"
)
dir.create(figure1h_output_dir, recursive = TRUE, showWarnings = FALSE)

figure1h_cluster_labels <- c(
  C0 = "Naïve CD8.1",
  C1 = "Exhausted CD8",
  C2 = "Naïve CD8.2",
  C3 = "Naïve CD8.3",
  C4 = "Th1",
  C5 = "CD8 Treg",
  C6 = "Naïve CD8.4",
  C7 = "Unconventional T"
)
figure1h_second_layer_cluster <- c(
  C0 = 3,
  C1 = 2,
  C2 = 3,
  C3 = 3,
  C4 = 2,
  C5 = 2,
  C6 = 3,
  C7 = 4
)
figure1h_sample_order <- c("M1", "M2", "L1", "L2")
figure1h_cluster_order <- names(figure1h_cluster_labels)

if (!file.exists(figure1h_object_path)) {
  stop("Figure 1h T-cell object not found: ", figure1h_object_path)
}

figure1h_obj <- readRDS(figure1h_object_path)
figure1h_md <- figure1h_obj@meta.data
figure1h_required_columns <- c("numeric_cluster", "sample_short")
figure1h_missing_columns <- setdiff(
  figure1h_required_columns,
  colnames(figure1h_md)
)
if (length(figure1h_missing_columns) > 0) {
  stop(
    "Figure 1h is missing required metadata columns: ",
    paste(figure1h_missing_columns, collapse = ", ")
  )
}

figure1h_cell_metadata <- figure1h_md |>
  transmute(
    cell = rownames(figure1h_md),
    sample = as.character(sample_short),
    cluster = as.character(numeric_cluster)
  ) |>
  filter(
    sample %in% figure1h_sample_order,
    cluster %in% figure1h_cluster_order
  )

if (nrow(figure1h_cell_metadata) == 0) {
  stop("No cells matched the Figure 1h samples and clusters.")
}
if (!setequal(
  unique(figure1h_cell_metadata$cluster),
  figure1h_cluster_order
)) {
  stop("Figure 1h observed clusters do not exactly match C0-C7.")
}

figure1h_mapping_check <- figure1h_cell_metadata |>
  count(cluster, name = "n_cells") |>
  mutate(
    label = unname(figure1h_cluster_labels[cluster]),
    second_layer_cluster = unname(
      figure1h_second_layer_cluster[cluster]
    )
  ) |>
  arrange(match(cluster, figure1h_cluster_order))
write.csv(
  figure1h_mapping_check,
  file.path(
    figure1h_output_dir,
    "tcell_cluster_label_mapping_check.csv"
  ),
  row.names = FALSE,
  quote = TRUE
)

if ("annotation_label" %in% colnames(figure1h_md)) {
  figure1h_existing_annotation_check <- figure1h_md |>
    transmute(
      cluster = as.character(numeric_cluster),
      existing_annotation_label = as.character(annotation_label)
    ) |>
    count(cluster, existing_annotation_label, name = "n_cells") |>
    filter(n_cells > 0) |>
    arrange(
      match(cluster, figure1h_cluster_order),
      existing_annotation_label
    )
  write.csv(
    figure1h_existing_annotation_check,
    file.path(
      figure1h_output_dir,
      "tcell_cluster_existing_annotation_check.csv"
    ),
    row.names = FALSE,
    quote = TRUE
  )
}

figure1h_sample_totals <- figure1h_cell_metadata |>
  count(sample, name = "total_t_cells") |>
  complete(
    sample = figure1h_sample_order,
    fill = list(total_t_cells = 0)
  )

figure1h_composition <- figure1h_cell_metadata |>
  count(sample, cluster, name = "n_cells") |>
  complete(
    sample = figure1h_sample_order,
    cluster = figure1h_cluster_order,
    fill = list(n_cells = 0)
  ) |>
  left_join(figure1h_sample_totals, by = "sample") |>
  mutate(
    cluster_label = unname(figure1h_cluster_labels[cluster]),
    proportion = if_else(
      total_t_cells > 0,
      n_cells / total_t_cells,
      NA_real_
    ),
    arithmetic_mean_proportion = proportion
  ) |>
  arrange(
    match(cluster, figure1h_cluster_order),
    match(sample, figure1h_sample_order)
  ) |>
  group_by(cluster) |>
  mutate(
    cluster_min_proportion = min(proportion, na.rm = TRUE),
    cluster_max_proportion = max(proportion, na.rm = TRUE),
    relative_proportion_scaled_0_1 = if_else(
      cluster_max_proportion > cluster_min_proportion,
      (proportion - cluster_min_proportion) /
        (cluster_max_proportion - cluster_min_proportion),
      0
    ),
    enriched_sample = sample[which.max(proportion)]
  ) |>
  ungroup()

write.csv(
  figure1h_composition,
  file.path(
    figure1h_output_dir,
    "tcell_cluster_composition_by_sample.csv"
  ),
  row.names = FALSE,
  quote = TRUE
)

figure1h_enrichment <- figure1h_composition |>
  group_by(cluster, cluster_label, enriched_sample) |>
  summarise(
    mean_proportion_across_samples = mean(proportion, na.rm = TRUE),
    max_proportion = max(proportion, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(match(cluster, figure1h_cluster_order))
write.csv(
  figure1h_enrichment,
  file.path(
    figure1h_output_dir,
    "tcell_cluster_enrichment_summary.csv"
  ),
  row.names = FALSE,
  quote = TRUE
)

figure1h_m1_m2_enrichment <- figure1h_composition |>
  group_by(cluster, cluster_label) |>
  summarise(
    mean_M1_M2_proportion = mean(
      proportion[sample %in% c("M1", "M2")],
      na.rm = TRUE
    ),
    mean_L1_L2_proportion = mean(
      proportion[sample %in% c("L1", "L2")],
      na.rm = TRUE
    ),
    M1_M2_enriched = mean_M1_M2_proportion >
      mean_L1_L2_proportion,
    .groups = "drop"
  ) |>
  arrange(match(cluster, figure1h_cluster_order))
write.csv(
  figure1h_m1_m2_enrichment,
  file.path(
    figure1h_output_dir,
    "tcell_M1_M2_enrichment_summary.csv"
  ),
  row.names = FALSE,
  quote = TRUE
)

figure1h_proportion_matrix <- figure1h_composition |>
  select(cluster, sample, proportion) |>
  pivot_wider(names_from = sample, values_from = proportion) |>
  arrange(match(cluster, figure1h_cluster_order))
figure1h_scaled_matrix <- figure1h_composition |>
  select(cluster, sample, relative_proportion_scaled_0_1) |>
  pivot_wider(
    names_from = sample,
    values_from = relative_proportion_scaled_0_1
  ) |>
  arrange(match(cluster, figure1h_cluster_order))

write.csv(
  figure1h_proportion_matrix,
  file.path(
    figure1h_output_dir,
    "tcell_cluster_proportion_matrix.csv"
  ),
  row.names = FALSE,
  quote = TRUE
)
write.csv(
  figure1h_scaled_matrix,
  file.path(
    figure1h_output_dir,
    "tcell_cluster_relative_proportion_scaled_matrix.csv"
  ),
  row.names = FALSE,
  quote = TRUE
)

figure1h_scaled_mat <- t(as.matrix(
  figure1h_scaled_matrix[
    ,
    figure1h_sample_order,
    drop = FALSE
  ]
))
rownames(figure1h_scaled_mat) <- figure1h_sample_order
colnames(figure1h_scaled_mat) <- unname(
  figure1h_cluster_labels[figure1h_scaled_matrix$cluster]
)
write.csv(
  as.data.frame(figure1h_scaled_mat),
  file.path(
    figure1h_output_dir,
    "tcell_cluster_heatmap_matrix_sample_rows.csv"
  ),
  quote = TRUE
)

figure1h_colors <- colorRamp2(
  c(0, 0.2, 0.4, 0.6, 0.8, 1),
  c(
    "#440154",
    "#3B528B",
    "#21918C",
    "#5EC962",
    "#B8DE29",
    "#FDE725"
  )
)
figure1h_column_branch_ids <- factor(
  unname(
    figure1h_second_layer_cluster[
      figure1h_scaled_matrix$cluster
    ]
  ),
  levels = c(2, 3, 4)
)
figure1h_row_blocks <- factor(
  c("M1-M2", "M1-M2", "L1-L2", "L1-L2"),
  levels = c("M1-M2", "L1-L2")
)

make_figure1h_heatmap <- function() {
  Heatmap(
    figure1h_scaled_mat,
    name = "Relative proportion",
    col = figure1h_colors,
    cluster_rows = FALSE,
    row_split = figure1h_row_blocks,
    cluster_columns = TRUE,
    column_split = figure1h_column_branch_ids,
    cluster_column_slices = TRUE,
    column_gap = unit(2, "mm"),
    row_gap = unit(2.5, "mm"),
    row_title = NULL,
    column_title = NULL,
    show_parent_dend_line = FALSE,
    column_dend_height = unit(25, "mm"),
    height = unit(55, "mm"),
    width = unit(150, "mm"),
    show_row_dend = FALSE,
    show_column_dend = TRUE,
    row_names_side = "right",
    row_names_gp = gpar(fontsize = 11, fontface = "bold"),
    column_names_side = "bottom",
    column_names_rot = 90,
    column_names_gp = gpar(fontsize = 10, fontface = "bold"),
    rect_gp = gpar(col = "white", lwd = 1),
    show_heatmap_legend = FALSE
  )
}

figure1h_legend <- Legend(
  col_fun = figure1h_colors,
  title = "Proportion of cluster\nin each sample",
  title_position = "topcenter",
  direction = "horizontal",
  legend_width = unit(48, "mm"),
  legend_height = unit(5, "mm"),
  at = c(0, 0.2, 0.4, 0.6, 0.8, 1),
  labels = c("0", "0.2", "0.4", "0.6", "0.8", "1"),
  labels_gp = gpar(fontsize = 10, fontface = "bold"),
  title_gp = gpar(fontsize = 13, fontface = "bold"),
  border = "black"
)

figure1h_tiff <- file.path(
  figure1_dir,
  "figure_1h_tcell_cluster_composition_heatmap.tiff"
)
tiff(
  figure1h_tiff,
  width = 2400,
  height = 1800,
  res = 300,
  compression = "lzw",
  bg = "white"
)
draw(
  make_figure1h_heatmap(),
  heatmap_legend_list = list(figure1h_legend),
  heatmap_legend_side = "bottom",
  padding = unit(c(5, 8, 5, 15), "mm")
)
invisible(dev.off())

################################################################################
# Figure 1h statistical follow-up
################################################################################

# Each two-sided Fisher's exact test compares the target sample (M1 or M2)
# against pooled L1+L2, using the target cluster versus all other T cells as
# the 2 x 2 contingency table.
figure1h_fisher_clusters <- c(
  C1 = "Exhausted CD8",
  C4 = "Th1",
  C5 = "CD8 Treg"
)
figure1h_fisher_samples <- c("M1", "M2")

figure1h_fisher_results <- do.call(
  rbind,
  lapply(names(figure1h_fisher_clusters), function(cluster_i) {
    do.call(
      rbind,
      lapply(figure1h_fisher_samples, function(sample_i) {
        target_cluster <- sum(
          figure1h_cell_metadata$sample == sample_i &
            figure1h_cell_metadata$cluster == cluster_i
        )
        target_other <- sum(
          figure1h_cell_metadata$sample == sample_i &
            figure1h_cell_metadata$cluster != cluster_i
        )
        reference_cluster <- sum(
          figure1h_cell_metadata$sample %in% c("L1", "L2") &
            figure1h_cell_metadata$cluster == cluster_i
        )
        reference_other <- sum(
          figure1h_cell_metadata$sample %in% c("L1", "L2") &
            figure1h_cell_metadata$cluster != cluster_i
        )
        fisher_result <- fisher.test(
          matrix(
            c(
              target_cluster,
              target_other,
              reference_cluster,
              reference_other
            ),
            nrow = 2,
            byrow = TRUE,
            dimnames = list(
              c(sample_i, "L1+L2"),
              c(cluster_i, "all_other_T_cells")
            )
          ),
          alternative = "two.sided"
        )
        data.frame(
          cluster = cluster_i,
          cluster_label = unname(
            figure1h_fisher_clusters[cluster_i]
          ),
          target_sample = sample_i,
          target_cluster_cells = target_cluster,
          target_other_T_cells = target_other,
          reference_cluster_cells = reference_cluster,
          reference_other_T_cells = reference_other,
          odds_ratio = unname(fisher_result$estimate),
          p_value = fisher_result$p.value,
          row.names = NULL
        )
      })
    )
  })
)
figure1h_fisher_results <- figure1h_fisher_results[
  order(
    match(
      figure1h_fisher_results$cluster,
      names(figure1h_fisher_clusters)
    ),
    figure1h_fisher_results$target_sample
  ),
]
write.csv(
  figure1h_fisher_results,
  file.path(
    figure1h_output_dir,
    "tcell_M1_M2_fisher_exact_results.csv"
  ),
  row.names = FALSE,
  quote = TRUE
)

# Manuscript finding from the frozen Figure 1h object:
# Th1 was enriched in M1 (p = 9.00 x 10^-146) and M2
# (p = 2.05 x 10^-138). Exhausted CD8 and CD8 Treg were enriched in M1
# (p = 8.77 x 10^-96 and p = 8.71 x 10^-9) and M2
# (p = 8.04 x 10^-75 and p = 4.59 x 10^-4), respectively; all tests were
# two-sided Fisher's exact tests.

cat(
  "Saved Figure 1h:\n",
  figure1h_tiff,
  "\nCells analyzed: ",
  nrow(figure1h_cell_metadata),
  "\n",
  sep = ""
)
rm(figure1h_obj, figure1h_md)
invisible(gc())

required_meta_cols <- c("celltype_major", "celltype_minor", "sample_id")
missing_meta_cols <- setdiff(required_meta_cols, colnames(obj@meta.data))
if (length(missing_meta_cols) > 0) {
  stop("Missing metadata columns: ", paste(missing_meta_cols, collapse = ", "))
}
if (!"umap" %in% Reductions(obj)) {
  stop("The object does not contain a UMAP reduction named 'umap'.")
}

################################################################################
# Figure 1b: cell-type marker heatmap
################################################################################

figure1b_output_dir <- file.path(
  base_dir,
  "outs",
  "output",
  "global_marker_heatmap_reference_style"
)
dir.create(figure1b_output_dir, recursive = TRUE, showWarnings = FALSE)

figure1b_tiff <- file.path(figure1_dir, "figure_1b_celltype_marker_heatmap.tiff")
figure1b_csv <- file.path(
  figure1b_output_dir,
  "figure_1b_global_celltype_marker_heatmap_values.csv"
)

figure1b_row_levels <- c(
  "B cell",
  "Th1",
  "Naive CD8 T cell",
  "Exhausted CD8 T cell",
  "CD8 Treg",
  "Unconventional T cell",
  "Macrophage",
  "DAM1",
  "DAM2",
  "Monocyte",
  "Migratory DC",
  "Non-migratory DC",
  "Neutrophil",
  "Cycling Neutrophil",
  "NK cell"
)

figure1b_genes <- c(
  "Ptprc", "Ms4a1", "Cd3d", "Cd4", "Ifng", "Cd8a", "Lef1", "Pdcd1",
  "Foxp3", "Ikzf2", "Cx3cr1", "Csf1r", "C1qa", "Marco", "Folr2",
  "Ccr7", "Fscn1", "Ccl22", "Xcr1", "Clec9a", "Irf8", "Batf3",
  "Sirpa", "Klf4", "S100a9", "Csf3r", "Mki67", "Top2a", "Klrb1c",
  "Ncr1"
)

figure1b_group <- rep(NA_character_, ncol(obj))
names(figure1b_group) <- colnames(obj)
figure1b_major <- as.character(obj$celltype_major)
figure1b_minor <- as.character(obj$celltype_minor)

figure1b_group[figure1b_major == "B_cell"] <- "B cell"
figure1b_group[figure1b_major == "CD4_T" & figure1b_minor == "Th1"] <- "Th1"
figure1b_group[
  figure1b_major == "CD8_T" & figure1b_minor == "naive_CD8_T_cell"
] <- "Naive CD8 T cell"
figure1b_group[
  figure1b_major == "CD8_T" & figure1b_minor == "exhausted_CD8_T_cell"
] <- "Exhausted CD8 T cell"
figure1b_group[
  figure1b_major == "CD8_T" & figure1b_minor == "CD8_Treg"
] <- "CD8 Treg"
figure1b_group[
  figure1b_major == "UNC_T" & figure1b_minor == "unconventional_T_cell"
] <- "Unconventional T cell"
figure1b_group[
  figure1b_major == "macrophage/microglia" & figure1b_minor == "macrophage"
] <- "Macrophage"
figure1b_group[
  figure1b_major == "macrophage/microglia" & figure1b_minor == "DAM1"
] <- "DAM1"
figure1b_group[
  figure1b_major == "macrophage/microglia" & figure1b_minor == "DAM2"
] <- "DAM2"
figure1b_group[
  figure1b_major == "monocyte" & figure1b_minor == "monocyte"
] <- "Monocyte"
figure1b_group[
  figure1b_major == "DC" & figure1b_minor == "migratory_DC"
] <- "Migratory DC"
figure1b_group[
  figure1b_major == "DC" & figure1b_minor == "non_migratory_DC"
] <- "Non-migratory DC"
figure1b_group[
  figure1b_major == "Neutrophil" & figure1b_minor == "Neutrophil"
] <- "Neutrophil"
figure1b_group[
  figure1b_major == "Neutrophil" & figure1b_minor == "cycling_neutrophil"
] <- "Cycling Neutrophil"
figure1b_group[
  figure1b_major == "NK_cell" & figure1b_minor == "NK_cell"
] <- "NK cell"

figure1b_cells <- names(figure1b_group)[!is.na(figure1b_group)]
figure1b_genes_present <- figure1b_genes[figure1b_genes %in% rownames(obj)]

if (length(figure1b_genes_present) != length(figure1b_genes)) {
  stop(
    "Figure 1b is missing marker genes: ",
    paste(setdiff(figure1b_genes, figure1b_genes_present), collapse = ", ")
  )
}
if (!setequal(unique(figure1b_group[figure1b_cells]), figure1b_row_levels)) {
  stop("Figure 1b annotated groups do not exactly match the expected 15 groups.")
}

figure1b_expression <- tryCatch(
  GetAssayData(obj, assay = "RNA", layer = "data"),
  error = function(e) SeuratObject::LayerData(obj, assay = "RNA", layer = "data")
)
figure1b_expression <- figure1b_expression[
  figure1b_genes_present,
  figure1b_cells,
  drop = FALSE
]
figure1b_cells_by_group <- split(
  figure1b_cells,
  factor(
    figure1b_group[figure1b_cells],
    levels = figure1b_row_levels
  )
)

figure1b_average_expression <- sapply(figure1b_row_levels, function(group_i) {
  Matrix::rowMeans(
    figure1b_expression[, figure1b_cells_by_group[[group_i]], drop = FALSE]
  )
})
rownames(figure1b_average_expression) <- figure1b_genes_present
colnames(figure1b_average_expression) <- figure1b_row_levels

figure1b_minmax01 <- function(x) {
  value_range <- range(x, na.rm = TRUE)
  if (!all(is.finite(value_range)) || diff(value_range) == 0) {
    return(rep(0, length(x)))
  }
  (x - value_range[1]) / diff(value_range)
}

figure1b_normalized_expression <- t(apply(
  figure1b_average_expression,
  1,
  figure1b_minmax01
))
rownames(figure1b_normalized_expression) <- figure1b_genes_present
colnames(figure1b_normalized_expression) <- figure1b_row_levels

figure1b_plot_df <- expand.grid(
  gene = figure1b_genes_present,
  celltype = figure1b_row_levels,
  stringsAsFactors = FALSE
) |>
  mutate(
    avg_expression = mapply(
      function(gene_i, celltype_i) {
        figure1b_average_expression[gene_i, celltype_i]
      },
      gene,
      celltype
    ),
    normalized_expression = mapply(
      function(gene_i, celltype_i) {
        figure1b_normalized_expression[gene_i, celltype_i]
      },
      gene,
      celltype
    ),
    gene = factor(gene, levels = figure1b_genes_present),
    celltype = factor(celltype, levels = rev(figure1b_row_levels))
  )

write.csv(figure1b_plot_df, figure1b_csv, row.names = FALSE)

figure1b_plot <- ggplot(
  figure1b_plot_df,
  aes(x = gene, y = celltype, fill = normalized_expression)
) +
  geom_tile(color = "#321143", linewidth = 0.20) +
  scale_fill_gradientn(
    colours = c(
      "#440154",
      "#414487",
      "#2A788E",
      "#22A884",
      "#7AD151",
      "#FDE725"
    ),
    values = scales::rescale(c(0, 0.2, 0.4, 0.6, 0.8, 1)),
    limits = c(0, 1),
    name = "Normalized\nexpression"
  ) +
  scale_x_discrete(position = "bottom", expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
    axis.text.x = element_text(
      angle = 55,
      hjust = 1,
      vjust = 1,
      face = "bold",
      color = "black",
      size = 13
    ),
    axis.text.y = element_text(
      face = "bold",
      color = "black",
      size = 13
    ),
    axis.ticks = element_line(color = "black", linewidth = 0.7),
    axis.ticks.length = unit(3.2, "pt"),
    plot.margin = margin(18, 22, 24, 22),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 15),
    legend.text = element_text(size = 13),
    legend.key.height = unit(1.8, "in"),
    legend.key.width = unit(0.30, "in")
  ) +
  guides(
    fill = guide_colorbar(
      barheight = unit(1.9, "in"),
      barwidth = unit(0.35, "in"),
      frame.colour = "black",
      ticks.colour = "black"
    )
  )

ggsave(
  filename = figure1b_tiff,
  plot = figure1b_plot,
  device = "tiff",
  width = 15,
  height = 7.3,
  units = "in",
  dpi = 300,
  compression = "lzw",
  bg = "white"
)

cat("Saved Figure 1b:\n", figure1b_tiff, "\n", figure1b_csv, "\n", sep = "")

rm(
  figure1b_expression,
  figure1b_average_expression,
  figure1b_normalized_expression
)
invisible(gc())

################################################################################
# Figure 1c: all-cell UMAP colored by major cell type
################################################################################

umap_df <- as.data.frame(Embeddings(obj, reduction = "umap"))
colnames(umap_df)[1:2] <- c("UMAP_1", "UMAP_2")
umap_df$cell <- rownames(umap_df)
umap_df$celltype_major <- as.character(obj$celltype_major)
umap_df$celltype_minor <- as.character(obj$celltype_minor)
umap_df$sample_id <- obj@meta.data[umap_df$cell, "sample_id", drop = TRUE]
umap_df$sample_panel <- dplyr::recode(
  umap_df$sample_id,
  "MGI0279_1_Wu2020_M1-lib1" = "M1",
  "MGI0279_1_Wu2020_M2-lib1" = "M2",
  "MGI0279_1_Wu2020_L1-lib1" = "L1",
  "MGI0279_1_Wu2020_L2-lib1" = "L2",
  .default = as.character(umap_df$sample_id)
)

umap_df$plot_celltype <- dplyr::case_when(
  umap_df$celltype_major %in% c("CD4_T", "CD8_T", "UNC_T", "T_cell") ~ "T cells",
  umap_df$celltype_major == "B_cell" ~ "B cells",
  umap_df$celltype_major == "Neutrophil" & umap_df$celltype_minor == "cycling_neutrophil" ~ "Cycling Neutrophil",
  umap_df$celltype_major == "Neutrophil" ~ "Neutrophil",
  umap_df$celltype_major == "macrophage/microglia" ~ "Macrophage/microglia",
  umap_df$celltype_major == "DC" ~ "DC",
  umap_df$celltype_major == "NK_cell" ~ "NK cell",
  umap_df$celltype_major == "monocyte" ~ "Monocyte",
  is.na(umap_df$celltype_major) ~ "NA",
  TRUE ~ gsub("_", " ", umap_df$celltype_major)
)

celltype_levels <- c(
  "B cells",
  "T cells",
  "Neutrophil",
  "Cycling Neutrophil",
  "Monocyte",
  "Macrophage/microglia",
  "DC",
  "NK cell",
  "NA"
)
celltype_levels <- c(
  celltype_levels[celltype_levels %in% unique(umap_df$plot_celltype)],
  sort(setdiff(unique(umap_df$plot_celltype), celltype_levels))
)
umap_df$plot_celltype <- factor(umap_df$plot_celltype, levels = celltype_levels)

celltype_colors <- c(
  "B cells" = "#159DC1",
  "T cells" = "#C8192E",
  "Neutrophil" = "#D18A00",
  "Cycling Neutrophil" = "#E6C229",
  "Macrophage/microglia" = "#5E49A8",
  "DC" = "#1F68B3",
  "NK cell" = "#3F9D43",
  "Monocyte" = "#9D4D1F",
  "NA" = "#D0D0D0"
)

missing_colors <- setdiff(levels(umap_df$plot_celltype), names(celltype_colors))
if (length(missing_colors) > 0) {
  fallback_colors <- setNames(
    grDevices::hcl.colors(length(missing_colors), palette = "Dark 3"),
    missing_colors
  )
  celltype_colors <- c(celltype_colors, fallback_colors)
}
celltype_colors <- celltype_colors[levels(umap_df$plot_celltype)]

count_df <- umap_df |>
  count(plot_celltype, name = "n_cells") |>
  mutate(color = celltype_colors[as.character(plot_celltype)])
write.csv(
  count_df,
  file.path(output_dir, "all_cells_celltype_major_tcells_merged_counts_and_palette.csv"),
  row.names = FALSE
)

plot_df <- umap_df |>
  left_join(count_df |> select(plot_celltype, n_cells), by = "plot_celltype") |>
  arrange(desc(n_cells))

x_range <- range(plot_df$UMAP_1, na.rm = TRUE)
y_range <- range(plot_df$UMAP_2, na.rm = TRUE)
x_span <- diff(x_range)
y_span <- diff(y_range)

axis_x0 <- x_range[1] - 0.180 * x_span
axis_y0 <- y_range[1] - 0.085 * y_span
axis_x1 <- axis_x0 + 0.285 * x_span
axis_y1 <- axis_y0 + 0.240 * y_span
axis_df <- data.frame(
  x = c(axis_x1, axis_x0, axis_x0),
  y = c(axis_y0, axis_y0, axis_y1)
)

p_umap <- ggplot(plot_df, aes(UMAP_1, UMAP_2, color = plot_celltype)) +
  geom_point(size = 0.42, alpha = 1, stroke = 0) +
  scale_color_manual(values = celltype_colors, drop = FALSE) +
  coord_fixed(
    xlim = c(axis_x0 - 0.050 * x_span, x_range[2] + 0.035 * x_span),
    ylim = c(axis_y0 - 0.075 * y_span, y_range[2] + 0.030 * y_span),
    ratio = 1,
    expand = FALSE,
    clip = "off"
  ) +
  theme_void(base_size = 12) +
  theme(
    legend.position = "none",
    plot.margin = margin(0, 4, 3, 0)
  )

legend_direct_df <- tibble(
  plot_celltype = levels(plot_df$plot_celltype),
  label = levels(plot_df$plot_celltype),
  x_dot = 0.595,
  x_label = 0.636,
  y = 0.691 - (seq_along(levels(plot_df$plot_celltype)) - 1) * 0.052
)

final_plot <- ggdraw() +
  draw_plot(p_umap, x = -0.071, y = 0.043, width = 0.698, height = 0.829) +
  draw_line(
    c(0.179, 0.248),
    c(0.911, 0.911),
    color = "black",
    linewidth = 1.55,
    lineend = "butt"
  ) +
  draw_label(
    "All cells",
    x = 0.329,
    y = 0.911,
    hjust = 0.5,
    vjust = 0.5,
    fontface = "bold",
    size = 23.5
  ) +
  draw_line(
    c(0.409, 0.479),
    c(0.911, 0.911),
    color = "black",
    linewidth = 1.55,
    lineend = "butt"
  ) +
  draw_line(
    c(0.059, 0.059, 0.147),
    c(0.231, 0.105, 0.105),
    color = "black",
    linewidth = 1.35,
    lineend = "butt",
    linejoin = "mitre"
  ) +
  draw_label(
    "UMAP1",
    x = 0.103,
    y = 0.078,
    hjust = 0.5,
    vjust = 0.5,
    fontface = "bold",
    size = 15.5
  ) +
  draw_label(
    "UMAP2",
    x = 0.043,
    y = 0.166,
    angle = 90,
    hjust = 0.5,
    vjust = 0.5,
    fontface = "bold",
    size = 15.5
  ) +
  draw_line(
    c(0.611, 0.551, 0.551, 0.941, 0.941, 0.881),
    c(0.768, 0.768, 0.268, 0.268, 0.768, 0.768),
    color = "black",
    linewidth = 1.75,
    lineend = "butt",
    linejoin = "mitre"
  ) +
  draw_label(
    "Legend c (All cells)",
    x = 0.747,
    y = 0.768,
    hjust = 0.5,
    vjust = 0.5,
    fontface = "bold",
    size = 17.5
  ) +
  geom_point(
    data = legend_direct_df,
    aes(x = x_dot, y = y, fill = plot_celltype),
    inherit.aes = FALSE,
    shape = 21,
    size = 5.3,
    stroke = 0.12,
    color = "white",
    show.legend = FALSE
  ) +
  geom_text(
    data = legend_direct_df,
    aes(x = x_label, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 0.5,
    size = 5.65,
    fontface = "bold",
    color = "black"
  ) +
  scale_fill_manual(values = celltype_colors, drop = FALSE)

out_tiff <- file.path(figure1_dir, "figure_1c_all_cells_umap.tiff")

ggsave(
  filename = out_tiff,
  plot = final_plot,
  width = 9.14,
  height = 5.50,
  dpi = 300,
  bg = "white",
  compression = "lzw"
)

cat("Cells plotted:", nrow(plot_df), "\n")
cat("Legend groups:", length(levels(plot_df$plot_celltype)), "\n")
cat("Saved:\n")
cat(out_tiff, "\n")
cat(file.path(output_dir, "all_cells_celltype_major_tcells_merged_counts_and_palette.csv"), "\n")

################################################################################
# Figure 1f: broad cell-type composition by sample
################################################################################

composition_samples <- c("M1", "M2", "L1", "L2")
composition_grid <- expand.grid(
  sample_panel = composition_samples,
  plot_celltype = levels(umap_df$plot_celltype),
  stringsAsFactors = FALSE
)

composition_counts <- umap_df |>
  filter(sample_panel %in% composition_samples) |>
  count(sample_panel, plot_celltype, name = "n_cells") |>
  mutate(plot_celltype = as.character(plot_celltype))

composition_df <- composition_grid |>
  left_join(
    composition_counts,
    by = c("sample_panel", "plot_celltype")
  ) |>
  mutate(n_cells = ifelse(is.na(n_cells), 0L, n_cells)) |>
  group_by(sample_panel) |>
  mutate(
    total_cells = sum(n_cells),
    proportion = n_cells / total_cells,
    tissue_panel = ifelse(
      sample_panel %in% c("M1", "M2"),
      "Meninges",
      "dCLN"
    ),
    plot_celltype = factor(
      plot_celltype,
      levels = levels(umap_df$plot_celltype)
    )
  ) |>
  ungroup()

figure1f_composition_csv <- file.path(
  output_dir,
  "figure_1f_celltype_composition_by_sample.csv"
)
write.csv(composition_df, figure1f_composition_csv, row.names = FALSE)

make_composition_panel <- function(
    df,
    sample_levels,
    panel_title,
    show_y_title = TRUE) {
  df <- df |>
    mutate(sample_panel = factor(sample_panel, levels = sample_levels))

  panel <- ggplot(
    df,
    aes(x = sample_panel, y = proportion, fill = plot_celltype)
  ) +
    geom_col(width = 0.84, color = "white", linewidth = 0.34) +
    scale_fill_manual(values = celltype_colors, drop = FALSE) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = c(0, 0.25, 0.50, 0.75, 1.00),
      labels = function(x) sprintf("%.2f", x),
      expand = c(0, 0)
    ) +
    scale_x_discrete(drop = FALSE, expand = expansion(add = 0.28)) +
    labs(
      title = panel_title,
      x = "Sample ID",
      y = if (show_y_title) "Proportion of nuclei" else NULL
    ) +
    coord_cartesian(clip = "off") +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "none",
      plot.title = element_text(
        face = "bold",
        size = 20,
        hjust = 0.5,
        margin = margin(b = 12)
      ),
      axis.title.x = element_text(
        face = "bold",
        size = 17,
        margin = margin(t = 6)
      ),
      axis.title.y = element_text(
        face = "bold",
        size = 17,
        margin = margin(r = 8)
      ),
      axis.text.x = element_text(
        face = "bold",
        size = 12.5,
        color = "black"
      ),
      axis.text.y = element_text(
        face = "bold",
        size = 12.5,
        color = "black"
      ),
      axis.line = element_line(color = "black", linewidth = 1.15),
      axis.ticks = element_line(color = "black", linewidth = 1.05),
      axis.ticks.length = unit(0.12, "cm"),
      panel.grid = element_blank(),
      plot.margin = margin(0, 4, 0, 4)
    )

  if (!show_y_title) {
    panel <- panel + theme(axis.title.y = element_blank())
  }

  panel
}

p_composition_meninges <- composition_df |>
  filter(sample_panel %in% c("M1", "M2")) |>
  make_composition_panel(
    sample_levels = c("M1", "M2"),
    panel_title = "Meninges",
    show_y_title = TRUE
  )

p_composition_dcln <- composition_df |>
  filter(sample_panel %in% c("L1", "L2")) |>
  make_composition_panel(
    sample_levels = c("L1", "L2"),
    panel_title = "dCLN",
    show_y_title = FALSE
  )

composition_panel_grid <- plot_grid(
  p_composition_meninges,
  p_composition_dcln,
  nrow = 1,
  rel_widths = c(1.04, 1.00),
  align = "hv",
  axis = "tb"
)

final_composition_plot <- ggdraw() +
  draw_plot(
    composition_panel_grid,
    x = 0.050,
    y = 0.045,
    width = 0.850,
    height = 0.885
  )

figure1f_tiff <- file.path(
  figure1_dir,
  "figure_1f_celltype_composition_stacked_plot.tiff"
)
ggsave(
  filename = figure1f_tiff,
  plot = final_composition_plot,
  width = 5.20,
  height = 5.45,
  dpi = 300,
  bg = "white",
  compression = "lzw"
)

cat("Figure 1f samples: ", paste(composition_samples, collapse = ", "), "\n", sep = "")
cat("Saved Figure 1f:\n", figure1f_tiff, "\n", figure1f_composition_csv, "\n", sep = "")

################################################################################
# Figure 1g: scCODA broad cell-type composition
################################################################################

figure1g_sccoda_dir <- file.path(
  base_dir,
  "outs",
  "output",
  "scCODA_meninges_vs_dcln"
)
dir.create(figure1g_sccoda_dir, recursive = TRUE, showWarnings = FALSE)

figure1g_python_path <- file.path(
  figure1g_sccoda_dir,
  "run_sccoda_broad_celltypes.py"
)
figure1g_summary_path <- file.path(
  figure1g_sccoda_dir,
  "sccoda_broad_celltypes_summary_table.csv"
)
figure1g_python_bin <- file.path(
  figure1g_sccoda_dir,
  ".venv_sccoda_latest",
  "bin",
  "python"
)

# The Python source is embedded so Figure 1g can be regenerated from the
# Figure 1f composition table when the cached scCODA summary is unavailable.
figure1g_python <- r"---{#!/usr/bin/env python3
"""Run scCODA for Figure 1g broad cell-type composition."""

from __future__ import annotations

import json
import os
import platform
import warnings
from pathlib import Path

import anndata
import arviz
import numpy as np
import pandas as pd
import scipy
import sccoda
import sccoda.util.cell_composition_data as cell_data
import sccoda.util.comp_ana as comp_ana
import tensorflow as tf
import tensorflow_probability as tfp


OUT_DIR = Path(
    "/Users/shingheimok/Desktop/Phd_Study/research/Wu_lab/b_cells/"
    "outs/output/scCODA_meninges_vs_dcln"
)
SOURCE_CSV = Path(
    "/Users/shingheimok/Desktop/Phd_Study/research/Wu_lab/b_cells/"
    "outs/output/all_cells_umap_celltype_major_template_fit/"
    "figure_1f_celltype_composition_by_sample.csv"
)

BROAD_ORDER = [
    "B cells",
    "T cells",
    "Neutrophil",
    "Cycling Neutrophil",
    "Monocyte",
    "Macrophage/microglia",
    "DC",
    "NK cell",
]


def build_counts() -> pd.DataFrame:
    long_df = pd.read_csv(SOURCE_CSV)
    wide = (
        long_df.pivot(
            index=["sample_panel", "tissue_panel"],
            columns="plot_celltype",
            values="n_cells",
        )
        .reset_index()
        .rename(
            columns={
                "sample_panel": "sample_short",
                "tissue_panel": "group",
            }
        )
    )
    wide = wide.loc[:, ["sample_short", "group", *BROAD_ORDER]]
    wide["sample_short"] = pd.Categorical(
        wide["sample_short"],
        ["M1", "M2", "L1", "L2"],
        ordered=True,
    )
    wide = wide.sort_values("sample_short").reset_index(drop=True)
    wide["sample_short"] = wide["sample_short"].astype(str)
    wide.to_csv(
        OUT_DIR / "broad_celltype_counts_by_sample_for_sccoda.csv",
        index=False,
    )
    return wide


def flatten_effects(result) -> pd.DataFrame:
    effects = result.effect_df.reset_index()
    effects = effects.rename(
        columns={
            "Cell Type": "cell_type",
            "Covariate": "covariate",
            "Final Parameter": "final_parameter",
            "SD": "sd",
            "Inclusion probability": "pi_inclusion_probability",
            "Expected Sample": "expected_sample",
            "log2-fold change": "log2_fold_change",
        }
    )
    for col in [c for c in effects.columns if c.startswith("HDI ")]:
        effects = effects.rename(
            columns={col: col.lower().replace(" ", "_").replace("%", "pct")}
        )
    effects["cell_type"] = pd.Categorical(
        effects["cell_type"],
        categories=BROAD_ORDER,
        ordered=True,
    )
    effects["significant_by_sccoda_threshold"] = (
        effects["final_parameter"].abs() > 1e-12
    )
    effects["sccoda_threshold_probability"] = result.model_specs.get(
        "threshold_prob",
        np.nan,
    )
    effects["reference_cell_type"] = BROAD_ORDER[result.model_specs["reference"]]
    effects["pairwise_comparison"] = "dCLN vs Meninges"
    return effects.sort_values("cell_type").reset_index(drop=True)


def main() -> None:
    os.environ.setdefault("TF_CPP_MIN_LOG_LEVEL", "2")
    warnings.filterwarnings("ignore")
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    counts = build_counts().set_index("sample_short")
    data = cell_data.from_pandas(counts, covariate_columns=["group"])
    model = comp_ana.CompositionalAnalysis(
        data,
        formula="C(group, Treatment('Meninges'))",
        reference_cell_type="automatic",
    )
    result = model.sample_hmc_da(
        num_results=20_000,
        num_burnin=5_000,
        verbose=False,
    )
    result.set_fdr(est_fdr=0.05)

    effects = flatten_effects(result)
    effects.to_csv(
        OUT_DIR / "sccoda_broad_celltypes_effects.csv",
        index=False,
    )
    result.to_netcdf(str(OUT_DIR / "sccoda_broad_celltypes_result.nc"))

    group_counts = (
        pd.read_csv(SOURCE_CSV)
        .query("plot_celltype in @BROAD_ORDER")
        .groupby(["tissue_panel", "plot_celltype"], as_index=False)
        .agg(n=("n_cells", "sum"))
    )
    totals = group_counts.groupby("tissue_panel")["n"].transform("sum")
    group_counts["proportion"] = group_counts["n"] / totals
    wide = group_counts.pivot(
        index="plot_celltype",
        columns="tissue_panel",
        values=["n", "proportion"],
    )
    wide.columns = [f"{a}_{b}" for a, b in wide.columns]
    wide = wide.reset_index().rename(
        columns={
            "plot_celltype": "cell_type",
            "n_Meninges": "meninges_count",
            "n_dCLN": "dcln_count",
            "proportion_Meninges": "meninges_proportion",
            "proportion_dCLN": "dcln_proportion",
        }
    )
    summary = effects.merge(wide, on="cell_type", how="left")
    summary["observed_delta_dcln_minus_meninges"] = (
        summary["dcln_proportion"] - summary["meninges_proportion"]
    )
    summary["cell_type"] = pd.Categorical(
        summary["cell_type"],
        categories=BROAD_ORDER,
        ordered=True,
    )
    summary = summary.sort_values("cell_type")
    summary.to_csv(
        OUT_DIR / "sccoda_broad_celltypes_summary_table.csv",
        index=False,
    )

    metadata = {
        "analysis": (
            "scCODA Figure 1g broad cell-type composition, "
            "dCLN vs Meninges"
        ),
        "source_csv": str(SOURCE_CSV),
        "formula": "C(group, Treatment('Meninges'))",
        "sampler": "sample_hmc_da",
        "num_results": 20_000,
        "num_burnin": 5_000,
        "automatic_reference_cell_type": (
            BROAD_ORDER[result.model_specs["reference"]]
        ),
        "sccoda_threshold_probability": result.model_specs.get(
            "threshold_prob"
        ),
        "sampling_stats": {
            "chain_length": int(result.sampling_stats["chain_length"]),
            "num_burnin": int(result.sampling_stats["num_burnin"]),
            "acc_rate": float(result.sampling_stats["acc_rate"]),
            "duration_seconds": float(result.sampling_stats["duration"]),
        },
        "versions": {
            "python": platform.python_version(),
            "sccoda": getattr(sccoda, "__version__", "unknown"),
            "tensorflow": tf.__version__,
            "tensorflow_probability": tfp.__version__,
            "anndata": anndata.__version__,
            "arviz": arviz.__version__,
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy.__version__,
        },
    }
    with open(
        OUT_DIR / "sccoda_broad_celltypes_run_metadata.json",
        "w",
        encoding="utf-8",
    ) as handle:
        json.dump(metadata, handle, indent=2)

    print(json.dumps(metadata, indent=2))
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
}---"

writeLines(figure1g_python, figure1g_python_path, useBytes = TRUE)

if (!file.exists(figure1g_summary_path)) {
  if (!file.exists(figure1g_python_bin)) {
    stop(
      "Missing Figure 1g scCODA summary and Python environment: ",
      figure1g_python_bin,
      "\nRun the generated Python script manually in a scCODA environment."
    )
  }
  system2(figure1g_python_bin, figure1g_python_path)
}

figure1g_celltype_order <- c(
  "B cells",
  "T cells",
  "Neutrophil",
  "Cycling Neutrophil",
  "Monocyte",
  "Macrophage/microglia",
  "DC",
  "NK cell"
)
figure1g_x_positions <- setNames(
  1 + (seq_along(figure1g_celltype_order) - 1) * 0.5,
  figure1g_celltype_order
)

figure1g_plot_df <- read.csv(
  figure1f_composition_csv,
  check.names = FALSE
) |>
  filter(plot_celltype %in% figure1g_celltype_order) |>
  mutate(
    tissue_panel = factor(tissue_panel, levels = c("Meninges", "dCLN")),
    plot_celltype = factor(
      plot_celltype,
      levels = figure1g_celltype_order
    ),
    plot_proportion = pmax(proportion, 5e-4),
    x_pos = figure1g_x_positions[as.character(plot_celltype)]
  )

figure1g_dot_df <- figure1g_plot_df |>
  mutate(
    group_offset = ifelse(tissue_panel == "Meninges", -0.115, 0.115),
    sample_offset = ifelse(
      sample_panel %in% c("M1", "L1"),
      -0.018,
      0.018
    ),
    x_dot = x_pos + group_offset + sample_offset
  )

figure1g_legend_df <- data.frame(
  tissue_panel = factor(
    c("Meninges", "dCLN"),
    levels = c("Meninges", "dCLN")
  ),
  x_legend = c(0.62, 0.62),
  y_legend = c(0.002, 0.002)
)

figure1g_summary_df <- read.csv(
  figure1g_summary_path,
  check.names = FALSE
) |>
  mutate(
    cell_type = factor(cell_type, levels = figure1g_celltype_order),
    significant_by_sccoda_threshold = tolower(
      as.character(significant_by_sccoda_threshold)
    ) %in% c("true", "t", "1"),
    pi_inclusion_probability = as.numeric(pi_inclusion_probability),
    x = figure1g_x_positions[as.character(cell_type)],
    stars = case_when(
      significant_by_sccoda_threshold &
        pi_inclusion_probability >= 0.95 ~ "**",
      significant_by_sccoda_threshold ~ "*",
      TRUE ~ NA_character_
    )
  ) |>
  filter(!is.na(stars))

figure1g_brackets <- figure1g_summary_df |>
  mutate(
    xmin = x - 0.11,
    xmax = x + 0.11,
    y = 1.03,
    yend = 0.96,
    text_y = 1.08
  )

figure1g_plot <- ggplot(
  figure1g_plot_df,
  aes(
    x = x_pos,
    y = plot_proportion,
    fill = tissue_panel,
    group = interaction(plot_celltype, tissue_panel)
  )
) +
  geom_boxplot(
    position = position_dodge(width = 0.26),
    width = 0.24,
    linewidth = 0.45,
    outlier.shape = NA,
    color = "black",
    show.legend = FALSE
  ) +
  geom_point(
    data = figure1g_legend_df,
    aes(x = x_legend, y = y_legend, fill = tissue_panel),
    inherit.aes = FALSE,
    shape = 22,
    size = 4.0,
    stroke = 0.7,
    color = "black",
    alpha = 0,
    show.legend = TRUE
  ) +
  geom_point(
    data = figure1g_dot_df,
    aes(x = x_dot, y = plot_proportion),
    inherit.aes = FALSE,
    shape = 21,
    size = 0.9,
    stroke = 0.25,
    fill = "black",
    color = "white",
    show.legend = FALSE
  ) +
  geom_segment(
    data = figure1g_brackets,
    aes(x = xmin, xend = xmax, y = y, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.55,
    color = "black"
  ) +
  geom_segment(
    data = figure1g_brackets,
    aes(x = xmin, xend = xmin, y = yend, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.55,
    color = "black"
  ) +
  geom_segment(
    data = figure1g_brackets,
    aes(x = xmax, xend = xmax, y = yend, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.55,
    color = "black"
  ) +
  geom_text(
    data = figure1g_brackets,
    aes(x = x, y = text_y, label = stars),
    inherit.aes = FALSE,
    size = 3.4,
    fontface = "bold",
    color = "black"
  ) +
  annotate(
    "segment",
    x = 0.68,
    xend = 1.60,
    y = 1.245,
    yend = 1.245,
    linewidth = 0.58,
    color = "black"
  ) +
  annotate(
    "segment",
    x = 3.98,
    xend = 4.87,
    y = 1.245,
    yend = 1.245,
    linewidth = 0.58,
    color = "black"
  ) +
  annotate(
    "text",
    x = 2.75,
    y = 1.245,
    label = "Broad cell types",
    fontface = "bold",
    size = 5.2,
    color = "black"
  ) +
  scale_fill_manual(
    values = c("Meninges" = "white", "dCLN" = "grey60"),
    breaks = c("Meninges", "dCLN")
  ) +
  scale_y_sqrt(
    breaks = c(0.01, 0.10, 0.50, 1.00),
    labels = c("0.01", "0.10", "0.50", "1.00"),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_x_continuous(
    limits = c(0.55, 4.95),
    breaks = figure1g_x_positions,
    labels = names(figure1g_x_positions),
    expand = c(0, 0)
  ) +
  labs(x = NULL, y = "Proportion of nuclei", fill = NULL) +
  coord_cartesian(ylim = c(0, 1.06), clip = "off") +
  theme_classic(base_size = 10) +
  theme(
    axis.line = element_line(linewidth = 0.55, color = "black"),
    axis.ticks = element_line(linewidth = 0.50, color = "black"),
    axis.ticks.length = unit(0.09, "in"),
    axis.title.y = element_text(
      size = 13,
      face = "bold",
      color = "black",
      margin = margin(r = 6)
    ),
    axis.text.y = element_text(
      size = 10,
      face = "bold",
      color = "black"
    ),
    axis.text.x = element_text(
      size = 11.4,
      face = "bold",
      color = "black",
      angle = 35,
      hjust = 1,
      vjust = 1
    ),
    legend.position = c(0.96, 0.81),
    legend.justification = c(0, 1),
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.key.size = unit(0.20, "in"),
    legend.text = element_text(
      size = 11.4,
      face = "bold",
      color = "black"
    ),
    plot.margin = margin(t = 30, r = 115, b = 11, l = 8)
  ) +
  guides(
    fill = guide_legend(
      override.aes = list(
        alpha = 1,
        shape = 22,
        size = 4.4,
        stroke = 0.7,
        color = "black"
      )
    )
  )

figure1g_tiff <- file.path(
  figure1_dir,
  "figure_1g_sccoda_broad_celltypes_composition_boxplot.tiff"
)
ggsave(
  filename = figure1g_tiff,
  plot = figure1g_plot,
  width = 5.9,
  height = 4.8,
  dpi = 300,
  compression = "lzw",
  bg = "white"
)

cat("Saved Figure 1g:\n", figure1g_tiff, "\n", sep = "")
cat("scCODA script:\n", figure1g_python_path, "\n", sep = "")

################################################################################
# Figures 1d1-1d4: sample-specific all-cell UMAPs
################################################################################

save_sample_umap_figure <- function(sample_name, figure_id) {
  sample_df <- umap_df |>
    filter(sample_panel == sample_name)

  if (nrow(sample_df) == 0) {
    stop("No cells found for sample panel: ", sample_name)
  }

  p_sample_umap <- ggplot(sample_df, aes(UMAP_1, UMAP_2)) +
    geom_point(color = "#1F1F1F", size = 0.50, alpha = 1, stroke = 0) +
    coord_fixed(
      xlim = c(axis_x0 - 0.052 * x_span, x_range[2] + 0.034 * x_span),
      ylim = c(axis_y0 - 0.073 * y_span, y_range[2] + 0.031 * y_span),
      ratio = 1,
      expand = FALSE,
      clip = "off"
    ) +
    theme_void(base_size = 12) +
    theme(
      plot.margin = margin(0, 0, 0, 0)
    )

  final_sample_plot <- ggdraw() +
    draw_plot(p_sample_umap, x = 0.092, y = 0.070, width = 0.817, height = 0.753) +
    draw_line(c(0.347, 0.520), c(0.878, 0.878), color = "black", linewidth = 1.55) +
    draw_label(
      sample_name,
      x = 0.573,
      y = 0.878,
      hjust = 0.5,
      vjust = 0.5,
      fontface = "bold",
      size = 23.5
    ) +
    draw_line(c(0.624, 0.790), c(0.878, 0.878), color = "black", linewidth = 1.55) +
    draw_line(
      c(0.168, 0.168, 0.325),
      c(0.277, 0.124, 0.124),
      color = "black",
      linewidth = 1.55,
      lineend = "butt",
      linejoin = "mitre"
    ) +
    draw_label(
      "UMAP1",
      x = 0.249,
      y = 0.102,
      hjust = 0.5,
      vjust = 0.5,
      fontface = "bold",
      size = 15.5
    ) +
    draw_label(
      "UMAP2",
      x = 0.138,
      y = 0.194,
      angle = 90,
      hjust = 0.5,
      vjust = 0.5,
      fontface = "bold",
      size = 15.5
    )

  out_tiff_sample <- file.path(figure1_dir, paste0(figure_id, ".tiff"))

  ggsave(
    filename = out_tiff_sample,
    plot = final_sample_plot,
    width = 5.36,
    height = 5.48,
    dpi = 300,
    bg = "white",
    compression = "lzw"
  )

  cat(figure_id, sample_name, "cells plotted:", nrow(sample_df), "\n")
  cat("Saved:\n")
  cat(out_tiff_sample, "\n")
}

sample_figures <- data.frame(
  sample_name = c("M1", "M2", "L1", "L2"),
  figure_id = c("figure_1d1", "figure_1d2", "figure_1d3", "figure_1d4"),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(sample_figures))) {
  save_sample_umap_figure(sample_figures$sample_name[i], sample_figures$figure_id[i])
}

################################################################################
# T-cell-only reclustering for Figure 1e1
################################################################################

if (!"celltype_minor" %in% colnames(obj@meta.data)) {
  stop("Missing metadata column: celltype_minor")
}

tcell_major_values <- c("CD4_T", "CD8_T", "UNC_T", "T_cell")
tcell_cells <- rownames(obj@meta.data)[obj$celltype_major %in% tcell_major_values]

if (length(tcell_cells) == 0) {
  stop("No T cells found using celltype_major values: ", paste(tcell_major_values, collapse = ", "))
}

tcell_obj <- subset(obj, cells = tcell_cells)
DefaultAssay(tcell_obj) <- "RNA"
tcell_obj <- tryCatch(
  JoinLayers(tcell_obj, assay = "RNA"),
  error = function(e) tcell_obj
)

set.seed(1234)
tcell_obj <- NormalizeData(tcell_obj, assay = "RNA", verbose = FALSE)

set.seed(1234)
tcell_obj <- FindVariableFeatures(
  tcell_obj,
  assay = "RNA",
  selection.method = "vst",
  nfeatures = 2000,
  verbose = FALSE
)

set.seed(1234)
tcell_obj <- ScaleData(
  tcell_obj,
  assay = "RNA",
  features = VariableFeatures(tcell_obj),
  verbose = FALSE
)

set.seed(1234)
tcell_obj <- RunPCA(
  tcell_obj,
  assay = "RNA",
  features = VariableFeatures(tcell_obj),
  npcs = 10,
  verbose = FALSE
)

set.seed(1234)
tcell_obj <- FindNeighbors(
  tcell_obj,
  dims = 1:10,
  k.param = 30,
  verbose = FALSE
)

set.seed(1234)
tcell_obj <- FindClusters(
  tcell_obj,
  resolution = 0.4,
  verbose = FALSE
)
tcell_obj$tcell_only_cluster_pc10_k30_res04 <- as.character(Idents(tcell_obj))

set.seed(1234)
tcell_obj <- RunUMAP(
  tcell_obj,
  dims = 1:10,
  n.neighbors = 30,
  seed.use = 1234,
  verbose = FALSE
)

tcell_umap_df <- as.data.frame(Embeddings(tcell_obj, reduction = "umap"))
colnames(tcell_umap_df)[1:2] <- c("UMAP_1", "UMAP_2")
tcell_umap_df$cell <- rownames(tcell_umap_df)
tcell_umap_df$celltype_minor <- as.character(tcell_obj$celltype_minor)

tcell_minor_levels <- c(
  "Th1",
  "naive_CD8_T_cell",
  "exhausted_CD8_T_cell",
  "NK_like_CD8_T_cell",
  "CD8_Treg",
  "unconventional_T_cell"
)
tcell_minor_levels <- c(
  tcell_minor_levels[tcell_minor_levels %in% unique(tcell_umap_df$celltype_minor)],
  sort(setdiff(unique(tcell_umap_df$celltype_minor), tcell_minor_levels))
)
tcell_umap_df$celltype_minor <- factor(tcell_umap_df$celltype_minor, levels = tcell_minor_levels)

tcell_minor_labels <- c(
  "Th1" = "Th1",
  "naive_CD8_T_cell" = "Naive CD8 T cell",
  "exhausted_CD8_T_cell" = "Exhausted CD8 T cell",
  "NK_like_CD8_T_cell" = "NK-like CD8 T cell",
  "CD8_Treg" = "CD8 Treg",
  "unconventional_T_cell" = "Unconventional T cell"
)
missing_tcell_labels <- setdiff(levels(tcell_umap_df$celltype_minor), names(tcell_minor_labels))
if (length(missing_tcell_labels) > 0) {
  tcell_minor_labels <- c(
    tcell_minor_labels,
    setNames(gsub("_", " ", missing_tcell_labels), missing_tcell_labels)
  )
}
tcell_minor_labels <- tcell_minor_labels[levels(tcell_umap_df$celltype_minor)]

tcell_minor_colors <- c(
  "Th1" = "#A323A6",
  "naive_CD8_T_cell" = "#4659B8",
  "exhausted_CD8_T_cell" = "#C91F33",
  "NK_like_CD8_T_cell" = "#1F78B4",
  "CD8_Treg" = "#16864F",
  "unconventional_T_cell" = "#43B5D8"
)
missing_tcell_colors <- setdiff(levels(tcell_umap_df$celltype_minor), names(tcell_minor_colors))
if (length(missing_tcell_colors) > 0) {
  tcell_minor_colors <- c(
    tcell_minor_colors,
    setNames(grDevices::hcl.colors(length(missing_tcell_colors), palette = "Dark 3"), missing_tcell_colors)
  )
}
tcell_minor_colors <- tcell_minor_colors[levels(tcell_umap_df$celltype_minor)]

tcell_count_df <- tcell_umap_df |>
  count(celltype_minor, name = "n_cells") |>
  mutate(
    label = tcell_minor_labels[as.character(celltype_minor)],
    color = tcell_minor_colors[as.character(celltype_minor)]
  )
write.csv(
  tcell_count_df,
  file.path(output_dir, "figure_1e1_tcell_celltype_minor_counts_and_palette.csv"),
  row.names = FALSE
)

tcell_plot_df <- tcell_umap_df |>
  left_join(tcell_count_df |> select(celltype_minor, n_cells), by = "celltype_minor") |>
  arrange(desc(n_cells))

tcell_x_range <- range(tcell_plot_df$UMAP_1, na.rm = TRUE)
tcell_y_range <- range(tcell_plot_df$UMAP_2, na.rm = TRUE)
tcell_x_span <- diff(tcell_x_range)
tcell_y_span <- diff(tcell_y_range)

################################################################################
# Figure 1e1: T-cell-only UMAP colored by minor cell type
################################################################################

p_tcell_umap_1b6 <- ggplot(tcell_plot_df, aes(UMAP_1, UMAP_2, color = celltype_minor)) +
  geom_point(size = 0.62, alpha = 0.98, stroke = 0) +
  scale_color_manual(values = tcell_minor_colors, drop = FALSE) +
  coord_fixed(
    xlim = c(tcell_x_range[1] - 0.055 * tcell_x_span, tcell_x_range[2] + 0.055 * tcell_x_span),
    ylim = c(tcell_y_range[1] - 0.055 * tcell_y_span, tcell_y_range[2] + 0.055 * tcell_y_span),
    ratio = 1,
    expand = FALSE,
    clip = "off"
  ) +
  theme_void(base_size = 12) +
  theme(
    legend.position = "none",
    plot.margin = margin(0, 4, 3, 0)
  )

tcell_legend_direct_df <- tibble(
  celltype_minor = levels(tcell_umap_df$celltype_minor),
  label = unname(tcell_minor_labels),
  x_dot = 0.595,
  x_label = 0.636,
  y = 0.691 - (seq_along(levels(tcell_umap_df$celltype_minor)) - 1) * 0.052
)

final_tcell_1b6_plot <- ggdraw() +
  draw_plot(p_tcell_umap_1b6, x = -0.071, y = 0.043, width = 0.698, height = 0.829) +
  draw_line(
    c(0.179, 0.268),
    c(0.911, 0.911),
    color = "black",
    linewidth = 1.55,
    lineend = "butt"
  ) +
  draw_label(
    "T cells",
    x = 0.329,
    y = 0.911,
    hjust = 0.5,
    vjust = 0.5,
    fontface = "bold",
    size = 23.5
  ) +
  draw_line(
    c(0.390, 0.479),
    c(0.911, 0.911),
    color = "black",
    linewidth = 1.55,
    lineend = "butt"
  ) +
  draw_line(
    c(0.059, 0.059, 0.147),
    c(0.231, 0.105, 0.105),
    color = "black",
    linewidth = 1.35,
    lineend = "butt",
    linejoin = "mitre"
  ) +
  draw_label(
    "UMAP1",
    x = 0.103,
    y = 0.078,
    hjust = 0.5,
    vjust = 0.5,
    fontface = "bold",
    size = 15.5
  ) +
  draw_label(
    "UMAP2",
    x = 0.043,
    y = 0.166,
    angle = 90,
    hjust = 0.5,
    vjust = 0.5,
    fontface = "bold",
    size = 15.5
  ) +
  geom_point(
    data = tcell_legend_direct_df,
    aes(x = x_dot, y = y, fill = celltype_minor),
    inherit.aes = FALSE,
    shape = 21,
    size = 5.3,
    stroke = 0.12,
    color = "white",
    show.legend = FALSE
  ) +
  geom_text(
    data = tcell_legend_direct_df,
    aes(x = x_label, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 0.5,
    size = 5.65,
    fontface = "bold",
    color = "black"
  ) +
  scale_fill_manual(values = tcell_minor_colors, drop = FALSE)

out_tcell_1b6_tiff <- file.path(figure1_dir, "figure_1e1_tcell_umap.tiff")

ggsave(
  filename = out_tcell_1b6_tiff,
  plot = final_tcell_1b6_plot,
  width = 9.14,
  height = 5.50,
  dpi = 300,
  bg = "white",
  compression = "lzw"
)

write.csv(
  tcell_count_df,
  file.path(output_dir, "figure_1e1_tcell_celltype_minor_counts_and_palette.csv"),
  row.names = FALSE
)

cat("figure_1e1 T cells plotted:", nrow(tcell_plot_df), "\n")
cat("T-cell clusters:", length(unique(tcell_obj$tcell_only_cluster_pc10_k30_res04)), "\n")
cat("Legend groups:", length(levels(tcell_umap_df$celltype_minor)), "\n")
cat("Saved:\n")
cat(out_tcell_1b6_tiff, "\n")
cat(file.path(output_dir, "figure_1e1_tcell_celltype_minor_counts_and_palette.csv"), "\n")

################################################################################
# Figure 1e2: macrophage/microglia-only UMAP colored by minor cell type
################################################################################

macrophage_major_value <- "macrophage/microglia"
macrophage_cells <- rownames(obj@meta.data)[obj$celltype_major == macrophage_major_value]

if (length(macrophage_cells) == 0) {
  stop("No macrophage/microglia cells found using celltype_major value: ", macrophage_major_value)
}

macrophage_obj <- subset(obj, cells = macrophage_cells)
DefaultAssay(macrophage_obj) <- "RNA"
macrophage_obj <- tryCatch(
  JoinLayers(macrophage_obj, assay = "RNA"),
  error = function(e) macrophage_obj
)

set.seed(1234)
macrophage_obj <- NormalizeData(macrophage_obj, assay = "RNA", verbose = FALSE)

set.seed(1234)
macrophage_obj <- FindVariableFeatures(
  macrophage_obj,
  assay = "RNA",
  selection.method = "vst",
  nfeatures = 2000,
  verbose = FALSE
)

set.seed(1234)
macrophage_obj <- ScaleData(
  macrophage_obj,
  assay = "RNA",
  features = VariableFeatures(macrophage_obj),
  verbose = FALSE
)

set.seed(1234)
macrophage_obj <- RunPCA(
  macrophage_obj,
  assay = "RNA",
  features = VariableFeatures(macrophage_obj),
  npcs = 10,
  verbose = FALSE
)

set.seed(1234)
macrophage_obj <- FindNeighbors(
  macrophage_obj,
  dims = 1:10,
  k.param = 30,
  verbose = FALSE
)

set.seed(1234)
macrophage_obj <- FindClusters(
  macrophage_obj,
  resolution = 0.6,
  verbose = FALSE
)
macrophage_obj$macrophage_microglia_cluster_pc10_k30_res06 <- as.character(Idents(macrophage_obj))

set.seed(1234)
macrophage_obj <- RunUMAP(
  macrophage_obj,
  dims = 1:10,
  n.neighbors = 30,
  seed.use = 1234,
  verbose = FALSE
)

macrophage_umap_df <- as.data.frame(Embeddings(macrophage_obj, reduction = "umap"))
colnames(macrophage_umap_df)[1:2] <- c("UMAP_1", "UMAP_2")
macrophage_umap_df$cell <- rownames(macrophage_umap_df)
macrophage_umap_df$celltype_minor <- as.character(macrophage_obj$celltype_minor)

macrophage_minor_levels <- c("macrophage", "DAM1", "DAM2")
macrophage_minor_levels <- c(
  macrophage_minor_levels[macrophage_minor_levels %in% unique(macrophage_umap_df$celltype_minor)],
  sort(setdiff(unique(macrophage_umap_df$celltype_minor), macrophage_minor_levels))
)
macrophage_umap_df$celltype_minor <- factor(
  macrophage_umap_df$celltype_minor,
  levels = macrophage_minor_levels
)

macrophage_minor_labels <- c(
  "macrophage" = "Macrophages",
  "DAM1" = "DAM1",
  "DAM2" = "DAM2"
)
missing_macrophage_labels <- setdiff(levels(macrophage_umap_df$celltype_minor), names(macrophage_minor_labels))
if (length(missing_macrophage_labels) > 0) {
  macrophage_minor_labels <- c(
    macrophage_minor_labels,
    setNames(gsub("_", " ", missing_macrophage_labels), missing_macrophage_labels)
  )
}
macrophage_minor_labels <- macrophage_minor_labels[levels(macrophage_umap_df$celltype_minor)]

macrophage_minor_colors <- c(
  "macrophage" = "#A95524",
  "DAM1" = "#5E49A8",
  "DAM2" = "#159DC1"
)
missing_macrophage_colors <- setdiff(levels(macrophage_umap_df$celltype_minor), names(macrophage_minor_colors))
if (length(missing_macrophage_colors) > 0) {
  macrophage_minor_colors <- c(
    macrophage_minor_colors,
    setNames(grDevices::hcl.colors(length(missing_macrophage_colors), palette = "Dark 3"), missing_macrophage_colors)
  )
}
macrophage_minor_colors <- macrophage_minor_colors[levels(macrophage_umap_df$celltype_minor)]

macrophage_count_df <- macrophage_umap_df |>
  count(celltype_minor, name = "n_cells") |>
  mutate(
    label = macrophage_minor_labels[as.character(celltype_minor)],
    color = macrophage_minor_colors[as.character(celltype_minor)]
  )

macrophage_plot_df <- macrophage_umap_df |>
  left_join(macrophage_count_df |> select(celltype_minor, n_cells), by = "celltype_minor") |>
  arrange(desc(n_cells))

macrophage_x_range <- range(macrophage_plot_df$UMAP_1, na.rm = TRUE)
macrophage_y_range <- range(macrophage_plot_df$UMAP_2, na.rm = TRUE)
macrophage_x_span <- diff(macrophage_x_range)
macrophage_y_span <- diff(macrophage_y_range)

p_macrophage_umap_1b7 <- ggplot(macrophage_plot_df, aes(UMAP_1, UMAP_2, color = celltype_minor)) +
  geom_point(size = 0.82, alpha = 0.98, stroke = 0) +
  scale_color_manual(values = macrophage_minor_colors, drop = FALSE) +
  coord_fixed(
    xlim = c(macrophage_x_range[1] - 0.060 * macrophage_x_span, macrophage_x_range[2] + 0.060 * macrophage_x_span),
    ylim = c(macrophage_y_range[1] - 0.060 * macrophage_y_span, macrophage_y_range[2] + 0.060 * macrophage_y_span),
    ratio = 1,
    expand = FALSE,
    clip = "off"
  ) +
  theme_void(base_size = 12) +
  theme(
    legend.position = "none",
    plot.margin = margin(0, 4, 3, 0)
  )

macrophage_legend_direct_df <- tibble(
  celltype_minor = levels(macrophage_umap_df$celltype_minor),
  label = unname(macrophage_minor_labels),
  x_dot = 0.595,
  x_label = 0.636,
  y = 0.660 - (seq_along(levels(macrophage_umap_df$celltype_minor)) - 1) * 0.065
)

final_macrophage_1b7_plot <- ggdraw() +
  draw_plot(p_macrophage_umap_1b7, x = 0.050, y = 0.155, width = 0.500, height = 0.700) +
  draw_line(
    c(0.130, 0.213),
    c(0.911, 0.911),
    color = "black",
    linewidth = 1.55,
    lineend = "butt"
  ) +
  draw_label(
    "Macrophages/\nMicroglia",
    x = 0.329,
    y = 0.907,
    hjust = 0.5,
    vjust = 0.5,
    fontface = "bold",
    size = 21.0,
    lineheight = 0.84
  ) +
  draw_line(
    c(0.445, 0.528),
    c(0.911, 0.911),
    color = "black",
    linewidth = 1.55,
    lineend = "butt"
  ) +
  draw_line(
    c(0.059, 0.059, 0.147),
    c(0.231, 0.105, 0.105),
    color = "black",
    linewidth = 1.35,
    lineend = "butt",
    linejoin = "mitre"
  ) +
  draw_label(
    "UMAP1",
    x = 0.103,
    y = 0.078,
    hjust = 0.5,
    vjust = 0.5,
    fontface = "bold",
    size = 15.5
  ) +
  draw_label(
    "UMAP2",
    x = 0.043,
    y = 0.166,
    angle = 90,
    hjust = 0.5,
    vjust = 0.5,
    fontface = "bold",
    size = 15.5
  ) +
  geom_point(
    data = macrophage_legend_direct_df,
    aes(x = x_dot, y = y, fill = celltype_minor),
    inherit.aes = FALSE,
    shape = 21,
    size = 5.3,
    stroke = 0.12,
    color = "white",
    show.legend = FALSE
  ) +
  geom_text(
    data = macrophage_legend_direct_df,
    aes(x = x_label, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 0.5,
    size = 5.65,
    fontface = "bold",
    color = "black"
  ) +
  scale_fill_manual(values = macrophage_minor_colors, drop = FALSE)

out_macrophage_1b7_tiff <- file.path(
  figure1_dir,
  "figure_1e2_macrophage_microglia_umap.tiff"
)

ggsave(
  filename = out_macrophage_1b7_tiff,
  plot = final_macrophage_1b7_plot,
  width = 9.14,
  height = 5.50,
  dpi = 300,
  bg = "white",
  compression = "lzw"
)

write.csv(
  macrophage_count_df,
  file.path(output_dir, "figure_1e2_macrophage_microglia_celltype_minor_counts_and_palette.csv"),
  row.names = FALSE
)

cat("figure_1e2 macrophage/microglia cells plotted:", nrow(macrophage_plot_df), "\n")
cat("Macrophage/microglia clusters:", length(unique(macrophage_obj$macrophage_microglia_cluster_pc10_k30_res06)), "\n")
cat("Legend groups:", length(levels(macrophage_umap_df$celltype_minor)), "\n")
cat("Saved:\n")
cat(out_macrophage_1b7_tiff, "\n")
cat(file.path(output_dir, "figure_1e2_macrophage_microglia_celltype_minor_counts_and_palette.csv"), "\n")

################################################################################
# Figure 1e3: DC-only UMAP colored by minor cell type
################################################################################

dc_major_value <- "DC"
dc_cells <- rownames(obj@meta.data)[obj$celltype_major == dc_major_value]

if (length(dc_cells) == 0) {
  stop("No DC cells found using celltype_major value: ", dc_major_value)
}

dc_obj <- subset(obj, cells = dc_cells)
DefaultAssay(dc_obj) <- "RNA"
dc_obj <- tryCatch(
  JoinLayers(dc_obj, assay = "RNA"),
  error = function(e) dc_obj
)

set.seed(1234)
dc_obj <- NormalizeData(dc_obj, assay = "RNA", verbose = FALSE)

set.seed(1234)
dc_obj <- FindVariableFeatures(
  dc_obj,
  assay = "RNA",
  selection.method = "vst",
  nfeatures = 2000,
  verbose = FALSE
)

set.seed(1234)
dc_obj <- ScaleData(
  dc_obj,
  assay = "RNA",
  features = VariableFeatures(dc_obj),
  verbose = FALSE
)

set.seed(1234)
dc_obj <- RunPCA(
  dc_obj,
  assay = "RNA",
  features = VariableFeatures(dc_obj),
  npcs = 15,
  verbose = FALSE
)

set.seed(1234)
dc_obj <- FindNeighbors(
  dc_obj,
  dims = 1:15,
  k.param = 30,
  verbose = FALSE
)

set.seed(1234)
dc_obj <- FindClusters(
  dc_obj,
  resolution = 0.8,
  verbose = FALSE
)
dc_obj$dc_cluster_pc15_k30_res08 <- as.character(Idents(dc_obj))

set.seed(1234)
dc_obj <- RunUMAP(
  dc_obj,
  dims = 1:15,
  n.neighbors = 30,
  seed.use = 1234,
  verbose = FALSE
)

dc_umap_df <- as.data.frame(Embeddings(dc_obj, reduction = "umap"))
colnames(dc_umap_df)[1:2] <- c("UMAP_1", "UMAP_2")
dc_umap_df$cell <- rownames(dc_umap_df)
dc_umap_df$celltype_minor <- as.character(dc_obj$celltype_minor)

dc_minor_levels <- c("migratory_DC", "non_migratory_DC")
dc_minor_levels <- c(
  dc_minor_levels[dc_minor_levels %in% unique(dc_umap_df$celltype_minor)],
  sort(setdiff(unique(dc_umap_df$celltype_minor), dc_minor_levels))
)
dc_umap_df$celltype_minor <- factor(dc_umap_df$celltype_minor, levels = dc_minor_levels)

dc_minor_labels <- c(
  "migratory_DC" = "Migratory DC",
  "non_migratory_DC" = "Non-migratory DC"
)
missing_dc_labels <- setdiff(levels(dc_umap_df$celltype_minor), names(dc_minor_labels))
if (length(missing_dc_labels) > 0) {
  dc_minor_labels <- c(
    dc_minor_labels,
    setNames(gsub("_", " ", missing_dc_labels), missing_dc_labels)
  )
}
dc_minor_labels <- dc_minor_labels[levels(dc_umap_df$celltype_minor)]

dc_minor_colors <- c(
  "migratory_DC" = "#1F68B3",
  "non_migratory_DC" = "#D18A00"
)
missing_dc_colors <- setdiff(levels(dc_umap_df$celltype_minor), names(dc_minor_colors))
if (length(missing_dc_colors) > 0) {
  dc_minor_colors <- c(
    dc_minor_colors,
    setNames(grDevices::hcl.colors(length(missing_dc_colors), palette = "Dark 3"), missing_dc_colors)
  )
}
dc_minor_colors <- dc_minor_colors[levels(dc_umap_df$celltype_minor)]

dc_count_df <- dc_umap_df |>
  count(celltype_minor, name = "n_cells") |>
  mutate(
    label = dc_minor_labels[as.character(celltype_minor)],
    color = dc_minor_colors[as.character(celltype_minor)]
  )

dc_plot_df <- dc_umap_df |>
  left_join(dc_count_df |> select(celltype_minor, n_cells), by = "celltype_minor") |>
  arrange(desc(n_cells))

dc_x_range <- range(dc_plot_df$UMAP_1, na.rm = TRUE)
dc_y_range <- range(dc_plot_df$UMAP_2, na.rm = TRUE)
dc_x_span <- diff(dc_x_range)
dc_y_span <- diff(dc_y_range)

p_dc_umap_1b8 <- ggplot(dc_plot_df, aes(UMAP_1, UMAP_2, color = celltype_minor)) +
  geom_point(size = 0.92, alpha = 0.98, stroke = 0) +
  scale_color_manual(values = dc_minor_colors, drop = FALSE) +
  coord_fixed(
    xlim = c(dc_x_range[1] - 0.060 * dc_x_span, dc_x_range[2] + 0.060 * dc_x_span),
    ylim = c(dc_y_range[1] - 0.060 * dc_y_span, dc_y_range[2] + 0.060 * dc_y_span),
    ratio = 1,
    expand = FALSE,
    clip = "off"
  ) +
  theme_void(base_size = 12) +
  theme(
    legend.position = "none",
    plot.margin = margin(0, 4, 3, 0)
  )

dc_legend_direct_df <- tibble(
  celltype_minor = levels(dc_umap_df$celltype_minor),
  label = unname(dc_minor_labels),
  x_dot = 0.595,
  x_label = 0.636,
  y = 0.660 - (seq_along(levels(dc_umap_df$celltype_minor)) - 1) * 0.065
)

final_dc_1b8_plot <- ggdraw() +
  draw_plot(p_dc_umap_1b8, x = 0.050, y = 0.155, width = 0.500, height = 0.700) +
  draw_line(
    c(0.179, 0.292),
    c(0.911, 0.911),
    color = "black",
    linewidth = 1.55,
    lineend = "butt"
  ) +
  draw_label(
    "DC",
    x = 0.329,
    y = 0.911,
    hjust = 0.5,
    vjust = 0.5,
    fontface = "bold",
    size = 23.5
  ) +
  draw_line(
    c(0.366, 0.479),
    c(0.911, 0.911),
    color = "black",
    linewidth = 1.55,
    lineend = "butt"
  ) +
  draw_line(
    c(0.059, 0.059, 0.147),
    c(0.231, 0.105, 0.105),
    color = "black",
    linewidth = 1.35,
    lineend = "butt",
    linejoin = "mitre"
  ) +
  draw_label(
    "UMAP1",
    x = 0.103,
    y = 0.078,
    hjust = 0.5,
    vjust = 0.5,
    fontface = "bold",
    size = 15.5
  ) +
  draw_label(
    "UMAP2",
    x = 0.043,
    y = 0.166,
    angle = 90,
    hjust = 0.5,
    vjust = 0.5,
    fontface = "bold",
    size = 15.5
  ) +
  geom_point(
    data = dc_legend_direct_df,
    aes(x = x_dot, y = y, fill = celltype_minor),
    inherit.aes = FALSE,
    shape = 21,
    size = 5.3,
    stroke = 0.12,
    color = "white",
    show.legend = FALSE
  ) +
  geom_text(
    data = dc_legend_direct_df,
    aes(x = x_label, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 0.5,
    size = 5.65,
    fontface = "bold",
    color = "black"
  ) +
  scale_fill_manual(values = dc_minor_colors, drop = FALSE)

out_dc_1b8_tiff <- file.path(figure1_dir, "figure_1e3_dc_umap.tiff")

ggsave(
  filename = out_dc_1b8_tiff,
  plot = final_dc_1b8_plot,
  width = 9.14,
  height = 5.50,
  dpi = 300,
  bg = "white",
  compression = "lzw"
)

write.csv(
  dc_count_df,
  file.path(output_dir, "figure_1e3_dc_celltype_minor_counts_and_palette.csv"),
  row.names = FALSE
)

cat("figure_1e3 DC cells plotted:", nrow(dc_plot_df), "\n")
cat("DC clusters:", length(unique(dc_obj$dc_cluster_pc15_k30_res08)), "\n")
cat("Legend groups:", length(levels(dc_umap_df$celltype_minor)), "\n")
cat("Saved:\n")
cat(out_dc_1b8_tiff, "\n")
cat(file.path(output_dir, "figure_1e3_dc_celltype_minor_counts_and_palette.csv"), "\n")

################################################################################
# Disabled legacy duplicate T-cell-only UMAP
################################################################################

if (FALSE) {

tcell_axis_x0 <- tcell_x_range[1] - 0.115 * tcell_x_span
tcell_axis_y0 <- tcell_y_range[1] - 0.070 * tcell_y_span
tcell_axis_x1 <- tcell_axis_x0 + 0.210 * tcell_x_span
tcell_axis_y1 <- tcell_axis_y0 + 0.215 * tcell_y_span
tcell_axis_df <- data.frame(
  x = c(tcell_axis_x1, tcell_axis_x0, tcell_axis_x0),
  y = c(tcell_axis_y0, tcell_axis_y0, tcell_axis_y1)
)

p_tcell_umap <- ggplot(tcell_plot_df, aes(UMAP_1, UMAP_2, color = celltype_minor)) +
  geom_point(size = 0.62, alpha = 0.98, stroke = 0) +
  geom_path(
    data = tcell_axis_df,
    aes(x = x, y = y),
    inherit.aes = FALSE,
    linewidth = 0.86,
    lineend = "butt",
    linejoin = "mitre",
    color = "black"
  ) +
  annotate(
    "text",
    x = tcell_axis_x0,
    y = tcell_axis_y0 - 0.024 * tcell_y_span,
    label = "UMAP1",
    hjust = 0,
    vjust = 1,
    size = 5.2,
    fontface = "bold"
  ) +
  annotate(
    "text",
    x = tcell_axis_x0 - 0.060 * tcell_x_span,
    y = tcell_axis_y0 + 0.010 * tcell_y_span,
    label = "UMAP2",
    angle = 90,
    hjust = 0,
    vjust = 1,
    size = 5.2,
    fontface = "bold"
  ) +
  scale_color_manual(values = tcell_minor_colors, drop = FALSE) +
  coord_fixed(
    xlim = c(tcell_axis_x0 - 0.055 * tcell_x_span, tcell_x_range[2] + 0.050 * tcell_x_span),
    ylim = c(tcell_axis_y0 - 0.075 * tcell_y_span, tcell_y_range[2] + 0.055 * tcell_y_span),
    ratio = 1,
    expand = FALSE,
    clip = "off"
  ) +
  theme_void(base_size = 12) +
  theme(
    legend.position = "none",
    plot.margin = margin(0, 0, 0, 0)
  )

tcell_legend_df <- tibble(
  celltype_minor = levels(tcell_umap_df$celltype_minor),
  label = unname(tcell_minor_labels),
  y = rev(seq(0.62, by = 0.82, length.out = length(levels(tcell_umap_df$celltype_minor))))
)

tcell_legend_panel <- ggplot(tcell_legend_df) +
  geom_point(
    aes(x = 0.08, y = y, fill = celltype_minor),
    shape = 21,
    size = 5.2,
    stroke = 0.12,
    color = "white",
    show.legend = FALSE
  ) +
  geom_text(
    aes(x = 0.19, y = y, label = label),
    hjust = 0,
    vjust = 0.5,
    size = 5.4,
    fontface = "bold",
    color = "black"
  ) +
  scale_fill_manual(values = tcell_minor_colors, drop = FALSE) +
  coord_cartesian(
    xlim = c(0, 1.35),
    ylim = c(0.28, max(tcell_legend_df$y) + 0.44),
    expand = FALSE,
    clip = "off"
  ) +
  theme_void(base_size = 12) +
  theme(
    plot.margin = margin(0, 0, 0, 0)
  )

final_tcell_plot <- ggdraw() +
  draw_plot(p_tcell_umap, x = 0.055, y = 0.072, width = 0.585, height = 0.825) +
  draw_line(c(0.145, 0.235), c(0.925, 0.925), color = "black", linewidth = 1.35) +
  draw_label(
    "Lymphoid",
    x = 0.370,
    y = 0.925,
    hjust = 0.5,
    vjust = 0.5,
    fontface = "bold",
    size = 28.0
  ) +
  draw_line(c(0.510, 0.600), c(0.925, 0.925), color = "black", linewidth = 1.35) +
  draw_plot(tcell_legend_panel, x = 0.665, y = 0.315, width = 0.315, height = 0.505)

out_tcell_tiff <- file.path(output_dir, "figure_1c1.tiff")

ggsave(
  filename = out_tcell_tiff,
  plot = final_tcell_plot,
  width = 9.46,
  height = 5.30,
  dpi = 300,
  bg = "white",
  compression = "lzw"
)

cat("figure_1c1 T cells plotted:", nrow(tcell_plot_df), "\n")
cat("T-cell clusters:", length(unique(tcell_obj$tcell_only_cluster_pc10_k30_res04)), "\n")
cat("Legend groups:", length(levels(tcell_umap_df$celltype_minor)), "\n")
cat("Saved:\n")
cat(out_tcell_tiff, "\n")
cat(file.path(output_dir, "figure_1c1_tcell_celltype_minor_counts_and_palette.csv"), "\n")

}
