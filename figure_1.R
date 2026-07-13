################################################################################
# Figure 1
################################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(cowplot)
  library(dplyr)
})

################################################################################
# Figure 1b1: all-cell UMAP colored by major cell type
################################################################################

rds_file <- "/Users/shingheimok/Desktop/Phd_Study/research/Wu_lab/b_cells/srt_fullannot_with_bcell_dbscan_eps017_min5_min31.rds"
qs_file <- sub("\\.rds$", ".qs", rds_file)
output_dir <- "/Users/shingheimok/Desktop/Phd_Study/research/Wu_lab/b_cells/outs/output/all_cells_umap_celltype_major_template_fit"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

obj <- if (file.exists(qs_file) && requireNamespace("qs", quietly = TRUE)) {
  qs::qread(qs_file)
} else {
  readRDS(rds_file)
}

required_meta_cols <- c("celltype_major", "celltype_minor", "sample_id")
missing_meta_cols <- setdiff(required_meta_cols, colnames(obj@meta.data))
if (length(missing_meta_cols) > 0) {
  stop("Missing metadata columns: ", paste(missing_meta_cols, collapse = ", "))
}
if (!"umap" %in% Reductions(obj)) {
  stop("The object does not contain a UMAP reduction named 'umap'.")
}

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
    "Legend b (All cells)",
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

out_tiff <- file.path(output_dir, "figure_1b1.tiff")

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
# Figures 1b2-1b5: sample-specific all-cell UMAPs
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

  out_tiff_sample <- file.path(output_dir, paste0(figure_id, ".tiff"))

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
  figure_id = c("figure_1b2", "figure_1b3", "figure_1b4", "figure_1b5"),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(sample_figures))) {
  save_sample_umap_figure(sample_figures$sample_name[i], sample_figures$figure_id[i])
}

################################################################################
# T-cell-only reclustering for Figures 1b6 and 1c1
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
  file.path(output_dir, "figure_1c1_tcell_celltype_minor_counts_and_palette.csv"),
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
# Figure 1b6: T-cell-only UMAP colored by minor cell type
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

out_tcell_1b6_tiff <- file.path(output_dir, "figure_1b6.tiff")

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
  file.path(output_dir, "figure_1b6_tcell_celltype_minor_counts_and_palette.csv"),
  row.names = FALSE
)

cat("figure_1b6 T cells plotted:", nrow(tcell_plot_df), "\n")
cat("T-cell clusters:", length(unique(tcell_obj$tcell_only_cluster_pc10_k30_res04)), "\n")
cat("Legend groups:", length(levels(tcell_umap_df$celltype_minor)), "\n")
cat("Saved:\n")
cat(out_tcell_1b6_tiff, "\n")
cat(file.path(output_dir, "figure_1b6_tcell_celltype_minor_counts_and_palette.csv"), "\n")

################################################################################
# Figure 1b7: macrophage/microglia-only UMAP colored by minor cell type
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

out_macrophage_1b7_tiff <- file.path(output_dir, "figure_1b7.tiff")

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
  file.path(output_dir, "figure_1b7_macrophage_microglia_celltype_minor_counts_and_palette.csv"),
  row.names = FALSE
)

cat("figure_1b7 macrophage/microglia cells plotted:", nrow(macrophage_plot_df), "\n")
cat("Macrophage/microglia clusters:", length(unique(macrophage_obj$macrophage_microglia_cluster_pc10_k30_res06)), "\n")
cat("Legend groups:", length(levels(macrophage_umap_df$celltype_minor)), "\n")
cat("Saved:\n")
cat(out_macrophage_1b7_tiff, "\n")
cat(file.path(output_dir, "figure_1b7_macrophage_microglia_celltype_minor_counts_and_palette.csv"), "\n")

################################################################################
# Figure 1b8: DC-only UMAP colored by minor cell type
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

out_dc_1b8_tiff <- file.path(output_dir, "figure_1b8.tiff")

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
  file.path(output_dir, "figure_1b8_dc_celltype_minor_counts_and_palette.csv"),
  row.names = FALSE
)

cat("figure_1b8 DC cells plotted:", nrow(dc_plot_df), "\n")
cat("DC clusters:", length(unique(dc_obj$dc_cluster_pc15_k30_res08)), "\n")
cat("Legend groups:", length(levels(dc_umap_df$celltype_minor)), "\n")
cat("Saved:\n")
cat(out_dc_1b8_tiff, "\n")
cat(file.path(output_dir, "figure_1b8_dc_celltype_minor_counts_and_palette.csv"), "\n")

################################################################################
# Figure 1c1: T-cell-only UMAP colored by minor cell type
################################################################################

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
