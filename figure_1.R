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

if (!"celltype_major" %in% colnames(obj@meta.data)) {
  stop("Missing metadata column: celltype_major")
}
if (!"umap" %in% Reductions(obj)) {
  stop("The object does not contain a UMAP reduction named 'umap'.")
}

umap_df <- as.data.frame(Embeddings(obj, reduction = "umap"))
colnames(umap_df)[1:2] <- c("UMAP_1", "UMAP_2")
umap_df$cell <- rownames(umap_df)
umap_df$celltype_major <- as.character(obj$celltype_major)

umap_df$plot_celltype <- dplyr::case_when(
  umap_df$celltype_major %in% c("CD4_T", "CD8_T", "UNC_T", "T_cell") ~ "T cells",
  umap_df$celltype_major == "B_cell" ~ "B cells",
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
axis_x1 <- axis_x0 + 0.240 * x_span
axis_y1 <- axis_y0 + 0.180 * y_span
axis_df <- data.frame(
  x = c(axis_x1, axis_x0, axis_x0),
  y = c(axis_y0, axis_y0, axis_y1)
)

p_umap <- ggplot(plot_df, aes(UMAP_1, UMAP_2, color = plot_celltype)) +
  geom_point(size = 0.42, alpha = 1, stroke = 0) +
  geom_path(
    data = axis_df,
    aes(x = x, y = y),
    inherit.aes = FALSE,
    linewidth = 0.78,
    lineend = "butt",
    linejoin = "mitre",
    color = "black"
  ) +
  annotate(
    "text",
    x = axis_x0,
    y = axis_y0 - 0.018 * y_span,
    label = "UMAP1",
    hjust = 0,
    vjust = 1,
    size = 4.8,
    fontface = "bold"
  ) +
  annotate(
    "text",
    x = axis_x0 - 0.060 * x_span,
    y = axis_y0 + 0.002 * y_span,
    label = "UMAP2",
    angle = 90,
    hjust = 0,
    vjust = 1,
    size = 4.8,
    fontface = "bold"
  ) +
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

legend_df <- tibble(
  plot_celltype = levels(plot_df$plot_celltype),
  label = levels(plot_df$plot_celltype),
  y = rev(seq(0.78, by = 0.96, length.out = length(levels(plot_df$plot_celltype))))
)

legend_panel <- ggplot(legend_df) +
  annotate("segment", x = 0, xend = 0, y = 0, yend = 7.90, linewidth = 0.82, color = "black") +
  annotate("segment", x = 1, xend = 1, y = 0, yend = 7.90, linewidth = 0.82, color = "black") +
  annotate("segment", x = 0, xend = 1, y = 0, yend = 0, linewidth = 0.82, color = "black") +
  annotate("segment", x = 0, xend = 0.155, y = 7.90, yend = 7.90, linewidth = 0.82, color = "black") +
  annotate("segment", x = 0.845, xend = 1, y = 7.90, yend = 7.90, linewidth = 0.82, color = "black") +
  annotate("rect", xmin = 0.155, xmax = 0.845, ymin = 7.62, ymax = 8.16, fill = "white", color = NA) +
  geom_point(
    aes(x = 0.110, y = y, fill = plot_celltype),
    shape = 21,
    size = 4.75,
    stroke = 0.12,
    color = "white",
    show.legend = FALSE
  ) +
  geom_text(
    aes(x = 0.215, y = y, label = label),
    hjust = 0,
    vjust = 0.5,
    size = 5.15,
    fontface = "bold",
    color = "black"
  ) +
  annotate(
    "text",
    x = 0.50,
    y = 7.90,
    label = "Legend b (All cells)",
    hjust = 0.5,
    vjust = 0.5,
    size = 5.70,
    fontface = "bold",
    color = "black"
  ) +
  scale_fill_manual(values = celltype_colors, drop = FALSE) +
  coord_cartesian(
    xlim = c(0, 1),
    ylim = c(-0.02, 8.18),
    expand = FALSE,
    clip = "off"
  ) +
  theme_void(base_size = 12) +
  theme(
    plot.margin = margin(0, 0, 0, 0)
  )

final_plot <- ggdraw() +
  draw_plot(p_umap, x = 0.060, y = 0.070, width = 0.525, height = 0.825) +
  draw_line(c(0.228, 0.294), c(0.932, 0.932), color = "black", linewidth = 1.12) +
  draw_label(
    "All cells",
    x = 0.371,
    y = 0.932,
    hjust = 0.5,
    vjust = 0.5,
    fontface = "bold",
    size = 20.0
  ) +
  draw_line(c(0.446, 0.513), c(0.932, 0.932), color = "black", linewidth = 1.12) +
  draw_plot(legend_panel, x = 0.582, y = 0.256, width = 0.368, height = 0.552)

out_png <- file.path(output_dir, "all_cells_umap_celltype_major_tcells_merged_template_fit.png")
out_tiff <- file.path(output_dir, "all_cells_umap_celltype_major_tcells_merged_template_fit.tiff")

ggsave(
  filename = out_png,
  plot = final_plot,
  width = 9.14,
  height = 5.24,
  dpi = 100,
  bg = "white",
  device = ragg::agg_png
)

ggsave(
  filename = out_tiff,
  plot = final_plot,
  width = 9.14,
  height = 5.24,
  dpi = 300,
  bg = "white",
  compression = "lzw"
)

cat("Cells plotted:", nrow(plot_df), "\n")
cat("Legend groups:", length(levels(plot_df$plot_celltype)), "\n")
cat("Saved:\n")
cat(out_png, "\n")
cat(out_tiff, "\n")
cat(file.path(output_dir, "all_cells_celltype_major_tcells_merged_counts_and_palette.csv"), "\n")
