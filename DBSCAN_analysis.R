################################################################################
# MGI0279 B-cell DBSCAN program analysis
################################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(ggplot2)
})

################################################################################
# 1. Paths and inputs
################################################################################

path.home <- "/storage3/fs1/gfwu/Active/David/mng_dcln_project"
output_dir <- file.path(path.home, "eae_scrna_vdj_preprocess")

dbscan_dir <- file.path(output_dir, "DBSCAN", "bcell_tissuewise_eps020_min10_min21")
obj_file <- file.path(dbscan_dir, "validated_bcell_obj_tissue_dbscan_eps020_min10_min21.rds")

if (!file.exists(obj_file)) {
  stop("Input RDS not found: ", obj_file, call. = FALSE)
}

out_dir <- file.path(output_dir, "DBSCAN", "bcell_program_signatures_and_overlap")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cluster_col <- "tissue_dbscan_cluster_eps020_min10_min21"
tissue_col <- "TissueGroup"
min_signature_genes <- 50

obj <- readRDS(obj_file)
DefaultAssay(obj) <- "RNA"

if (!tissue_col %in% colnames(obj@meta.data)) {
  obj[[tissue_col]] <- dplyr::case_when(
    grepl("_L", obj$sample_id) ~ "dCLN",
    grepl("_M", obj$sample_id) ~ "MNG",
    TRUE ~ NA_character_
  )
}

if (!cluster_col %in% colnames(obj@meta.data)) {
  stop("Cluster column not found: ", cluster_col, call. = FALSE)
}

################################################################################
# 2. DBSCAN programs to test
################################################################################

job_df <- obj@meta.data |>
  dplyr::mutate(
    TissueGroup = .data[[tissue_col]],
    cluster = .data[[cluster_col]]
  ) |>
  dplyr::filter(
    !is.na(TissueGroup),
    !is.na(cluster),
    !grepl("DBSCAN_noise$", cluster)
  ) |>
  dplyr::count(TissueGroup, cluster, name = "n_cells") |>
  dplyr::arrange(TissueGroup, dplyr::desc(n_cells))

if (nrow(job_df) == 0) {
  stop("No retained DBSCAN clusters found in ", cluster_col, call. = FALSE)
}

cat("\n====================\n")
cat("DBSCAN programs to test\n")
cat("====================\n")
print(job_df, row.names = FALSE)

write.csv(
  job_df,
  file.path(out_dir, "dbscan_programs_tested.csv"),
  row.names = FALSE
)

################################################################################
# 3. Top-50 signatures within tissue
################################################################################

run_one_signature <- function(i, obj, job_df, tissue_col, cluster_col) {
  grp <- job_df$TissueGroup[i]
  cl <- job_df$cluster[i]
  n_cl <- job_df$n_cells[i]

  cells_grp <- colnames(obj)[obj[[tissue_col, drop = TRUE]] == grp]
  obj_grp <- subset(obj, cells = cells_grp)
  Idents(obj_grp) <- obj_grp[[cluster_col, drop = TRUE]]

  deg <- FindMarkers(
    obj_grp,
    ident.1 = cl,
    assay = "RNA",
    test.use = "t",
    only.pos = TRUE,
    logfc.threshold = 0,
    min.pct = 0,
    verbose = FALSE
  )

  if (nrow(deg) == 0) return(data.frame())

  deg$gene <- rownames(deg)
  fc_col <- intersect(c("avg_log2FC", "avg_logFC"), colnames(deg))[1]
  if (is.na(fc_col)) stop("Could not find logFC column in FindMarkers output.")

  deg |>
    dplyr::filter(.data[[fc_col]] >= 1, p_val < 0.001) |>
    dplyr::arrange(dplyr::desc(.data[[fc_col]])) |>
    dplyr::slice_head(n = min_signature_genes) |>
    dplyr::mutate(
      TissueGroup = grp,
      cluster = cl,
      n_cells_cluster = n_cl,
      program_id = cl,
      rank = dplyr::row_number(),
      .before = 1
    )
}

signature_list <- lapply(seq_len(nrow(job_df)), function(i) {
  message("Running signature test for ", job_df$cluster[i])
  run_one_signature(i, obj, job_df, tissue_col, cluster_col)
})

signature_df <- dplyr::bind_rows(signature_list)

if (nrow(signature_df) == 0) {
  stop("No program signatures passed avg_log2FC >= 1 and raw P < 0.001.", call. = FALSE)
}

write.csv(
  signature_df,
  file.path(out_dir, "dbscan_program_signatures_top50_ttest.csv"),
  row.names = FALSE
)

program_meta <- signature_df |>
  dplyr::group_by(program_id) |>
  dplyr::summarise(
    TissueGroup = dplyr::first(TissueGroup),
    cluster = dplyr::first(cluster),
    n_cells_cluster = dplyr::first(n_cells_cluster),
    n_signature_genes = dplyr::n(),
    .groups = "drop"
  ) |>
  dplyr::arrange(TissueGroup, cluster)

write.csv(
  program_meta,
  file.path(out_dir, "dbscan_program_signature_sizes.csv"),
  row.names = FALSE
)

program_meta_use <- program_meta |>
  dplyr::filter(n_signature_genes == min_signature_genes)

excluded_programs <- program_meta |>
  dplyr::filter(n_signature_genes != min_signature_genes)

write.csv(
  excluded_programs,
  file.path(out_dir, "dbscan_programs_excluded_from_overlap_clustering.csv"),
  row.names = FALSE
)

if (nrow(program_meta_use) < 2) {
  stop("Fewer than 2 programs have full 50-gene signatures.", call. = FALSE)
}

signature_use <- signature_df |>
  dplyr::filter(program_id %in% program_meta_use$program_id)

program_gene_list <- split(signature_use$gene, signature_use$program_id)
program_gene_list <- lapply(program_gene_list, unique)
program_ids <- program_meta_use$program_id

################################################################################
# 4. Pairwise overlap, Fisher tests, and hierarchical clustering
################################################################################

pair_index <- utils::combn(program_ids, 2, simplify = FALSE)
universe_genes <- rownames(obj)
universe_size <- length(universe_genes)

pair_res <- dplyr::bind_rows(lapply(pair_index, function(pair_ids) {
  p1 <- pair_ids[1]
  p2 <- pair_ids[2]
  g1 <- program_gene_list[[p1]]
  g2 <- program_gene_list[[p2]]
  shared <- sort(intersect(g1, g2))

  data.frame(
    program_1 = p1,
    program_2 = p2,
    tissue_1 = program_meta_use$TissueGroup[match(p1, program_meta_use$program_id)],
    tissue_2 = program_meta_use$TissueGroup[match(p2, program_meta_use$program_id)],
    overlap_n = length(shared),
    distance_50_minus_overlap = min_signature_genes - length(shared),
    shared_genes = paste(shared, collapse = "/"),
    stringsAsFactors = FALSE
  )
})) |>
  dplyr::arrange(distance_50_minus_overlap, dplyr::desc(overlap_n))

write.csv(
  pair_res,
  file.path(out_dir, "dbscan_program_pairwise_overlap_50minusdistance.csv"),
  row.names = FALSE
)

fisher_res <- dplyr::bind_rows(lapply(pair_index, function(pair_ids) {
  p1 <- pair_ids[1]
  p2 <- pair_ids[2]
  g1 <- intersect(program_gene_list[[p1]], universe_genes)
  g2 <- intersect(program_gene_list[[p2]], universe_genes)

  a <- length(intersect(g1, g2))
  b <- length(setdiff(g1, g2))
  c <- length(setdiff(g2, g1))
  d <- universe_size - a - b - c

  ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "greater")

  data.frame(
    program_1 = p1,
    program_2 = p2,
    tissue_1 = program_meta_use$TissueGroup[match(p1, program_meta_use$program_id)],
    tissue_2 = program_meta_use$TissueGroup[match(p2, program_meta_use$program_id)],
    overlap_top50 = a,
    fisher_p_value = ft$p.value,
    cross_tissue = program_meta_use$TissueGroup[match(p1, program_meta_use$program_id)] !=
      program_meta_use$TissueGroup[match(p2, program_meta_use$program_id)],
    stringsAsFactors = FALSE
  )
})) |>
  dplyr::mutate(
    fdr_bh = p.adjust(fisher_p_value, method = "BH"),
    significant_p_lt_0_001 = fisher_p_value < 0.001,
    significant_fdr_lt_0_05 = fdr_bh < 0.05
  ) |>
  dplyr::arrange(fisher_p_value, dplyr::desc(overlap_top50))

write.csv(
  fisher_res,
  file.path(out_dir, "dbscan_program_top50_fisher_overlap_full.csv"),
  row.names = FALSE
)

write.csv(
  fisher_res |>
    dplyr::filter(significant_p_lt_0_001),
  file.path(out_dir, "dbscan_program_top50_fisher_overlap_significant_p_lt_0.001.csv"),
  row.names = FALSE
)

write.csv(
  fisher_res |>
    dplyr::filter(cross_tissue, significant_fdr_lt_0_05),
  file.path(out_dir, "dbscan_program_top50_fisher_overlap_cross_tissue_significant_fdr_lt_0.05.csv"),
  row.names = FALSE
)

overlap_mat <- matrix(0, length(program_ids), length(program_ids), dimnames = list(program_ids, program_ids))
diag(overlap_mat) <- min_signature_genes

dist_mat <- matrix(0, length(program_ids), length(program_ids), dimnames = list(program_ids, program_ids))
pval_mat <- matrix(1, length(program_ids), length(program_ids), dimnames = list(program_ids, program_ids))
fdr_mat <- matrix(1, length(program_ids), length(program_ids), dimnames = list(program_ids, program_ids))
diag(pval_mat) <- 0
diag(fdr_mat) <- 0

for (i in seq_len(nrow(pair_res))) {
  p1 <- pair_res$program_1[i]
  p2 <- pair_res$program_2[i]
  overlap_mat[p1, p2] <- pair_res$overlap_n[i]
  overlap_mat[p2, p1] <- pair_res$overlap_n[i]
  dist_mat[p1, p2] <- pair_res$distance_50_minus_overlap[i]
  dist_mat[p2, p1] <- pair_res$distance_50_minus_overlap[i]
}

for (i in seq_len(nrow(fisher_res))) {
  p1 <- fisher_res$program_1[i]
  p2 <- fisher_res$program_2[i]
  pval_mat[p1, p2] <- fisher_res$fisher_p_value[i]
  pval_mat[p2, p1] <- fisher_res$fisher_p_value[i]
  fdr_mat[p1, p2] <- fisher_res$fdr_bh[i]
  fdr_mat[p2, p1] <- fisher_res$fdr_bh[i]
}

hc <- hclust(as.dist(dist_mat), method = "average")
ordered_ids <- hc$labels[hc$order]

write.csv(overlap_mat[ordered_ids, ordered_ids], file.path(out_dir, "dbscan_program_overlap_matrix.csv"))
write.csv(dist_mat[ordered_ids, ordered_ids], file.path(out_dir, "dbscan_program_distance_matrix_50minusoverlap.csv"))
write.csv(pval_mat[ordered_ids, ordered_ids], file.path(out_dir, "dbscan_program_fisher_pvalue_matrix.csv"))
write.csv(fdr_mat[ordered_ids, ordered_ids], file.path(out_dir, "dbscan_program_fisher_fdr_matrix.csv"))
saveRDS(hc, file.path(out_dir, "dbscan_program_hclust_50minusoverlap.rds"))

program_map <- program_meta_use |>
  dplyr::mutate(program_id = factor(program_id, levels = ordered_ids)) |>
  dplyr::arrange(program_id) |>
  dplyr::mutate(program_number = dplyr::row_number()) |>
  dplyr::select(program_number, program_id, TissueGroup, cluster, n_cells_cluster, n_signature_genes)

write.csv(
  program_map,
  file.path(out_dir, "dbscan_program_order_mapping.csv"),
  row.names = FALSE
)

################################################################################
# 5. Shared versus tissue-restricted program calls
################################################################################

nearest_cross_tissue <- dplyr::bind_rows(lapply(program_ids, function(pid) {
  partners <- pair_res |>
    dplyr::filter(program_1 == pid | program_2 == pid) |>
    dplyr::mutate(
      partner = ifelse(program_1 == pid, program_2, program_1),
      partner_tissue = program_meta_use$TissueGroup[match(partner, program_meta_use$program_id)]
    ) |>
    dplyr::filter(partner_tissue != program_meta_use$TissueGroup[match(pid, program_meta_use$program_id)]) |>
    dplyr::arrange(distance_50_minus_overlap, dplyr::desc(overlap_n))

  if (nrow(partners) == 0) {
    return(data.frame(program_id = pid, nearest_cross_tissue_program = NA_character_,
                      nearest_cross_tissue_overlap = NA_integer_,
                      nearest_cross_tissue_distance = NA_integer_))
  }

  data.frame(
    program_id = pid,
    nearest_cross_tissue_program = partners$partner[1],
    nearest_cross_tissue_overlap = partners$overlap_n[1],
    nearest_cross_tissue_distance = partners$distance_50_minus_overlap[1],
    stringsAsFactors = FALSE
  )
}))

shared_partner_df <- fisher_res |>
  dplyr::filter(cross_tissue, significant_fdr_lt_0_05) |>
  dplyr::select(program_1, program_2, overlap_top50, fisher_p_value, fdr_bh)

program_classification <- program_meta_use |>
  dplyr::left_join(nearest_cross_tissue, by = "program_id") |>
  dplyr::rowwise() |>
  dplyr::mutate(
    significant_cross_tissue_partners = paste(unique(c(
      shared_partner_df$program_2[shared_partner_df$program_1 == program_id],
      shared_partner_df$program_1[shared_partner_df$program_2 == program_id]
    )), collapse = "/"),
    program_class = ifelse(
      significant_cross_tissue_partners == "",
      "tissue_restricted",
      "shared_cross_tissue"
    )
  ) |>
  dplyr::ungroup()

write.csv(
  program_classification,
  file.path(out_dir, "dbscan_program_shared_vs_tissue_restricted.csv"),
  row.names = FALSE
)

################################################################################
# 6. Plots
################################################################################

overlap_df <- as.data.frame(as.table(overlap_mat[ordered_ids, ordered_ids]))
colnames(overlap_df) <- c("program_1", "program_2", "overlap_n")
overlap_df$program_1 <- factor(overlap_df$program_1, levels = ordered_ids)
overlap_df$program_2 <- factor(overlap_df$program_2, levels = rev(ordered_ids))

heatmap_plot <- ggplot(overlap_df, aes(program_1, program_2, fill = overlap_n)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_gradientn(
    colors = c("white", "#FFF4A3", "#FDB366", "#F46D43", "#9C179E", "#1B1B5F"),
    limits = c(0, min_signature_genes)
  ) +
  theme_classic() +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8)
  ) +
  labs(fill = "Shared genes")

tiff(
  file.path(out_dir, "dbscan_program_overlap_heatmap.tiff"),
  width = 8,
  height = 7,
  units = "in",
  res = 300,
  compression = "lzw"
)
print(heatmap_plot)
dev.off()

tiff(
  file.path(out_dir, "dbscan_program_hclust_dendrogram.tiff"),
  width = 8,
  height = 6,
  units = "in",
  res = 300,
  compression = "lzw"
)
plot(
  hc,
  main = "Hierarchical clustering of DBSCAN programs",
  xlab = "",
  sub = "",
  ylab = "Distance (50 - overlapping genes)"
)
dev.off()

cat("\n====================\n")
cat("Saved DBSCAN program analysis outputs to:\n")
cat(out_dir, "\n")
