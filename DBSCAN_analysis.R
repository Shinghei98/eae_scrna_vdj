################################################################################
# MGI0279 B-cell DBSCAN and program analysis
################################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(ggplot2)
  library(dbscan)
})

################################################################################
# 1. Paths and inputs
################################################################################

path.home <- "/storage3/fs1/gfwu/Active/David/mng_dcln_project"
output_dir <- file.path(path.home, "eae_scrna_vdj_preprocess")

annotated_rds <- file.path(output_dir, "srt_fullannot.rds")
if (!file.exists(annotated_rds)) {
  stop("Annotated object not found: ", annotated_rds, call. = FALSE)
}

dbscan_dir <- file.path(output_dir, "DBSCAN", "bcell_tissuewise_eps017_min5_min31")
program_out_dir <- file.path(dbscan_dir, "bcell_program_signatures_and_overlap")
if (!dir.exists(dbscan_dir)) dir.create(dbscan_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(program_out_dir)) dir.create(program_out_dir, recursive = TRUE, showWarnings = FALSE)

full_with_dbscan_rds <- file.path(output_dir, "srt_fullannot_with_bcell_dbscan_eps017_min5_min31.rds")
cluster_col <- "tissue_dbscan_cluster_eps017_min5_min31"
cluster_simple_col <- "tissue_dbscan_cluster_simple_eps017_min5_min31"
tissue_col <- "TissueGroup"
bcell_minor_backup_col <- "celltype_minor_pre_bcell_dbscan"
min_signature_genes <- 50
eps_use <- 0.17
minPts_use <- 5
min_cluster_size <- 30

################################################################################
# 2. Helpers
################################################################################

strict_ig_regex <- paste0(
  "^(Ighv|Ighd($|[0-9-])|Ighj($|[0-9-])|",
  "Igha$|Ighe$|Ighm$|Ighg[0-9a-z]*$|",
  "Igkv|Igkj($|[0-9-])|Igkc$|",
  "Iglv|Iglj($|[0-9-])|Iglc($|[0-9-])|",
  "Igll)"
)

save_tiff_plot <- function(plot_obj,
                           out_file,
                           width = 7,
                           height = 6,
                           units = "in",
                           res = 300) {
  if (!dir.exists(dirname(out_file))) {
    dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  }
  tiff(out_file, width = width, height = height, units = units, res = res, compression = "lzw")
  print(plot_obj)
  dev.off()
  out_file
}

build_noig_umap <- function(obj,
                            npcs = 15,
                            dims_use = 1:15,
                            nfeatures = 2000,
                            umap_name = "umap_noig_fixed") {
  DefaultAssay(obj) <- "RNA"

  ig_genes_strict <- rownames(obj)[grepl(strict_ig_regex, rownames(obj))]

  set.seed(1234)
  obj <- FindVariableFeatures(
    obj,
    selection.method = "vst",
    nfeatures = nfeatures,
    verbose = FALSE
  )

  vf_no_ig <- setdiff(VariableFeatures(obj), ig_genes_strict)
  VariableFeatures(obj) <- vf_no_ig

  set.seed(1234)
  obj <- ScaleData(obj, features = vf_no_ig, verbose = FALSE)

  set.seed(1234)
  obj <- RunPCA(
    obj,
    features = vf_no_ig,
    npcs = npcs,
    verbose = FALSE
  )

  set.seed(1234)
  obj <- RunUMAP(
    obj,
    reduction = "pca",
    dims = dims_use,
    reduction.name = umap_name,
    reduction.key = "UMAPNOIG_",
    verbose = FALSE
  )

  obj
}

renumber_dbscan_by_tissue <- function(cluster_chr) {
  out <- cluster_chr
  keep <- sort(unique(cluster_chr[cluster_chr != "DBSCAN_0"]))
  new <- paste0("DBSCAN_", seq_along(keep))
  names(new) <- keep
  out[cluster_chr != "DBSCAN_0"] <- unname(new[cluster_chr[cluster_chr != "DBSCAN_0"]])
  out
}

################################################################################
# 3. Validated B-cell object and tissue-wise DBSCAN
################################################################################

full_obj <- readRDS(annotated_rds)
DefaultAssay(full_obj) <- "RNA"

if (!tissue_col %in% colnames(full_obj@meta.data)) {
  full_obj[[tissue_col]] <- dplyr::case_when(
    grepl("_L", full_obj$sample_id) ~ "dCLN",
    grepl("_M", full_obj$sample_id) ~ "MNG",
    TRUE ~ NA_character_
  )
}

validated_bcell_obj <- subset(
  full_obj,
  subset = celltype_major == "B_cell" & celltype_minor == "B_cell"
)
validated_bcell_obj <- subset(
  validated_bcell_obj,
  cells = colnames(validated_bcell_obj)[!is.na(validated_bcell_obj[[tissue_col, drop = TRUE]])]
)

cat("\n====================\n")
cat("Validated B cells for DBSCAN\n")
cat("====================\n")
print(table(validated_bcell_obj[[tissue_col, drop = TRUE]], useNA = "ifany"))

tissue_objs <- lapply(c("MNG", "dCLN"), function(tissue_name) {
  obj_tissue <- subset(
    validated_bcell_obj,
    cells = colnames(validated_bcell_obj)[validated_bcell_obj[[tissue_col, drop = TRUE]] == tissue_name]
  )
  build_noig_umap(
    obj_tissue,
    npcs = 15,
    dims_use = 1:15,
    nfeatures = 2000,
    umap_name = "umap_noig_fixed"
  )
})
names(tissue_objs) <- c("MNG", "dCLN")

eps_grid <- sort(unique(c(
  0.05, 0.08, 0.10, 0.12, 0.15, 0.18, 0.20, 0.25, 0.30,
  0.16, 0.17, 0.175
)))

summarize_dbscan_result <- function(tissue_name, eps_value, coords) {
  set.seed(1234)
  db <- dbscan::dbscan(coords, eps = eps_value, minPts = minPts_use)

  raw_cluster <- paste0("DBSCAN_", db$cluster)
  raw_size_tab <- table(raw_cluster)
  raw_non_noise <- raw_size_tab[names(raw_size_tab) != "DBSCAN_0"]

  small_clusters <- names(raw_size_tab)[
    names(raw_size_tab) != "DBSCAN_0" & raw_size_tab <= min_cluster_size
  ]
  posthoc_cluster <- ifelse(
    raw_cluster %in% small_clusters,
    "DBSCAN_0",
    raw_cluster
  )

  posthoc_size_tab <- table(posthoc_cluster)
  posthoc_non_noise <- posthoc_size_tab[names(posthoc_size_tab) != "DBSCAN_0"]

  summary_df <- data.frame(
    TissueGroup = tissue_name,
    eps = eps_value,
    minPts = minPts_use,
    n_cells = length(raw_cluster),
    raw_n_noise = unname(ifelse("DBSCAN_0" %in% names(raw_size_tab), raw_size_tab[["DBSCAN_0"]], 0)),
    raw_noise_pct = round(100 * mean(raw_cluster == "DBSCAN_0"), 1),
    raw_n_clusters = length(raw_non_noise),
    raw_min_cluster_n = ifelse(length(raw_non_noise) > 0, min(raw_non_noise), NA),
    raw_median_cluster_n = ifelse(length(raw_non_noise) > 0, median(raw_non_noise), NA),
    raw_max_cluster_n = ifelse(length(raw_non_noise) > 0, max(raw_non_noise), NA),
    n_small_clusters_collapsed = length(small_clusters),
    posthoc_n_noise = unname(ifelse("DBSCAN_0" %in% names(posthoc_size_tab), posthoc_size_tab[["DBSCAN_0"]], 0)),
    posthoc_noise_pct = round(100 * mean(posthoc_cluster == "DBSCAN_0"), 1),
    posthoc_n_clusters = length(posthoc_non_noise),
    posthoc_min_cluster_n = ifelse(length(posthoc_non_noise) > 0, min(posthoc_non_noise), NA),
    posthoc_median_cluster_n = ifelse(length(posthoc_non_noise) > 0, median(posthoc_non_noise), NA),
    posthoc_max_cluster_n = ifelse(length(posthoc_non_noise) > 0, max(posthoc_non_noise), NA),
    stringsAsFactors = FALSE
  )

  cluster_size_df <- data.frame(
    TissueGroup = tissue_name,
    eps = eps_value,
    minPts = minPts_use,
    dbscan_raw = names(raw_size_tab),
    raw_n_cells = as.integer(raw_size_tab),
    stringsAsFactors = FALSE
  ) |>
    dplyr::left_join(
      data.frame(
        dbscan_posthoc = names(posthoc_size_tab),
        posthoc_n_cells = as.integer(posthoc_size_tab),
        stringsAsFactors = FALSE
      ),
      by = c("dbscan_raw" = "dbscan_posthoc")
    ) |>
    dplyr::mutate(
      collapsed_to_noise = dbscan_raw %in% small_clusters
    )

  list(summary = summary_df, cluster_sizes = cluster_size_df)
}

scan_results <- lapply(names(tissue_objs), function(tissue_name) {
  coords <- Embeddings(tissue_objs[[tissue_name]], "umap_noig_fixed")[, 1:2, drop = FALSE]
  lapply(eps_grid, function(eps_value) {
    summarize_dbscan_result(tissue_name, eps_value, coords)
  })
}) |>
  unlist(recursive = FALSE)

dbscan_eps_scan_summary <- dplyr::bind_rows(lapply(scan_results, `[[`, "summary")) |>
  dplyr::arrange(TissueGroup, eps)
dbscan_eps_scan_cluster_sizes <- dplyr::bind_rows(lapply(scan_results, `[[`, "cluster_sizes")) |>
  dplyr::arrange(TissueGroup, eps, dbscan_raw)

cat("\n====================\n")
cat("DBSCAN eps sensitivity scan\n")
cat("====================\n")
print(dbscan_eps_scan_summary, row.names = FALSE)

write.csv(
  dbscan_eps_scan_summary,
  file.path(dbscan_dir, "bcell_dbscan_eps_scan_summary.csv"),
  row.names = FALSE
)
write.csv(
  dbscan_eps_scan_cluster_sizes,
  file.path(dbscan_dir, "bcell_dbscan_eps_scan_cluster_sizes.csv"),
  row.names = FALSE
)

# Selected after scanning DBSCAN parameters for overfragmentation/noise at
# smaller eps values, collapse into dominant clusters at larger eps values, and
# recovery of the prior MNG_DBSCAN_2 B-cell program.

assignment_list <- lapply(names(tissue_objs), function(tissue_name) {
  obj_tissue <- tissue_objs[[tissue_name]]
  coords <- Embeddings(obj_tissue, "umap_noig_fixed")[, 1:2, drop = FALSE]

  set.seed(1234)
  db <- dbscan::dbscan(coords, eps = eps_use, minPts = minPts_use)

  raw_cluster <- paste0("DBSCAN_", db$cluster)
  size_tab <- table(raw_cluster)
  small_clusters <- names(size_tab)[
    names(size_tab) != "DBSCAN_0" & size_tab <= min_cluster_size
  ]
  posthoc_cluster <- ifelse(
    raw_cluster %in% small_clusters,
    "DBSCAN_0",
    raw_cluster
  )

  data.frame(
    cell = colnames(obj_tissue),
    TissueGroup = tissue_name,
    UMAP_1 = coords[, 1],
    UMAP_2 = coords[, 2],
    dbscan_raw = raw_cluster,
    dbscan_posthoc = posthoc_cluster,
    stringsAsFactors = FALSE
  )
})

dbscan_assignments <- dplyr::bind_rows(assignment_list) |>
  dplyr::group_by(TissueGroup) |>
  dplyr::mutate(
    dbscan_posthoc_renumbered = renumber_dbscan_by_tissue(dbscan_posthoc),
    tissue_dbscan_cluster_simple = dplyr::if_else(
      dbscan_posthoc_renumbered == "DBSCAN_0",
      "DBSCAN_noise",
      dbscan_posthoc_renumbered
    ),
    tissue_dbscan_cluster = dplyr::if_else(
      dbscan_posthoc_renumbered == "DBSCAN_0",
      paste0(TissueGroup, "_DBSCAN_noise"),
      paste0(TissueGroup, "_", dbscan_posthoc_renumbered)
    )
  ) |>
  dplyr::ungroup() |>
  dplyr::rename(
    !!cluster_simple_col := tissue_dbscan_cluster_simple,
    !!cluster_col := tissue_dbscan_cluster
  )

dbscan_summary <- dbscan_assignments |>
  dplyr::count(TissueGroup, .data[[cluster_col]], name = "n_cells") |>
  dplyr::group_by(TissueGroup) |>
  dplyr::mutate(cluster_fraction = n_cells / sum(n_cells)) |>
  dplyr::ungroup() |>
  dplyr::arrange(TissueGroup, dplyr::desc(n_cells))

validated_bcell_obj[[cluster_simple_col]] <- NA_character_
validated_bcell_obj[[cluster_col]] <- NA_character_
validated_bcell_obj@meta.data[dbscan_assignments$cell, cluster_simple_col] <-
  dbscan_assignments[[cluster_simple_col]]
validated_bcell_obj@meta.data[dbscan_assignments$cell, cluster_col] <-
  dbscan_assignments[[cluster_col]]

full_obj[[cluster_simple_col]] <- NA_character_
full_obj[[cluster_col]] <- NA_character_
full_obj@meta.data[dbscan_assignments$cell, cluster_simple_col] <-
  dbscan_assignments[[cluster_simple_col]]
full_obj@meta.data[dbscan_assignments$cell, cluster_col] <-
  dbscan_assignments[[cluster_col]]

################################################################################
# 4. Map DBSCAN labels back to B-cell celltype_minor in the global object
################################################################################

if (anyDuplicated(dbscan_assignments$cell) > 0) {
  stop("Duplicated cell barcodes found in DBSCAN assignments.", call. = FALSE)
}

missing_from_full <- setdiff(dbscan_assignments$cell, colnames(full_obj))
missing_from_bcell <- setdiff(dbscan_assignments$cell, colnames(validated_bcell_obj))
if (length(missing_from_full) > 0 || length(missing_from_bcell) > 0) {
  stop(
    "DBSCAN assignments do not align with saved objects. Missing from full object: ",
    length(missing_from_full),
    "; missing from B-cell object: ",
    length(missing_from_bcell),
    call. = FALSE
  )
}

if (!bcell_minor_backup_col %in% colnames(full_obj@meta.data)) {
  full_obj[[bcell_minor_backup_col]] <- full_obj$celltype_minor
}
if (!bcell_minor_backup_col %in% colnames(validated_bcell_obj@meta.data)) {
  validated_bcell_obj[[bcell_minor_backup_col]] <- validated_bcell_obj$celltype_minor
}

full_obj@meta.data[dbscan_assignments$cell, "celltype_minor"] <-
  dbscan_assignments[[cluster_col]]
validated_bcell_obj@meta.data[dbscan_assignments$cell, "celltype_minor"] <-
  dbscan_assignments[[cluster_col]]

if (any(is.na(full_obj@meta.data[dbscan_assignments$cell, "celltype_minor"]))) {
  stop("Some B cells still have missing celltype_minor after DBSCAN mapping.", call. = FALSE)
}
if (!all(full_obj@meta.data[dbscan_assignments$cell, "celltype_minor"] ==
         full_obj@meta.data[dbscan_assignments$cell, cluster_col])) {
  stop("B-cell celltype_minor does not match DBSCAN labels after mapping.", call. = FALSE)
}

cat("\n====================\n")
cat("B-cell celltype_minor after DBSCAN mapping\n")
cat("====================\n")
print(table(full_obj@meta.data[dbscan_assignments$cell, "celltype_minor"], useNA = "ifany"))

bcell_metadata_with_dbscan <- full_obj@meta.data[dbscan_assignments$cell, , drop = FALSE]
bcell_metadata_with_dbscan$cell <- rownames(bcell_metadata_with_dbscan)

write.csv(
  dbscan_assignments,
  file.path(dbscan_dir, "bcell_dbscan_eps017_min5_min31_assignments.csv"),
  row.names = FALSE
)
write.csv(
  dbscan_summary,
  file.path(dbscan_dir, "bcell_dbscan_eps017_min5_min31_summary.csv"),
  row.names = FALSE
)
write.csv(
  bcell_metadata_with_dbscan,
  file.path(dbscan_dir, "validated_bcell_metadata_with_dbscan_eps017_min5_min31.csv"),
  row.names = FALSE
)
saveRDS(
  validated_bcell_obj,
  file.path(dbscan_dir, "validated_bcell_obj_tissue_dbscan_eps017_min5_min31.rds")
)
saveRDS(
  full_obj,
  file.path(dbscan_dir, "srt_fullannot_with_bcell_dbscan_eps017_min5_min31.rds")
)
saveRDS(
  full_obj,
  full_with_dbscan_rds
)
saveRDS(
  tissue_objs,
  file.path(dbscan_dir, "tissuewise_noig_umap_objects_eps017_min5_min31.rds")
)

dbscan_plot_df <- dbscan_assignments |>
  dplyr::mutate(
    draw_order = ifelse(grepl("DBSCAN_noise$", .data[[cluster_col]]), 1, 2)
  ) |>
  dplyr::arrange(draw_order)

dbscan_plot <- ggplot(
  dbscan_plot_df,
  aes(UMAP_1, UMAP_2, color = .data[[cluster_col]])
) +
  geom_point(size = 0.2) +
  facet_wrap(~ TissueGroup, scales = "free") +
  theme_classic() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  ) +
  ggtitle(NULL)

save_tiff_plot(
  dbscan_plot,
  file.path(dbscan_dir, "bcell_tissuewise_dbscan_eps017_min5_min31.tiff"),
  width = 8,
  height = 4
)

################################################################################
# 5. DBSCAN programs to test
################################################################################

obj <- tryCatch(
  JoinLayers(validated_bcell_obj, assay = "RNA"),
  error = function(e) validated_bcell_obj
)
DefaultAssay(obj) <- "RNA"

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
  file.path(program_out_dir, "dbscan_programs_tested.csv"),
  row.names = FALSE
)

################################################################################
# 6. Top-50 signatures within tissue
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
  file.path(program_out_dir, "dbscan_program_signatures_top50_ttest.csv"),
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
  file.path(program_out_dir, "dbscan_program_signature_sizes.csv"),
  row.names = FALSE
)

program_meta_use <- program_meta |>
  dplyr::filter(n_signature_genes == min_signature_genes)

excluded_programs <- program_meta |>
  dplyr::filter(n_signature_genes != min_signature_genes)

write.csv(
  excluded_programs,
  file.path(program_out_dir, "dbscan_programs_excluded_from_overlap_clustering.csv"),
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
# 7. Pairwise overlap, Fisher tests, and hierarchical clustering
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
  file.path(program_out_dir, "dbscan_program_pairwise_overlap_50minusdistance.csv"),
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
  file.path(program_out_dir, "dbscan_program_top50_fisher_overlap_full.csv"),
  row.names = FALSE
)

write.csv(
  fisher_res |>
    dplyr::filter(significant_p_lt_0_001),
  file.path(program_out_dir, "dbscan_program_top50_fisher_overlap_significant_p_lt_0.001.csv"),
  row.names = FALSE
)

write.csv(
  fisher_res |>
    dplyr::filter(cross_tissue, significant_fdr_lt_0_05),
  file.path(program_out_dir, "dbscan_program_top50_fisher_overlap_cross_tissue_significant_fdr_lt_0.05.csv"),
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

write.csv(overlap_mat[ordered_ids, ordered_ids], file.path(program_out_dir, "dbscan_program_overlap_matrix.csv"))
write.csv(dist_mat[ordered_ids, ordered_ids], file.path(program_out_dir, "dbscan_program_distance_matrix_50minusoverlap.csv"))
write.csv(pval_mat[ordered_ids, ordered_ids], file.path(program_out_dir, "dbscan_program_fisher_pvalue_matrix.csv"))
write.csv(fdr_mat[ordered_ids, ordered_ids], file.path(program_out_dir, "dbscan_program_fisher_fdr_matrix.csv"))
saveRDS(hc, file.path(program_out_dir, "dbscan_program_hclust_50minusoverlap.rds"))

program_map <- program_meta_use |>
  dplyr::mutate(program_id = factor(program_id, levels = ordered_ids)) |>
  dplyr::arrange(program_id) |>
  dplyr::mutate(program_number = dplyr::row_number()) |>
  dplyr::select(program_number, program_id, TissueGroup, cluster, n_cells_cluster, n_signature_genes)

write.csv(
  program_map,
  file.path(program_out_dir, "dbscan_program_order_mapping.csv"),
  row.names = FALSE
)

################################################################################
# 8. Shared versus tissue-restricted program calls
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
  file.path(program_out_dir, "dbscan_program_shared_vs_tissue_restricted.csv"),
  row.names = FALSE
)

################################################################################
# 9. Plots
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

save_tiff_plot(
  heatmap_plot,
  file.path(program_out_dir, "dbscan_program_overlap_heatmap.tiff"),
  width = 8,
  height = 7
)

tiff(
  file.path(program_out_dir, "dbscan_program_hclust_dendrogram.tiff"),
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
cat("Saved B-cell DBSCAN outputs to:\n")
cat(dbscan_dir, "\n")
cat("Saved full annotated object with B-cell DBSCAN labels to:\n")
cat(full_with_dbscan_rds, "\n")
cat("Saved DBSCAN program analysis outputs to:\n")
cat(program_out_dir, "\n")
