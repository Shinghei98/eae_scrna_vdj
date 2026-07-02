################################################################################
# MGI0279 cell type annotation
################################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(ggplot2)
  library(dplyr)
})

################################################################################
# 1. Paths and input object
################################################################################

path.home <- "/storage3/fs1/gfwu/Active/David/mng_dcln_project"
output_dir <- file.path(path.home, "eae_scrna_vdj_preprocess")

input_rds <- file.path(output_dir, "srt_integrated_clustered_pruned_with_vdj_metadata.rds")
if (!file.exists(input_rds)) {
  input_rds <- file.path(output_dir, "MGI0279_L1_L2_M1_M2_integrated_seurat.rds")
}
if (!file.exists(input_rds)) {
  stop("Input RDS not found in ", output_dir, call. = FALSE)
}

srt_integrated_all_cells <- readRDS(input_rds)

################################################################################
# 2. Helpers
################################################################################

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
  path
}

save_tiff_plot <- function(plot_obj,
                           out_file,
                           width = 7,
                           height = 6,
                           units = "in",
                           res = 300) {
  ensure_dir(dirname(out_file))
  tiff(out_file, width = width, height = height, units = units, res = res)
  print(plot_obj)
  dev.off()
  out_file
}

plot_umap <- function(obj,
                      output_dir,
                      prefix,
                      group = "seurat_clusters",
                      label = TRUE) {
  out_file <- file.path(output_dir, "QC", "umap", paste0(prefix, "_umap.tiff"))
  p <- DimPlot(obj, reduction = "umap", group.by = group, label = label) +
    labs(x = NULL, y = NULL) +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank()
    ) +
    ggtitle(NULL)
  save_tiff_plot(p, out_file)
}

recluster_rna_subset <- function(obj,
                                 cluster_col,
                                 npcs,
                                 dims_use,
                                 k_param,
                                 resolution,
                                 nfeatures = 2000) {
  DefaultAssay(obj) <- "RNA"

  obj <- tryCatch(
    JoinLayers(obj, assay = "RNA"),
    error = function(e) obj
  )

  set.seed(1234)
  obj <- NormalizeData(obj, assay = "RNA", verbose = FALSE)

  set.seed(1234)
  obj <- FindVariableFeatures(
    obj,
    assay = "RNA",
    selection.method = "vst",
    nfeatures = nfeatures,
    verbose = FALSE
  )

  set.seed(1234)
  obj <- ScaleData(
    obj,
    assay = "RNA",
    features = VariableFeatures(obj),
    verbose = FALSE
  )

  set.seed(1234)
  obj <- RunPCA(
    obj,
    assay = "RNA",
    features = VariableFeatures(obj),
    npcs = npcs,
    verbose = FALSE
  )

  set.seed(1234)
  obj <- FindNeighbors(obj, dims = dims_use, k.param = k_param, verbose = FALSE)

  set.seed(1234)
  obj <- FindClusters(obj, resolution = resolution, verbose = FALSE)

  set.seed(1234)
  obj <- RunUMAP(
    obj,
    dims = dims_use,
    n.neighbors = k_param,
    seed.use = 1234,
    verbose = FALSE
  )

  obj[[cluster_col]] <- as.character(Idents(obj))
  obj
}

recluster_integrated_global <- function(obj,
                                        dims_use = 1:15,
                                        k_param = 20,
                                        resolution = 0.6) {
  DefaultAssay(obj) <- "integrated"

  if (!"pca" %in% names(obj@reductions)) {
    set.seed(1234)
    obj <- RunPCA(obj, npcs = max(dims_use), verbose = FALSE)
  }

  set.seed(1234)
  obj <- FindNeighbors(obj, dims = dims_use, k.param = k_param, verbose = FALSE)

  set.seed(1234)
  obj <- FindClusters(obj, resolution = resolution, verbose = FALSE)

  set.seed(1234)
  obj <- RunUMAP(obj, dims = dims_use, verbose = FALSE)

  obj
}

init_annotation_cols <- function(obj) {
  obj$celltype_major <- NA_character_
  obj$celltype_minor <- NA_character_
  obj$annotation_source <- NA_character_
  obj$annotation_cluster <- NA_character_
  obj$annotation_label <- NA_character_
  obj
}

set_cells_label <- function(obj, cells, major, minor) {
  cells <- intersect(cells, colnames(obj))
  obj@meta.data[cells, "celltype_major"] <- major
  obj@meta.data[cells, "celltype_minor"] <- minor
  obj
}

cells_in_clusters <- function(obj, cluster_col, clusters) {
  rownames(obj@meta.data)[as.character(obj@meta.data[[cluster_col]]) %in% clusters]
}

map_refined_cells <- function(obj,
                              cells,
                              major,
                              minor,
                              source,
                              cluster_name) {
  cells <- intersect(cells, colnames(obj))
  obj@meta.data[cells, "celltype_major"] <- major
  obj@meta.data[cells, "celltype_minor"] <- minor
  obj@meta.data[cells, "annotation_source"] <- source
  obj@meta.data[cells, "annotation_cluster"] <- cluster_name
  obj@meta.data[cells, "annotation_label"] <- paste(major, minor, sep = " / ")

  audit_row <- data.frame(
    source = source,
    cluster = cluster_name,
    celltype_major = major,
    celltype_minor = minor,
    n_cells_mapped = length(cells),
    stringsAsFactors = FALSE
  )
  annotation_audit_list[[length(annotation_audit_list) + 1]] <<- audit_row
  obj
}

################################################################################
# 3. Global clustering and broad annotation
################################################################################

annotation_out_dir <- output_dir
tcell_out_dir <- ensure_dir(file.path(annotation_out_dir, "tcell_annotation"))
myeloid_out_dir <- ensure_dir(file.path(annotation_out_dir, "myeloid_decision_tree"))
bcell_out_dir <- ensure_dir(file.path(annotation_out_dir, "bcell_annotation"))

srt_integrated_all_cells <- recluster_integrated_global(
  srt_integrated_all_cells,
  dims_use = 1:15,
  k_param = 20,
  resolution = 0.6
)
srt_integrated_all_cells <- init_annotation_cols(srt_integrated_all_cells)
srt_integrated_all_cells$TissueGroup <- dplyr::case_when(
  grepl("_M", srt_integrated_all_cells$sample_id) ~ "MNG",
  grepl("_L", srt_integrated_all_cells$sample_id) ~ "dCLN",
  TRUE ~ NA_character_
)

plot_umap(
  srt_integrated_all_cells,
  output_dir = annotation_out_dir,
  prefix = "annotation_global_pc15_k20_res06",
  group = "seurat_clusters",
  label = TRUE
)

annotation_audit_list <- list()

# All cells are first clustered into 18 clusters with PC=15, k_param=20,
# res=0.6. This yields B cells (C0, C2-4, C10), T cells (C1, C7-9,
# C11), neutrophils (C6, C15), cycling neutrophils (C14),
# macrophages/microglia (C5), monocytes (C16), NK cells (C12), DC
# (C13), and doublets with unknown identity (C17).
srt_integrated_all_cells <- set_cells_label(
  srt_integrated_all_cells,
  cells_in_clusters(srt_integrated_all_cells, "seurat_clusters", c("0", "2", "3", "4", "10")),
  "B_cell",
  "B_cell"
)
srt_integrated_all_cells <- set_cells_label(
  srt_integrated_all_cells,
  cells_in_clusters(srt_integrated_all_cells, "seurat_clusters", c("1", "7", "8", "9", "11")),
  "T_cell",
  "T_cell"
)
srt_integrated_all_cells <- set_cells_label(
  srt_integrated_all_cells,
  cells_in_clusters(srt_integrated_all_cells, "seurat_clusters", c("6", "15")),
  "Neutrophil",
  "Neutrophil"
)
srt_integrated_all_cells <- set_cells_label(
  srt_integrated_all_cells,
  cells_in_clusters(srt_integrated_all_cells, "seurat_clusters", "14"),
  "Neutrophil",
  "cycling_neutrophil"
)
srt_integrated_all_cells <- set_cells_label(
  srt_integrated_all_cells,
  cells_in_clusters(srt_integrated_all_cells, "seurat_clusters", "5"),
  "macrophage/microglia",
  "macrophage/microglia"
)
srt_integrated_all_cells <- set_cells_label(
  srt_integrated_all_cells,
  cells_in_clusters(srt_integrated_all_cells, "seurat_clusters", "16"),
  "monocyte",
  "monocyte"
)
srt_integrated_all_cells <- set_cells_label(
  srt_integrated_all_cells,
  cells_in_clusters(srt_integrated_all_cells, "seurat_clusters", "12"),
  "NK_cell",
  "NK_cell"
)
srt_integrated_all_cells <- set_cells_label(
  srt_integrated_all_cells,
  cells_in_clusters(srt_integrated_all_cells, "seurat_clusters", "13"),
  "DC",
  "DC"
)
srt_integrated_all_cells <- set_cells_label(
  srt_integrated_all_cells,
  cells_in_clusters(srt_integrated_all_cells, "seurat_clusters", "17"),
  "doublet",
  "unknown_identity"
)

################################################################################
# 4. T-cell annotation
################################################################################

# First cluster T cells (global C1, C7-9, C11) with PC=10, k_param=30,
# res=0.5. This validates T cells (C0-2, C4-5, C7, C9), B/T cell
# doublets (C3, C6), and neutrophil/T cell doublets (C8).
tcell_obj <- subset(
  srt_integrated_all_cells,
  cells = cells_in_clusters(srt_integrated_all_cells, "seurat_clusters", c("1", "7", "8", "9", "11"))
)

tcell_obj <- recluster_rna_subset(
  tcell_obj,
  cluster_col = "tcell_first_cluster_pc10_res05",
  npcs = 10,
  dims_use = 1:10,
  k_param = 30,
  resolution = 0.5,
  nfeatures = 2000
)

plot_umap(
  tcell_obj,
  output_dir = annotation_out_dir,
  prefix = "tcell_subset_pc10_k30_res05",
  group = "seurat_clusters",
  label = TRUE
)

saveRDS(tcell_obj, file.path(tcell_out_dir, "tcell_obj_pc10_k30_res05.rds"))

srt_integrated_all_cells <- map_refined_cells(
  srt_integrated_all_cells,
  cells_in_clusters(tcell_obj, "tcell_first_cluster_pc10_res05", c("3", "6")),
  "doublet",
  "doublet_B_T",
  "tcell_obj",
  "C3_C6"
)
srt_integrated_all_cells <- map_refined_cells(
  srt_integrated_all_cells,
  cells_in_clusters(tcell_obj, "tcell_first_cluster_pc10_res05", "8"),
  "doublet",
  "doublet_T_neutrophil",
  "tcell_obj",
  "C8"
)

tcell_validated_obj <- subset(
  tcell_obj,
  cells = cells_in_clusters(
    tcell_obj,
    "tcell_first_cluster_pc10_res05",
    c("0", "1", "2", "4", "5", "7", "9")
  )
)

# Recluster validated T cells with PC=10, k_param=30, res=0.4. This yields
# Th1 (C4), naive CD8 T cells (C0, C3, C6), exhausted CD8 T cells (C1),
# NK-like CD8 T cells (C2), CD8+ Tregs (C5), and unconventional T cells
# (C7).
tcell_validated_obj <- recluster_rna_subset(
  tcell_validated_obj,
  cluster_col = "tcell_validated_cluster_pc10_res04",
  npcs = 10,
  dims_use = 1:10,
  k_param = 30,
  resolution = 0.4,
  nfeatures = 2000
)

plot_umap(
  tcell_validated_obj,
  output_dir = annotation_out_dir,
  prefix = "tcell_validated_pc10_k30_res04",
  group = "seurat_clusters",
  label = TRUE
)

saveRDS(
  tcell_validated_obj,
  file.path(tcell_out_dir, "tcell_validated_obj_pc10_k30_res04.rds")
)

srt_integrated_all_cells <- map_refined_cells(
  srt_integrated_all_cells,
  cells_in_clusters(tcell_validated_obj, "tcell_validated_cluster_pc10_res04", c("0", "3", "6")),
  "CD8_T",
  "naive_CD8_T_cell",
  "tcell_validated_obj",
  "C0_C3_C6"
)
srt_integrated_all_cells <- map_refined_cells(
  srt_integrated_all_cells,
  cells_in_clusters(tcell_validated_obj, "tcell_validated_cluster_pc10_res04", "1"),
  "CD8_T",
  "exhausted_CD8_T_cell",
  "tcell_validated_obj",
  "C1"
)
srt_integrated_all_cells <- map_refined_cells(
  srt_integrated_all_cells,
  cells_in_clusters(tcell_validated_obj, "tcell_validated_cluster_pc10_res04", "2"),
  "CD8_T",
  "NK_like_CD8_T_cell",
  "tcell_validated_obj",
  "C2"
)
srt_integrated_all_cells <- map_refined_cells(
  srt_integrated_all_cells,
  cells_in_clusters(tcell_validated_obj, "tcell_validated_cluster_pc10_res04", "4"),
  "CD4_T",
  "Th1",
  "tcell_validated_obj",
  "C4"
)
srt_integrated_all_cells <- map_refined_cells(
  srt_integrated_all_cells,
  cells_in_clusters(tcell_validated_obj, "tcell_validated_cluster_pc10_res04", "5"),
  "CD8_T",
  "CD8_Treg",
  "tcell_validated_obj",
  "C5"
)
srt_integrated_all_cells <- map_refined_cells(
  srt_integrated_all_cells,
  cells_in_clusters(tcell_validated_obj, "tcell_validated_cluster_pc10_res04", "7"),
  "UNC_T",
  "unconventional_T_cell",
  "tcell_validated_obj",
  "C7"
)

################################################################################
# 5. Myeloid, macrophage/microglia, and DC annotation
################################################################################

# First cluster macrophages/microglia (global C5), monocytes (global C16),
# and DC (global C13) with PC=15, k_param=30, res=0.6. This yields
# macrophages/microglia (C0-2, C4), DC (C3, C6-7), monocytes (C5),
# and residual B cells (C8).
sct.myeloid <- subset(
  srt_integrated_all_cells,
  cells = cells_in_clusters(srt_integrated_all_cells, "seurat_clusters", c("5", "16", "13"))
)

sct.myeloid <- recluster_rna_subset(
  sct.myeloid,
  cluster_col = "myeloid_cluster_pc15_res06",
  npcs = 15,
  dims_use = 1:15,
  k_param = 30,
  resolution = 0.6,
  nfeatures = 2000
)

plot_umap(
  sct.myeloid,
  output_dir = annotation_out_dir,
  prefix = "myeloid_subset_pc15_k30_res06",
  group = "seurat_clusters",
  label = TRUE
)

saveRDS(sct.myeloid, file.path(myeloid_out_dir, "sct.myeloid_pc15_res06.rds"))

srt_integrated_all_cells <- map_refined_cells(
  srt_integrated_all_cells,
  cells_in_clusters(sct.myeloid, "myeloid_cluster_pc15_res06", "5"),
  "monocyte",
  "monocyte",
  "sct.myeloid",
  "C5"
)
srt_integrated_all_cells <- map_refined_cells(
  srt_integrated_all_cells,
  cells_in_clusters(sct.myeloid, "myeloid_cluster_pc15_res06", "8"),
  "B_cell",
  "residual_B_cell",
  "sct.myeloid",
  "C8"
)

sct.macrophage <- subset(
  sct.myeloid,
  cells = cells_in_clusters(sct.myeloid, "myeloid_cluster_pc15_res06", c("0", "1", "2", "4"))
)

# Recluster macrophages/microglia with PC=10, k_param=30, res=0.6. This
# yields DAM2 (C0), DAM1 (C1), macrophages (C2), and low-quality
# macrophage/microglia (C3).
sct.macrophage <- recluster_rna_subset(
  sct.macrophage,
  cluster_col = "macrophage_cluster_pc10_res06",
  npcs = 10,
  dims_use = 1:10,
  k_param = 30,
  resolution = 0.6,
  nfeatures = 2000
)

plot_umap(
  sct.macrophage,
  output_dir = annotation_out_dir,
  prefix = "macrophage_microglia_subset_pc10_k30_res06",
  group = "seurat_clusters",
  label = TRUE
)

saveRDS(sct.macrophage, file.path(myeloid_out_dir, "sct.macrophage_pc10_res06.rds"))

srt_integrated_all_cells <- map_refined_cells(
  srt_integrated_all_cells,
  cells_in_clusters(sct.macrophage, "macrophage_cluster_pc10_res06", "0"),
  "macrophage/microglia",
  "DAM2",
  "sct.macrophage",
  "C0"
)
srt_integrated_all_cells <- map_refined_cells(
  srt_integrated_all_cells,
  cells_in_clusters(sct.macrophage, "macrophage_cluster_pc10_res06", "1"),
  "macrophage/microglia",
  "DAM1",
  "sct.macrophage",
  "C1"
)
srt_integrated_all_cells <- map_refined_cells(
  srt_integrated_all_cells,
  cells_in_clusters(sct.macrophage, "macrophage_cluster_pc10_res06", "2"),
  "macrophage/microglia",
  "macrophage",
  "sct.macrophage",
  "C2"
)
srt_integrated_all_cells <- map_refined_cells(
  srt_integrated_all_cells,
  cells_in_clusters(sct.macrophage, "macrophage_cluster_pc10_res06", "3"),
  "low_quality",
  "macrophage_microglia_low_quality",
  "sct.macrophage",
  "C3"
)

sct.dc <- subset(
  sct.myeloid,
  cells = cells_in_clusters(sct.myeloid, "myeloid_cluster_pc15_res06", c("3", "6", "7"))
)

# Recluster DC with PC=15, k_param=30, res=0.8. This yields non-migratory
# DC (C1-2) and migratory DC (C0, C3).
sct.dc <- recluster_rna_subset(
  sct.dc,
  cluster_col = "dc_cluster_pc15_res08",
  npcs = 15,
  dims_use = 1:15,
  k_param = 30,
  resolution = 0.8,
  nfeatures = 2000
)

plot_umap(
  sct.dc,
  output_dir = annotation_out_dir,
  prefix = "dc_subset_pc15_k30_res08",
  group = "seurat_clusters",
  label = TRUE
)

saveRDS(sct.dc, file.path(myeloid_out_dir, "sct.dc_pc15_res08.rds"))

srt_integrated_all_cells <- map_refined_cells(
  srt_integrated_all_cells,
  cells_in_clusters(sct.dc, "dc_cluster_pc15_res08", c("0", "3")),
  "DC",
  "migratory_DC",
  "sct.dc",
  "C0_C3"
)
srt_integrated_all_cells <- map_refined_cells(
  srt_integrated_all_cells,
  cells_in_clusters(sct.dc, "dc_cluster_pc15_res08", c("1", "2")),
  "DC",
  "non_migratory_DC",
  "sct.dc",
  "C1_C2"
)

################################################################################
# 6. B-cell candidate reclustering
################################################################################

# Collect B-lineage candidates from global B-cell clusters (C0, C2-4, C10),
# residual B cells from myeloid reclustering (myeloid C8), and B/T
# doublet-like residual B cells from T-cell reclustering (T-cell C3, C6).
# Recluster with PC=15, k_param=30, res=0.2. This yields B cells
# (CD79bHIGH, C0-1, C3), doublets with T cells (C2), and doublets with
# neutrophils (C4).
global_b_cells <- cells_in_clusters(
  srt_integrated_all_cells,
  "seurat_clusters",
  c("0", "2", "3", "4", "10")
)
myeloid_residual_b_cells <- cells_in_clusters(
  sct.myeloid,
  "myeloid_cluster_pc15_res06",
  "8"
)
tcell_residual_b_cells <- cells_in_clusters(
  tcell_obj,
  "tcell_first_cluster_pc10_res05",
  c("3", "6")
)

bcell_pool_cells <- unique(c(
  global_b_cells,
  myeloid_residual_b_cells,
  tcell_residual_b_cells
))

sct.bcell <- subset(srt_integrated_all_cells, cells = bcell_pool_cells)

sct.bcell <- recluster_rna_subset(
  sct.bcell,
  cluster_col = "bcell_cluster_pc15_res02",
  npcs = 15,
  dims_use = 1:15,
  k_param = 30,
  resolution = 0.2,
  nfeatures = 2000
)

plot_umap(
  sct.bcell,
  output_dir = annotation_out_dir,
  prefix = "bcell_pc15_k30_res02",
  group = "seurat_clusters",
  label = TRUE
)

saveRDS(sct.bcell, file.path(bcell_out_dir, "sct.bcell_pc15_k30_res02.rds"))

srt_integrated_all_cells <- map_refined_cells(
  srt_integrated_all_cells,
  cells_in_clusters(sct.bcell, "bcell_cluster_pc15_res02", c("0", "1", "3")),
  "B_cell",
  "B_cell",
  "sct.bcell",
  "C0_C1_C3"
)
srt_integrated_all_cells <- map_refined_cells(
  srt_integrated_all_cells,
  cells_in_clusters(sct.bcell, "bcell_cluster_pc15_res02", "2"),
  "doublet",
  "doublet_B_T",
  "sct.bcell",
  "C2"
)
srt_integrated_all_cells <- map_refined_cells(
  srt_integrated_all_cells,
  cells_in_clusters(sct.bcell, "bcell_cluster_pc15_res02", "4"),
  "doublet",
  "doublet_B_neutrophil",
  "sct.bcell",
  "C4"
)

################################################################################
# 7. Save annotated object and summaries
################################################################################

annotation_audit_df <- dplyr::bind_rows(annotation_audit_list)
write.csv(
  annotation_audit_df,
  file.path(annotation_out_dir, "annotation_mapping_audit.csv"),
  row.names = FALSE
)

annotation_summary <- srt_integrated_all_cells@meta.data |>
  dplyr::count(celltype_major, celltype_minor, name = "n_cells")
write.csv(
  annotation_summary,
  file.path(annotation_out_dir, "annotation_celltype_summary.csv"),
  row.names = FALSE
)

saveRDS(
  srt_integrated_all_cells,
  file.path(annotation_out_dir, "srt_fullannot.rds")
)

message("Saved annotated object: ", file.path(annotation_out_dir, "srt_fullannot.rds"))
message("Saved annotation summary: ", file.path(annotation_out_dir, "annotation_celltype_summary.csv"))
