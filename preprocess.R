################################################################################
# MGI0279 scRNA + V(D)J preprocessing
################################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(SoupX)
  library(scDblFinder)
  library(SingleCellExperiment)
  library(ggplot2)
  library(dplyr)
})

################################################################################
# 1. Paths and sample metadata
################################################################################

path.home <- "/storage3/fs1/gfwu/Active/david/mng_dcln_project"
path.data <- "/storage3/fs1/gfwu/Active/MGI_SEQUENCE_DATA/Sravanthi Data"

path.data.scrna <- file.path(path.data, "Wu_MGI0279_1_10X")
path.data.vdj <- file.path(path.home, "Wu_vdj_MGI0279_1_10X_cr_v7")
path.data.vdj.fallback <- file.path(path.data, "Wu_vdj_MGI0279_1_10X_cr_v7")

if (!dir.exists(path.data.vdj) && dir.exists(path.data.vdj.fallback)) {
  message("Using fallback VDJ path: ", path.data.vdj.fallback)
  path.data.vdj <- path.data.vdj.fallback
}

output_dir <- file.path(path.home, "eae_scrna_vdj_preprocess")

samples <- data.frame(
  sample_short = c("L1", "L2", "M1", "M2"),
  sample_id = c(
    "MGI0279_1_Wu2020_L1-lib1",
    "MGI0279_1_Wu2020_L2-lib1",
    "MGI0279_1_Wu2020_M1-lib1",
    "MGI0279_1_Wu2020_M2-lib1"
  ),
  vdj_id = c(
    "MGI0279_WUAB-Wu2020-Wu2020_L1",
    "MGI0279_WUAB-Wu2020-Wu2020_L2",
    "MGI0279_WUAB-Wu2020-Wu2020_M1",
    "MGI0279_WUAB-Wu2020-Wu2020_M2"
  ),
  tissue = c("lymph_node", "lymph_node", "meninges", "meninges"),
  stringsAsFactors = FALSE
)

################################################################################
# 2. General helpers
################################################################################

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
  path
}

require_files <- function(paths, label) {
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0) {
    stop(
      "Missing required ", label, " file(s):\n",
      paste(missing, collapse = "\n"),
      call. = FALSE
    )
  }
  invisible(paths)
}

save_tiff_plot <- function(plot_obj,
                           out_file,
                           width = 10,
                           height = 6,
                           units = "in",
                           res = 300) {
  ensure_dir(dirname(out_file))
  tiff(out_file, width = width, height = height, units = units, res = res)
  print(plot_obj)
  dev.off()
  out_file
}

read_cellranger_h5_matrix <- function(h5_file) {
  mat <- Read10X_h5(h5_file)

  if (is.list(mat)) {
    if ("Gene Expression" %in% names(mat)) {
      mat <- mat[["Gene Expression"]]
    } else {
      mat <- mat[[1]]
    }
  }

  mat
}

restrict_to_called_barcodes <- function(mat, barcode_file = NULL) {
  if (is.null(barcode_file) || !file.exists(barcode_file)) return(mat)

  barcode_df <- read.csv(barcode_file, header = TRUE, stringsAsFactors = FALSE)
  barcode_vec <- unique(as.character(barcode_df[[1]]))
  barcode_keep <- intersect(colnames(mat), barcode_vec)

  if (length(barcode_keep) == 0) {
    warning("No overlapping barcodes found in ", barcode_file, "; returning original matrix.")
    return(mat)
  }

  mat[, barcode_keep, drop = FALSE]
}

get_scrna_paths <- function(path.data.scrna, sample_id) {
  outs_dir <- file.path(path.data.scrna, sample_id, "outs")

  list(
    outs_dir = outs_dir,
    filtered_h5 = file.path(outs_dir, "filtered_feature_bc_matrix.h5"),
    raw_h5 = file.path(outs_dir, "raw_feature_bc_matrix.h5"),
    barcode_csv = file.path(outs_dir, "final_cell_barcodes.csv")
  )
}

read_sample_scrna <- function(path.data.scrna, sample_id) {
  sample_paths <- get_scrna_paths(path.data.scrna, sample_id)
  require_files(
    c(sample_paths$filtered_h5, sample_paths$raw_h5, sample_paths$barcode_csv),
    paste0("scRNA input for ", sample_id)
  )

  raw_counts <- read_cellranger_h5_matrix(sample_paths$raw_h5)
  filtered_counts <- read_cellranger_h5_matrix(sample_paths$filtered_h5)
  filtered_counts <- restrict_to_called_barcodes(
    filtered_counts,
    barcode_file = sample_paths$barcode_csv
  )

  list(
    raw_counts = raw_counts,
    filtered_counts = filtered_counts,
    sample_paths = sample_paths
  )
}

################################################################################
# 3. V(D)J helpers
################################################################################

resolve_existing_file <- function(dir, candidates) {
  hits <- file.path(dir, candidates)
  hits <- hits[file.exists(hits)]
  if (length(hits) == 0) return(NA_character_)
  hits[[1]]
}

resolve_clonotype_file <- function(dir) {
  direct <- resolve_existing_file(
    dir,
    c("clonotypes.csv", "clonotype.csv", "clonotype", "clonotypes")
  )
  if (!is.na(direct) && file.info(direct)$isdir) {
    nested <- list.files(direct, pattern = "\\.csv$", full.names = TRUE)
    if (length(nested) > 0) return(nested[[1]])
  }
  direct
}

get_vdj_dir_candidates <- function(path.data.vdj, sample_id, vdj_id, chain_dir) {
  c(
    file.path(path.data.vdj, vdj_id, "outs", "per_sample_outs", vdj_id, chain_dir),
    file.path(path.data.vdj, sample_id, "outs", "per_sample_outs", sample_id, chain_dir),
    file.path(path.data.vdj, vdj_id, "outs", chain_dir),
    file.path(path.data.vdj, sample_id, "outs", chain_dir)
  )
}

get_vdj_paths <- function(path.data.vdj, sample_id, vdj_id, chain_dir) {
  dir_candidates <- get_vdj_dir_candidates(path.data.vdj, sample_id, vdj_id, chain_dir)

  path_table <- lapply(dir_candidates, function(vdj_dir) {
    list(
      vdj_dir = vdj_dir,
      filtered_contigs = resolve_existing_file(
        vdj_dir,
        c("filtered_contig_annotations.csv", "filtered_contig_annotation.csv")
      ),
      all_contigs = resolve_existing_file(
        vdj_dir,
        c("all_contig_annotations.csv", "all_contig_annotation.csv")
      ),
      clonotypes = resolve_clonotype_file(vdj_dir)
    )
  })

  complete <- vapply(
    path_table,
    function(x) all(!is.na(c(x$filtered_contigs, x$all_contigs, x$clonotypes))),
    logical(1)
  )

  if (any(complete)) {
    return(path_table[[which(complete)[[1]]]])
  }

  existing_dir <- vapply(path_table, function(x) dir.exists(x$vdj_dir), logical(1))
  if (any(existing_dir)) {
    return(path_table[[which(existing_dir)[[1]]]])
  }

  warning(
    "No VDJ directory found for ", sample_id, " / ", chain_dir, ". Tried:\n",
    paste(dir_candidates, collapse = "\n")
  )
  path_table[[1]]
}

read_optional_csv <- function(path) {
  if (is.na(path) || !file.exists(path) || file.info(path)$isdir) return(NULL)
  read.csv(path, stringsAsFactors = FALSE)
}

warn_missing_vdj_files <- function(paths, sample_id, chain_label) {
  named_paths <- c(
    filtered_contigs = paths$filtered_contigs,
    all_contigs = paths$all_contigs,
    clonotypes = paths$clonotypes
  )
  missing <- names(named_paths)[is.na(named_paths) | !file.exists(named_paths)]

  if (length(missing) > 0) {
    warning(
      "Missing ", chain_label, " VDJ file(s) for ", sample_id, ": ",
      paste(missing, collapse = ", "),
      ". Directory used: ", paths$vdj_dir
    )
  }
}

summarize_vdj_by_barcode <- function(filtered_contigs,
                                     sample_short,
                                     chain_label) {
  if (is.null(filtered_contigs) || nrow(filtered_contigs) == 0) {
    return(data.frame(cell = character(), stringsAsFactors = FALSE))
  }

  if (!"barcode" %in% colnames(filtered_contigs)) {
    warning("VDJ table has no barcode column for ", sample_short, " ", chain_label)
    return(data.frame(cell = character(), stringsAsFactors = FALSE))
  }

  clonotype_col <- intersect(
    c("raw_clonotype_id", "clonotype_id", "clone_id"),
    colnames(filtered_contigs)
  )
  clonotype_col <- if (length(clonotype_col) > 0) clonotype_col[[1]] else NA_character_

  cdr3_col <- intersect(
    c("cdr3", "cdr3_aa"),
    colnames(filtered_contigs)
  )
  cdr3_col <- if (length(cdr3_col) > 0) cdr3_col[[1]] else NA_character_

  contigs <- filtered_contigs
  contigs$cell <- paste(sample_short, contigs$barcode, sep = "_")

  split_contigs <- split(contigs, contigs$cell)
  out <- lapply(split_contigs, function(x) {
    data.frame(
      cell = x$cell[[1]],
      contig_count = nrow(x),
      clonotype_id = if (is.na(clonotype_col)) NA_character_ else paste(unique(x[[clonotype_col]]), collapse = ";"),
      cdr3 = if (is.na(cdr3_col)) NA_character_ else paste(unique(x[[cdr3_col]]), collapse = ";"),
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, out)
  colnames(out) <- c(
    "cell",
    paste0(chain_label, "_contig_count"),
    paste0(chain_label, "_clonotype_id"),
    paste0(chain_label, "_cdr3")
  )
  rownames(out) <- NULL
  out
}

read_sample_vdj <- function(path.data.vdj, sample_id, vdj_id, sample_short) {
  bcr_paths <- get_vdj_paths(path.data.vdj, sample_id, vdj_id, "vdj_b")
  tcr_paths <- get_vdj_paths(path.data.vdj, sample_id, vdj_id, "vdj_t")
  warn_missing_vdj_files(bcr_paths, sample_id, "BCR")
  warn_missing_vdj_files(tcr_paths, sample_id, "TCR")

  bcr <- list(
    paths = bcr_paths,
    filtered_contigs = read_optional_csv(bcr_paths$filtered_contigs),
    all_contigs = read_optional_csv(bcr_paths$all_contigs),
    clonotypes = read_optional_csv(bcr_paths$clonotypes)
  )
  tcr <- list(
    paths = tcr_paths,
    filtered_contigs = read_optional_csv(tcr_paths$filtered_contigs),
    all_contigs = read_optional_csv(tcr_paths$all_contigs),
    clonotypes = read_optional_csv(tcr_paths$clonotypes)
  )

  bcr_meta <- summarize_vdj_by_barcode(bcr$filtered_contigs, sample_short, "bcr")
  tcr_meta <- summarize_vdj_by_barcode(tcr$filtered_contigs, sample_short, "tcr")

  meta <- merge(bcr_meta, tcr_meta, by = "cell", all = TRUE)
  if (nrow(meta) > 0) {
    required_cols <- c(
      "bcr_contig_count", "bcr_clonotype_id", "bcr_cdr3",
      "tcr_contig_count", "tcr_clonotype_id", "tcr_cdr3"
    )
    missing_cols <- setdiff(required_cols, colnames(meta))
    for (col in missing_cols) meta[[col]] <- NA

    meta$has_bcr <- !is.na(meta$bcr_contig_count) & meta$bcr_contig_count > 0
    meta$has_tcr <- !is.na(meta$tcr_contig_count) & meta$tcr_contig_count > 0
    rownames(meta) <- meta$cell
    meta$cell <- NULL
  }

  list(
    bcr = bcr,
    tcr = tcr,
    metadata = meta
  )
}

add_vdj_metadata <- function(obj, vdj) {
  if (!is.null(vdj$metadata) && nrow(vdj$metadata) > 0) {
    aligned <- vdj$metadata[colnames(obj), , drop = FALSE]
    obj <- AddMetaData(obj, aligned)
  }

  if (!"has_bcr" %in% colnames(obj@meta.data)) obj$has_bcr <- FALSE
  if (!"has_tcr" %in% colnames(obj@meta.data)) obj$has_tcr <- FALSE

  obj$has_bcr[is.na(obj$has_bcr)] <- FALSE
  obj$has_tcr[is.na(obj$has_tcr)] <- FALSE
  obj@misc$vdj <- vdj

  obj
}

################################################################################
# 4. Plotting and preprocessing helpers
################################################################################

plot_qc_violin <- function(obj, output_dir, prefix) {
  out_file <- file.path(output_dir, "QC", "violin", paste0(prefix, "_violin.tiff"))

  p <- VlnPlot(
    obj,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    ncol = 3,
    pt.size = 0.05
  ) &
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  save_tiff_plot(p, out_file, width = 8, height = 5)
}

plot_elbow <- function(obj, output_dir, prefix = "integrated", ndims = 50) {
  out_file <- file.path(output_dir, "QC", "elbow", paste0(prefix, "_elbow.tiff"))
  p <- ElbowPlot(obj, ndims = ndims)
  save_tiff_plot(p, out_file, width = 6, height = 5)
}

plot_umap <- function(obj,
                      output_dir,
                      prefix = "integrated",
                      group = "seurat_clusters",
                      label = TRUE) {
  out_file <- file.path(output_dir, "QC", "umap", paste0(prefix, "_umap.tiff"))

  p <- DimPlot(
    obj,
    reduction = "umap",
    group.by = group,
    label = label
  ) +
    labs(x = NULL, y = NULL) +
    ggtitle(NULL)

  save_tiff_plot(p, out_file, width = 6, height = 5)
}

run_soupx <- function(raw_mat,
                      filtered_mat,
                      n_pcs = 30,
                      dims_use = 1:30,
                      resolution = 0.2,
                      out = output_dir,
                      plot_file = "soupx_autoEstCont.pdf") {
  temp_obj <- CreateSeuratObject(counts = filtered_mat, project = "SoupXTemp")
  soupx_obj <- SoupChannel(tod = raw_mat, toc = filtered_mat)

  DefaultAssay(temp_obj) <- "RNA"
  temp_obj <- SCTransform(temp_obj, verbose = FALSE)
  temp_obj <- RunPCA(temp_obj, npcs = n_pcs, verbose = FALSE)
  temp_obj <- RunUMAP(temp_obj, dims = dims_use, verbose = FALSE)
  temp_obj <- FindNeighbors(temp_obj, dims = dims_use, verbose = FALSE)
  temp_obj <- FindClusters(temp_obj, resolution = resolution, verbose = TRUE)

  meta <- temp_obj@meta.data
  umap <- temp_obj@reductions$umap@cell.embeddings

  soupx_obj <- setClusters(soupx_obj, setNames(meta$seurat_clusters, rownames(meta)))
  soupx_obj <- setDR(soupx_obj, umap)

  ensure_dir(out)
  pdf(file.path(out, plot_file), width = 7, height = 7)
  soupx_obj <- autoEstCont(soupx_obj)
  dev.off()

  adjusted_counts <- adjustCounts(soupx_obj, roundToInt = TRUE)

  list(
    soupx_channel = soupx_obj,
    soupx_seurat_seed = temp_obj,
    adjusted_counts = adjusted_counts
  )
}

apply_qc_cutoffs <- function(obj,
                             min_features = 200,
                             max_features = 2500,
                             max_percent_mt = 5) {
  subset(
    obj,
    subset =
      nFeature_RNA > min_features &
      nFeature_RNA < max_features &
      percent.mt < max_percent_mt
  )
}

run_doublet_filter <- function(obj) {
  sce_tmp <- as.SingleCellExperiment(obj)
  sce_tmp <- scDblFinder(sce_tmp, dbr = NULL)
  new_obj <- as.Seurat(sce_tmp)
  subset(new_obj, subset = scDblFinder.class == "singlet")
}

preprocess_one_sample <- function(sample_row) {
  sample_short <- sample_row$sample_short
  sample_id <- sample_row$sample_id
  vdj_id <- sample_row$vdj_id

  message("Reading scRNA data for ", sample_short, " (", sample_id, ")")
  scrna <- read_sample_scrna(path.data.scrna, sample_id)

  message("Running SoupX for ", sample_short)
  soupx_results <- run_soupx(
    raw_mat = scrna$raw_counts,
    filtered_mat = scrna$filtered_counts,
    n_pcs = 30,
    dims_use = 1:30,
    resolution = 0.2,
    out = output_dir,
    plot_file = paste0(sample_short, "_soupx_autoEstCont.pdf")
  )

  obj <- CreateSeuratObject(
    counts = soupx_results$adjusted_counts,
    project = sample_short
  )
  obj <- RenameCells(obj, add.cell.id = sample_short)

  obj$sample_short <- sample_short
  obj$sample_id <- sample_id
  obj$vdj_id <- vdj_id
  obj$tissue <- sample_row$tissue
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-|^mt-")

  message("Reading BCR/TCR data for ", sample_short)
  vdj <- read_sample_vdj(path.data.vdj, sample_id, vdj_id, sample_short)
  obj <- add_vdj_metadata(obj, vdj)
  obj@misc$scrna_paths <- scrna$sample_paths

  plot_qc_violin(obj, output_dir, paste0(sample_short, "_preQC"))

  message("Filtering cells and doublets for ", sample_short)
  obj <- apply_qc_cutoffs(obj)
  obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
  obj <- run_doublet_filter(obj)

  DefaultAssay(obj) <- "RNA"
  obj <- SCTransform(
    obj,
    vars.to.regress = "nCount_RNA",
    method = "glmGamPoi",
    verbose = FALSE
  )

  plot_qc_violin(obj, output_dir, paste0(sample_short, "_postQC"))

  sample_out <- file.path(output_dir, "individual_samples", paste0(sample_short, "_preprocessed.rds"))
  ensure_dir(dirname(sample_out))
  saveRDS(obj, sample_out)

  obj
}

integrate_list <- function(seurat_list) {
  features <- SelectIntegrationFeatures(object.list = seurat_list, nfeatures = 3000)
  seurat_list <- PrepSCTIntegration(object.list = seurat_list, anchor.features = features)
  anchors <- FindIntegrationAnchors(
    object.list = seurat_list,
    normalization.method = "SCT",
    anchor.features = features
  )
  IntegrateData(anchorset = anchors, normalization.method = "SCT")
}

process_integrated <- function(obj,
                               output_dir,
                               prefix = "MGI0279_integrated",
                               npcs = 50,
                               dims_use = 1:30,
                               k_param = 20,
                               resolution = 0.2,
                               ndims = 50) {
  DefaultAssay(obj) <- "integrated"
  obj <- RunPCA(obj, npcs = npcs, verbose = FALSE)
  plot_elbow(obj, output_dir = output_dir, prefix = prefix, ndims = ndims)
  obj <- FindNeighbors(obj, dims = dims_use, k.param = k_param, verbose = FALSE)
  obj <- FindClusters(obj, resolution = resolution, verbose = FALSE)
  obj <- RunUMAP(obj, dims = dims_use, verbose = FALSE)
  plot_umap(obj, output_dir = output_dir, prefix = prefix, group = "seurat_clusters", label = TRUE)
  obj
}

export_umap_metadata <- function(obj, output_dir, prefix = "MGI0279_integrated") {
  if (!"umap" %in% names(obj@reductions)) return(NULL)

  out_file <- file.path(output_dir, "tables", paste0(prefix, "_umap_metadata.csv"))
  ensure_dir(dirname(out_file))

  umap_df <- as.data.frame(Embeddings(obj, reduction = "umap"))
  umap_df$cell <- rownames(umap_df)

  meta_df <- obj@meta.data
  meta_df$cell <- rownames(meta_df)

  out_df <- dplyr::left_join(umap_df, meta_df, by = "cell")
  write.csv(out_df, out_file, row.names = FALSE)
  out_file
}

################################################################################
# 5. Run preprocessing, integration, and save outputs
################################################################################

ensure_dir(output_dir)

message("Output directory: ", output_dir)
message("scRNA input directory: ", path.data.scrna)
message("VDJ input directory: ", path.data.vdj)

seurat_list <- setNames(vector("list", nrow(samples)), samples$sample_short)

for (i in seq_len(nrow(samples))) {
  seurat_list[[samples$sample_short[[i]]]] <- preprocess_one_sample(samples[i, , drop = FALSE])
}

message("Integrating samples")
srt_integrated <- integrate_list(seurat_list)
srt_integrated@misc$sample_table <- samples
srt_integrated@misc$path.home <- path.home
srt_integrated@misc$path.data <- path.data
srt_integrated@misc$path.data.scrna <- path.data.scrna
srt_integrated@misc$path.data.vdj <- path.data.vdj

message("Running PCA, clustering, and UMAP")
srt_integrated <- process_integrated(
  srt_integrated,
  output_dir = output_dir,
  prefix = "MGI0279_L1_L2_M1_M2_integrated",
  npcs = 50,
  dims_use = 1:30,
  k_param = 20,
  resolution = 0.2,
  ndims = 50
)

umap_csv_file <- export_umap_metadata(
  srt_integrated,
  output_dir = output_dir,
  prefix = "MGI0279_L1_L2_M1_M2_integrated"
)

integrated_rds <- file.path(output_dir, "MGI0279_L1_L2_M1_M2_integrated_seurat.rds")
saveRDS(srt_integrated, integrated_rds)

message("Saved integrated object: ", integrated_rds)
if (!is.null(umap_csv_file)) message("Saved UMAP metadata: ", umap_csv_file)
