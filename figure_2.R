################################################################################
# Figure 2
################################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(ggplot2)
  library(cowplot)
  library(Matrix)
})

################################################################################
# Figure 2a: independent no-immunoglobulin UMAPs of four B-cell datasets
################################################################################

# This panel recreates the UMAP settings used in DBSCAN_analysis.R.  All four
# panels are independently embedded, so their axes describe within-dataset
# geometry and should not be compared quantitatively between datasets.

seed_use <- 1234L
set.seed(seed_use)

local_analysis_root <- Sys.getenv(
  "EAE_LOCAL_ANALYSIS_ROOT",
  unset = "/Users/shingheimok/Desktop/Phd_Study/research/Wu_lab/b_cells"
)
storage_project_root <- Sys.getenv("EAE_PROJECT_ROOT", unset = "")
if (!nzchar(storage_project_root)) {
  storage_project_root <- if (dir.exists("/Volumes/gfwu/Active/David/mng_dcln_project")) {
    "/Volumes/gfwu/Active/David/mng_dcln_project"
  } else {
    "/storage3/fs1/gfwu/Active/David/mng_dcln_project"
  }
}
output_dir <- Sys.getenv(
  "EAE_FIGURE2_OUTPUT_DIR",
  unset = file.path(local_analysis_root, "figures", "figure_2")
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

first_existing <- function(paths, description) {
  paths <- unique(paths[nzchar(paths)])
  existing <- paths[file.exists(paths)]
  if (length(existing) == 0L) {
    stop(
      "Could not find ", description, ". Checked:\n",
      paste(paths, collapse = "\n"),
      call. = FALSE
    )
  }
  existing[[1]]
}

annotated_rds <- first_existing(
  c(
    file.path(storage_project_root, "eae_scrna_vdj_preprocess", "srt_fullannot.rds"),
    file.path(local_analysis_root, "validated_bcell_obj.rds"),
    file.path(local_analysis_root, "srt_fullannot.rds"),
    file.path(local_analysis_root, "srt_fullannot_with_bcell_dbscan_eps017_min5_min31.rds")
  ),
  "annotated EAE object"
)

kolz_rds <- first_existing(
  c(
    file.path(local_analysis_root, "external", "GSE279684", "GSE279684_2024.04.18_Th1-Th17.rds"),
    file.path(storage_project_root, "kolz dataset", "GSE279684_2024.04.18_Th1-Th17.rds.gz"),
    file.path(storage_project_root, "koltz dataset", "GSE279684_2024.04.18_Th1-Th17.rds.gz")
  ),
  "Kolz Th1/Th17 object"
)

dbscan_metadata_csv <- first_existing(
  c(
    file.path(
      local_analysis_root,
      "outs", "output", "DBSCAN", "four_datasets_unbiased_dbscan_merge",
      "four_dataset_dbscan_cell_metadata.csv"
    ),
    file.path(
      storage_project_root, "eae_scrna_vdj_preprocess", "DBSCAN",
      "four_datasets_unbiased_dbscan_merge", "four_dataset_dbscan_cell_metadata.csv"
    )
  ),
  "four-dataset stable DBSCAN metadata"
)

# The raw Storage3 Kolz file can be gzip-wrapped twice.  The local .rds copy is
# read directly; this helper keeps the script deployable on Storage3.
is_gzip_file <- function(path) {
  con <- file(path, open = "rb")
  on.exit(close(con), add = TRUE)
  bytes <- readBin(con, what = "raw", n = 2L)
  length(bytes) == 2L && identical(as.integer(bytes), c(31L, 139L))
}

decompress_once <- function(src, dst) {
  input_con <- gzfile(src, open = "rb")
  output_con <- file(dst, open = "wb")
  on.exit({
    close(input_con)
    close(output_con)
  }, add = TRUE)
  repeat {
    block <- readBin(input_con, what = "raw", n = 1024L * 1024L)
    if (length(block) == 0L) break
    writeBin(block, output_con)
  }
}

read_rds_auto <- function(path) {
  if (!grepl("\\.gz$", path, ignore.case = TRUE)) return(readRDS(path))
  outer_tmp <- tempfile(fileext = ".rds-or-gz")
  on.exit(unlink(outer_tmp), add = TRUE)
  decompress_once(path, outer_tmp)
  if (is_gzip_file(outer_tmp)) return(readRDS(gzfile(outer_tmp, open = "rb")))
  readRDS(outer_tmp)
}

strict_ig_regex <- paste0(
  "^(Ighv|Ighd($|[0-9-])|Ighj($|[0-9-])|",
  "Igha$|Ighe$|Ighm$|Ighg[0-9a-z]*$|",
  "Igkv|Igkj($|[0-9-])|Igkc$|",
  "Iglv|Iglj($|[0-9-])|Iglc($|[0-9-])|",
  "Igll)"
)

build_noig_umap <- function(obj, dataset_name) {
  DefaultAssay(obj) <- "RNA"

  set.seed(seed_use)
  obj <- NormalizeData(
    obj,
    normalization.method = "LogNormalize",
    scale.factor = 10000,
    verbose = FALSE
  )
  ig_genes <- rownames(obj)[grepl(strict_ig_regex, rownames(obj))]

  set.seed(seed_use)
  obj <- FindVariableFeatures(
    obj,
    selection.method = "vst",
    nfeatures = 2000,
    verbose = FALSE
  )
  vf_no_ig <- setdiff(VariableFeatures(obj), ig_genes)
  if (length(vf_no_ig) < 50L) {
    stop("Too few non-immunoglobulin variable genes in ", dataset_name, call. = FALSE)
  }
  VariableFeatures(obj) <- vf_no_ig

  set.seed(seed_use)
  obj <- ScaleData(obj, features = vf_no_ig, verbose = FALSE)
  set.seed(seed_use)
  obj <- RunPCA(obj, features = vf_no_ig, npcs = 15, verbose = FALSE)
  set.seed(seed_use)
  obj <- RunUMAP(
    obj,
    reduction = "pca",
    dims = 1:15,
    reduction.name = "umap_noig",
    reduction.key = "UMAPNOIG_",
    n.neighbors = 30,
    min.dist = 0.3,
    metric = "cosine",
    seed.use = seed_use,
    verbose = FALSE
  )
  Embeddings(obj, "umap_noig")[, 1:2, drop = FALSE]
}

make_tissue_group <- function(obj) {
  md <- obj@meta.data
  if ("TissueGroup" %in% colnames(md)) return(obj)
  sample_col <- intersect(c("sample_id", "orig.ident", "sample"), colnames(md))[1]
  if (is.na(sample_col)) stop("EAE object lacks TissueGroup and a sample identifier.")
  obj$TissueGroup <- dplyr::case_when(
    grepl("_M", as.character(md[[sample_col]])) ~ "MNG",
    grepl("_L", as.character(md[[sample_col]])) ~ "dCLN",
    TRUE ~ NA_character_
  )
  obj
}

load_datasets <- function() {
  our_full <- make_tissue_group(readRDS(annotated_rds))
  DefaultAssay(our_full) <- "RNA"
  our_md <- our_full@meta.data
  bcell_pre_col <- if ("celltype_minor_pre_bcell_dbscan" %in% colnames(our_md)) {
    "celltype_minor_pre_bcell_dbscan"
  } else {
    "celltype_minor"
  }
  if (!all(c("celltype_major", bcell_pre_col, "TissueGroup") %in% colnames(our_md))) {
    stop("Annotated EAE object is missing the validated B-cell metadata columns.")
  }
  is_validated_bcell <-
    as.character(our_md$celltype_major) == "B_cell" &
    as.character(our_md[[bcell_pre_col]]) == "B_cell"

  kolz_full <- read_rds_auto(kolz_rds)
  DefaultAssay(kolz_full) <- "RNA"
  kolz_md <- kolz_full@meta.data
  if (!all(c("sample", "orig.ident", "compartment") %in% colnames(kolz_md))) {
    stop("Kolz object lacks sample, orig.ident, or compartment metadata.")
  }
  kolz_counts <- tryCatch(
    GetAssayData(kolz_full, assay = "RNA", layer = "counts"),
    error = function(e) GetAssayData(kolz_full, assay = "RNA", slot = "counts")
  )
  nonzero_rna <- Matrix::colSums(kolz_counts) > 0

  datasets <- list(
    our_MNG = subset(
      our_full,
      cells = rownames(our_md)[is_validated_bcell & our_md$TissueGroup == "MNG"]
    ),
    our_dCLN = subset(
      our_full,
      cells = rownames(our_md)[is_validated_bcell & our_md$TissueGroup == "dCLN"]
    ),
    kolz_Th1 = subset(
      kolz_full,
      cells = rownames(kolz_md)[kolz_md$sample == "Th1" & nonzero_rna[rownames(kolz_md)]]
    ),
    kolz_Th17 = subset(
      kolz_full,
      cells = rownames(kolz_md)[
        kolz_md$orig.ident %in% c("M1", "M2", "M3") &
          toupper(as.character(kolz_md$compartment)) == "CNS"
      ]
    )
  )
  datasets
}

dataset_titles <- c(
  our_MNG = "Meninges",
  our_dCLN = "dCLN",
  kolz_Th1 = "Kolz et al. Th1",
  kolz_Th17 = "Kolz et al. Th17"
)
dataset_order <- names(dataset_titles)

# Panel-specific colors make each independent DBSCAN partition readable.  Gray
# indicates cells outside the two selected stable clusters.
panel_palettes <- list(
  our_MNG = c(unstable = "#C9C9C9", stable_cluster_2 = "#4DBBD5", stable_cluster_1 = "#E64B35"),
  our_dCLN = c(unstable = "#C9C9C9", stable_cluster_2 = "#F39B7F", stable_cluster_1 = "#00A087"),
  kolz_Th1 = c(unstable = "#C9C9C9", stable_cluster_2 = "#4DBBD5", stable_cluster_1 = "#E64B35"),
  kolz_Th17 = c(unstable = "#C9C9C9", stable_cluster_2 = "#F39B7F", stable_cluster_1 = "#00A087")
)

dbscan_metadata <- read.csv(dbscan_metadata_csv, stringsAsFactors = FALSE)
required_dbscan_cols <- c("dataset", "cell", "stable_label")
if (!all(required_dbscan_cols %in% colnames(dbscan_metadata))) {
  stop("DBSCAN metadata is missing: ", paste(setdiff(required_dbscan_cols, colnames(dbscan_metadata)), collapse = ", "))
}

datasets <- load_datasets()

make_panel <- function(dataset, show_x_title, show_y_title) {
  coords <- build_noig_umap(datasets[[dataset]], dataset)
  umap_df <- data.frame(
    UMAP_1 = coords[, 1],
    UMAP_2 = coords[, 2],
    cell = rownames(coords),
    stringsAsFactors = FALSE
  ) |>
    left_join(
      dbscan_metadata |>
        filter(dataset == !!dataset) |>
        select(cell, stable_label),
      by = "cell"
    ) |>
    mutate(
      dbscan_display = case_when(
        grepl("stable_cluster_1$", stable_label) ~ "stable_cluster_1",
        grepl("stable_cluster_2$", stable_label) ~ "stable_cluster_2",
        TRUE ~ "unstable"
      ),
      dbscan_display = factor(
        dbscan_display,
        levels = c("unstable", "stable_cluster_2", "stable_cluster_1")
      )
    )

  missing_labels <- sum(is.na(umap_df$stable_label))
  if (missing_labels > 0L) {
    warning(dataset, ": ", missing_labels, " cells have no stable DBSCAN label and are shown in gray.")
  }

  ggplot(umap_df, aes(UMAP_1, UMAP_2, color = dbscan_display)) +
    geom_point(size = 0.48, alpha = 1, stroke = 0) +
    scale_color_manual(values = panel_palettes[[dataset]], drop = FALSE) +
    coord_fixed(expand = expansion(mult = 0.05)) +
    labs(
      title = dataset_titles[[dataset]],
      x = if (show_x_title) "UMAP 1" else NULL,
      y = if (show_y_title) "UMAP 2" else NULL
    ) +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 15, margin = margin(b = 5)),
      axis.title = element_text(size = 13),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.75),
      legend.position = "none",
      plot.margin = margin(5, 5, 5, 5)
    )
}

p_mng <- make_panel("our_MNG", show_x_title = FALSE, show_y_title = TRUE)
p_dcln <- make_panel("our_dCLN", show_x_title = FALSE, show_y_title = FALSE)
p_th1 <- make_panel("kolz_Th1", show_x_title = TRUE, show_y_title = TRUE)
p_th17 <- make_panel("kolz_Th17", show_x_title = TRUE, show_y_title = FALSE)

figure_2a <- cowplot::plot_grid(
  p_mng, p_dcln,
  p_th1, p_th17,
  ncol = 2,
  align = "hv",
  axis = "tblr",
  rel_widths = c(1, 1),
  rel_heights = c(1, 1)
)

tiff_file <- file.path(output_dir, "figure_2a_stable_dbscan_umaps.tiff")
png_file <- file.path(output_dir, "figure_2a_stable_dbscan_umaps.png")

ggsave(
  tiff_file,
  figure_2a,
  width = 9.6,
  height = 9.4,
  dpi = 300,
  compression = "lzw",
  bg = "white"
)
ggsave(
  png_file,
  figure_2a,
  width = 9.6,
  height = 9.4,
  dpi = 300,
  bg = "white"
)

cat("Saved Figure 2a:\n", tiff_file, "\n", png_file, "\n", sep = "")
