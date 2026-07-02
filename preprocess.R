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

path.home <- "/storage3/fs1/gfwu/Active/David/mng_dcln_project"
path.data <- "/storage3/fs1/gfwu/Active/MGI_SEQUENCE_DATA/Sravanthi Data"

path.data.scrna <- file.path(path.data, "Wu_MGI0279_1_10X")
path.data.vdj <- file.path(path.data, "Wu_vdj_MGI0279_1_10X_cr_v7")

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
samples$cell_prefix <- samples$sample_id

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

parse_10x_bool <- function(x) {
  tolower(as.character(x)) %in% c("true", "t", "1")
}

add_percent_ribo <- function(obj) {
  obj[["percent.ribo"]] <- PercentageFeatureSet(
    obj,
    pattern = "^(Rpl|Rps|RPL|RPS|Mrpl|Mrps|MRPL|MRPS)"
  )
  obj
}

extract_10x_barcode <- function(x) {
  out <- sub(
    ".*?([ACGTN]+-[0-9]+)(?:_[0-9]+)?$",
    "\\1",
    x,
    perl = TRUE
  )

  unchanged <- identical(length(out), length(x)) && all(out == x)
  if (unchanged) {
    out <- sub("(_[0-9]+)$", "", x, perl = TRUE)
  }

  out
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
    c(sample_paths$filtered_h5, sample_paths$raw_h5),
    paste0("scRNA input for ", sample_id)
  )

  raw_counts <- read_cellranger_h5_matrix(sample_paths$raw_h5)
  filtered_counts <- read_cellranger_h5_matrix(sample_paths$filtered_h5)

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

summarize_vdj_metadata <- function(path.data.vdj,
                                   sample_id,
                                   vdj_id,
                                   cell_prefix = sample_id,
                                   receptor = c("BCR", "TCR")) {
  receptor <- match.arg(receptor)
  prefix <- tolower(receptor)

  if (receptor == "BCR") {
    vdj_paths <- get_vdj_paths(path.data.vdj, sample_id, vdj_id, "vdj_b")
    chain_pattern <- "^IG"
  } else {
    vdj_paths <- get_vdj_paths(path.data.vdj, sample_id, vdj_id, "vdj_t")
    chain_pattern <- "^TR"
  }

  contig_file <- NULL
  if (!is.na(vdj_paths$filtered_contigs) && file.exists(vdj_paths$filtered_contigs)) {
    contig_file <- vdj_paths$filtered_contigs
  } else if (!is.na(vdj_paths$all_contigs) && file.exists(vdj_paths$all_contigs)) {
    contig_file <- vdj_paths$all_contigs
  }

  if (is.null(contig_file)) {
    warning("Missing ", receptor, " contig annotation CSV for ", sample_id)
    return(data.frame())
  }

  contig_df <- read.csv(contig_file, stringsAsFactors = FALSE)
  if (nrow(contig_df) == 0) return(data.frame())

  required_cols <- c(
    "barcode", "is_cell", "high_confidence", "productive", "chain",
    "raw_clonotype_id", "raw_consensus_id", "cdr3_nt", "cdr3",
    "v_gene", "d_gene", "j_gene", "c_gene"
  )
  missing_cols <- setdiff(required_cols, colnames(contig_df))
  if (length(missing_cols) > 0) {
    warning(
      "Missing expected ", receptor, " column(s) for ", sample_id, ": ",
      paste(missing_cols, collapse = ", ")
    )
    for (col in missing_cols) contig_df[[col]] <- NA
  }

  contig_df <- contig_df |>
    dplyr::mutate(
      barcode_raw = barcode,
      barcode_prefixed = paste(cell_prefix, barcode, sep = "_"),
      is_cell = parse_10x_bool(is_cell),
      high_confidence = parse_10x_bool(high_confidence),
      productive = parse_10x_bool(productive),
      raw_clonotype_id = dplyr::na_if(raw_clonotype_id, ""),
      raw_consensus_id = dplyr::na_if(raw_consensus_id, "")
    ) |>
    dplyr::filter(
      is_cell,
      high_confidence,
      productive,
      grepl(chain_pattern, chain),
      !is.na(raw_clonotype_id)
    ) |>
    dplyr::mutate(
      sequence_signature = ifelse(
        !is.na(raw_consensus_id),
        raw_consensus_id,
        paste(chain, cdr3_nt, sep = ":")
      )
    )

  if (nrow(contig_df) == 0) return(data.frame())

  barcode_meta <- contig_df |>
    dplyr::arrange(barcode_raw, chain) |>
    dplyr::group_by(barcode_raw) |>
    dplyr::summarise(
      sample_id = sample_id,
      barcode_prefixed = dplyr::first(barcode_prefixed),
      chain_count = dplyr::n(),
      chains = paste(sort(unique(chain)), collapse = ";"),
      clonotype_id = paste(sort(unique(raw_clonotype_id)), collapse = "|"),
      multiple_clonotypes = dplyr::n_distinct(raw_clonotype_id) > 1,
      subclonotype_id = paste(sort(unique(sequence_signature)), collapse = "|"),
      v_genes = paste(sort(unique(v_gene[v_gene != "None"])), collapse = ";"),
      d_genes = paste(sort(unique(d_gene[d_gene != "None"])), collapse = ";"),
      j_genes = paste(sort(unique(j_gene[j_gene != "None"])), collapse = ";"),
      isotypes = paste(sort(unique(c_gene[c_gene != "None"])), collapse = ";"),
      cdr3s_aa = paste(sort(unique(paste(chain, cdr3, sep = ":"))), collapse = ";"),
      cdr3s_nt = paste(sort(unique(paste(chain, cdr3_nt, sep = ":"))), collapse = ";"),
      .groups = "drop"
    )

  receptor_fields <- c(
    "chain_count", "chains", "clonotype_id", "multiple_clonotypes",
    "subclonotype_id", "v_genes", "d_genes", "j_genes", "isotypes",
    "cdr3s_aa", "cdr3s_nt"
  )
  rename_idx <- match(receptor_fields, names(barcode_meta))
  names(barcode_meta)[rename_idx] <- paste0(prefix, "_", receptor_fields)

  if (!is.na(vdj_paths$clonotypes) && file.exists(vdj_paths$clonotypes)) {
    clonotype_df <- read.csv(vdj_paths$clonotypes, stringsAsFactors = FALSE)
    expected_clonotype_cols <- c("clonotype_id", "frequency", "proportion", "cdr3s_aa", "cdr3s_nt")
    if (all(expected_clonotype_cols %in% colnames(clonotype_df))) {
      clonotype_df <- clonotype_df |>
        dplyr::rename(
          !!paste0(prefix, "_clonotype_id") := clonotype_id,
          !!paste0(prefix, "_clonotype_frequency_cellranger") := frequency,
          !!paste0(prefix, "_clonotype_proportion_cellranger") := proportion,
          !!paste0(prefix, "_clonotype_cdr3s_aa_cellranger") := cdr3s_aa,
          !!paste0(prefix, "_clonotype_cdr3s_nt_cellranger") := cdr3s_nt
        )

      barcode_meta <- barcode_meta |>
        dplyr::left_join(clonotype_df, by = paste0(prefix, "_clonotype_id"))
    }
  }

  barcode_meta
}

summarize_bcr_metadata <- function(path.data.vdj, sample_id, vdj_id, cell_prefix = sample_id) {
  summarize_vdj_metadata(
    path.data.vdj = path.data.vdj,
    sample_id = sample_id,
    vdj_id = vdj_id,
    cell_prefix = cell_prefix,
    receptor = "BCR"
  )
}

summarize_tcr_metadata <- function(path.data.vdj, sample_id, vdj_id, cell_prefix = sample_id) {
  summarize_vdj_metadata(
    path.data.vdj = path.data.vdj,
    sample_id = sample_id,
    vdj_id = vdj_id,
    cell_prefix = cell_prefix,
    receptor = "TCR"
  )
}

ensure_vdj_metadata_columns <- function(meta, prefix) {
  expected <- c(
    "barcode_raw", "sample_id", "barcode_prefixed",
    paste0(prefix, "_chain_count"),
    paste0(prefix, "_chains"),
    paste0(prefix, "_clonotype_id"),
    paste0(prefix, "_multiple_clonotypes"),
    paste0(prefix, "_subclonotype_id"),
    paste0(prefix, "_v_genes"),
    paste0(prefix, "_d_genes"),
    paste0(prefix, "_j_genes"),
    paste0(prefix, "_isotypes"),
    paste0(prefix, "_cdr3s_aa"),
    paste0(prefix, "_cdr3s_nt"),
    paste0(prefix, "_clonotype_frequency_cellranger"),
    paste0(prefix, "_clonotype_proportion_cellranger"),
    paste0(prefix, "_clonotype_cdr3s_aa_cellranger"),
    paste0(prefix, "_clonotype_cdr3s_nt_cellranger")
  )

  if (nrow(meta) == 0 && ncol(meta) == 0) {
    meta <- data.frame(barcode_raw = character(), sample_id = character())
  }

  missing <- setdiff(expected, colnames(meta))
  for (col in missing) meta[[col]] <- NA
  meta[, expected, drop = FALSE]
}

add_bcr_tcr_metadata_to_object <- function(obj,
                                           sample_metadata,
                                           path.data.vdj,
                                           remove_missing_bcr = TRUE,
                                           remove_missing_tcr = FALSE) {
  bcr_meta_list <- lapply(seq_len(nrow(sample_metadata)), function(i) {
    summarize_bcr_metadata(
      path.data.vdj = path.data.vdj,
      sample_id = sample_metadata$sample_id[i],
      vdj_id = sample_metadata$vdj_id[i],
      cell_prefix = sample_metadata$cell_prefix[i]
    )
  })

  tcr_meta_list <- lapply(seq_len(nrow(sample_metadata)), function(i) {
    summarize_tcr_metadata(
      path.data.vdj = path.data.vdj,
      sample_id = sample_metadata$sample_id[i],
      vdj_id = sample_metadata$vdj_id[i],
      cell_prefix = sample_metadata$cell_prefix[i]
    )
  })

  bcr_meta <- dplyr::bind_rows(bcr_meta_list)
  tcr_meta <- dplyr::bind_rows(tcr_meta_list)
  bcr_meta <- ensure_vdj_metadata_columns(bcr_meta, "bcr")
  tcr_meta <- ensure_vdj_metadata_columns(tcr_meta, "tcr")

  meta_df <- obj@meta.data
  meta_df$barcode <- rownames(meta_df)
  meta_df$barcode_raw <- extract_10x_barcode(meta_df$barcode)

  meta_df <- meta_df |>
    dplyr::left_join(bcr_meta, by = c("barcode_raw", "sample_id")) |>
    dplyr::left_join(tcr_meta, by = c("barcode_raw", "sample_id")) |>
    dplyr::mutate(
      has_bcr = !is.na(bcr_clonotype_id),
      has_tcr = !is.na(tcr_clonotype_id)
    )

  bcr_clone_size_df <- meta_df |>
    dplyr::filter(has_bcr) |>
    dplyr::count(sample_id, bcr_clonotype_id, name = "bcr_clonotype_size")

  bcr_subclone_size_df <- meta_df |>
    dplyr::filter(has_bcr) |>
    dplyr::count(
      sample_id,
      bcr_clonotype_id,
      bcr_subclonotype_id,
      name = "bcr_subclonotype_size"
    )

  tcr_clone_size_df <- meta_df |>
    dplyr::filter(has_tcr) |>
    dplyr::count(sample_id, tcr_clonotype_id, name = "tcr_clonotype_size")

  tcr_subclone_size_df <- meta_df |>
    dplyr::filter(has_tcr) |>
    dplyr::count(
      sample_id,
      tcr_clonotype_id,
      tcr_subclonotype_id,
      name = "tcr_subclonotype_size"
    )

  meta_df <- meta_df |>
    dplyr::left_join(
      bcr_clone_size_df,
      by = c("sample_id", "bcr_clonotype_id")
    ) |>
    dplyr::left_join(
      bcr_subclone_size_df,
      by = c("sample_id", "bcr_clonotype_id", "bcr_subclonotype_id")
    ) |>
    dplyr::left_join(
      tcr_clone_size_df,
      by = c("sample_id", "tcr_clonotype_id")
    ) |>
    dplyr::left_join(
      tcr_subclone_size_df,
      by = c("sample_id", "tcr_clonotype_id", "tcr_subclonotype_id")
    ) |>
    dplyr::mutate(
      bcr_expansion_status = dplyr::case_when(
        !has_bcr ~ "missing_bcr",
        bcr_clonotype_size >= 2 ~ "expanded",
        TRUE ~ "unexpanded"
      ),
      bcr_subclonotype_uid = dplyr::if_else(
        has_bcr,
        paste(sample_id, bcr_clonotype_id, bcr_subclonotype_id, sep = "::"),
        NA_character_
      ),
      bcr_clonotype_uid = dplyr::if_else(
        has_bcr,
        paste(sample_id, bcr_clonotype_id, sep = "::"),
        NA_character_
      ),
      tcr_expansion_status = dplyr::case_when(
        !has_tcr ~ "missing_tcr",
        tcr_clonotype_size >= 2 ~ "expanded",
        TRUE ~ "unexpanded"
      ),
      tcr_subclonotype_uid = dplyr::if_else(
        has_tcr,
        paste(sample_id, tcr_clonotype_id, tcr_subclonotype_id, sep = "::"),
        NA_character_
      ),
      tcr_clonotype_uid = dplyr::if_else(
        has_tcr,
        paste(sample_id, tcr_clonotype_id, sep = "::"),
        NA_character_
      )
    )

  rownames(meta_df) <- meta_df$barcode
  meta_df$barcode_raw <- NULL
  obj@meta.data <- meta_df[colnames(obj), , drop = FALSE]

  obj_with_vdj_metadata <- obj

  bcr_summary <- obj_with_vdj_metadata@meta.data |>
    dplyr::count(sample_id, has_bcr, bcr_expansion_status, name = "n_cells")

  tcr_summary <- obj_with_vdj_metadata@meta.data |>
    dplyr::count(sample_id, has_tcr, tcr_expansion_status, name = "n_cells")

  obj_bcr_only <- obj
  obj_tcr_only <- obj

  if (remove_missing_bcr) {
    obj_bcr_only <- subset(obj_bcr_only, subset = has_bcr)
  }

  if (remove_missing_tcr) {
    obj_tcr_only <- subset(obj_tcr_only, subset = has_tcr)
  }

  list(
    obj_with_vdj_metadata = obj_with_vdj_metadata,
    obj_bcr_only = obj_bcr_only,
    obj_tcr_only = obj_tcr_only,
    barcode_bcr_metadata = bcr_meta,
    barcode_tcr_metadata = tcr_meta,
    bcr_summary = bcr_summary,
    tcr_summary = tcr_summary
  )
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
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank()
    ) +
    ggtitle(NULL)

  save_tiff_plot(p, out_file, width = 6, height = 5)
}

run_soupx <- function(raw_mat,
                      filtered_mat,
                      n_pcs = 30,
                      dims_use = 1:30,
                      resolution = 0.2,
                      out = output_dir,
                      plot_file = "soupx_autoEstCont.tiff") {
  temp_obj <- CreateSeuratObject(counts = filtered_mat, project = "SoupXTemp")
  soupx_obj <- SoupChannel(tod = raw_mat, toc = filtered_mat)

  DefaultAssay(temp_obj) <- "RNA"
  temp_obj <- SCTransform(temp_obj, verbose = FALSE)
  set.seed(1234)
  temp_obj <- RunPCA(temp_obj, npcs = n_pcs, verbose = FALSE)
  set.seed(1234)
  temp_obj <- RunUMAP(temp_obj, dims = dims_use, verbose = FALSE)
  set.seed(1234)
  temp_obj <- FindNeighbors(temp_obj, dims = dims_use, verbose = FALSE)
  set.seed(1234)
  temp_obj <- FindClusters(temp_obj, resolution = resolution, verbose = TRUE)

  meta <- temp_obj@meta.data
  umap <- temp_obj@reductions$umap@cell.embeddings

  soupx_obj <- setClusters(soupx_obj, setNames(as.character(meta$seurat_clusters), rownames(meta)))
  soupx_obj <- setDR(soupx_obj, umap)

  ensure_dir(out)
  tiff(file.path(out, plot_file), width = 7, height = 7, units = "in", res = 300)
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
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-|^mt-")
  obj <- add_percent_ribo(obj)
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
  set.seed(1234)
  sce_tmp <- scDblFinder(
    sce_tmp,
    samples = "sample_id",
    clusters = TRUE,
    dbr = 0.08
  )
  new_obj <- as.Seurat(sce_tmp)
  subset(new_obj, subset = scDblFinder.class == "singlet" & scDblFinder.score < 0.15)
}

preprocess_one_sample <- function(sample_row) {
  sample_short <- sample_row$sample_short
  sample_id <- sample_row$sample_id

  message("Reading scRNA data for ", sample_short, " (", sample_id, ")")
  scrna <- read_sample_scrna(path.data.scrna, sample_id)

  message("Running SoupX for ", sample_short)
  soupx_results <- run_soupx(
    raw_mat = scrna$raw_counts,
    filtered_mat = scrna$filtered_counts,
    n_pcs = 30,
    dims_use = 1:30,
    resolution = 0.2,
    out = file.path(output_dir, "QC", "soupx"),
    plot_file = paste0(sample_id, "_soupx_autoEstCont.tiff")
  )

  obj <- CreateSeuratObject(
    counts = soupx_results$adjusted_counts,
    project = sample_id
  )
  obj <- RenameCells(obj, add.cell.id = sample_id)

  obj$sample_short <- sample_short
  obj$sample_id <- sample_id
  obj$vdj_id <- sample_row$vdj_id
  obj$tissue <- sample_row$tissue
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-|^mt-")
  obj <- add_percent_ribo(obj)
  obj@misc$scrna_paths <- scrna$sample_paths

  plot_qc_violin(obj, output_dir, paste0(sample_id, "_preQC"))

  message("Filtering cells and doublets for ", sample_short)
  obj <- apply_qc_cutoffs(obj)
  obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
  obj <- run_doublet_filter(obj)

  DefaultAssay(obj) <- "RNA"
  obj <- SCTransform(
    obj,
    vars.to.regress = c("nCount_RNA", "percent.mt", "percent.ribo"),
    method = "glmGamPoi",
    verbose = FALSE
  )

  plot_qc_violin(obj, output_dir, paste0(sample_id, "_postQC"))

  sample_out <- file.path(output_dir, "individual_samples", paste0(sample_id, "_preprocessed.rds"))
  ensure_dir(dirname(sample_out))
  saveRDS(obj, sample_out)

  obj
}

integrate_list <- function(seurat_list) {
  set.seed(1234)
  features <- SelectIntegrationFeatures(object.list = seurat_list, nfeatures = 3000)
  seurat_list <- PrepSCTIntegration(object.list = seurat_list, anchor.features = features)
  set.seed(1234)
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
                               k_param = 50,
                               resolution = 0.2,
                               ndims = 50) {
  DefaultAssay(obj) <- "integrated"
  set.seed(1234)
  obj <- RunPCA(obj, npcs = npcs, verbose = FALSE)
  plot_elbow(obj, output_dir = output_dir, prefix = prefix, ndims = ndims)
  set.seed(1234)
  obj <- FindNeighbors(obj, dims = dims_use, k.param = k_param, verbose = FALSE)
  set.seed(1234)
  obj <- FindClusters(obj, resolution = resolution, verbose = FALSE)
  set.seed(1234)
  obj <- RunUMAP(obj, dims = dims_use, verbose = FALSE)
  plot_umap(obj, output_dir = output_dir, prefix = prefix, group = "seurat_clusters", label = TRUE)
  obj
}

rerun_umap_with_granularity <- function(obj,
                                        dims_use = 1:15,
                                        k_param = 50,
                                        resolution = 0.5) {
  DefaultAssay(obj) <- "integrated"
  set.seed(1234)
  obj <- FindNeighbors(obj, dims = dims_use, k.param = k_param, verbose = FALSE)
  set.seed(1234)
  obj <- FindClusters(obj, resolution = resolution, verbose = FALSE)
  set.seed(1234)
  obj <- RunUMAP(obj, dims = dims_use, verbose = FALSE)
  obj
}

export_umap_metadata <- function(obj, output_dir, prefix = "MGI0279_integrated") {
  if (!"umap" %in% names(obj@reductions)) return(NULL)

  out_file <- file.path(output_dir, "tables", paste0(prefix, "_umap_metadata.csv"))
  ensure_dir(dirname(out_file))

  umap_df <- as.data.frame(Embeddings(obj, reduction = "umap"))
  umap_df$barcode <- rownames(umap_df)

  meta_df <- obj@meta.data
  meta_df$barcode <- rownames(meta_df)

  out_df <- dplyr::left_join(umap_df, meta_df, by = "barcode")
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

seurat_list <- setNames(vector("list", nrow(samples)), samples$sample_id)

for (i in seq_len(nrow(samples))) {
  seurat_list[[samples$sample_id[[i]]]] <- preprocess_one_sample(samples[i, , drop = FALSE])
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
  k_param = 50,
  resolution = 0.2,
  ndims = 50
)

srt_integrated_clustered_rds <- file.path(output_dir, "srt_integrated_clustered.rds")
saveRDS(srt_integrated, srt_integrated_clustered_rds)

message("Re-running UMAP and clustering with reference granularity")
srt_integrated <- rerun_umap_with_granularity(
  srt_integrated,
  dims_use = 1:15,
  k_param = 20,
  resolution = 0.6
)

plot_umap(
  srt_integrated,
  output_dir = output_dir,
  prefix = "MGI0279_integrated_15_20_04",
  group = "seurat_clusters",
  label = TRUE
)

srt_integrated_pruned_rds <- file.path(output_dir, "srt_integrated_clustered_pruned.rds")
saveRDS(srt_integrated, srt_integrated_pruned_rds)

message("Adding BCR/TCR metadata")
vdj_results <- add_bcr_tcr_metadata_to_object(
  srt_integrated,
  sample_metadata = samples,
  path.data.vdj = path.data.vdj,
  remove_missing_bcr = TRUE,
  remove_missing_tcr = FALSE
)

srt_integrated_all_cells <- vdj_results$obj_with_vdj_metadata
srt_integrated_bcr_only <- vdj_results$obj_bcr_only
srt_integrated_tcr_only <- vdj_results$obj_tcr_only
srt_integrated <- srt_integrated_all_cells

vdj_table_dir <- ensure_dir(file.path(output_dir, "tables"))
vdj_metadata_csv <- file.path(vdj_table_dir, "MGI0279_integrated_vdj_metadata_all_cells.csv")
bcr_summary_csv <- file.path(vdj_table_dir, "MGI0279_integrated_bcr_summary.csv")
tcr_summary_csv <- file.path(vdj_table_dir, "MGI0279_integrated_tcr_summary.csv")

all_cell_vdj_metadata <- srt_integrated_all_cells@meta.data
all_cell_vdj_metadata$barcode <- rownames(all_cell_vdj_metadata)
write.csv(all_cell_vdj_metadata, vdj_metadata_csv, row.names = FALSE)
write.csv(vdj_results$bcr_summary, bcr_summary_csv, row.names = FALSE)
write.csv(vdj_results$tcr_summary, tcr_summary_csv, row.names = FALSE)

vdj_rds <- file.path(output_dir, "srt_integrated_clustered_pruned_with_vdj_metadata.rds")
saveRDS(srt_integrated_all_cells, vdj_rds)

umap_csv_file <- export_umap_metadata(
  srt_integrated,
  output_dir = output_dir,
  prefix = "MGI0279_L1_L2_M1_M2_integrated"
)

integrated_rds <- file.path(output_dir, "MGI0279_L1_L2_M1_M2_integrated_seurat.rds")
saveRDS(srt_integrated, integrated_rds)

message("Saved integrated object: ", integrated_rds)
message("Saved VDJ-annotated object: ", vdj_rds)
if (!is.null(umap_csv_file)) message("Saved UMAP metadata: ", umap_csv_file)
