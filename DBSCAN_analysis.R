################################################################################
# Four-dataset unbiased DBSCAN stability pipeline
#
# Purpose
# -------
# Starting from the product of annotation.R and the raw Kolz Th1/Th17 object:
#
#   1. Select B cells from our MNG and dCLN data.
#   2. Select Th1-CNS and Th17-CNS cells from the Kolz object.
#   3. Build an independent no-immunoglobulin UMAP for each dataset.
#   4. Scan DBSCAN across minPts = 5, 10, 15 and 30 adaptive eps values per
#      minPts. Raw DBSCAN labels are retained; small clusters are not collapsed
#      during the scan.
#   5. Match clusters across neighboring eps values by cell-set Jaccard
#      similarity and define stable trajectories.
#   6. Retain stable trajectories with at least 5% of the dataset, merge
#      equivalent trajectories across minPts, and choose a deterministic
#      representative state for stable cluster 1 and stable cluster 2.
#   7. Label all four datasets, write labeled objects and metadata, and compare
#      selected clusters across datasets by top-50 marker overlap with
#      one-sided Fisher exact tests.
#
# This supersedes the earlier fixed-eps, MNG_DBSCAN_2-benchmark pipeline. The
# MNG_DBSCAN_13/MNG_DBSCAN_2 labels are not used to select the final clusters.
################################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(dbscan)
  library(Matrix)
})

################################################################################
# 1. Reproducibility, paths, and hyperparameters
################################################################################

seed_use <- 1234L
set.seed(seed_use)

project_root <- Sys.getenv(
  "EAE_PROJECT_ROOT",
  unset = "/storage3/fs1/gfwu/Active/David/mng_dcln_project"
)
preprocess_dir <- file.path(project_root, "eae_scrna_vdj_preprocess")

annotated_rds <- file.path(preprocess_dir, "srt_fullannot.rds")
kolz_candidates <- c(
  file.path(project_root, "kolz dataset", "GSE279684_2024.04.18_Th1-Th17.rds.gz"),
  file.path(project_root, "koltz dataset", "GSE279684_2024.04.18_Th1-Th17.rds.gz"),
  file.path(preprocess_dir, "kolz dataset", "GSE279684_2024.04.18_Th1-Th17.rds.gz"),
  file.path(preprocess_dir, "koltz dataset", "GSE279684_2024.04.18_Th1-Th17.rds.gz")
)
kolz_rds <- kolz_candidates[file.exists(kolz_candidates)][1]

if (!file.exists(annotated_rds)) {
  stop(
    "Annotated object not found: ", annotated_rds,
    "\nRun annotation.R first, or set EAE_PROJECT_ROOT.",
    call. = FALSE
  )
}
if (length(kolz_rds) == 0L || is.na(kolz_rds)) {
  stop(
    "Kolz Th1/Th17 object not found. Checked:\n",
    paste(kolz_candidates, collapse = "\n"),
    call. = FALSE
  )
}

output_dir <- file.path(
  preprocess_dir,
  "DBSCAN",
  "four_datasets_unbiased_dbscan_merge"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# DBSCAN scan parameters.
minPts_grid <- c(5L, 10L, 15L)
eps_quantiles <- seq(0.01, 0.99, length.out = 30L)

# Clusters larger than this are retained for stability testing. This is not a
# post-hoc collapse rule: every raw DBSCAN state remains in the assignments.
min_cluster_size_for_testing <- 30L

# A trajectory must contain at least three neighboring eps states to be called
# stable. The final report applies the stronger 5% size rule below.
jaccard_threshold <- 0.80
minimum_stable_links <- 2L
minimum_stable_states_for_final_selection <- 5L
minimum_stable_fraction <- 0.05

# Cross-minPts consensus matching is based on representative cell sets. The
# cross-dataset merge rule is deliberately more stringent and uses marker
# overlap plus Fisher FDR.
consensus_jaccard_threshold <- 0.50
merge_min_overlap <- 20L
merge_min_jaccard <- 0.25
merge_fdr_threshold <- 0.05
signature_size <- 50L

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

is_gzip_file <- function(path) {
  con <- file(path, open = "rb")
  on.exit(close(con), add = TRUE)
  bytes <- readBin(con, what = "raw", n = 2L)
  length(bytes) == 2L && identical(as.integer(bytes), c(31L, 139L))
}

decompress_once <- function(src, dst) {
  in_con <- gzfile(src, open = "rb")
  out_con <- file(dst, open = "wb")
  on.exit({
    close(in_con)
    close(out_con)
  }, add = TRUE)

  repeat {
    block <- readBin(in_con, what = "raw", n = 1024L * 1024L)
    if (length(block) == 0L) break
    writeBin(block, out_con)
  }
}

read_rds_auto <- function(path) {
  if (!grepl("\\.gz$", path, ignore.case = TRUE)) {
    return(readRDS(path))
  }

  # The Storage3 Kolz file is gzip-wrapped twice. Decompress one layer to a
  # temporary file, detect whether the result is another gzip stream, and then
  # let readRDS consume the inner stream without loading it into a raw vector.
  outer_tmp <- tempfile(fileext = ".rds-or-gz")
  on.exit(unlink(outer_tmp), add = TRUE)
  decompress_once(path, outer_tmp)

  if (is_gzip_file(outer_tmp)) {
    message("Detected nested gzip compression; reading the inner RDS stream.")
    return(readRDS(gzfile(outer_tmp, open = "rb")))
  }
  readRDS(outer_tmp)
}

state_key <- function(dataset, minPts, eps) {
  paste(
    dataset,
    as.integer(minPts),
    formatC(eps, digits = 12, format = "fg", flag = "#"),
    sep = "::"
  )
}

node_key <- function(dataset, minPts, eps, raw_cluster) {
  paste(
    state_key(dataset, minPts, eps),
    raw_cluster,
    sep = "::"
  )
}

jaccard <- function(x, y) {
  x <- unique(as.character(x))
  y <- unique(as.character(y))
  if (length(x) == 0L && length(y) == 0L) return(1)
  u <- union(x, y)
  if (length(u) == 0L) return(0)
  length(intersect(x, y)) / length(u)
}

fisher_overlap <- function(candidate_genes, reference_genes, universe_genes) {
  candidate_genes <- unique(intersect(candidate_genes, universe_genes))
  reference_genes <- unique(intersect(reference_genes, universe_genes))
  shared <- sort(intersect(candidate_genes, reference_genes))

  a <- length(shared)
  b <- length(setdiff(candidate_genes, reference_genes))
  c <- length(setdiff(reference_genes, candidate_genes))
  d <- length(universe_genes) - a - b - c
  if (d < 0L) stop("Invalid Fisher table: negative fourth cell.", call. = FALSE)

  ft <- fisher.test(
    matrix(c(a, b, c, d), nrow = 2L, byrow = TRUE),
    alternative = "greater"
  )

  data.frame(
    overlap_n = a,
    candidate_n = length(candidate_genes),
    reference_n = length(reference_genes),
    universe_n = length(universe_genes),
    overlap_jaccard = ifelse(
      length(union(candidate_genes, reference_genes)) == 0L,
      0,
      a / length(union(candidate_genes, reference_genes))
    ),
    overlap_genes = paste(shared, collapse = "/"),
    fisher_odds_ratio = unname(ft$estimate),
    fisher_p_value = ft$p.value,
    stringsAsFactors = FALSE
  )
}

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
    stop("Too few non-immunoglobulin variable genes in ", dataset_name, ".", call. = FALSE)
  }
  VariableFeatures(obj) <- vf_no_ig

  set.seed(seed_use)
  obj <- ScaleData(obj, features = vf_no_ig, verbose = FALSE)

  set.seed(seed_use)
  obj <- RunPCA(
    obj,
    features = vf_no_ig,
    npcs = 15,
    verbose = FALSE
  )

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

  coords <- Embeddings(obj, "umap_noig")[, 1:2, drop = FALSE]
  rownames(coords) <- colnames(obj)
  list(obj = obj, coords = coords)
}

find_root <- function(parent, x) {
  while (!identical(parent[[x]], x)) {
    parent[[x]] <- parent[[parent[[x]]]]
    x <- parent[[x]]
  }
  x
}

top_marker_signature <- function(obj, labels, cluster_name, signature_size = 50L) {
  labels <- labels[colnames(obj)]
  if (anyNA(labels)) stop("Missing DBSCAN labels in marker calculation.", call. = FALSE)
  if (!cluster_name %in% labels) {
    return(data.frame())
  }

  marker_obj <- obj
  Idents(marker_obj) <- factor(labels, levels = unique(c(cluster_name, setdiff(labels, cluster_name))))

  deg <- tryCatch(
    FindMarkers(
      marker_obj,
      ident.1 = cluster_name,
      assay = "RNA",
      test.use = "t",
      only.pos = TRUE,
      logfc.threshold = 0,
      min.pct = 0,
      verbose = FALSE
    ),
    error = function(e) {
      warning(
        "FindMarkers failed for ", cluster_name, ": ", conditionMessage(e),
        call. = FALSE
      )
      data.frame()
    }
  )
  if (nrow(deg) == 0L) return(data.frame())

  deg$gene <- rownames(deg)
  fc_col <- intersect(c("avg_log2FC", "avg_logFC"), colnames(deg))[1]
  if (is.na(fc_col)) stop("FindMarkers returned no log-fold-change column.", call. = FALSE)

  deg |>
    dplyr::filter(.data[[fc_col]] >= 1, p_val < 0.001) |>
    dplyr::arrange(dplyr::desc(.data[[fc_col]])) |>
    dplyr::slice_head(n = signature_size) |>
    dplyr::transmute(
      gene = gene,
      avg_log2FC = .data[[fc_col]],
      p_val = p_val
    )
}

################################################################################
# 3. Read annotation.R output and the raw Kolz object
################################################################################

cat("Loading annotated object:\n", annotated_rds, "\n", sep = "")
our_full <- readRDS(annotated_rds)
DefaultAssay(our_full) <- "RNA"
our_md <- our_full@meta.data

if (!"TissueGroup" %in% colnames(our_md)) {
  sample_col <- intersect(c("sample_id", "orig.ident", "sample"), colnames(our_md))[1]
  if (is.na(sample_col)) {
    stop("Our object has no TissueGroup or sample identifier.", call. = FALSE)
  }
  our_full$TissueGroup <- dplyr::case_when(
    grepl("_M", as.character(our_md[[sample_col]])) ~ "MNG",
    grepl("_L", as.character(our_md[[sample_col]])) ~ "dCLN",
    TRUE ~ NA_character_
  )
  our_md <- our_full@meta.data
}

bcell_pre_col <- if ("celltype_minor_pre_bcell_dbscan" %in% colnames(our_md)) {
  "celltype_minor_pre_bcell_dbscan"
} else {
  "celltype_minor"
}

is_validated_bcell <- as.character(our_md$celltype_major) == "B_cell" &
  as.character(our_md[[bcell_pre_col]]) == "B_cell"

our_mng_cells <- rownames(our_md)[
  is_validated_bcell & as.character(our_md$TissueGroup) == "MNG"
]
our_dcln_cells <- rownames(our_md)[
  is_validated_bcell & as.character(our_md$TissueGroup) == "dCLN"
]

cat("Loading Kolz object:\n", kolz_rds, "\n", sep = "")
kolz_full <- read_rds_auto(kolz_rds)
DefaultAssay(kolz_full) <- "RNA"
kolz_md <- kolz_full@meta.data

required_kolz_metadata <- c("sample", "orig.ident", "compartment")
missing_kolz_metadata <- setdiff(required_kolz_metadata, colnames(kolz_md))
if (length(missing_kolz_metadata) > 0L) {
  stop(
    "Kolz object lacks metadata: ",
    paste(missing_kolz_metadata, collapse = ", "),
    call. = FALSE
  )
}

# These are the exact CNS definitions used for the final table:
#   Th1  = sample == "Th1" and compartment == "CNS"
#   Th17 = orig.ident in M1/M2/M3 and compartment == "CNS"
kolz_th1_cells <- rownames(kolz_md)[
  as.character(kolz_md$sample) == "Th1" &
    toupper(as.character(kolz_md$compartment)) == "CNS"
]
kolz_th17_cells <- rownames(kolz_md)[
  as.character(kolz_md$orig.ident) %in% c("M1", "M2", "M3") &
    toupper(as.character(kolz_md$compartment)) == "CNS"
]

dataset_cells <- list(
  our_MNG = unique(our_mng_cells),
  our_dCLN = unique(our_dcln_cells),
  kolz_Th1 = unique(kolz_th1_cells),
  kolz_Th17 = unique(kolz_th17_cells)
)

expected_counts <- c(
  our_MNG = 4408L,
  our_dCLN = 8534L,
  kolz_Th1 = 4164L,
  kolz_Th17 = 6039L
)
observed_counts <- vapply(dataset_cells, length, integer(1))
cat("\nSelected cells by dataset:\n")
print(observed_counts)
if (!identical(as.integer(observed_counts), as.integer(expected_counts[names(observed_counts)]))) {
  stop(
    "Selected cell counts do not match the final pipeline inputs.\nExpected:\n",
    paste(names(expected_counts), expected_counts, sep = " = ", collapse = "\n"),
    "\nObserved:\n",
    paste(names(observed_counts), observed_counts, sep = " = ", collapse = "\n"),
    "\nCheck the annotation.R object, Kolz raw object, and CNS metadata filters.",
    call. = FALSE
  )
}

dataset_objects <- list(
  our_MNG = subset(our_full, cells = dataset_cells$our_MNG),
  our_dCLN = subset(our_full, cells = dataset_cells$our_dCLN),
  kolz_Th1 = subset(kolz_full, cells = dataset_cells$kolz_Th1),
  kolz_Th17 = subset(kolz_full, cells = dataset_cells$kolz_Th17)
)

dataset_order <- c("our_MNG", "our_dCLN", "kolz_Th1", "kolz_Th17")
dataset_objects <- dataset_objects[dataset_order]
dataset_cells <- dataset_cells[dataset_order]

if (anyDuplicated(unlist(dataset_cells))) {
  stop("Cell barcodes are duplicated across the four datasets.", call. = FALSE)
}

################################################################################
# 4. Build one no-Ig UMAP and one adaptive eps grid per dataset
################################################################################

embedding_cache <- list()
eps_grid_by_dataset <- list()
eps_grid_rows <- list()
cluster_state_rows <- list()
state_scan_rows <- list()
assignments_cache <- list()

for (dataset in dataset_order) {
  cat("\nBuilding embedding for ", dataset, " (", length(dataset_cells[[dataset]]), " cells)\n", sep = "")
  embedding_cache[[dataset]] <- build_noig_umap(dataset_objects[[dataset]], dataset)
  coords <- embedding_cache[[dataset]]$coords

  eps_grid_by_dataset[[dataset]] <- lapply(minPts_grid, function(minPts_use) {
    kdist <- dbscan::kNNdist(coords, k = minPts_use)
    sort(unique(round(
      quantile(
        kdist,
        probs = eps_quantiles,
        na.rm = TRUE,
        names = FALSE
      ),
      digits = 6
    )))
  })
  names(eps_grid_by_dataset[[dataset]]) <- as.character(minPts_grid)

  eps_grid_rows[[dataset]] <- dplyr::bind_rows(lapply(minPts_grid, function(minPts_use) {
    eps_values <- eps_grid_by_dataset[[dataset]][[as.character(minPts_use)]]
    data.frame(
      dataset = dataset,
      minPts = minPts_use,
      eps_rank = seq_along(eps_values),
      eps = eps_values,
      stringsAsFactors = FALSE
    )
  }))

  coords_cells <- rownames(coords)
  for (minPts_use in minPts_grid) {
    eps_values <- eps_grid_by_dataset[[dataset]][[as.character(minPts_use)]]
    for (eps_use in eps_values) {
      set.seed(seed_use)
      db <- dbscan::dbscan(
        coords,
        eps = eps_use,
        minPts = minPts_use
      )

      raw_labels <- paste0("DBSCAN_", db$cluster)
      names(raw_labels) <- coords_cells
      assignments_cache[[state_key(dataset, minPts_use, eps_use)]] <- raw_labels

      state_tab <- sort(table(raw_labels), decreasing = TRUE)
      non_noise <- setdiff(names(state_tab), "DBSCAN_0")
      eligible <- non_noise[
        as.integer(state_tab[non_noise]) > min_cluster_size_for_testing
      ]

      noise_n <- if ("DBSCAN_0" %in% names(state_tab)) {
        as.integer(state_tab[["DBSCAN_0"]])
      } else {
        0L
      }

      state_scan_rows[[length(state_scan_rows) + 1L]] <- data.frame(
        dataset = dataset,
        minPts = minPts_use,
        eps = eps_use,
        n_cells = length(raw_labels),
        n_raw_clusters = length(non_noise),
        n_eligible_clusters = length(eligible),
        n_noise = noise_n,
        noise_fraction = noise_n / length(raw_labels),
        stringsAsFactors = FALSE
      )

      if (length(eligible) == 0L) next
      for (raw_cluster in eligible) {
        n_cells <- as.integer(state_tab[[raw_cluster]])
        cluster_state_rows[[length(cluster_state_rows) + 1L]] <- data.frame(
          dataset = dataset,
          minPts = minPts_use,
          eps = eps_use,
          state_id = state_key(dataset, minPts_use, eps_use),
          node = node_key(dataset, minPts_use, eps_use, raw_cluster),
          raw_cluster = raw_cluster,
          n_cells = n_cells,
          fraction = n_cells / length(raw_labels),
          stringsAsFactors = FALSE
        )
      }
    }
  }
}

eps_grid_df <- dplyr::bind_rows(eps_grid_rows)
state_scan_df <- dplyr::bind_rows(state_scan_rows)
cluster_state_df <- dplyr::bind_rows(cluster_state_rows)

if (nrow(cluster_state_df) == 0L) {
  stop("No non-noise DBSCAN clusters exceeded the testing threshold.", call. = FALSE)
}

write.csv(
  eps_grid_df,
  file.path(output_dir, "all_dataset_eps_grid.csv"),
  row.names = FALSE
)
write.csv(
  state_scan_df,
  file.path(output_dir, "all_dataset_state_scan.csv"),
  row.names = FALSE
)
write.csv(
  dplyr::bind_rows(lapply(names(assignments_cache), function(key) {
    labels <- assignments_cache[[key]]
    key_parts <- strsplit(key, "::", fixed = TRUE)[[1]]
    data.frame(
      dataset = key_parts[1],
      minPts = as.integer(key_parts[2]),
      eps = as.numeric(key_parts[3]),
      cell = names(labels),
      dbscan_raw = unname(labels),
      stringsAsFactors = FALSE
    )
  })),
  file.path(output_dir, "all_dataset_dbscan_assignments.csv"),
  row.names = FALSE
)

threshold_df <- data.frame(
  dataset = dataset_order,
  n_cells = unname(observed_counts[dataset_order]),
  minimum_cluster_size = ceiling(unname(observed_counts[dataset_order]) * minimum_stable_fraction),
  minimum_cluster_fraction = minimum_stable_fraction,
  stringsAsFactors = FALSE
)
write.csv(
  threshold_df,
  file.path(output_dir, "dataset_cluster_size_thresholds.csv"),
  row.names = FALSE
)

################################################################################
# 5. Match adjacent eps states and build stable trajectories
################################################################################

stability_rows <- list()

for (dataset in dataset_order) {
  for (minPts_use in minPts_grid) {
    eps_values <- eps_grid_by_dataset[[dataset]][[as.character(minPts_use)]]
    if (length(eps_values) < 2L) next

    for (i in seq_len(length(eps_values) - 1L)) {
      eps_left <- eps_values[i]
      eps_right <- eps_values[i + 1L]
      labels_left <- assignments_cache[[state_key(dataset, minPts_use, eps_left)]]
      labels_right <- assignments_cache[[state_key(dataset, minPts_use, eps_right)]]

      sizes_left <- table(labels_left)
      sizes_right <- table(labels_right)
      clusters_left <- setdiff(
        names(sizes_left)[as.integer(sizes_left) > min_cluster_size_for_testing],
        "DBSCAN_0"
      )
      clusters_right <- setdiff(
        names(sizes_right)[as.integer(sizes_right) > min_cluster_size_for_testing],
        "DBSCAN_0"
      )
      if (length(clusters_left) == 0L || length(clusters_right) == 0L) next

      for (cluster_left in clusters_left) {
        scores <- vapply(clusters_right, function(cluster_right) {
          jaccard(
            names(labels_left)[labels_left == cluster_left],
            names(labels_right)[labels_right == cluster_right]
          )
        }, numeric(1))
        best_index <- which.max(scores)
        cluster_right <- clusters_right[best_index]

        stability_rows[[length(stability_rows) + 1L]] <- data.frame(
          dataset = dataset,
          minPts = minPts_use,
          eps_left = eps_left,
          cluster_left = cluster_left,
          eps_right = eps_right,
          cluster_right = cluster_right,
          best_jaccard = scores[best_index],
          stable_link = scores[best_index] >= jaccard_threshold,
          stringsAsFactors = FALSE
        )
      }
    }
  }
}

stability_df <- dplyr::bind_rows(stability_rows)
write.csv(
  stability_df,
  file.path(output_dir, "adjacent_eps_stability_all_four_datasets.csv"),
  row.names = FALSE
)

trajectory_map_rows <- list()
trajectory_states_rows <- list()
trajectory_summary_rows <- list()

for (dataset in dataset_order) {
  for (minPts_use in minPts_grid) {
    dataset_nodes <- cluster_state_df |>
      dplyr::filter(dataset == !!dataset, minPts == !!minPts_use) |>
      dplyr::pull(node) |>
      unique()
    if (length(dataset_nodes) == 0L) next

    parent <- setNames(dataset_nodes, dataset_nodes)
    stable_links <- stability_df |>
      dplyr::filter(
        dataset == !!dataset,
        minPts == !!minPts_use,
        stable_link
      )

    if (nrow(stable_links) > 0L) {
      for (link_i in seq_len(nrow(stable_links))) {
        left_node <- node_key(
          dataset,
          minPts_use,
          stable_links$eps_left[link_i],
          stable_links$cluster_left[link_i]
        )
        right_node <- node_key(
          dataset,
          minPts_use,
          stable_links$eps_right[link_i],
          stable_links$cluster_right[link_i]
        )
        if (left_node %in% dataset_nodes && right_node %in% dataset_nodes) {
          left_root <- find_root(parent, left_node)
          right_root <- find_root(parent, right_node)
          if (!identical(left_root, right_root)) {
            parent[[right_root]] <- left_root
          }
        }
      }
    }

    map_df <- data.frame(
      dataset = dataset,
      minPts = minPts_use,
      node = dataset_nodes,
      trajectory_id = vapply(dataset_nodes, function(x) find_root(parent, x), character(1)),
      stringsAsFactors = FALSE
    )
    trajectory_map_rows[[length(trajectory_map_rows) + 1L]] <- map_df

    states_df <- cluster_state_df |>
      dplyr::filter(dataset == !!dataset, minPts == !!minPts_use) |>
      dplyr::left_join(map_df, by = c("dataset", "minPts", "node"))
    trajectory_states_rows[[length(trajectory_states_rows) + 1L]] <- states_df

    summary_df <- states_df |>
      dplyr::group_by(dataset, minPts, trajectory_id) |>
      dplyr::summarise(
        eps_min = min(eps),
        eps_max = max(eps),
        n_eps_states = dplyr::n(),
        representative_eps = median(eps),
        mean_cluster_fraction = mean(fraction),
        stable_trajectory = n_eps_states >= (minimum_stable_links + 1L),
        .groups = "drop"
      )
    trajectory_summary_rows[[length(trajectory_summary_rows) + 1L]] <- summary_df
  }
}

trajectory_map_df <- dplyr::bind_rows(trajectory_map_rows)
trajectory_states_df <- dplyr::bind_rows(trajectory_states_rows)
trajectory_summary_df <- dplyr::bind_rows(trajectory_summary_rows)

trajectory_representatives_df <- trajectory_states_df |>
  dplyr::left_join(
    trajectory_summary_df,
    by = c("dataset", "minPts", "trajectory_id")
  ) |>
  dplyr::filter(stable_trajectory) |>
  dplyr::mutate(distance_to_representative = abs(eps - representative_eps)) |>
  dplyr::group_by(dataset, minPts, trajectory_id) |>
  dplyr::arrange(distance_to_representative, dplyr::desc(n_cells), eps, .by_group = TRUE) |>
  dplyr::slice_head(n = 1L) |>
  dplyr::ungroup()

write.csv(
  trajectory_map_df,
  file.path(output_dir, "eps_cluster_trajectory_map_all_four_datasets.csv"),
  row.names = FALSE
)
write.csv(
  trajectory_states_df,
  file.path(output_dir, "stable_cluster_trajectory_states_all_four_datasets.csv"),
  row.names = FALSE
)
write.csv(
  trajectory_summary_df |>
    dplyr::filter(stable_trajectory),
  file.path(output_dir, "stable_cluster_trajectories_all_four_datasets.csv"),
  row.names = FALSE
)

################################################################################
# 6. Merge stable trajectories across minPts and select final representatives
################################################################################

# Use only identified stable trajectories with a representative cluster >=5%.
# This excludes small stable trajectories from the final stable1/stable2 labels
# while preserving them in the complete trajectory output.
candidate_trajectories <- trajectory_representatives_df |>
  dplyr::filter(
    stable_trajectory,
    n_eps_states >= minimum_stable_links + 1L,
    fraction >= minimum_stable_fraction
  ) |>
  dplyr::mutate(
    consensus_node = paste(dataset, trajectory_id, sep = "||")
  )

if (nrow(candidate_trajectories) == 0L) {
  stop("No stable trajectory reached the 5% final-selection threshold.", call. = FALSE)
}

consensus_edge_rows <- list()
consensus_parent <- setNames(
  candidate_trajectories$consensus_node,
  candidate_trajectories$consensus_node
)

for (dataset in dataset_order) {
  candidates_dataset <- candidate_trajectories |>
    dplyr::filter(dataset == !!dataset)
  if (nrow(candidates_dataset) < 2L) next

  for (i in seq_len(nrow(candidates_dataset) - 1L)) {
    for (j in (i + 1L):nrow(candidates_dataset)) {
      if (candidates_dataset$minPts[i] == candidates_dataset$minPts[j]) next

      cells_i <- names(
        assignments_cache[[
          state_key(
            dataset,
            candidates_dataset$minPts[i],
            candidates_dataset$eps[i]
          )
        ]]
      )
      labels_i <- assignments_cache[[
        state_key(
          dataset,
          candidates_dataset$minPts[i],
          candidates_dataset$eps[i]
        )
      ]]
      labels_j <- assignments_cache[[
        state_key(
          dataset,
          candidates_dataset$minPts[j],
          candidates_dataset$eps[j]
        )
      ]]
      cells_i <- cells_i[labels_i[cells_i] == candidates_dataset$raw_cluster[i]]
      cells_j <- names(labels_j)[labels_j == candidates_dataset$raw_cluster[j]]
      score <- jaccard(cells_i, cells_j)

      consensus_edge_rows[[length(consensus_edge_rows) + 1L]] <- data.frame(
        dataset = dataset,
        trajectory_1 = candidates_dataset$trajectory_id[i],
        trajectory_2 = candidates_dataset$trajectory_id[j],
        minPts_1 = candidates_dataset$minPts[i],
        minPts_2 = candidates_dataset$minPts[j],
        jaccard = score,
        merged = score >= consensus_jaccard_threshold,
        stringsAsFactors = FALSE
      )

      if (score >= consensus_jaccard_threshold) {
        node_i <- candidates_dataset$consensus_node[i]
        node_j <- candidates_dataset$consensus_node[j]
        root_i <- find_root(consensus_parent, node_i)
        root_j <- find_root(consensus_parent, node_j)
        if (!identical(root_i, root_j)) consensus_parent[[root_j]] <- root_i
      }
    }
  }
}

consensus_edges_df <- dplyr::bind_rows(consensus_edge_rows)
write.csv(
  consensus_edges_df,
  file.path(output_dir, "consensus_trajectory_edges_all_four_datasets.csv"),
  row.names = FALSE
)

candidate_trajectories$consensus_component <- vapply(
  candidate_trajectories$consensus_node,
  function(x) find_root(consensus_parent, x),
  character(1)
)

consensus_component_order <- candidate_trajectories |>
  dplyr::group_by(dataset, consensus_component) |>
  dplyr::summarise(
    max_n_cells = max(n_cells),
    max_fraction = max(fraction),
    .groups = "drop"
  ) |>
  dplyr::arrange(dataset, dplyr::desc(max_n_cells), dplyr::desc(max_fraction), consensus_component) |>
  dplyr::group_by(dataset) |>
  dplyr::mutate(
    consensus_label = paste0(dataset, "_stable_cluster_", dplyr::row_number())
  ) |>
  dplyr::ungroup()

candidate_trajectories <- candidate_trajectories |>
  dplyr::left_join(
    consensus_component_order |>
      dplyr::select(dataset, consensus_component, consensus_label),
    by = c("dataset", "consensus_component")
  )

consensus_summary_df <- candidate_trajectories |>
  dplyr::group_by(dataset, consensus_label) |>
  dplyr::summarise(
    n_trajectories_merged = dplyr::n(),
    minPts_values = paste(sort(unique(minPts)), collapse = "/"),
    eps_min = min(eps_min),
    eps_max = max(eps_max),
    representative_eps = median(representative_eps),
    n_cells = max(n_cells),
    fraction = max(fraction),
    .groups = "drop"
  ) |>
  dplyr::arrange(dataset, consensus_label)

# Pick the longest stable trajectory in each consensus group. This is the
# deterministic rule that produced the final table: ties go to the lower
# minPts, then to the larger representative cluster.
selected_trajectories <- candidate_trajectories |>
  dplyr::filter(
    n_eps_states >= minimum_stable_states_for_final_selection
  ) |>
  dplyr::group_by(dataset, consensus_label) |>
  dplyr::arrange(
    dplyr::desc(n_eps_states),
    minPts,
    dplyr::desc(n_cells),
    trajectory_id,
    .by_group = TRUE
  ) |>
  dplyr::slice_head(n = 1L) |>
  dplyr::ungroup()

missing_selected_consensus <- setdiff(
  unique(candidate_trajectories$consensus_label),
  unique(selected_trajectories$consensus_label)
)
if (length(missing_selected_consensus) > 0L) {
  stop(
    "A 5%-eligible consensus cluster has no representative with at least ",
    minimum_stable_states_for_final_selection,
    " eps states: ", paste(missing_selected_consensus, collapse = ", "),
    call. = FALSE
  )
}

selected_states <- trajectory_states_df |>
  dplyr::inner_join(
    selected_trajectories |>
      dplyr::select(
        dataset, trajectory_id, consensus_label,
        trajectory_representative_eps = representative_eps
      ),
    by = c("dataset", "trajectory_id")
  ) |>
  dplyr::mutate(
    distance_to_representative = abs(eps - trajectory_representative_eps)
  ) |>
  dplyr::group_by(dataset, consensus_label) |>
  dplyr::arrange(
    distance_to_representative,
    dplyr::desc(n_cells),
    eps,
    .by_group = TRUE
  ) |>
  dplyr::slice_head(n = 1L) |>
  dplyr::ungroup()

selected_representatives_df <- selected_states |>
  dplyr::left_join(
    trajectory_summary_df |>
      dplyr::select(
        dataset, minPts, trajectory_id,
        n_eps_states, eps_min, eps_max,
        trajectory_representative_eps = representative_eps
      ),
    by = c("dataset", "minPts", "trajectory_id")
  ) |>
  dplyr::select(
    dataset, consensus_label, minPts, eps, raw_cluster,
    n_cells, fraction, node, trajectory_id, n_eps_states,
    eps_min, eps_max, trajectory_representative_eps
  )

all_representatives_df <- candidate_trajectories |>
  dplyr::select(
    dataset, consensus_label, minPts, eps, raw_cluster,
    n_cells, fraction, node, trajectory_id, n_eps_states,
    eps_min, eps_max, representative_eps
  ) |>
  dplyr::mutate(
    selected_for_final_label = node %in% selected_representatives_df$node
  ) |>
  dplyr::arrange(dataset, consensus_label, dplyr::desc(selected_for_final_label), minPts)

write.csv(
  consensus_summary_df,
  file.path(output_dir, "stable_consensus_clusters_all_four_datasets.csv"),
  row.names = FALSE
)
write.csv(
  all_representatives_df,
  file.path(output_dir, "stable_cluster_representatives_all_four_datasets.csv"),
  row.names = FALSE
)

################################################################################
# 7. Assign stable1/stable2/unstable labels to all four datasets
################################################################################

cell_metadata_rows <- list()

for (dataset in dataset_order) {
  dataset_cells_now <- rownames(embedding_cache[[dataset]]$coords)
  labels_out <- data.frame(
    dataset = dataset,
    cell = dataset_cells_now,
    raw_dbscan_label = paste0(dataset, "_DBSCAN_unstable"),
    stable_label = paste0(dataset, "_DBSCAN_unstable"),
    selected_minPts = NA_integer_,
    selected_eps = NA_real_,
    selected_raw_cluster = NA_character_,
    stringsAsFactors = FALSE
  )

  reps_dataset <- selected_representatives_df |>
    dplyr::filter(dataset == !!dataset)
  if (nrow(reps_dataset) > 0L) {
    for (rep_i in seq_len(nrow(reps_dataset))) {
      rep_row <- reps_dataset[rep_i, , drop = FALSE]
      labels_state <- assignments_cache[[
        state_key(dataset, rep_row$minPts, rep_row$eps)
      ]]
      cluster_cells <- names(labels_state)[
        labels_state == rep_row$raw_cluster
      ]
      cluster_cells <- intersect(cluster_cells, labels_out$cell)
      if (length(cluster_cells) == 0L) {
        stop(
          "Selected representative has no cells: ",
          rep_row$node,
          call. = FALSE
        )
      }
      hits <- match(cluster_cells, labels_out$cell)
      if (any(!is.na(labels_out$stable_label[hits]) &
              labels_out$stable_label[hits] != paste0(dataset, "_DBSCAN_unstable"))) {
        stop("Selected stable representatives overlap in ", dataset, ".", call. = FALSE)
      }
      labels_out$raw_dbscan_label[hits] <- paste0(dataset, "_", rep_row$raw_cluster)
      labels_out$stable_label[hits] <- rep_row$consensus_label
      labels_out$selected_minPts[hits] <- rep_row$minPts
      labels_out$selected_eps[hits] <- rep_row$eps
      labels_out$selected_raw_cluster[hits] <- rep_row$raw_cluster
    }
  }
  cell_metadata_rows[[dataset]] <- labels_out
}

four_dataset_metadata <- dplyr::bind_rows(cell_metadata_rows) |>
  dplyr::arrange(factor(dataset, levels = dataset_order), cell)

write.csv(
  four_dataset_metadata,
  file.path(output_dir, "four_dataset_dbscan_cell_metadata.csv"),
  row.names = FALSE
)

raw_proportions <- four_dataset_metadata |>
  dplyr::count(dataset, raw_dbscan_label, name = "n_cells") |>
  dplyr::group_by(dataset) |>
  dplyr::mutate(
    proportion = n_cells / sum(n_cells),
    percentage = 100 * proportion
  ) |>
  dplyr::ungroup() |>
  dplyr::arrange(dataset, dplyr::desc(n_cells))

write.csv(
  raw_proportions,
  file.path(output_dir, "cluster_proportions_by_dataset_raw_dbscan_label.csv"),
  row.names = FALSE
)

################################################################################
# 8. Marker signatures and cross-dataset Fisher comparisons
################################################################################

common_genes <- Reduce(
  intersect,
  lapply(embedding_cache[dataset_order], function(x) rownames(x$obj))
)
if (length(common_genes) < 100L) {
  stop("Fewer than 100 genes are shared across all four datasets.", call. = FALSE)
}

signature_rows <- list()
signature_genes <- list()

for (rep_i in seq_len(nrow(selected_representatives_df))) {
  rep_row <- selected_representatives_df[rep_i, , drop = FALSE]
  dataset <- rep_row$dataset
  labels_state <- assignments_cache[[
    state_key(dataset, rep_row$minPts, rep_row$eps)
  ]]
  markers <- top_marker_signature(
    embedding_cache[[dataset]]$obj,
    labels_state,
    rep_row$raw_cluster,
    signature_size = signature_size
  )
  if (nrow(markers) == 0L) {
    warning(
      "No marker signature passed thresholds for ",
      rep_row$consensus_label,
      call. = FALSE
    )
    signature_genes[[rep_row$consensus_label]] <- character()
    next
  }

  markers$dataset <- dataset
  markers$cluster <- rep_row$raw_cluster
  markers$consensus_label <- rep_row$consensus_label
  markers$rank <- seq_len(nrow(markers))
  signature_rows[[length(signature_rows) + 1L]] <- markers |>
    dplyr::select(dataset, consensus_label, cluster, rank, gene, avg_log2FC, p_val)
  signature_genes[[rep_row$consensus_label]] <- unique(markers$gene)
}

signature_df <- dplyr::bind_rows(signature_rows)
write.csv(
  signature_df,
  file.path(output_dir, "selected_stable_cluster_top50_marker_signatures.csv"),
  row.names = FALSE
)

selected_labels <- selected_representatives_df$consensus_label
selected_labels <- unique(selected_labels)
signature_genes <- signature_genes[selected_labels]
signature_genes <- lapply(signature_genes, function(x) intersect(x, common_genes))

pair_index <- utils::combn(selected_labels, 2L, simplify = FALSE)
fisher_rows <- lapply(pair_index, function(pair_labels) {
  label_1 <- pair_labels[1]
  label_2 <- pair_labels[2]
  dataset_1 <- selected_representatives_df$dataset[
    match(label_1, selected_representatives_df$consensus_label)
  ]
  dataset_2 <- selected_representatives_df$dataset[
    match(label_2, selected_representatives_df$consensus_label)
  ]
  overlap <- fisher_overlap(
    signature_genes[[label_1]],
    signature_genes[[label_2]],
    common_genes
  )
  cbind(
    data.frame(
      dataset_1 = dataset_1,
      cluster_1 = selected_representatives_df$raw_cluster[
        match(label_1, selected_representatives_df$consensus_label)
      ],
      stable_label_1 = label_1,
      dataset_2 = dataset_2,
      cluster_2 = selected_representatives_df$raw_cluster[
        match(label_2, selected_representatives_df$consensus_label)
      ],
      stable_label_2 = label_2,
      stringsAsFactors = FALSE
    ),
    overlap
  )
})

fisher_df <- dplyr::bind_rows(fisher_rows) |>
  dplyr::mutate(
    fdr_bh = p.adjust(fisher_p_value, method = "BH"),
    cross_dataset = dataset_1 != dataset_2
  ) |>
  dplyr::arrange(fdr_bh, fisher_p_value)

write.csv(
  fisher_df,
  file.path(output_dir, "all_cross_dataset_cluster_fisher_comparisons.csv"),
  row.names = FALSE
)

################################################################################
# 9. Merge marker-equivalent clusters and compute proportions
################################################################################

merge_candidates <- fisher_df |>
  dplyr::filter(
    cross_dataset,
    overlap_n >= merge_min_overlap,
    overlap_jaccard >= merge_min_jaccard,
    fdr_bh < merge_fdr_threshold
  )

merge_nodes <- selected_representatives_df |>
  dplyr::transmute(
    node = paste(dataset, raw_cluster, sep = "::"),
    dataset,
    cluster = raw_cluster,
    stable_label = consensus_label
  )
merge_parent <- setNames(merge_nodes$node, merge_nodes$node)

if (nrow(merge_candidates) > 0L) {
  for (merge_i in seq_len(nrow(merge_candidates))) {
    node_1 <- paste(
      merge_candidates$dataset_1[merge_i],
      merge_candidates$cluster_1[merge_i],
      sep = "::"
    )
    node_2 <- paste(
      merge_candidates$dataset_2[merge_i],
      merge_candidates$cluster_2[merge_i],
      sep = "::"
    )
    root_1 <- find_root(merge_parent, node_1)
    root_2 <- find_root(merge_parent, node_2)
    if (!identical(root_1, root_2)) merge_parent[[root_2]] <- root_1
  }
}

merge_nodes$merge_component <- vapply(
  merge_nodes$node,
  function(x) find_root(merge_parent, x),
  character(1)
)
merge_component_labels <- merge_nodes |>
  dplyr::distinct(merge_component) |>
  dplyr::arrange(merge_component) |>
  dplyr::mutate(merged_label = paste0("DBSCAN_merged_", dplyr::row_number()))
merge_nodes <- merge_nodes |>
  dplyr::left_join(merge_component_labels, by = "merge_component")

merge_edges_df <- merge_candidates |>
  dplyr::mutate(
    node_1 = paste(dataset_1, cluster_1, sep = "::"),
    node_2 = paste(dataset_2, cluster_2, sep = "::")
  ) |>
  dplyr::select(
    node_1, dataset_1, cluster_1,
    node_2, dataset_2, cluster_2,
    overlap_n, overlap_jaccard,
    fisher_odds_ratio, fisher_p_value, fdr_bh,
    overlap_genes
  )

write.csv(
  merge_edges_df,
  file.path(output_dir, "cross_dataset_cluster_merge_edges.csv"),
  row.names = FALSE
)

merged_label_by_stable <- merge_nodes |>
  dplyr::select(stable_label, merged_label)
four_dataset_metadata <- four_dataset_metadata |>
  dplyr::left_join(merged_label_by_stable, by = "stable_label") |>
  dplyr::mutate(
    merged_label = ifelse(
      is.na(merged_label) | grepl("_DBSCAN_unstable$", stable_label),
      stable_label,
      merged_label
    )
  )
write.csv(
  four_dataset_metadata,
  file.path(output_dir, "four_dataset_dbscan_cell_metadata.csv"),
  row.names = FALSE
)

merged_proportions <- four_dataset_metadata |>
  dplyr::count(dataset, merged_label, name = "n_cells") |>
  dplyr::group_by(dataset) |>
  dplyr::mutate(
    proportion = n_cells / sum(n_cells),
    percentage = 100 * proportion
  ) |>
  dplyr::ungroup() |>
  dplyr::arrange(dataset, dplyr::desc(n_cells))

write.csv(
  merged_proportions,
  file.path(output_dir, "cluster_proportions_by_dataset_merged_label.csv"),
  row.names = FALSE
)

################################################################################
# 10. Map labels back to our object and the raw Kolz object
################################################################################

add_metadata_columns <- function(obj, metadata_df) {
  new_cols <- c(
    "four_dataset_dbscan_dataset",
    "four_dataset_dbscan_raw_label",
    "four_dataset_dbscan_stable_label",
    "four_dataset_dbscan_merged_label",
    "four_dataset_dbscan_minPts",
    "four_dataset_dbscan_eps"
  )
  for (col in new_cols) {
    if (col %in% c(
      "four_dataset_dbscan_minPts",
      "four_dataset_dbscan_eps"
    )) {
      obj[[col]] <- NA_real_
    } else {
      obj[[col]] <- NA_character_
    }
  }

  keep <- intersect(colnames(obj), metadata_df$cell)
  if (length(keep) == 0L) return(obj)
  md <- metadata_df[match(keep, metadata_df$cell), , drop = FALSE]
  rownames(md) <- keep
  obj@meta.data[keep, "four_dataset_dbscan_dataset"] <- md$dataset
  obj@meta.data[keep, "four_dataset_dbscan_raw_label"] <- md$raw_dbscan_label
  obj@meta.data[keep, "four_dataset_dbscan_stable_label"] <- md$stable_label
  obj@meta.data[keep, "four_dataset_dbscan_merged_label"] <- md$merged_label
  obj@meta.data[keep, "four_dataset_dbscan_minPts"] <- md$selected_minPts
  obj@meta.data[keep, "four_dataset_dbscan_eps"] <- md$selected_eps
  obj
}

our_labeled <- add_metadata_columns(
  our_full,
  four_dataset_metadata |>
    dplyr::filter(dataset %in% c("our_MNG", "our_dCLN"))
)
kolz_labeled <- add_metadata_columns(
  kolz_full,
  four_dataset_metadata |>
    dplyr::filter(dataset %in% c("kolz_Th1", "kolz_Th17"))
)

saveRDS(
  our_labeled,
  file.path(output_dir, "srt_fullannot_with_four_dataset_stable_dbscan_labels.rds")
)
saveRDS(
  subset(our_labeled, cells = unique(c(our_mng_cells, our_dcln_cells))),
  file.path(output_dir, "validated_bcell_obj_four_dataset_stable_dbscan_labels.rds")
)
saveRDS(
  kolz_labeled,
  file.path(output_dir, "GSE279684_Th1_Th17_with_four_dataset_stable_dbscan_labels.rds")
)

################################################################################
# 11. Final report
################################################################################

final_long <- four_dataset_metadata |>
  dplyr::count(dataset, stable_label, name = "n_cells") |>
  dplyr::group_by(dataset) |>
  dplyr::mutate(
    proportion = n_cells / sum(n_cells),
    percentage = 100 * proportion
  ) |>
  dplyr::ungroup()

stable1_label <- setNames(
  paste0(dataset_order, "_stable_cluster_1"),
  dataset_order
)
stable2_label <- setNames(
  paste0(dataset_order, "_stable_cluster_2"),
  dataset_order
)
unstable_label <- setNames(
  paste0(dataset_order, "_DBSCAN_unstable"),
  dataset_order
)

report_df <- dplyr::bind_rows(lapply(dataset_order, function(dataset) {
  get_pct <- function(label) {
    value <- final_long$percentage[
      final_long$dataset == dataset & final_long$stable_label == label
    ]
    if (length(value) == 0L) 0 else value[1]
  }
  data.frame(
    Dataset = dataset,
    Stable_cluster_1 = get_pct(stable1_label[[dataset]]),
    Stable_cluster_2 = get_pct(stable2_label[[dataset]]),
    Unstable = get_pct(unstable_label[[dataset]]),
    stringsAsFactors = FALSE
  )
}))

write.csv(
  report_df,
  file.path(output_dir, "stable_cluster_composition_summary.csv"),
  row.names = FALSE
)

cat("\n====================\n")
cat("Final stable-cluster composition\n")
cat("====================\n")
print(report_df, row.names = FALSE)

cat("\nOutput directory:\n", output_dir, "\n", sep = "")
cat("Stable trajectories: ", sum(trajectory_summary_df$stable_trajectory), "\n", sep = "")
cat("Final selected clusters: ", nrow(selected_representatives_df), "\n", sep = "")
cat("Cross-dataset Fisher comparisons: ", nrow(fisher_df), "\n", sep = "")
cat("Cross-dataset merge edges: ", nrow(merge_edges_df), "\n", sep = "")
