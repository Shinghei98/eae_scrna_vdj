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
  library(msigdbr)
  library(gprofiler2)
  library(fgsea)
  library(AnnotationDbi)
  library(org.Mm.eg.db)
  library(KEGGREST)
  library(jsonlite)
  library(babelgene)
  library(grid)
})

################################################################################
# Figure 2a: MNG versus dCLN B-cell MAST DEG, enrichment, and volcano plot
################################################################################

# This section is intentionally self-contained. It starts from annotation.R's
# srt_fullannot.rds product, runs cell-level MAST on validated B cells, performs
# directional Hallmark and GO:BP over-representation analysis, and generates
# the final ontology-balanced volcano plot.
#
# Direction:
#   positive avg_log2FC = higher in MNG
#   negative avg_log2FC = higher in dCLN
#
# DEG/ORA and volcano thresholds:
#   abs(avg_log2FC) > 0.2
#   BH-adjusted MAST q_hurdle < 0.05
#
# Statistical limitation:
#   MAST treats cells as observations. M1/M2 and L1/L2 do not provide enough
#   independent replicates for a definitive sample-level tissue comparison.

figure_2a <- local({
  seed_use <- 1234L
  set.seed(seed_use)

  project_root <- Sys.getenv(
    "EAE_PROJECT_ROOT",
    unset = "/storage3/fs1/gfwu/Active/David/mng_dcln_project"
  )
  if (!dir.exists(project_root) &&
      dir.exists("/Volumes/gfwu/Active/David/mng_dcln_project")) {
    project_root <- "/Volumes/gfwu/Active/David/mng_dcln_project"
  }
  local_analysis_root <- Sys.getenv(
    "EAE_LOCAL_ANALYSIS_ROOT",
    unset = "/Users/shingheimok/Desktop/Phd_Study/research/Wu_lab/b_cells"
  )
  preprocess_dir <- file.path(project_root, "eae_scrna_vdj_preprocess")
  annotated_rds <- file.path(preprocess_dir, "srt_fullannot.rds")
  if (!file.exists(annotated_rds)) {
    stop("Annotated object not found: ", annotated_rds, call. = FALSE)
  }

  deg_output_dir <- file.path(
    preprocess_dir,
    "DEG",
    "validated_bcells_MNG_vs_dCLN_MAST_msigdbr"
  )
  volcano_output_dir <- file.path(deg_output_dir, "volcano")
  figure_output_dir <- Sys.getenv(
    "EAE_FIGURE2_OUTPUT_DIR",
    unset = "/Users/shingheimok/Desktop/Phd_Study/research/Wu_lab/b_cells/figures/figure_2"
  )
  dir.create(deg_output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(volcano_output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figure_output_dir, recursive = TRUE, showWarnings = FALSE)

  sample_order <- c("M1", "M2", "L1", "L2")
  deg_fdr_threshold <- 0.05
  minimum_abs_log2fc <- 0.2
  minimum_detection_fraction <- 0.01
  minimum_gene_set_size_in_universe <- 5L
  y_display_cap <- as.numeric(Sys.getenv("VOLCANO_Y_LIMIT", unset = "200"))

  if (!requireNamespace("MAST", quietly = TRUE)) {
    stop("The MAST package is required for Figure 2a.", call. = FALSE)
  }
  if (!is.finite(y_display_cap) || y_display_cap <= 0) {
    stop("VOLCANO_Y_LIMIT must be positive.", call. = FALSE)
  }

  extract_sample_id <- function(values) {
    values <- as.character(values)
    hits <- regmatches(
      values,
      regexpr(
        "(?<![A-Za-z0-9])(M1|M2|L1|L2)(?![A-Za-z0-9])",
        values,
        perl = TRUE
      )
    )
    hits[hits == ""] <- NA_character_
    hits
  }

  choose_sample_column <- function(metadata) {
    candidates <- intersect(
      c("sample_id", "orig.ident", "sample", "Sample", "replicate"),
      colnames(metadata)
    )
    if (length(candidates) == 0L) {
      stop("No sample identifier column found.", call. = FALSE)
    }
    parsed <- lapply(candidates, function(column) {
      extract_sample_id(metadata[[column]])
    })
    names(parsed) <- candidates
    scores <- vapply(parsed, function(x) sum(!is.na(x)), integer(1))
    selected <- candidates[which.max(scores)]
    if (scores[selected] == 0L) {
      stop("Could not parse M1/M2/L1/L2 from metadata.", call. = FALSE)
    }
    list(column = selected, parsed = parsed[[selected]])
  }

  load_msigdb_sets <- function(collection, subcollections = NULL) {
    mouse_native_collection <- dplyr::recode(
      collection,
      H = "MH",
      C5 = "M5",
      .default = collection
    )
    call_msigdbr <- function(use_mouse_database) {
      available_arguments <- names(formals(msigdbr::msigdbr))
      args <- list(species = "Mus musculus")
      collection_for_call <- if (use_mouse_database) {
        mouse_native_collection
      } else {
        collection
      }
      if ("collection" %in% available_arguments) {
        args$collection <- collection_for_call
      } else if ("category" %in% available_arguments) {
        args$category <- collection_for_call
      }
      if (use_mouse_database && "db_species" %in% available_arguments) {
        args$db_species <- "MM"
      }
      tryCatch(do.call(msigdbr::msigdbr, args), error = function(e) NULL)
    }

    sets <- call_msigdbr(TRUE)
    if (is.null(sets) || nrow(sets) == 0L) sets <- call_msigdbr(FALSE)
    if (is.null(sets) || nrow(sets) == 0L) {
      sets <- tryCatch(
        msigdbr::msigdbr(species = "Mus musculus"),
        error = function(e) NULL
      )
    }
    if (is.null(sets) || nrow(sets) == 0L) {
      stop("msigdbr returned no gene sets.", call. = FALSE)
    }

    collection_column <- intersect(c("gs_collection", "gs_cat"), colnames(sets))[1]
    subcollection_column <- intersect(
      c("gs_subcollection", "gs_subcat"),
      colnames(sets)
    )[1]
    if (!is.na(collection_column)) {
      sets <- sets[
        sets[[collection_column]] %in% c(collection, mouse_native_collection),
        ,
        drop = FALSE
      ]
    }
    if (!is.null(subcollections) && !is.na(subcollection_column)) {
      sets <- sets[
        sets[[subcollection_column]] %in% subcollections,
        ,
        drop = FALSE
      ]
    }
    if (!all(c("gs_name", "gene_symbol") %in% colnames(sets))) {
      stop("Unexpected msigdbr output columns.", call. = FALSE)
    }
    sets |>
      dplyr::select(.data$gs_name, .data$gene_symbol) |>
      dplyr::filter(
        !is.na(.data$gs_name),
        !is.na(.data$gene_symbol),
        nzchar(.data$gene_symbol)
      ) |>
      dplyr::distinct()
  }

  run_overrepresentation <- function(
      selected_genes,
      universe_genes,
      set_table,
      collection_name,
      direction_name) {
    selected_genes <- unique(intersect(selected_genes, universe_genes))
    if (length(selected_genes) == 0L) return(data.frame())
    set_list <- lapply(split(set_table$gene_symbol, set_table$gs_name), unique)

    result_rows <- lapply(names(set_list), function(term) {
      set_genes <- intersect(set_list[[term]], universe_genes)
      if (length(set_genes) < minimum_gene_set_size_in_universe) return(NULL)
      selected_in <- intersect(selected_genes, set_genes)
      selected_out <- setdiff(selected_genes, set_genes)
      background_in <- setdiff(set_genes, selected_genes)
      background_out <- setdiff(
        setdiff(universe_genes, selected_genes),
        set_genes
      )
      test <- fisher.test(
        matrix(
          c(
            length(selected_in), length(selected_out),
            length(background_in), length(background_out)
          ),
          nrow = 2L,
          byrow = TRUE
        ),
        alternative = "greater"
      )
      data.frame(
        collection = collection_name,
        direction = direction_name,
        term = term,
        n_set_genes_in_universe = length(set_genes),
        n_selected_genes_in_set = length(selected_in),
        n_selected_genes = length(selected_genes),
        n_universe_genes = length(universe_genes),
        odds_ratio = unname(test$estimate),
        p_value = test$p.value,
        selected_genes_in_set = paste(sort(selected_in), collapse = ";"),
        stringsAsFactors = FALSE
      )
    })
    result <- dplyr::bind_rows(result_rows)
    if (nrow(result) == 0L) return(result)
    result |>
      dplyr::mutate(fdr_bh = p.adjust(.data$p_value, method = "BH")) |>
      dplyr::arrange(.data$fdr_bh, .data$p_value)
  }

  write_enrichment_outputs <- function(result, prefix) {
    write.csv(
      result,
      file.path(deg_output_dir, paste0(prefix, "_all_tested_terms.csv")),
      row.names = FALSE
    )
    significant <- if (nrow(result) > 0L) {
      result[result$fdr_bh < 0.05, , drop = FALSE]
    } else {
      result
    }
    write.csv(
      significant,
      file.path(deg_output_dir, paste0(prefix, "_FDR_lt_0.05.csv")),
      row.names = FALSE
    )
    invisible(significant)
  }

  cat("Loading annotation.R output:\n", annotated_rds, "\n", sep = "")
  obj <- readRDS(annotated_rds)
  DefaultAssay(obj) <- "RNA"
  metadata <- obj@meta.data

  if (!"celltype_major" %in% colnames(metadata)) {
    stop("Missing celltype_major in annotation object.", call. = FALSE)
  }
  bcell_minor_column <- if (
    "celltype_minor_pre_bcell_dbscan" %in% colnames(metadata)
  ) {
    "celltype_minor_pre_bcell_dbscan"
  } else if ("celltype_minor" %in% colnames(metadata)) {
    "celltype_minor"
  } else {
    stop("No B-cell minor annotation column found.", call. = FALSE)
  }

  is_validated_bcell <-
    as.character(metadata$celltype_major) == "B_cell" &
    as.character(metadata[[bcell_minor_column]]) == "B_cell"
  sample_parse <- choose_sample_column(metadata)
  metadata$MAST_sample <- sample_parse$parsed
  metadata$MAST_tissue <- ifelse(
    metadata$MAST_sample %in% c("M1", "M2"),
    "MNG",
    ifelse(metadata$MAST_sample %in% c("L1", "L2"), "dCLN", NA_character_)
  )

  selected_cells <- rownames(metadata)[
    is_validated_bcell &
      metadata$MAST_sample %in% sample_order &
      !is.na(metadata$MAST_tissue)
  ]
  selected_cells <- intersect(selected_cells, colnames(obj))
  if (length(selected_cells) == 0L) {
    stop("No validated M1/M2/L1/L2 B cells were selected.", call. = FALSE)
  }

  sample_counts <- table(
    factor(metadata[selected_cells, "MAST_sample"], levels = sample_order)
  )
  if (any(sample_counts < 100L)) {
    stop(
      "A sample has fewer than 100 validated B cells: ",
      paste(names(sample_counts), sample_counts, sep = "=", collapse = ", "),
      call. = FALSE
    )
  }

  bcell_obj <- subset(obj, cells = selected_cells)
  bcell_obj$MAST_sample <- metadata[colnames(bcell_obj), "MAST_sample"]
  bcell_obj$MAST_tissue <- metadata[colnames(bcell_obj), "MAST_tissue"]
  DefaultAssay(bcell_obj) <- "RNA"
  bcell_obj <- tryCatch(
    SeuratObject::JoinLayers(bcell_obj, assay = "RNA"),
    error = function(e) bcell_obj
  )
  bcell_obj <- NormalizeData(
    bcell_obj,
    normalization.method = "LogNormalize",
    scale.factor = 10000,
    verbose = FALSE
  )

  sample_audit <- data.frame(
    sample = sample_order,
    tissue = c("MNG", "MNG", "dCLN", "dCLN"),
    n_validated_bcells = as.integer(sample_counts),
    sample_metadata_column = sample_parse$column,
    bcell_minor_column = bcell_minor_column,
    stringsAsFactors = FALSE
  )
  write.csv(
    sample_audit,
    file.path(deg_output_dir, "selected_cell_audit.csv"),
    row.names = FALSE
  )

  technical_covariates <- intersect(
    c("nCount_RNA", "percent.mt", "percent_mt"),
    colnames(bcell_obj@meta.data)
  )
  if (all(c("percent.mt", "percent_mt") %in% technical_covariates)) {
    technical_covariates <- setdiff(technical_covariates, "percent_mt")
  }
  Idents(bcell_obj) <- factor(bcell_obj$MAST_tissue, levels = c("dCLN", "MNG"))

  cat("Running cell-level MAST: MNG versus dCLN\n")
  set.seed(seed_use)
  mast_result <- FindMarkers(
    bcell_obj,
    ident.1 = "MNG",
    ident.2 = "dCLN",
    assay = "RNA",
    test.use = "MAST",
    logfc.threshold = 0,
    min.pct = minimum_detection_fraction,
    only.pos = FALSE,
    latent.vars = if (length(technical_covariates) > 0L) {
      technical_covariates
    } else {
      NULL
    },
    verbose = TRUE
  )
  if (nrow(mast_result) == 0L) stop("MAST returned no genes.", call. = FALSE)

  mast_result$gene <- rownames(mast_result)
  fold_change_column <- intersect(
    c("avg_log2FC", "avg_logFC"),
    colnames(mast_result)
  )[1]
  if (is.na(fold_change_column) || !"p_val" %in% colnames(mast_result)) {
    stop("MAST output lacks its fold-change or p-value column.", call. = FALSE)
  }
  if ("p_val_adj" %in% colnames(mast_result)) {
    mast_result$p_val_adj_Seurat <- mast_result$p_val_adj
  }
  mast_result$avg_log2FC_MNG_vs_dCLN <- mast_result[[fold_change_column]]
  mast_result$q_hurdle <- p.adjust(mast_result$p_val, method = "BH")
  mast_result$direction <- dplyr::case_when(
    mast_result$avg_log2FC_MNG_vs_dCLN > 0 ~ "MNG_higher",
    mast_result$avg_log2FC_MNG_vs_dCLN < 0 ~ "dCLN_higher",
    TRUE ~ "no_change"
  )
  mast_result$significant_fdr_0.05 <-
    !is.na(mast_result$q_hurdle) &
    mast_result$q_hurdle < deg_fdr_threshold
  mast_result <- mast_result |>
    dplyr::select(
      .data$gene,
      .data$avg_log2FC_MNG_vs_dCLN,
      .data$p_val,
      .data$q_hurdle,
      .data$direction,
      .data$significant_fdr_0.05,
      dplyr::everything()
    ) |>
    dplyr::arrange(
      .data$q_hurdle,
      dplyr::desc(abs(.data$avg_log2FC_MNG_vs_dCLN))
    )
  write.csv(
    mast_result,
    file.path(deg_output_dir, "validated_bcells_MNG_vs_dCLN_MAST_all_genes.csv"),
    row.names = FALSE
  )

  universe_genes <- unique(mast_result$gene[!is.na(mast_result$gene)])
  mng_higher_genes <- unique(mast_result$gene[
    mast_result$q_hurdle < deg_fdr_threshold &
      mast_result$avg_log2FC_MNG_vs_dCLN > minimum_abs_log2fc
  ])
  dcln_higher_genes <- unique(mast_result$gene[
    mast_result$q_hurdle < deg_fdr_threshold &
      mast_result$avg_log2FC_MNG_vs_dCLN < -minimum_abs_log2fc
  ])
  write.csv(
    mast_result[mast_result$gene %in% mng_higher_genes, , drop = FALSE],
    file.path(deg_output_dir, "validated_bcells_MNG_higher_MAST_DEGs.csv"),
    row.names = FALSE
  )
  write.csv(
    mast_result[mast_result$gene %in% dcln_higher_genes, , drop = FALSE],
    file.path(deg_output_dir, "validated_bcells_dCLN_higher_MAST_DEGs.csv"),
    row.names = FALSE
  )

  cat("Loading mouse MSigDB Hallmark and GO:BP gene sets\n")
  hallmark_sets <- load_msigdb_sets("H")
  go_bp_sets <- load_msigdb_sets("C5", c("GO:BP", "BP", "GOBP"))
  if (nrow(hallmark_sets) == 0L || nrow(go_bp_sets) == 0L) {
    stop("No Hallmark or GO:BP sets were returned.", call. = FALSE)
  }

  hallmark_mng <- run_overrepresentation(
    mng_higher_genes, universe_genes, hallmark_sets, "Hallmark", "MNG_higher"
  )
  hallmark_dcln <- run_overrepresentation(
    dcln_higher_genes, universe_genes, hallmark_sets, "Hallmark", "dCLN_higher"
  )
  go_mng <- run_overrepresentation(
    mng_higher_genes, universe_genes, go_bp_sets, "GO_BP", "MNG_higher"
  )
  go_dcln <- run_overrepresentation(
    dcln_higher_genes, universe_genes, go_bp_sets, "GO_BP", "dCLN_higher"
  )
  write_enrichment_outputs(hallmark_mng, "MNG_higher_Hallmark_enrichment")
  write_enrichment_outputs(hallmark_dcln, "dCLN_higher_Hallmark_enrichment")
  write_enrichment_outputs(go_mng, "MNG_higher_GO_BP_enrichment")
  write_enrichment_outputs(go_dcln, "dCLN_higher_GO_BP_enrichment")
  write.csv(
    dplyr::bind_rows(hallmark_mng, hallmark_dcln),
    file.path(deg_output_dir, "combined_Hallmark_enrichment_all_directions.csv"),
    row.names = FALSE
  )
  write.csv(
    dplyr::bind_rows(go_mng, go_dcln),
    file.path(deg_output_dir, "combined_GO_BP_enrichment_all_directions.csv"),
    row.names = FALSE
  )

  analysis_settings <- data.frame(
    parameter = c(
      "seed", "comparison", "positive_logFC_direction", "selected_samples",
      "selected_cells", "MAST_test_level", "MAST_min_pct",
      "MAST_technical_covariates", "DEG_FDR_threshold",
      "DEG_minimum_abs_log2FC", "enrichment_method", "enrichment_universe",
      "minimum_gene_set_size_in_universe"
    ),
    value = c(
      seed_use,
      "MNG versus dCLN validated B cells",
      "MNG_higher",
      paste(sample_order, collapse = "/"),
      length(selected_cells),
      "cell-level MAST; not sample-level pseudobulk",
      minimum_detection_fraction,
      ifelse(
        length(technical_covariates) == 0L,
        "none",
        paste(technical_covariates, collapse = "/")
      ),
      deg_fdr_threshold,
      minimum_abs_log2fc,
      "one-sided Fisher over-representation using msigdbr",
      "all genes tested by MAST",
      minimum_gene_set_size_in_universe
    ),
    stringsAsFactors = FALSE
  )
  write.csv(
    analysis_settings,
    file.path(deg_output_dir, "analysis_settings.csv"),
    row.names = FALSE
  )

  deg <- mast_result |>
    dplyr::mutate(
      log2FC = as.numeric(.data$avg_log2FC_MNG_vs_dCLN),
      q_value = as.numeric(.data$q_hurdle)
    ) |>
    dplyr::filter(
      is.finite(.data$log2FC),
      !is.na(.data$q_value),
      .data$q_value >= 0
    )
  smallest_nonzero_q <- suppressWarnings(
    min(deg$q_value[deg$q_value > 0], na.rm = TRUE)
  )
  if (!is.finite(smallest_nonzero_q)) smallest_nonzero_q <- 1e-300
  q_floor <- max(smallest_nonzero_q / 10, 1e-300)
  deg <- deg |>
    dplyr::mutate(
      q_for_plot = pmax(.data$q_value, q_floor),
      neg_log10_q_uncapped = -log10(.data$q_for_plot),
      neg_log10_q = pmin(.data$neg_log10_q_uncapped, y_display_cap),
      volcano_group = dplyr::case_when(
        abs(.data$log2FC) > minimum_abs_log2fc &
          .data$q_value < deg_fdr_threshold & .data$log2FC < 0 ~ "dCLN higher",
        abs(.data$log2FC) > minimum_abs_log2fc &
          .data$q_value < deg_fdr_threshold & .data$log2FC > 0 ~ "MNG higher",
        TRUE ~ "Not significant"
      )
    )

  ontology_label_panel <- data.frame(
    gene = c(
      "Ndufa6", "Cox8a", "Socs1", "Psmb9", "B2m",
      "H2-Q7", "Txn1", "Gpx4", "Tpi1", "Aldoa",
      "Fcer2a", "Cr2", "Cd22", "Bcl6", "Setd1b",
      "Zfp36l1", "Zfp36l2", "Fgd2", "Stap1", "Ralgps2"
    ),
    expected_group = c(rep("MNG higher", 10), rep("dCLN higher", 10)),
    annotation_ontology = c(
      "OXPHOS/electron transport chain complex I",
      "OXPHOS/electron transport chain complex IV",
      "Interferon-gamma response",
      "MHC-I antigen processing/immunoproteasome",
      "MHC-I antigen processing",
      "MHC-I antigen processing",
      "Reactive oxygen species pathway",
      "Reactive oxygen species/OXPHOS",
      "Glucose catabolism",
      "Glucose catabolism",
      "Vesicle-mediated transport/endocytosis",
      "Mature B-cell/complement signaling",
      "B-cell activation",
      "Lymphocyte activation/chromatin regulation",
      "Chromatin organization",
      "mRNA stability and metabolism",
      "mRNA stability and B-cell differentiation",
      "Cytoskeleton organization",
      "Vesicle transport/B-cell signaling",
      "Small-GTPase signaling"
    ),
    stringsAsFactors = FALSE
  )
  label_table <- ontology_label_panel |>
    dplyr::left_join(deg, by = "gene")
  missing_label_genes <- label_table$gene[is.na(label_table$volcano_group)]
  invalid_label_genes <- label_table$gene[
    !is.na(label_table$volcano_group) &
      as.character(label_table$volcano_group) != label_table$expected_group
  ]
  if (length(missing_label_genes) > 0L || length(invalid_label_genes) > 0L) {
    stop(
      "Figure 2a ontology labels missing or outside required thresholds: ",
      paste(unique(c(missing_label_genes, invalid_label_genes)), collapse = ", "),
      call. = FALSE
    )
  }

  dcln_color <- "#5B9DF5"
  mng_color <- "#FF6B6B"
  neutral_color <- "#C8C8C8"
  group_levels <- c("dCLN higher", "Not significant", "MNG higher")
  deg$volcano_group <- factor(deg$volcano_group, levels = group_levels)
  label_table$volcano_group <- factor(
    label_table$volcano_group,
    levels = group_levels
  )
  x_quantile <- as.numeric(
    stats::quantile(abs(deg$log2FC), 0.995, na.rm = TRUE, names = FALSE)
  )
  x_limit <- max(
    2,
    ceiling(2 * max(x_quantile, abs(label_table$log2FC), na.rm = TRUE)) / 2
  )

  arrow_grob <- grid::grobTree(
    grid::polygonGrob(
      x = grid::unit(c(0.49, 0.17, 0.17, 0.01, 0.17, 0.17, 0.49), "npc"),
      y = grid::unit(c(0.32, 0.32, 0.08, 0.50, 0.92, 0.68, 0.68), "npc"),
      gp = grid::gpar(fill = dcln_color, col = NA)
    ),
    grid::polygonGrob(
      x = grid::unit(c(0.51, 0.83, 0.83, 0.99, 0.83, 0.83, 0.51), "npc"),
      y = grid::unit(c(0.32, 0.32, 0.08, 0.50, 0.92, 0.68, 0.68), "npc"),
      gp = grid::gpar(fill = mng_color, col = NA)
    ),
    grid::textGrob(
      "dCLN", x = grid::unit(0.32, "npc"), y = grid::unit(0.50, "npc"),
      gp = grid::gpar(col = "white", fontsize = 15, fontface = "bold")
    ),
    grid::textGrob(
      "MNG", x = grid::unit(0.68, "npc"), y = grid::unit(0.50, "npc"),
      gp = grid::gpar(col = "white", fontsize = 15, fontface = "bold")
    )
  )

  p <- ggplot(deg, aes(x = .data$log2FC, y = .data$neg_log10_q)) +
    geom_point(
      aes(color = .data$volcano_group),
      size = 1.45,
      alpha = 0.72,
      stroke = 0
    ) +
    geom_vline(
      xintercept = c(-minimum_abs_log2fc, minimum_abs_log2fc),
      color = "#8A8A8A",
      linewidth = 0.45,
      linetype = "dotted"
    ) +
    geom_vline(
      xintercept = 0,
      color = "#B8B8B8",
      linewidth = 0.35,
      linetype = "dashed"
    ) +
    geom_hline(
      yintercept = -log10(deg_fdr_threshold),
      color = "#8A8A8A",
      linewidth = 0.45,
      linetype = "dotted"
    ) +
    scale_color_manual(
      values = c(
        "dCLN higher" = dcln_color,
        "Not significant" = neutral_color,
        "MNG higher" = mng_color
      ),
      drop = FALSE
    ) +
    scale_x_continuous(
      breaks = scales::breaks_pretty(n = 8),
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    scale_y_continuous(
      breaks = scales::breaks_pretty(n = 6),
      expand = expansion(mult = c(0, 0.03))
    ) +
    labs(
      x = expression(log[2] * "(Fold change: MNG / dCLN)"),
      y = expression(-log[10] * "(MAST BH-adjusted q-value)")
    ) +
    annotation_custom(
      grob = arrow_grob,
      xmin = -x_limit,
      xmax = x_limit,
      ymin = -0.21 * y_display_cap,
      ymax = -0.055 * y_display_cap
    ) +
    coord_cartesian(
      xlim = c(-x_limit, x_limit),
      ylim = c(0, y_display_cap),
      clip = "off"
    ) +
    theme_classic(base_size = 15, base_family = "sans") +
    theme(
      legend.position = "none",
      axis.title = element_text(face = "bold", color = "black"),
      axis.text = element_text(color = "black"),
      axis.line = element_line(color = "black", linewidth = 0.7),
      axis.ticks = element_line(color = "black", linewidth = 0.6),
      plot.margin = margin(t = 14, r = 28, b = 104, l = 22)
    )

  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p <- p + ggrepel::geom_text_repel(
      data = label_table,
      aes(label = .data$gene, color = .data$volcano_group),
      size = 5,
      fontface = "bold",
      box.padding = 0.42,
      point.padding = 0.22,
      min.segment.length = 0,
      segment.size = 0.45,
      max.overlaps = Inf,
      seed = seed_use,
      direction = "both",
      force = 1.2,
      show.legend = FALSE
    )
  } else {
    p <- p + geom_text(
      data = label_table,
      aes(label = .data$gene, color = .data$volcano_group),
      size = 4.7,
      fontface = "bold",
      vjust = -0.5,
      show.legend = FALSE
    )
  }

  figure_tiff <- file.path(
    figure_output_dir,
    "figure_2a_MNG_vs_dCLN_volcano.tiff"
  )
  ggsave(
    figure_tiff,
    p,
    width = 8.2,
    height = 8.2,
    units = "in",
    dpi = 300,
    compression = "lzw",
    bg = "white"
  )
  write.csv(
    deg,
    file.path(volcano_output_dir, "validated_bcells_MNG_vs_dCLN_volcano_plot_data.csv"),
    row.names = FALSE
  )
  write.csv(
    label_table,
    file.path(volcano_output_dir, "validated_bcells_MNG_vs_dCLN_volcano_labeled_genes.csv"),
    row.names = FALSE
  )
  volcano_audit <- data.frame(
    parameter = c(
      "comparison", "positive_log2FC_direction", "negative_log2FC_direction",
      "absolute_log2FC_threshold", "q_value_column", "q_value_threshold",
      "MNG_higher_genes", "dCLN_higher_genes", "labels_per_direction",
      "label_selection_method", "y_values_capped_for_display_at"
    ),
    value = c(
      "validated B cells: MNG versus dCLN",
      "MNG higher",
      "dCLN higher",
      minimum_abs_log2fc,
      "q_hurdle",
      deg_fdr_threshold,
      sum(deg$volcano_group == "MNG higher"),
      sum(deg$volcano_group == "dCLN higher"),
      10,
      "ontology-balanced hard-coded panel; all labels pass FC and q thresholds",
      y_display_cap
    ),
    stringsAsFactors = FALSE
  )
  write.csv(
    volcano_audit,
    file.path(volcano_output_dir, "validated_bcells_MNG_vs_dCLN_volcano_audit.csv"),
    row.names = FALSE
  )

  cat("Saved Figure 2a:\n", figure_tiff, "\n", sep = "")
  p
})

################################################################################
# OLD Figure 2b: gprofiler2 ORA of MNG-higher and dCLN-higher DEGs
################################################################################

if (FALSE) {

# This section follows the referenced unordered multi-query workflow. MNG- and
# dCLN-higher genes are submitted as separate unordered lists, all genes tested
# by MAST are used as the custom background, and GO, Reactome, WikiPathways,
# KEGG, and CORUM are tested. Only g:SCS-adjusted P < 0.05 terms containing
# fewer than 500 genes are retained.

figure_2b_gprofiler2_result <- local({
  local_analysis_root <- Sys.getenv(
    "EAE_LOCAL_ANALYSIS_ROOT",
    unset = "/Users/shingheimok/Desktop/Phd_Study/research/Wu_lab/b_cells"
  )
  deg_dir <- file.path(
    local_analysis_root,
    "outs", "output", "DEG",
    "validated_bcells_MNG_vs_dCLN_MAST_msigdbr"
  )
  deg_csv <- file.path(
    deg_dir,
    "validated_bcells_MNG_vs_dCLN_MAST_all_genes.csv"
  )
  output_dir <- file.path(deg_dir, "gprofiler2_ORA")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  figure_output_dir <- file.path(
    local_analysis_root,
    "figures", "figure_2"
  )
  dir.create(figure_output_dir, recursive = TRUE, showWarnings = FALSE)

  output_csv <- file.path(
    output_dir,
    "figure_2b_MNG_vs_dCLN_gprofiler2_all_passing_pathways.csv"
  )
  audit_csv <- file.path(
    output_dir,
    "figure_2b_MNG_vs_dCLN_gprofiler2_analysis_audit.csv"
  )

  q_threshold <- 0.05
  minimum_abs_log2fc <- 0.2
  maximum_retained_term_size <- 500L
  organism_use <- "mmusculus"
  correction_method_use <- "g_SCS"
  sources_use <- c(
    "GO:BP", "GO:MF", "GO:CC",
    "REAC", "WP", "KEGG", "CORUM"
  )

  if (!file.exists(deg_csv)) {
    stop("Local Figure 2a DEG table not found: ", deg_csv, call. = FALSE)
  }

  deg <- read.csv(
    deg_csv,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  required_columns <- c(
    "gene",
    "avg_log2FC_MNG_vs_dCLN",
    "q_hurdle"
  )
  missing_columns <- setdiff(required_columns, colnames(deg))
  if (length(missing_columns) > 0L) {
    stop(
      "DEG table is missing required columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  deg$gene <- as.character(deg$gene)
  deg$avg_log2FC_MNG_vs_dCLN <- as.numeric(
    deg$avg_log2FC_MNG_vs_dCLN
  )
  deg$q_hurdle <- as.numeric(deg$q_hurdle)
  deg <- deg[
    !is.na(deg$gene) &
      nzchar(deg$gene) &
      is.finite(deg$avg_log2FC_MNG_vs_dCLN) &
      !is.na(deg$q_hurdle),
    ,
    drop = FALSE
  ]
  if (anyDuplicated(deg$gene)) {
    stop("DEG table contains duplicate gene symbols.", call. = FALSE)
  }

  universe_genes <- unique(deg$gene)
  mng_higher_genes <- unique(deg$gene[
    deg$q_hurdle < q_threshold &
      deg$avg_log2FC_MNG_vs_dCLN > minimum_abs_log2fc
  ])
  dcln_higher_genes <- unique(deg$gene[
    deg$q_hurdle < q_threshold &
      deg$avg_log2FC_MNG_vs_dCLN < -minimum_abs_log2fc
  ])
  query_list <- list(
    MNG_higher = sort(mng_higher_genes),
    dCLN_higher = sort(dcln_higher_genes)
  )

  if (any(lengths(query_list) == 0L)) {
    stop("At least one Figure 2b directional DEG query is empty.", call. = FALSE)
  }
  if (length(intersect(query_list$MNG_higher, query_list$dCLN_higher)) > 0L) {
    stop("MNG-higher and dCLN-higher query lists overlap.", call. = FALSE)
  }

  cat("Figure 2b gprofiler2 ORA input\n")
  cat("  MAST-tested background genes: ", length(universe_genes), "\n", sep = "")
  cat("  MNG-higher genes: ", length(query_list$MNG_higher), "\n", sep = "")
  cat("  dCLN-higher genes: ", length(query_list$dCLN_higher), "\n", sep = "")

  gost_result <- gprofiler2::gost(
    query = query_list,
    organism = organism_use,
    ordered_query = FALSE,
    multi_query = TRUE,
    significant = TRUE,
    exclude_iea = FALSE,
    measure_underrepresentation = FALSE,
    evcodes = FALSE,
    user_threshold = q_threshold,
    correction_method = correction_method_use,
    domain_scope = "custom",
    custom_bg = universe_genes,
    sources = sources_use
  )

  if (is.null(gost_result) || is.null(gost_result$result)) {
    stop("g:Profiler returned no Figure 2b result object.", call. = FALSE)
  }
  multi_result <- gost_result$result
  if (nrow(multi_result) == 0L) {
    stop("No g:Profiler term passed the significance threshold.", call. = FALSE)
  }

  required_result_columns <- c(
    "term_id", "term_name", "source", "term_size",
    "p_values", "significant", "query_sizes", "intersection_sizes",
    "effective_domain_size", "source_order"
  )
  missing_result_columns <- setdiff(
    required_result_columns,
    colnames(multi_result)
  )
  if (length(missing_result_columns) > 0L) {
    stop(
      "Unexpected g:Profiler multi-query output; missing: ",
      paste(missing_result_columns, collapse = ", "),
      call. = FALSE
    )
  }

  query_names <- names(query_list)
  long_rows <- lapply(seq_len(nrow(multi_result)), function(row_index) {
    row <- multi_result[row_index, , drop = FALSE]
    p_values <- as.numeric(unlist(row$p_values[[1]], use.names = FALSE))
    significant <- as.logical(
      unlist(row$significant[[1]], use.names = FALSE)
    )
    query_sizes <- as.numeric(
      unlist(row$query_sizes[[1]], use.names = FALSE)
    )
    intersection_sizes <- as.numeric(
      unlist(row$intersection_sizes[[1]], use.names = FALSE)
    )

    if (!all(
      lengths(list(p_values, significant, query_sizes, intersection_sizes)) ==
        length(query_names)
    )) {
      stop(
        "Unexpected multi-query vector length for ",
        row$term_id,
        call. = FALSE
      )
    }

    parents <- if ("parents" %in% colnames(row)) {
      paste(unlist(row$parents[[1]], use.names = FALSE), collapse = ";")
    } else {
      NA_character_
    }

    data.frame(
      query = query_names,
      direction = ifelse(
        query_names == "MNG_higher",
        "higher in MNG",
        "higher in dCLN"
      ),
      source = as.character(row$source),
      term_id = as.character(row$term_id),
      term_name = as.character(row$term_name),
      adjusted_p_value_gSCS = p_values,
      significant_gSCS_0.05 = significant,
      term_size = as.integer(row$term_size),
      query_size = query_sizes,
      intersection_size = intersection_sizes,
      effective_domain_size = as.integer(row$effective_domain_size),
      source_order = as.integer(row$source_order),
      parents = parents,
      stringsAsFactors = FALSE
    )
  })

  long_result <- do.call(rbind, long_rows)
  passing_result <- long_result[
    long_result$significant_gSCS_0.05 %in% TRUE &
      is.finite(long_result$adjusted_p_value_gSCS) &
      long_result$adjusted_p_value_gSCS < q_threshold &
      long_result$term_size < maximum_retained_term_size,
    ,
    drop = FALSE
  ]
  if (nrow(passing_result) == 0L) {
    stop(
      "No Figure 2b term remained after P < 0.05 and term_size < 500.",
      call. = FALSE
    )
  }

  passing_result <- passing_result[
    order(
      match(passing_result$query, query_names),
      passing_result$adjusted_p_value_gSCS,
      match(passing_result$source, sources_use),
      passing_result$term_name
    ),
    ,
    drop = FALSE
  ]
  rownames(passing_result) <- NULL
  write.csv(passing_result, output_csv, row.names = FALSE)

  database_metadata <- NA_character_
  if (!is.null(gost_result$meta)) {
    database_metadata <- paste(
      capture.output(str(gost_result$meta, max.level = 3L)),
      collapse = " | "
    )
  }
  audit <- data.frame(
    parameter = c(
      "analysis", "input_DEG_file", "organism", "ordered_query",
      "multi_query", "MAST_tested_background_genes",
      "MNG_higher_query_genes", "dCLN_higher_query_genes",
      "q_hurdle_threshold", "minimum_absolute_log2FC",
      "correction_method", "annotation_sources",
      "maximum_retained_term_size_exclusive",
      "passing_query_term_pairs", "gprofiler2_version",
      "gProfiler_metadata"
    ),
    value = c(
      "unordered directional DEG over-representation analysis",
      deg_csv,
      organism_use,
      FALSE,
      TRUE,
      length(universe_genes),
      length(query_list$MNG_higher),
      length(query_list$dCLN_higher),
      q_threshold,
      minimum_abs_log2fc,
      correction_method_use,
      paste(sources_use, collapse = "/"),
      maximum_retained_term_size,
      nrow(passing_result),
      as.character(utils::packageVersion("gprofiler2")),
      database_metadata
    ),
    stringsAsFactors = FALSE
  )
  write.csv(audit, audit_csv, row.names = FALSE)

  cat("Figure 2b passing terms by direction and source:\n")
  print(with(passing_result, addmargins(table(query, source))))
  cat("Saved Figure 2b pathway table:\n", output_csv, "\n", sep = "")
  cat("Saved Figure 2b audit:\n", audit_csv, "\n", sep = "")

  # Curated, nonredundant pathways for the displayed panel. Selection was
  # restricted to terms that passed the analysis above, then balanced across
  # B-cell biology, antigen handling, signaling, metabolism, translation, and
  # transcriptional/epigenetic regulation. No DNA-methylation or specific
  # cytokine-secretion term passed the stated filters, so none is implied here.
  selected_pathways <- data.frame(
    query = c(rep("MNG_higher", 10L), rep("dCLN_higher", 10L)),
    source = c(
      "GO:BP", "REAC", "GO:BP", "REAC", "REAC",
      "WP", "GO:BP", "REAC", "KEGG", "REAC",
      "GO:MF", "GO:BP", "GO:MF", "WP", "WP",
      "GO:BP", "WP", "GO:CC", "WP", "GO:BP"
    ),
    term_id = c(
      "GO:0006412", "REAC:R-MMU-611105", "GO:0006119",
      "REAC:R-MMU-1236977", "REAC:R-MMU-1236974", "WP:WP4466",
      "GO:0002474", "REAC:R-MMU-1236975", "KEGG:01230",
      "REAC:R-MMU-198933",
      "GO:0140993", "GO:0000902", "GO:0001067", "WP:WP258",
      "WP:WP88", "GO:0043488", "WP:WP493", "GO:0070161",
      "WP:WP274", "GO:0007264"
    ),
    display_term = c(
      "Translation",
      "Respiratory electron transport",
      "Oxidative phosphorylation",
      "Endosomal/vacuolar pathway",
      "ER-phagosome pathway",
      "Oxidative stress and redox",
      "MHC-I antigen processing",
      "Antigen cross-presentation",
      "Amino-acid biosynthesis",
      "Lymphoid-non-lymphoid immunoregulation",
      "Histone-modifying activity",
      "Cell morphogenesis",
      "Transcription-regulatory region binding",
      "TGF-beta receptor signaling",
      "Toll-like receptor signaling",
      "Regulation of mRNA stability",
      "MAPK signaling",
      "Anchoring junction",
      "B-cell receptor signaling",
      "Small-GTPase-mediated signaling"
    ),
    biological_theme = c(
      "translation", "oxidative metabolism", "oxidative metabolism",
      "antigen handling", "antigen handling", "redox metabolism",
      "antigen processing", "antigen processing", "anabolism",
      "immune interaction",
      "epigenetic regulation", "cell organization",
      "transcriptional regulation", "immune signaling",
      "innate immune signaling", "post-transcriptional regulation",
      "signal transduction", "adhesion/migration",
      "B-cell activation", "migration/signaling"
    ),
    stringsAsFactors = FALSE
  )

  selected_result <- merge(
    selected_pathways,
    passing_result,
    by = c("query", "source", "term_id"),
    all.x = TRUE,
    sort = FALSE
  )
  missing_selected <- selected_result$term_id[
    is.na(selected_result$adjusted_p_value_gSCS)
  ]
  if (length(missing_selected) > 0L) {
    stop(
      "Selected Figure 2b terms were not found among passing results: ",
      paste(missing_selected, collapse = ", "),
      call. = FALSE
    )
  }

  selected_result$enrichment_direction <- ifelse(
    selected_result$query == "MNG_higher", "MNG", "dCLN"
  )
  selected_result$minus_log2_adjusted_p <- -log2(
    pmax(selected_result$adjusted_p_value_gSCS, .Machine$double.xmin)
  )

  dcln_scores <- selected_result$minus_log2_adjusted_p[
    selected_result$enrichment_direction == "dCLN"
  ]
  dcln_terms <- selected_result$display_term[
    selected_result$enrichment_direction == "dCLN"
  ]
  mng_scores <- selected_result$minus_log2_adjusted_p[
    selected_result$enrichment_direction == "MNG"
  ]
  mng_terms <- selected_result$display_term[
    selected_result$enrichment_direction == "MNG"
  ]
  spacer_label <- " "
  term_levels <- c(
    dcln_terms[order(dcln_scores)],
    spacer_label,
    mng_terms[order(mng_scores)]
  )
  selected_result$display_term <- factor(
    selected_result$display_term,
    levels = term_levels
  )

  spacer_row <- selected_result[1L, , drop = FALSE]
  spacer_row$display_term <- factor(spacer_label, levels = term_levels)
  spacer_row$enrichment_direction <- NA_character_
  spacer_row$minus_log2_adjusted_p <- 0
  plot_data <- rbind(selected_result, spacer_row)

  mng_color <- "#FF6B6B"
  dcln_color <- "#5B9DF5"
  x_axis_max <- 5 * ceiling(
    max(plot_data$minus_log2_adjusted_p, na.rm = TRUE) / 5
  )
  pathway_plot <- ggplot(
    plot_data,
    aes(
      x = display_term,
      y = minus_log2_adjusted_p,
      fill = enrichment_direction
    )
  ) +
    geom_col(width = 0.82, na.rm = TRUE) +
    coord_flip(clip = "off") +
    scale_fill_manual(
      values = c(MNG = mng_color, dCLN = dcln_color),
      breaks = c("MNG", "dCLN"),
      drop = FALSE
    ) +
    scale_y_continuous(
      breaks = seq(0, x_axis_max, by = 5),
      limits = c(0, x_axis_max),
      expand = expansion(mult = c(0, 0))
    ) +
    labs(
      x = "Functional ontology term",
      y = "−log₂(g:SCS-adjusted P)",
      fill = "Enrichment\nin genes for"
    ) +
    theme_classic(base_size = 18) +
    theme(
      text = element_text(face = "bold", color = "black"),
      axis.title.x = element_text(
        size = 19,
        face = "bold",
        margin = margin(t = 9)
      ),
      axis.title.y = element_text(
        size = 19,
        face = "bold",
        margin = margin(r = 15)
      ),
      axis.text.x = element_text(size = 15, face = "bold", color = "black"),
      axis.text.y = element_text(size = 17, face = "bold", color = "black"),
      axis.line = element_line(linewidth = 1.4, color = "black"),
      axis.ticks = element_line(linewidth = 1.3, color = "black"),
      axis.ticks.length = unit(0.18, "cm"),
      legend.position = "none",
      plot.margin = margin(12, 18, 5, 18)
    )

  # Draw the source-style legend manually so each category is preceded by a
  # large true square and the two entries are stacked vertically.
  pathway_legend <- grid::grobTree(
    grid::textGrob(
      "Enrichment\nin genes for",
      x = unit(0.02, "npc"),
      y = unit(0.50, "npc"),
      hjust = 0,
      vjust = 0.5,
      gp = grid::gpar(fontsize = 17, fontface = "bold", lineheight = 0.95)
    ),
    grid::rectGrob(
      x = unit(1.58, "in"),
      y = unit(0.70, "npc"),
      width = unit(0.25, "in"),
      height = unit(0.25, "in"),
      gp = grid::gpar(fill = mng_color, col = NA)
    ),
    grid::textGrob(
      "MNG",
      x = unit(1.80, "in"),
      y = unit(0.70, "npc"),
      hjust = 0,
      gp = grid::gpar(fontsize = 17, fontface = "bold")
    ),
    grid::rectGrob(
      x = unit(1.58, "in"),
      y = unit(0.28, "npc"),
      width = unit(0.25, "in"),
      height = unit(0.25, "in"),
      gp = grid::gpar(fill = dcln_color, col = NA)
    ),
    grid::textGrob(
      "dCLN",
      x = unit(1.80, "in"),
      y = unit(0.28, "npc"),
      hjust = 0,
      gp = grid::gpar(fontsize = 17, fontface = "bold")
    )
  )
  figure_2b <- cowplot::ggdraw() +
    cowplot::draw_plot(
      pathway_plot,
      x = 0,
      y = 0.075,
      width = 1,
      height = 0.895
    ) +
    cowplot::draw_grob(
      pathway_legend,
      x = 0.235,
      y = 0.015,
      width = 0.40,
      height = 0.08
    )

  figure_tiff <- file.path(
    figure_output_dir,
    "figure_2b_MNG_vs_dCLN_pathway_enrichment.tiff"
  )
  selected_csv <- file.path(
    output_dir,
    "figure_2b_MNG_vs_dCLN_displayed_pathways.csv"
  )

  ggsave(
    figure_tiff,
    figure_2b,
    width = 10.2,
    height = 9.72,
    units = "in",
    dpi = 300,
    compression = "lzw",
    bg = "white"
  )
  write.csv(
    selected_result[order(
      selected_result$enrichment_direction,
      selected_result$minus_log2_adjusted_p,
      decreasing = TRUE
    ), ],
    selected_csv,
    row.names = FALSE
  )

  cat("Saved Figure 2b enrichment panel:\n", figure_tiff, "\n", sep = "")
  cat("Saved displayed-pathway audit:\n", selected_csv, "\n", sep = "")

  list(
    all_passing_terms = passing_result,
    displayed_terms = selected_result,
    plot = figure_2b
  )
})

}

################################################################################
# Figure 2b: preranked GSEA of MNG versus dCLN B-cell differential expression
################################################################################

# All genes tested by MAST are ranked by signed avg_log2FC (positive = MNG;
# negative = dCLN). fgseaMultilevel tests GO, Reactome, WikiPathways, KEGG and
# CORUM gene sets together. Pathways must contain 10-499 genes in the ranked
# universe, and BH FDR is controlled across all seven sources. Positive NES
# denotes MNG enrichment and negative NES denotes dCLN enrichment.

figure_2b_gsea_result <- local({
  set.seed(1234L)

  local_analysis_root <- Sys.getenv(
    "EAE_LOCAL_ANALYSIS_ROOT",
    unset = "/Users/shingheimok/Desktop/Phd_Study/research/Wu_lab/b_cells"
  )
  deg_dir <- file.path(
    local_analysis_root,
    "outs", "output", "DEG",
    "validated_bcells_MNG_vs_dCLN_MAST_msigdbr"
  )
  deg_csv <- file.path(
    deg_dir,
    "validated_bcells_MNG_vs_dCLN_MAST_all_genes.csv"
  )
  output_dir <- file.path(deg_dir, "fgsea_GSEA")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  figure_output_dir <- file.path(local_analysis_root, "figures", "figure_2")
  dir.create(figure_output_dir, recursive = TRUE, showWarnings = FALSE)

  all_result_csv <- file.path(
    output_dir,
    "figure_2b_MNG_vs_dCLN_GSEA_all_pathways.csv"
  )
  passing_result_csv <- file.path(
    output_dir,
    "figure_2b_MNG_vs_dCLN_GSEA_passing_FDR_lt_0.05_term_size_10_499.csv"
  )
  selected_csv <- file.path(
    output_dir,
    "figure_2b_MNG_vs_dCLN_GSEA_displayed_pathways.csv"
  )
  pathway_catalog_csv <- file.path(
    output_dir,
    "figure_2b_MNG_vs_dCLN_GSEA_pathway_catalog.csv"
  )
  audit_csv <- file.path(
    output_dir,
    "figure_2b_MNG_vs_dCLN_GSEA_analysis_audit.csv"
  )
  corum_cache <- file.path(output_dir, "CORUM_allComplexes_current.json")

  fdr_threshold <- 0.05
  minimum_tested_term_size <- 10L
  maximum_retained_term_size <- 500L
  sources_use <- c(
    "GO:BP", "GO:MF", "GO:CC", "REAC", "WP", "KEGG", "CORUM"
  )

  if (!file.exists(deg_csv)) {
    stop("Local Figure 2a DEG table not found: ", deg_csv, call. = FALSE)
  }
  deg <- read.csv(deg_csv, check.names = FALSE, stringsAsFactors = FALSE)
  required_columns <- c("gene", "avg_log2FC_MNG_vs_dCLN")
  missing_columns <- setdiff(required_columns, colnames(deg))
  if (length(missing_columns) > 0L) {
    stop(
      "DEG table is missing required columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  deg$gene <- as.character(deg$gene)
  deg$avg_log2FC_MNG_vs_dCLN <- as.numeric(deg$avg_log2FC_MNG_vs_dCLN)
  deg <- deg[
    !is.na(deg$gene) & nzchar(deg$gene) &
      is.finite(deg$avg_log2FC_MNG_vs_dCLN),
    c("gene", "avg_log2FC_MNG_vs_dCLN"),
    drop = FALSE
  ]
  deg <- deg[order(deg$gene, -abs(deg$avg_log2FC_MNG_vs_dCLN)), , drop = FALSE]
  deg <- deg[!duplicated(deg$gene), , drop = FALSE]
  ranks <- deg$avg_log2FC_MNG_vs_dCLN
  names(ranks) <- deg$gene
  ranks <- sort(ranks, decreasing = TRUE)

  clean_term_name <- function(x, prefixes = character()) {
    for (prefix in prefixes) x <- sub(paste0("^", prefix), "", x)
    gsub("_", " ", x, fixed = TRUE)
  }

  load_msigdb_source <- function(collection, subcollection, source, prefixes) {
    x <- msigdbr::msigdbr(
      db_species = "MM",
      species = "Mus musculus",
      collection = collection,
      subcollection = subcollection
    )
    if (nrow(x) == 0L) {
      stop("No gene sets returned for ", source, ".", call. = FALSE)
    }
    x |>
      transmute(
        source = .env$source,
        term_id = ifelse(
          is.na(.data$gs_exact_source) | .data$gs_exact_source == "",
          .data$gs_id,
          .data$gs_exact_source
        ),
        term_name = clean_term_name(.data$gs_name, prefixes),
        description = .data$gs_description,
        gene = .data$gene_symbol
      ) |>
      filter(!is.na(.data$gene), .data$gene != "") |>
      distinct()
  }

  cat("Loading Figure 2b mouse pathway libraries\n")
  go_bp <- load_msigdb_source("M5", "GO:BP", "GO:BP", c("GOBP_"))
  go_mf <- load_msigdb_source("M5", "GO:MF", "GO:MF", c("GOMF_"))
  go_cc <- load_msigdb_source("M5", "GO:CC", "GO:CC", c("GOCC_"))
  reactome <- load_msigdb_source(
    "M2", "CP:REACTOME", "REAC", c("REACTOME_")
  )
  wikipathways <- load_msigdb_source(
    "M2", "CP:WIKIPATHWAYS", "WP", c("WP_")
  )

  kegg_links <- KEGGREST::keggLink("pathway", "mmu")
  kegg_names <- KEGGREST::keggList("pathway", "mmu")
  kegg_link_df <- data.frame(
    entrez = sub("^mmu:", "", names(kegg_links)),
    pathway = sub("^path:", "", unname(kegg_links)),
    stringsAsFactors = FALSE
  )
  kegg_symbols <- AnnotationDbi::mapIds(
    org.Mm.eg.db,
    keys = unique(kegg_link_df$entrez),
    keytype = "ENTREZID",
    column = "SYMBOL",
    multiVals = "first"
  )
  kegg_link_df$gene <- unname(kegg_symbols[kegg_link_df$entrez])
  kegg_name_df <- data.frame(
    pathway = names(kegg_names),
    term_name = sub(
      " - Mus musculus \\(house mouse\\)$",
      "",
      unname(kegg_names)
    ),
    stringsAsFactors = FALSE
  )
  kegg <- kegg_link_df |>
    left_join(kegg_name_df, by = "pathway") |>
    transmute(
      source = "KEGG",
      term_id = paste0("KEGG:", sub("^mmu", "", .data$pathway)),
      term_name = .data$term_name,
      description = .data$term_name,
      gene = .data$gene
    ) |>
    filter(!is.na(.data$gene), .data$gene != "") |>
    distinct()

  corum_url <- paste0(
    "https://mips.helmholtz-muenchen.de/fastapi-corum/",
    "public/file/download_current_file?file_id=complete&file_format=json"
  )
  if (!file.exists(corum_cache) || file.info(corum_cache)$size < 10000L) {
    download.file(corum_url, corum_cache, mode = "wb", quiet = TRUE)
  }
  corum_json <- jsonlite::fromJSON(corum_cache, simplifyVector = FALSE)
  corum_json <- Filter(
    function(x) !is.null(x$organism) && x$organism %in% c("Mouse", "Human"),
    corum_json
  )
  corum_raw <- bind_rows(lapply(corum_json, function(complex) {
    if (length(complex$subunits) == 0L) return(NULL)
    bind_rows(lapply(complex$subunits, function(subunit) {
      swissprot <- subunit$swissprot
      if (is.null(swissprot)) return(NULL)
      data.frame(
        organism = complex$organism,
        complex_id = as.character(complex$complex_id),
        complex_name = as.character(complex$complex_name),
        gene_name = if (is.null(swissprot$gene_name)) {
          NA_character_
        } else {
          as.character(swissprot$gene_name)
        },
        entrez_id = if (is.null(swissprot$entrez_id)) {
          NA_character_
        } else {
          as.character(swissprot$entrez_id)
        },
        stringsAsFactors = FALSE
      )
    }))
  }))
  mouse_corum <- corum_raw |>
    filter(.data$organism == "Mouse", !is.na(.data$entrez_id))
  mouse_corum_symbols <- AnnotationDbi::mapIds(
    org.Mm.eg.db,
    keys = unique(mouse_corum$entrez_id),
    keytype = "ENTREZID",
    column = "SYMBOL",
    multiVals = "first"
  )
  mouse_corum$gene <- unname(mouse_corum_symbols[mouse_corum$entrez_id])
  human_corum <- corum_raw |>
    filter(.data$organism == "Human", !is.na(.data$gene_name))
  human_map <- babelgene::orthologs(
    unique(human_corum$gene_name),
    species = "mouse",
    human = TRUE,
    min_support = 3,
    top = TRUE
  ) |>
    transmute(gene_name = .data$human_symbol, gene = .data$symbol) |>
    distinct()
  human_corum <- human_corum |>
    left_join(human_map, by = "gene_name", relationship = "many-to-many")
  corum <- bind_rows(mouse_corum, human_corum) |>
    transmute(
      source = "CORUM",
      term_id = paste0("CORUM:", .data$complex_id),
      term_name = .data$complex_name,
      description = .data$complex_name,
      gene = .data$gene
    ) |>
    filter(!is.na(.data$gene), .data$gene != "") |>
    distinct()

  catalog_full <- bind_rows(
    go_bp, go_mf, go_cc, reactome, wikipathways, kegg, corum
  ) |>
    mutate(
      source = factor(.data$source, levels = sources_use),
      pathway_key = paste(.data$source, .data$term_id, sep = "::")
    ) |>
    distinct(
      .data$pathway_key, .data$source, .data$term_id,
      .data$term_name, .data$description, .data$gene
    )
  term_metadata <- catalog_full |>
    group_by(
      .data$pathway_key, .data$source, .data$term_id,
      .data$term_name, .data$description
    ) |>
    summarise(database_term_size = n_distinct(.data$gene), .groups = "drop") |>
    filter(
      .data$database_term_size >= minimum_tested_term_size,
      .data$database_term_size < maximum_retained_term_size
    )
  catalog_ranked <- catalog_full |>
    semi_join(term_metadata, by = "pathway_key") |>
    filter(.data$gene %in% names(ranks)) |>
    distinct(.data$pathway_key, .data$gene)
  pathways <- split(catalog_ranked$gene, catalog_ranked$pathway_key)
  pathways <- lapply(pathways, unique)
  pathways <- pathways[lengths(pathways) >= minimum_tested_term_size]
  pathway_metadata <- term_metadata |>
    filter(.data$pathway_key %in% names(pathways)) |>
    mutate(ranked_universe_overlap_size = lengths(pathways[.data$pathway_key])) |>
    mutate(source = as.character(.data$source))
  write.csv(pathway_metadata, pathway_catalog_csv, row.names = FALSE)

  cat(
    "Running Figure 2b preranked GSEA on ", length(ranks),
    " MAST-tested genes and ", length(pathways), " pathways\n",
    sep = ""
  )
  set.seed(1234L)
  gsea <- fgsea::fgseaMultilevel(
    pathways = pathways,
    stats = ranks,
    minSize = minimum_tested_term_size,
    maxSize = Inf,
    eps = 1e-10,
    scoreType = "std",
    nPermSimple = 10000,
    nproc = 1
  ) |>
    as.data.frame()
  gsea$leading_edge_genes <- vapply(
    gsea$leadingEdge,
    function(x) paste(x, collapse = ";"),
    character(1)
  )
  gsea$leading_edge_n <- lengths(gsea$leadingEdge)
  gsea$leadingEdge <- NULL
  all_results <- gsea |>
    dplyr::rename(
      pathway_key = .data$pathway,
      p_value = .data$pval,
      fgsea_padj = .data$padj,
      enrichment_score = .data$ES,
      normalized_enrichment_score = .data$NES,
      ranked_universe_overlap_size_fgsea = .data$size
    ) |>
    left_join(pathway_metadata, by = "pathway_key") |>
    mutate(
      fdr_bh_all_sources = p.adjust(.data$p_value, method = "BH"),
      enrichment_direction = ifelse(
        .data$normalized_enrichment_score > 0,
        "MNG",
        "dCLN"
      )
    ) |>
    dplyr::select(
      .data$enrichment_direction, .data$source, .data$term_id,
      .data$term_name, .data$description, .data$database_term_size,
      .data$ranked_universe_overlap_size,
      .data$enrichment_score, .data$normalized_enrichment_score,
      .data$p_value, .data$fdr_bh_all_sources, .data$fgsea_padj,
      .data$log2err, .data$leading_edge_n, .data$leading_edge_genes,
      .data$pathway_key
    ) |>
    arrange(
      .data$fdr_bh_all_sources,
      desc(abs(.data$normalized_enrichment_score)),
      factor(.data$source, levels = sources_use),
      .data$term_name
    )
  write.csv(all_results, all_result_csv, row.names = FALSE)

  passing_result <- all_results |>
    filter(
      is.finite(.data$fdr_bh_all_sources),
      .data$fdr_bh_all_sources < fdr_threshold,
      .data$database_term_size < maximum_retained_term_size,
      .data$ranked_universe_overlap_size >= minimum_tested_term_size
    )
  if (nrow(passing_result) == 0L) {
    stop("No Figure 2b GSEA pathway passed FDR < 0.05.", call. = FALSE)
  }
  write.csv(passing_result, passing_result_csv, row.names = FALSE)

  # Curated pathway themes preserve the established Figure 2b design. For each
  # theme, the most significant passing GSEA term is selected; if a theme has
  # no passing term, the remaining slot is filled by the best nonduplicated
  # pathway in that direction. The exact displayed terms are saved to CSV.
  mng_preferences <- data.frame(
    display_term = c(
      "Translation", "Respiratory electron transport",
      "Oxidative phosphorylation", "Endosomal/vacuolar pathway",
      "ER-phagosome pathway", "Oxidative stress and redox",
      "MHC-I antigen processing", "Antigen cross-presentation",
      "Amino-acid biosynthesis", "Lymphoid-non-lymphoid immunoregulation"
    ),
    pattern = c(
      "^TRANSLATION$", "RESPIRATORY ELECTRON TRANSPORT|ELECTRON TRANSPORT CHAIN",
      "OXIDATIVE PHOSPHORYLATION", "ENDOSOM|VACUOLAR",
      "ER PHAGOSOME|ENDOPLASMIC RETICULUM.*PHAGOSOME",
      "OXIDATIVE STRESS|REDOX", "ANTIGEN PROCESSING.*MHC CLASS I",
      "CROSS PRESENTATION", "AMINO ACID BIOSYNTH",
      "LYMPHOID.*NON LYMPHOID|IMMUNOREGUL"
    ),
    stringsAsFactors = FALSE
  )
  dcln_preferences <- data.frame(
    display_term = c(
      "Histone-modifying activity", "Cell morphogenesis",
      "TGF-beta receptor signaling", "Transcription-regulatory region binding",
      "Toll-like receptor signaling", "Regulation of mRNA stability",
      "MAPK signaling", "Anchoring junction",
      "B-cell receptor signaling", "Small-GTPase-mediated signaling"
    ),
    pattern = c(
      "HISTONE MODIFY|HISTONE MODIFICATION", "CELL MORPHOGENESIS",
      "TGF.*RECEPTOR SIGNAL", "TRANSCRIPTION.*REGULATORY REGION BINDING",
      "TOLL LIKE RECEPTOR", "MRNA STABILITY", "MAPK",
      "ANCHORING JUNCTION", "B CELL RECEPTOR SIGNAL", "SMALL GTPASE"
    ),
    stringsAsFactors = FALSE
  )

  select_display_terms <- function(direction, preferences, number_to_keep = 10L) {
    pool <- passing_result |>
      filter(.data$enrichment_direction == direction) |>
      arrange(
        .data$fdr_bh_all_sources,
        desc(abs(.data$normalized_enrichment_score))
      )
    selected_rows <- list()
    used_keys <- character()
    for (preference_index in seq_len(nrow(preferences))) {
      candidates <- pool[
        !pool$pathway_key %in% used_keys &
          grepl(
            preferences$pattern[[preference_index]],
            pool$term_name,
            ignore.case = TRUE
          ),
        ,
        drop = FALSE
      ]
      if (nrow(candidates) > 0L) {
        chosen <- candidates[1L, , drop = FALSE]
        chosen$display_term <- preferences$display_term[[preference_index]]
        selected_rows[[length(selected_rows) + 1L]] <- chosen
        used_keys <- c(used_keys, chosen$pathway_key)
      }
    }
    if (length(selected_rows) < number_to_keep) {
      fillers <- pool[!pool$pathway_key %in% used_keys, , drop = FALSE]
      fillers <- head(fillers, number_to_keep - length(selected_rows))
      if (nrow(fillers) > 0L) {
        fillers$display_term <- tools::toTitleCase(
          tolower(fillers$term_name)
        )
        selected_rows[[length(selected_rows) + 1L]] <- fillers
      }
    }
    selected <- bind_rows(selected_rows)
    if (nrow(selected) < number_to_keep) {
      stop(
        "Fewer than ", number_to_keep, " passing GSEA pathways for ",
        direction, ".",
        call. = FALSE
      )
    }
    head(selected, number_to_keep)
  }

  # Lock the final panel to nonredundant, biologically interpretable pathways
  # that actually pass this GSEA. The previous ORA-only BCR/TLR/MAPK labels are
  # not carried forward because those terms do not pass the GSEA FDR cutoff.
  selected_pathways <- data.frame(
    enrichment_direction = c(rep("MNG", 10L), rep("dCLN", 6L)),
    source = c(
      "REAC", "REAC", "KEGG", "REAC", "REAC",
      "GO:BP", "GO:BP", "KEGG", "GO:MF", "CORUM",
      "GO:BP", "GO:CC", "REAC", "GO:BP", "GO:BP", "REAC"
    ),
    term_id = c(
      "R-MMU-72613", "R-MMU-611105", "KEGG:00190",
      "R-MMU-1236974", "R-MMU-9912633", "GO:0019883",
      "GO:0035456", "KEGG:01230", "GO:0016209", "CORUM:39",
      "GO:0007015", "GO:0031252", "R-MMU-9012999",
      "GO:0010810", "GO:0040029", "R-MMU-9958863"
    ),
    display_term = c(
      "Translation", "Respiratory electron transport",
      "Oxidative phosphorylation", "ER-phagosome pathway",
      "Proteasomal antigen processing", "MHC-I antigen processing",
      "Response to interferon-beta", "Amino-acid biosynthesis",
      "Antioxidant activity", "Immunoproteasome",
      "Actin-cytoskeleton remodeling", "Cell polarity and migration",
      "CDC42/RAC1 signaling", "Cell adhesion",
      "Epigenetic regulation", "Amino-acid transport"
    ),
    stringsAsFactors = FALSE
  )
  selected_result <- passing_result |>
    inner_join(
      selected_pathways,
      by = c("enrichment_direction", "source", "term_id"),
      relationship = "one-to-one"
    )
  if (nrow(selected_result) != nrow(selected_pathways)) {
    missing_selected <- anti_join(
      selected_pathways,
      passing_result,
      by = c("enrichment_direction", "source", "term_id")
    )
    stop(
      "Selected Figure 2b GSEA pathways are missing: ",
      paste(missing_selected$term_id, collapse = ", "),
      call. = FALSE
    )
  }
  selected_result <- selected_result |>
    mutate(
      minus_log2_fdr = -log2(
        pmax(.data$fdr_bh_all_sources, .Machine$double.xmin)
      )
    )

  dcln_scores <- selected_result$minus_log2_fdr[
    selected_result$enrichment_direction == "dCLN"
  ]
  dcln_terms <- selected_result$display_term[
    selected_result$enrichment_direction == "dCLN"
  ]
  mng_scores <- selected_result$minus_log2_fdr[
    selected_result$enrichment_direction == "MNG"
  ]
  mng_terms <- selected_result$display_term[
    selected_result$enrichment_direction == "MNG"
  ]
  spacer_label <- " "
  term_levels <- c(
    dcln_terms[order(dcln_scores)],
    spacer_label,
    mng_terms[order(mng_scores)]
  )
  selected_result$display_term <- factor(
    selected_result$display_term,
    levels = term_levels
  )
  spacer_row <- selected_result[1L, , drop = FALSE]
  spacer_row$display_term <- factor(spacer_label, levels = term_levels)
  spacer_row$enrichment_direction <- NA_character_
  spacer_row$minus_log2_fdr <- 0
  plot_data <- bind_rows(selected_result, spacer_row)

  mng_color <- "#FF6B6B"
  dcln_color <- "#5B9DF5"
  x_axis_max <- 5 * ceiling(max(plot_data$minus_log2_fdr, na.rm = TRUE) / 5)
  pathway_plot <- ggplot(
    plot_data,
    aes(
      x = display_term,
      y = minus_log2_fdr,
      fill = enrichment_direction
    )
  ) +
    geom_col(width = 0.82, na.rm = TRUE) +
    coord_flip(clip = "off") +
    scale_fill_manual(
      values = c(MNG = mng_color, dCLN = dcln_color),
      breaks = c("MNG", "dCLN"),
      drop = FALSE
    ) +
    scale_y_continuous(
      breaks = seq(0, x_axis_max, by = 5),
      limits = c(0, x_axis_max),
      expand = expansion(mult = c(0, 0))
    ) +
    labs(
      x = "Functional ontology term",
      y = "−log₂(BH-adjusted P)",
      fill = "Enrichment\nin genes for"
    ) +
    theme_classic(base_size = 18) +
    theme(
      text = element_text(face = "bold", color = "black"),
      axis.title.x = element_text(
        size = 19, face = "bold", margin = margin(t = 9)
      ),
      axis.title.y = element_text(
        size = 19, face = "bold", margin = margin(r = 15)
      ),
      axis.text.x = element_text(size = 15, face = "bold", color = "black"),
      axis.text.y = element_text(size = 17, face = "bold", color = "black"),
      axis.line = element_line(linewidth = 1.4, color = "black"),
      axis.ticks = element_line(linewidth = 1.3, color = "black"),
      axis.ticks.length = unit(0.18, "cm"),
      legend.position = "none",
      plot.margin = margin(12, 18, 5, 18)
    )

  pathway_legend <- grid::grobTree(
    grid::textGrob(
      "Enrichment\nin genes for",
      x = unit(0.02, "npc"), y = unit(0.50, "npc"),
      hjust = 0, vjust = 0.5,
      gp = grid::gpar(fontsize = 17, fontface = "bold", lineheight = 0.95)
    ),
    grid::rectGrob(
      x = unit(1.58, "in"), y = unit(0.70, "npc"),
      width = unit(0.25, "in"), height = unit(0.25, "in"),
      gp = grid::gpar(fill = mng_color, col = NA)
    ),
    grid::textGrob(
      "MNG", x = unit(1.80, "in"), y = unit(0.70, "npc"),
      hjust = 0, gp = grid::gpar(fontsize = 17, fontface = "bold")
    ),
    grid::rectGrob(
      x = unit(1.58, "in"), y = unit(0.28, "npc"),
      width = unit(0.25, "in"), height = unit(0.25, "in"),
      gp = grid::gpar(fill = dcln_color, col = NA)
    ),
    grid::textGrob(
      "dCLN", x = unit(1.80, "in"), y = unit(0.28, "npc"),
      hjust = 0, gp = grid::gpar(fontsize = 17, fontface = "bold")
    )
  )
  figure_2b <- cowplot::ggdraw() +
    cowplot::draw_plot(
      pathway_plot, x = 0, y = 0.075, width = 1, height = 0.895
    ) +
    cowplot::draw_grob(
      pathway_legend, x = 0.235, y = 0.015, width = 0.40, height = 0.08
    )

  figure_tiff <- file.path(
    figure_output_dir,
    "figure_2b_MNG_vs_dCLN_pathway_enrichment.tiff"
  )
  ggsave(
    figure_tiff,
    figure_2b,
    width = 10.2,
    height = 9.72,
    units = "in",
    dpi = 300,
    compression = "lzw",
    bg = "white"
  )
  write.csv(
    selected_result |>
      arrange(
        .data$enrichment_direction,
        desc(.data$minus_log2_fdr)
      ),
    selected_csv,
    row.names = FALSE
  )

  corum_release <- tryCatch(
    jsonlite::fromJSON(
      "https://mips.helmholtz-muenchen.de/fastapi-corum/public/releases/current"
    ),
    error = function(e) NULL
  )
  audit <- data.frame(
    parameter = c(
      "analysis", "input_DEG_file", "ranking_metric",
      "ranking_direction", "ranked_gene_count", "GSEA_engine", "seed",
      "annotation_sources", "minimum_tested_term_size",
      "maximum_source_term_size_exclusive", "multiple_testing",
      "FDR_threshold", "tested_pathway_count", "passing_pathway_count",
      "MNG_passing_pathways", "dCLN_passing_pathways",
      "MSigDB_version", "CORUM_release", "ORA_performed"
    ),
    value = c(
      "preranked GSEA of MNG versus dCLN MAST differential expression",
      deg_csv,
      "avg_log2FC_MNG_vs_dCLN",
      "positive = MNG; negative = dCLN",
      length(ranks),
      "fgseaMultilevel",
      1234,
      paste(sources_use, collapse = "/"),
      minimum_tested_term_size,
      maximum_retained_term_size,
      "Benjamini-Hochberg across all seven sources",
      fdr_threshold,
      length(pathways),
      nrow(passing_result),
      sum(passing_result$enrichment_direction == "MNG"),
      sum(passing_result$enrichment_direction == "dCLN"),
      as.character(packageVersion("msigdbr")),
      if (is.null(corum_release)) {
        "current release endpoint unavailable"
      } else {
        paste0(corum_release$version, " (", corum_release$date, ")")
      },
      FALSE
    ),
    stringsAsFactors = FALSE
  )
  write.csv(audit, audit_csv, row.names = FALSE)

  cat("Figure 2b GSEA passing pathways by direction and source:\n")
  print(with(passing_result, addmargins(table(enrichment_direction, source))))
  cat("Saved Figure 2b GSEA panel:\n", figure_tiff, "\n", sep = "")
  cat("Saved Figure 2b GSEA pathways:\n", passing_result_csv, "\n", sep = "")

  list(
    all_pathways = all_results,
    passing_pathways = passing_result,
    displayed_pathways = selected_result,
    plot = figure_2b
  )
})

################################################################################
# Figure 2c: Kolz Th17 CNS versus LN B-cell MAST DEG and volcano plot
################################################################################

# This sensitivity analysis applies the Figure 2a pipeline to the Kolz Th17
# adoptive-transfer data. Th17 B cells are defined as orig.ident M1/M2/M3;
# only CNS and LN are compared, and spleen is excluded. Positive log2 fold
# change means higher expression in CNS. MAST includes mouse, nCount_RNA, and
# percent.mito as covariates. Directional g:Profiler ORA uses all genes tested
# by MAST as the custom background.

figure_2c_kolz_th17 <- local({
  seed_use <- 1234L
  set.seed(seed_use)

  local_analysis_root <- Sys.getenv(
    "EAE_LOCAL_ANALYSIS_ROOT",
    unset = "/Users/shingheimok/Desktop/Phd_Study/research/Wu_lab/b_cells"
  )
  object_candidates <- c(
    Sys.getenv("KOLZ_TH1_TH17_RDS", unset = ""),
    file.path(
      local_analysis_root,
      "external", "GSE279684", "GSE279684_2024.04.18_Th1-Th17.rds"
    )
  )
  object_candidates <- object_candidates[nzchar(object_candidates)]
  object_path <- object_candidates[file.exists(object_candidates)][1]
  if (is.na(object_path)) {
    stop(
      "Kolz Th1/Th17 RDS was not found. Set KOLZ_TH1_TH17_RDS or place the ",
      "uncompressed RDS under external/GSE279684/.",
      call. = FALSE
    )
  }

  output_dir <- file.path(
    local_analysis_root,
    "outs", "output", "DEG", "Kolz_Th17_CNS_vs_LN_MAST_gprofiler2"
  )
  figure_output_dir <- file.path(
    local_analysis_root, "figures", "figure_2"
  )
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figure_output_dir, recursive = TRUE, showWarnings = FALSE)

  deg_csv <- file.path(output_dir, "Kolz_Th17_CNS_vs_LN_MAST_all_genes.csv")
  pathway_csv <- file.path(
    output_dir,
    "Kolz_Th17_CNS_vs_LN_gprofiler2_all_passing_pathways.csv"
  )
  audit_csv <- file.path(output_dir, "Kolz_Th17_CNS_vs_LN_analysis_audit.csv")

  q_threshold <- 0.05
  minimum_abs_log2fc <- 0.2
  minimum_detection_fraction <- 0.01
  maximum_retained_term_size <- 500L
  y_display_cap <- as.numeric(Sys.getenv("KOLZ_VOLCANO_Y_LIMIT", unset = "200"))
  sources_use <- c("GO:BP", "GO:MF", "GO:CC", "REAC", "WP", "KEGG", "CORUM")
  force_rerun <- tolower(Sys.getenv("RERUN_KOLZ_TH17_MAST", unset = "false")) %in%
    c("1", "true", "yes")

  if (force_rerun || !file.exists(deg_csv) || !file.exists(pathway_csv)) {
    if (!requireNamespace("MAST", quietly = TRUE)) {
      stop("The MAST package is required to rerun Figure 2c.", call. = FALSE)
    }

    cat("Loading Kolz object:\n", object_path, "\n", sep = "")
    kolz <- readRDS(object_path)
    DefaultAssay(kolz) <- "RNA"
    kolz_md <- kolz@meta.data
    required_md <- c("orig.ident", "compartment", "nCount_RNA", "percent.mito")
    missing_md <- setdiff(required_md, colnames(kolz_md))
    if (length(missing_md) > 0L) {
      stop(
        "Kolz metadata is missing: ", paste(missing_md, collapse = ", "),
        call. = FALSE
      )
    }

    selected_cells <- rownames(kolz_md)[
      as.character(kolz_md$orig.ident) %in% c("M1", "M2", "M3") &
        as.character(kolz_md$compartment) %in% c("CNS", "LN")
    ]
    th17 <- subset(kolz, cells = selected_cells)
    th17$MAST_compartment <- as.character(th17$compartment)
    th17$MAST_mouse <- factor(
      as.character(th17$orig.ident), levels = c("M1", "M2", "M3")
    )
    th17 <- tryCatch(
      SeuratObject::JoinLayers(th17, assay = "RNA"),
      error = function(e) th17
    )
    th17 <- NormalizeData(
      th17,
      assay = "RNA",
      normalization.method = "LogNormalize",
      scale.factor = 10000,
      verbose = FALSE
    )

    cell_audit <- as.data.frame.matrix(table(
      mouse = th17$MAST_mouse,
      compartment = th17$MAST_compartment
    ))
    cell_audit$mouse <- rownames(cell_audit)
    rownames(cell_audit) <- NULL
    write.csv(
      cell_audit,
      file.path(output_dir, "selected_cell_audit.csv"),
      row.names = FALSE
    )

    Idents(th17) <- factor(th17$MAST_compartment, levels = c("LN", "CNS"))
    latent_vars <- c("nCount_RNA", "percent.mito", "MAST_mouse")
    set.seed(seed_use)
    mast_result <- FindMarkers(
      th17,
      ident.1 = "CNS",
      ident.2 = "LN",
      assay = "RNA",
      test.use = "MAST",
      logfc.threshold = 0,
      min.pct = minimum_detection_fraction,
      only.pos = FALSE,
      latent.vars = latent_vars,
      verbose = TRUE
    )
    if (nrow(mast_result) == 0L) stop("MAST returned no genes.", call. = FALSE)

    mast_result$gene <- rownames(mast_result)
    fc_col <- intersect(c("avg_log2FC", "avg_logFC"), colnames(mast_result))[1]
    if (is.na(fc_col) || !"p_val" %in% colnames(mast_result)) {
      stop("Unexpected MAST output columns.", call. = FALSE)
    }
    if ("p_val_adj" %in% colnames(mast_result)) {
      mast_result$p_val_adj_Seurat <- mast_result$p_val_adj
    }
    mast_result$avg_log2FC_CNS_vs_LN <- as.numeric(mast_result[[fc_col]])
    mast_result$q_hurdle <- p.adjust(mast_result$p_val, method = "BH")
    mast_result$direction <- dplyr::case_when(
      mast_result$avg_log2FC_CNS_vs_LN > 0 ~ "CNS_higher",
      mast_result$avg_log2FC_CNS_vs_LN < 0 ~ "LN_higher",
      TRUE ~ "no_change"
    )
    mast_result$significant_fdr_0.05 <-
      !is.na(mast_result$q_hurdle) & mast_result$q_hurdle < q_threshold
    mast_result <- mast_result |>
      dplyr::select(
        gene, avg_log2FC_CNS_vs_LN, p_val, q_hurdle,
        direction, significant_fdr_0.05, dplyr::everything()
      ) |>
      dplyr::arrange(q_hurdle, dplyr::desc(abs(avg_log2FC_CNS_vs_LN)))
    write.csv(mast_result, deg_csv, row.names = FALSE)

    universe_genes <- unique(mast_result$gene[
      !is.na(mast_result$gene) & nzchar(mast_result$gene)
    ])
    query_list <- list(
      CNS_higher = sort(unique(mast_result$gene[
        mast_result$q_hurdle < q_threshold &
          mast_result$avg_log2FC_CNS_vs_LN > minimum_abs_log2fc
      ])),
      LN_higher = sort(unique(mast_result$gene[
        mast_result$q_hurdle < q_threshold &
          mast_result$avg_log2FC_CNS_vs_LN < -minimum_abs_log2fc
      ]))
    )
    if (any(lengths(query_list) == 0L)) {
      stop("A Kolz directional DEG query is empty.", call. = FALSE)
    }

    gost_result <- gprofiler2::gost(
      query = query_list,
      organism = "mmusculus",
      ordered_query = FALSE,
      multi_query = TRUE,
      significant = TRUE,
      exclude_iea = FALSE,
      measure_underrepresentation = FALSE,
      evcodes = FALSE,
      user_threshold = q_threshold,
      correction_method = "g_SCS",
      domain_scope = "custom",
      custom_bg = universe_genes,
      sources = sources_use
    )
    if (is.null(gost_result) || is.null(gost_result$result) ||
        nrow(gost_result$result) == 0L) {
      stop("g:Profiler returned no significant terms.", call. = FALSE)
    }

    multi_result <- gost_result$result
    query_names <- names(query_list)
    long_rows <- lapply(seq_len(nrow(multi_result)), function(i) {
      row <- multi_result[i, , drop = FALSE]
      p_values <- as.numeric(unlist(row$p_values[[1]], use.names = FALSE))
      significant <- as.logical(unlist(row$significant[[1]], use.names = FALSE))
      query_sizes <- as.numeric(unlist(row$query_sizes[[1]], use.names = FALSE))
      intersection_sizes <- as.numeric(
        unlist(row$intersection_sizes[[1]], use.names = FALSE)
      )
      data.frame(
        query = query_names,
        direction = ifelse(
          query_names == "CNS_higher", "higher in CNS", "higher in LN"
        ),
        source = as.character(row$source),
        term_id = as.character(row$term_id),
        term_name = as.character(row$term_name),
        adjusted_p_value_gSCS = p_values,
        significant_gSCS_0.05 = significant,
        term_size = as.integer(row$term_size),
        query_size = query_sizes,
        intersection_size = intersection_sizes,
        effective_domain_size = as.integer(row$effective_domain_size),
        source_order = as.integer(row$source_order),
        stringsAsFactors = FALSE
      )
    })
    pathways <- dplyr::bind_rows(long_rows) |>
      dplyr::filter(
        significant_gSCS_0.05 %in% TRUE,
        is.finite(adjusted_p_value_gSCS),
        adjusted_p_value_gSCS < q_threshold,
        term_size < maximum_retained_term_size
      ) |>
      dplyr::arrange(
        match(query, query_names), adjusted_p_value_gSCS,
        match(source, sources_use), term_name
      )
    write.csv(pathways, pathway_csv, row.names = FALSE)

    audit <- data.frame(
      parameter = c(
        "comparison", "CNS_cells", "LN_cells", "Spleen_cells_excluded",
        "mouse_block", "MAST_latent_variables", "MAST_tested_genes",
        "CNS_higher_DEGs", "LN_higher_DEGs", "q_hurdle_threshold",
        "minimum_absolute_log2FC", "gprofiler_sources",
        "maximum_retained_term_size_exclusive", "passing_query_term_pairs"
      ),
      value = c(
        "Kolz Th17 B cells: CNS versus LN; Spleen excluded",
        sum(th17$MAST_compartment == "CNS"),
        sum(th17$MAST_compartment == "LN"),
        sum(
          kolz_md$orig.ident %in% c("M1", "M2", "M3") &
            kolz_md$compartment == "Spleen"
        ),
        "orig.ident M1/M2/M3", paste(latent_vars, collapse = "/"),
        length(universe_genes), length(query_list$CNS_higher),
        length(query_list$LN_higher), q_threshold, minimum_abs_log2fc,
        paste(sources_use, collapse = "/"), maximum_retained_term_size,
        nrow(pathways)
      ),
      stringsAsFactors = FALSE
    )
    write.csv(audit, audit_csv, row.names = FALSE)
    rm(kolz, th17)
    invisible(gc())
  }

  deg <- read.csv(deg_csv, check.names = FALSE, stringsAsFactors = FALSE)
  required_deg <- c("gene", "avg_log2FC_CNS_vs_LN", "q_hurdle")
  missing_deg <- setdiff(required_deg, colnames(deg))
  if (length(missing_deg) > 0L) {
    stop("Kolz DEG table is missing: ", paste(missing_deg, collapse = ", "))
  }
  deg$log2FC <- as.numeric(deg$avg_log2FC_CNS_vs_LN)
  deg$q_value <- as.numeric(deg$q_hurdle)
  deg <- deg[
    is.finite(deg$log2FC) & !is.na(deg$q_value) & deg$q_value >= 0,
    ,
    drop = FALSE
  ]
  smallest_nonzero_q <- suppressWarnings(min(deg$q_value[deg$q_value > 0]))
  if (!is.finite(smallest_nonzero_q)) smallest_nonzero_q <- 1e-300
  q_floor <- max(smallest_nonzero_q / 10, 1e-300)
  deg$q_for_plot <- pmax(deg$q_value, q_floor)
  deg$neg_log10_q_uncapped <- -log10(deg$q_for_plot)
  deg$neg_log10_q <- pmin(deg$neg_log10_q_uncapped, y_display_cap)
  deg$volcano_group <- dplyr::case_when(
    abs(deg$log2FC) > minimum_abs_log2fc & deg$q_value < q_threshold &
      deg$log2FC < 0 ~ "LN higher",
    abs(deg$log2FC) > minimum_abs_log2fc & deg$q_value < q_threshold &
      deg$log2FC > 0 ~ "CNS higher",
    TRUE ~ "Not significant"
  )

  label_panel <- data.frame(
    gene = c(
      "Nr4a1", "Ccr5", "Fam46c", "H2-Q7", "Sdhaf1",
      "Tpi1", "Prdx6", "Hspa5", "Atf4", "Il2ra",
      "S1pr4", "Clec2i", "Cr2", "Fcer2a", "Notch2",
      "Ralgps2", "Fgd2", "Sell", "Bach2", "Cxcr5"
    ),
    expected_group = c(rep("CNS higher", 10L), rep("LN higher", 10L)),
    annotation_theme = c(
      "activation", "chemokine signaling", "antibody secretion",
      "MHC-I antigen processing", "oxidative metabolism",
      "glycolysis", "redox", "ER stress", "amino-acid/stress response",
      "cytokine receptor/activation",
      "lymphocyte trafficking", "immune signaling", "mature B-cell identity",
      "mature B-cell signaling", "Notch signaling", "small-GTPase signaling",
      "cytoskeleton organization", "lymphocyte trafficking",
      "B-cell differentiation", "lymphoid-follicle trafficking"
    ),
    stringsAsFactors = FALSE
  )
  label_table <- dplyr::left_join(label_panel, deg, by = "gene")
  invalid_labels <- label_table$gene[
    is.na(label_table$volcano_group) |
      label_table$volcano_group != label_table$expected_group
  ]
  if (length(invalid_labels) > 0L) {
    stop(
      "Kolz volcano labels do not pass the stated thresholds/direction: ",
      paste(invalid_labels, collapse = ", "),
      call. = FALSE
    )
  }

  cns_color <- "#FF6B6B"
  ln_color <- "#5B9DF5"
  neutral_color <- "#C8C8C8"
  group_levels <- c("LN higher", "Not significant", "CNS higher")
  deg$volcano_group <- factor(deg$volcano_group, levels = group_levels)
  label_table$volcano_group <- factor(
    label_table$volcano_group, levels = group_levels
  )
  x_quantile <- as.numeric(
    stats::quantile(abs(deg$log2FC), 0.995, na.rm = TRUE, names = FALSE)
  )
  x_limit <- max(
    2,
    ceiling(2 * max(x_quantile, abs(label_table$log2FC), na.rm = TRUE)) / 2
  )

  arrow_grob <- grid::grobTree(
    grid::polygonGrob(
      x = unit(c(0.49, 0.17, 0.17, 0.01, 0.17, 0.17, 0.49), "npc"),
      y = unit(c(0.32, 0.32, 0.08, 0.50, 0.92, 0.68, 0.68), "npc"),
      gp = grid::gpar(fill = ln_color, col = NA)
    ),
    grid::polygonGrob(
      x = unit(c(0.51, 0.83, 0.83, 0.99, 0.83, 0.83, 0.51), "npc"),
      y = unit(c(0.32, 0.32, 0.08, 0.50, 0.92, 0.68, 0.68), "npc"),
      gp = grid::gpar(fill = cns_color, col = NA)
    ),
    grid::textGrob(
      "LN", x = unit(0.32, "npc"), y = unit(0.50, "npc"),
      gp = grid::gpar(col = "white", fontsize = 15, fontface = "bold")
    ),
    grid::textGrob(
      "CNS", x = unit(0.68, "npc"), y = unit(0.50, "npc"),
      gp = grid::gpar(col = "white", fontsize = 15, fontface = "bold")
    )
  )

  volcano_plot <- ggplot(deg, aes(x = log2FC, y = neg_log10_q)) +
    geom_point(aes(color = volcano_group), size = 1.45, alpha = 0.72, stroke = 0) +
    geom_vline(
      xintercept = c(-minimum_abs_log2fc, minimum_abs_log2fc),
      color = "#8A8A8A", linewidth = 0.45, linetype = "dotted"
    ) +
    geom_vline(
      xintercept = 0, color = "#B8B8B8", linewidth = 0.35, linetype = "dashed"
    ) +
    geom_hline(
      yintercept = -log10(q_threshold), color = "#8A8A8A",
      linewidth = 0.45, linetype = "dotted"
    ) +
    scale_color_manual(
      values = c(
        "LN higher" = ln_color,
        "Not significant" = neutral_color,
        "CNS higher" = cns_color
      ),
      drop = FALSE
    ) +
    scale_x_continuous(
      breaks = scales::breaks_pretty(n = 8),
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    scale_y_continuous(
      breaks = scales::breaks_pretty(n = 6),
      expand = expansion(mult = c(0, 0.03))
    ) +
    labs(
      x = expression(log[2] * "(Fold change: CNS / LN)"),
      y = expression(-log[10] * "(MAST BH-adjusted q-value)")
    ) +
    annotation_custom(
      grob = arrow_grob,
      xmin = -x_limit, xmax = x_limit,
      ymin = -0.21 * y_display_cap, ymax = -0.055 * y_display_cap
    ) +
    coord_cartesian(
      xlim = c(-x_limit, x_limit), ylim = c(0, y_display_cap), clip = "off"
    ) +
    theme_classic(base_size = 15, base_family = "sans") +
    theme(
      legend.position = "none",
      axis.title = element_text(face = "bold", color = "black"),
      axis.text = element_text(color = "black"),
      axis.line = element_line(color = "black", linewidth = 0.7),
      axis.ticks = element_line(color = "black", linewidth = 0.6),
      plot.margin = margin(t = 14, r = 28, b = 104, l = 22)
    )
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    volcano_plot <- volcano_plot + ggrepel::geom_text_repel(
      data = label_table,
      aes(label = gene, color = volcano_group),
      size = 5,
      fontface = "bold",
      box.padding = 0.42,
      point.padding = 0.22,
      min.segment.length = 0,
      segment.size = 0.45,
      max.overlaps = Inf,
      seed = seed_use,
      direction = "both",
      force = 1.2,
      show.legend = FALSE
    )
  } else {
    volcano_plot <- volcano_plot + geom_text(
      data = label_table,
      aes(label = gene, color = volcano_group),
      size = 4.7, fontface = "bold", vjust = -0.5, show.legend = FALSE
    )
  }

  figure_tiff <- file.path(
    figure_output_dir, "figure_2c_Kolz_Th17_CNS_vs_LN_volcano.tiff"
  )
  ggsave(
    figure_tiff,
    volcano_plot,
    width = 8.2,
    height = 8.2,
    units = "in",
    dpi = 300,
    compression = "lzw",
    bg = "white"
  )
  write.csv(
    label_table,
    file.path(output_dir, "Kolz_Th17_CNS_vs_LN_volcano_labeled_genes.csv"),
    row.names = FALSE
  )
  cat("Saved Figure 2c:\n", figure_tiff, "\n", sep = "")

  list(
    deg = deg,
    pathways = read.csv(pathway_csv, check.names = FALSE),
    labels = label_table,
    plot = volcano_plot
  )
})

################################################################################
# Figure 2d: Kolz Th17 CNS versus LN directional pathway enrichment
################################################################################

figure_2d_kolz_th17 <- local({
  local_analysis_root <- Sys.getenv(
    "EAE_LOCAL_ANALYSIS_ROOT",
    unset = "/Users/shingheimok/Desktop/Phd_Study/research/Wu_lab/b_cells"
  )
  output_dir <- file.path(
    local_analysis_root,
    "outs", "output", "DEG", "Kolz_Th17_CNS_vs_LN_MAST_gprofiler2"
  )
  figure_output_dir <- file.path(local_analysis_root, "figures", "figure_2")
  pathway_csv <- file.path(
    output_dir,
    "Kolz_Th17_CNS_vs_LN_gprofiler2_all_passing_pathways.csv"
  )
  if (!file.exists(pathway_csv)) {
    stop("Kolz pathway table not found: ", pathway_csv, call. = FALSE)
  }
  pathways <- read.csv(pathway_csv, check.names = FALSE, stringsAsFactors = FALSE)

  # All LN-higher terms are shown because only two passed the prespecified
  # filters. Ten nonredundant CNS-higher terms summarize translation/anabolism,
  # activation, cytokine/chemokine signaling, glycolysis, and adhesion.
  selected_pathways <- data.frame(
    query = c(rep("CNS_higher", 10L), rep("LN_higher", 2L)),
    source = c(
      "GO:BP", "KEGG", "WP", "WP", "REAC",
      "KEGG", "GO:BP", "WP", "GO:BP", "GO:BP",
      "CORUM", "WP"
    ),
    term_id = c(
      "GO:0002181", "KEGG:01230", "WP:WP5242", "WP:WP157",
      "REAC:R-MMU-380108", "KEGG:04060", "GO:0006954", "WP:WP493",
      "GO:0007155", "GO:0043066", "CORUM:6887", "WP:WP29"
    ),
    display_term = c(
      "Cytoplasmic translation",
      "Amino-acid biosynthesis",
      "IL-17A signaling",
      "Glycolysis/gluconeogenesis",
      "Chemokine receptor signaling",
      "Cytokine-receptor interaction",
      "Inflammatory response",
      "MAPK signaling",
      "Cell adhesion",
      "Negative regulation of apoptosis",
      "Ragulator-AXIN/LKB1-AMPK complex",
      "Notch signaling"
    ),
    stringsAsFactors = FALSE
  )
  selected_result <- merge(
    selected_pathways,
    pathways,
    by = c("query", "source", "term_id"),
    all.x = TRUE,
    sort = FALSE
  )
  missing_selected <- selected_result$term_id[
    is.na(selected_result$adjusted_p_value_gSCS)
  ]
  if (length(missing_selected) > 0L) {
    stop(
      "Selected Figure 2d terms were not found among passing results: ",
      paste(missing_selected, collapse = ", "),
      call. = FALSE
    )
  }
  selected_result$enrichment_direction <- ifelse(
    selected_result$query == "CNS_higher", "CNS", "LN"
  )
  selected_result$minus_log2_adjusted_p <- -log2(pmax(
    selected_result$adjusted_p_value_gSCS, .Machine$double.xmin
  ))

  ln_scores <- selected_result$minus_log2_adjusted_p[
    selected_result$enrichment_direction == "LN"
  ]
  ln_terms <- selected_result$display_term[
    selected_result$enrichment_direction == "LN"
  ]
  cns_scores <- selected_result$minus_log2_adjusted_p[
    selected_result$enrichment_direction == "CNS"
  ]
  cns_terms <- selected_result$display_term[
    selected_result$enrichment_direction == "CNS"
  ]
  spacer_label <- " "
  term_levels <- c(
    ln_terms[order(ln_scores)], spacer_label, cns_terms[order(cns_scores)]
  )
  selected_result$display_term <- factor(
    selected_result$display_term, levels = term_levels
  )
  spacer_row <- selected_result[1L, , drop = FALSE]
  spacer_row$display_term <- factor(spacer_label, levels = term_levels)
  spacer_row$enrichment_direction <- NA_character_
  spacer_row$minus_log2_adjusted_p <- 0
  plot_data <- rbind(selected_result, spacer_row)

  cns_color <- "#FF6B6B"
  ln_color <- "#5B9DF5"
  x_axis_max <- 5 * ceiling(
    max(plot_data$minus_log2_adjusted_p, na.rm = TRUE) / 5
  )
  pathway_plot <- ggplot(
    plot_data,
    aes(x = display_term, y = minus_log2_adjusted_p, fill = enrichment_direction)
  ) +
    geom_col(width = 0.82, na.rm = TRUE) +
    coord_flip(clip = "off") +
    scale_fill_manual(
      values = c(CNS = cns_color, LN = ln_color),
      breaks = c("CNS", "LN"),
      drop = FALSE
    ) +
    scale_y_continuous(
      breaks = seq(0, x_axis_max, by = 5),
      limits = c(0, x_axis_max),
      expand = expansion(mult = c(0, 0))
    ) +
    labs(
      x = "Functional ontology term",
      y = "−log₂(g:SCS-adjusted P)",
      fill = "Enrichment\nin genes for"
    ) +
    theme_classic(base_size = 18) +
    theme(
      text = element_text(face = "bold", color = "black"),
      axis.title.x = element_text(size = 19, face = "bold", margin = margin(t = 9)),
      axis.title.y = element_text(size = 19, face = "bold", margin = margin(r = 15)),
      axis.text.x = element_text(size = 15, face = "bold", color = "black"),
      axis.text.y = element_text(size = 17, face = "bold", color = "black"),
      axis.line = element_line(linewidth = 1.4, color = "black"),
      axis.ticks = element_line(linewidth = 1.3, color = "black"),
      axis.ticks.length = unit(0.18, "cm"),
      legend.position = "none",
      plot.margin = margin(12, 18, 5, 18)
    )

  pathway_legend <- grid::grobTree(
    grid::textGrob(
      "Enrichment\nin genes for",
      x = unit(0.02, "npc"), y = unit(0.50, "npc"),
      hjust = 0, vjust = 0.5,
      gp = grid::gpar(fontsize = 17, fontface = "bold", lineheight = 0.95)
    ),
    grid::rectGrob(
      x = unit(1.58, "in"), y = unit(0.70, "npc"),
      width = unit(0.25, "in"), height = unit(0.25, "in"),
      gp = grid::gpar(fill = cns_color, col = NA)
    ),
    grid::textGrob(
      "CNS", x = unit(1.80, "in"), y = unit(0.70, "npc"), hjust = 0,
      gp = grid::gpar(fontsize = 17, fontface = "bold")
    ),
    grid::rectGrob(
      x = unit(1.58, "in"), y = unit(0.28, "npc"),
      width = unit(0.25, "in"), height = unit(0.25, "in"),
      gp = grid::gpar(fill = ln_color, col = NA)
    ),
    grid::textGrob(
      "LN", x = unit(1.80, "in"), y = unit(0.28, "npc"), hjust = 0,
      gp = grid::gpar(fontsize = 17, fontface = "bold")
    )
  )
  figure_2d <- cowplot::ggdraw() +
    cowplot::draw_plot(pathway_plot, x = 0, y = 0.095, width = 1, height = 0.865) +
    cowplot::draw_grob(pathway_legend, x = 0.235, y = 0.015, width = 0.40, height = 0.09)

  figure_tiff <- file.path(
    figure_output_dir,
    "figure_2d_Kolz_Th17_CNS_vs_LN_pathway_enrichment.tiff"
  )
  selected_csv <- file.path(
    output_dir, "Kolz_Th17_CNS_vs_LN_displayed_pathways.csv"
  )
  ggsave(
    figure_tiff,
    figure_2d,
    width = 10.2,
    height = 7.4,
    units = "in",
    dpi = 300,
    compression = "lzw",
    bg = "white"
  )
  write.csv(selected_result, selected_csv, row.names = FALSE)

  # Exact-ID sensitivity comparison with MNG-higher terms from Figure 2b.
  mng_pathway_csv <- file.path(
    local_analysis_root,
    "outs", "output", "DEG", "validated_bcells_MNG_vs_dCLN_MAST_msigdbr",
    "gprofiler2_ORA",
    "figure_2b_MNG_vs_dCLN_gprofiler2_all_passing_pathways.csv"
  )
  overlap_summary <- NULL
  if (file.exists(mng_pathway_csv)) {
    mng_pathways <- read.csv(
      mng_pathway_csv, check.names = FALSE, stringsAsFactors = FALSE
    )
    mng_terms <- unique(mng_pathways$term_id[mng_pathways$query == "MNG_higher"])
    cns_terms_all <- unique(pathways$term_id[pathways$query == "CNS_higher"])
    overlap_terms <- intersect(mng_terms, cns_terms_all)
    overlap_summary <- data.frame(
      comparison = "MNG-higher versus Kolz Th17 CNS-higher significant terms",
      MNG_term_count = length(mng_terms),
      CNS_term_count = length(cns_terms_all),
      exact_term_overlap = length(overlap_terms),
      jaccard = length(overlap_terms) /
        length(union(mng_terms, cns_terms_all)),
      overlap_coefficient = length(overlap_terms) /
        min(length(mng_terms), length(cns_terms_all)),
      stringsAsFactors = FALSE
    )
    write.csv(
      overlap_summary,
      file.path(output_dir, "Kolz_CNS_vs_EAE_MNG_pathway_overlap_summary.csv"),
      row.names = FALSE
    )
    write.csv(
      merge(
        mng_pathways[
          mng_pathways$query == "MNG_higher" &
            mng_pathways$term_id %in% overlap_terms,
          c("source", "term_id", "term_name", "adjusted_p_value_gSCS")
        ],
        pathways[
          pathways$query == "CNS_higher" & pathways$term_id %in% overlap_terms,
          c("source", "term_id", "term_name", "adjusted_p_value_gSCS")
        ],
        by = c("source", "term_id"),
        suffixes = c("_MNG", "_Kolz_CNS")
      ),
      file.path(output_dir, "Kolz_CNS_vs_EAE_MNG_exact_overlapping_pathways.csv"),
      row.names = FALSE
    )
  }

  cat("Saved Figure 2d:\n", figure_tiff, "\n", sep = "")
  list(
    all_passing_terms = pathways,
    displayed_terms = selected_result,
    overlap_summary = overlap_summary,
    plot = figure_2d
  )
})

################################################################################
# Figure 2d1: validated B-cell reclustering and tissue-colored UMAP
################################################################################

# Exact input cell set: the frozen 12,942-cell annotation.R-derived validated
# B-cell export used by the cNMF analysis (4,408 MNG and 8,534 dCLN cells).
# Cells are reprocessed using 2,000 VST HVGs, 10 PCs, k.param = 30,
# resolution = 0.2, and seed = 1234. The square UMAP is colored by tissue using
# the Figure 2a/2b palette (MNG red; dCLN blue); points are fully opaque and the
# legend is intentionally hidden.

figure_2d1_bcell_tissue_umap <- local({
  seed_use <- 1234L
  set.seed(seed_use)

  local_analysis_root <- Sys.getenv(
    "EAE_LOCAL_ANALYSIS_ROOT",
    unset = "/Users/shingheimok/Desktop/Phd_Study/research/Wu_lab/b_cells"
  )
  input_dir <- file.path(
    local_analysis_root, "outs", "output", "cnmf_eps017_min5_min31"
  )
  output_dir <- file.path(
    local_analysis_root, "outs", "output",
    "bcell_recluster_pc10_k30_res02_seed1234"
  )
  figure_output_dir <- file.path(
    local_analysis_root, "figures", "figure_2"
  )
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figure_output_dir, recursive = TRUE, showWarnings = FALSE)

  matrix_file <- file.path(
    input_dir, "validated_bcell_counts_cell_by_gene.mtx"
  )
  cell_file <- file.path(input_dir, "validated_bcell_cells.tsv")
  gene_file <- file.path(input_dir, "validated_bcell_genes.tsv")
  metadata_file <- file.path(input_dir, "validated_bcell_metadata.csv")
  required_files <- c(matrix_file, cell_file, gene_file, metadata_file)
  missing_files <- required_files[!file.exists(required_files)]
  if (length(missing_files) > 0L) {
    stop(
      "Missing required Figure 2d1 input files:\n",
      paste(missing_files, collapse = "\n"),
      call. = FALSE
    )
  }

  cat("Reading the frozen validated-B-cell count matrix for Figure 2d1...\n")
  cell_by_gene <- Matrix::readMM(matrix_file)
  cells <- read.delim(
    cell_file, header = FALSE, stringsAsFactors = FALSE
  )[[1]]
  genes <- read.delim(
    gene_file, header = FALSE, stringsAsFactors = FALSE
  )[[1]]
  metadata <- read.csv(
    metadata_file, stringsAsFactors = FALSE, check.names = FALSE
  )

  if (!identical(dim(cell_by_gene), c(length(cells), length(genes)))) {
    stop(
      "Figure 2d1 matrix dimensions do not match the manifests: matrix = ",
      paste(dim(cell_by_gene), collapse = " x "),
      "; manifests = ", length(cells), " x ", length(genes),
      call. = FALSE
    )
  }
  if (anyDuplicated(cells) || anyDuplicated(genes)) {
    stop("Figure 2d1 cell and gene identifiers must be unique.", call. = FALSE)
  }
  if (!"cell" %in% colnames(metadata) || anyDuplicated(metadata$cell)) {
    stop(
      "Figure 2d1 metadata must contain one unique row per cell.",
      call. = FALSE
    )
  }

  metadata <- metadata[match(cells, metadata$cell), , drop = FALSE]
  if (
    anyNA(metadata$cell) ||
      !identical(as.character(metadata$cell), as.character(cells))
  ) {
    stop(
      "Figure 2d1 metadata does not align to the count-matrix cells.",
      call. = FALSE
    )
  }

  required_metadata <- c(
    "celltype_major", "celltype_minor", "TissueGroup", "sample_id"
  )
  missing_metadata <- setdiff(required_metadata, colnames(metadata))
  if (length(missing_metadata) > 0L) {
    stop(
      "Missing Figure 2d1 metadata: ",
      paste(missing_metadata, collapse = ", "),
      call. = FALSE
    )
  }
  if (
    any(metadata$celltype_major != "B_cell") ||
      any(metadata$celltype_minor != "B_cell")
  ) {
    stop(
      "The Figure 2d1 input contains cells outside the validated B-cell set.",
      call. = FALSE
    )
  }

  metadata$TissueGroup <- factor(
    metadata$TissueGroup, levels = c("MNG", "dCLN")
  )
  if (anyNA(metadata$TissueGroup)) {
    stop(
      "Figure 2d1 TissueGroup must contain only MNG and dCLN.",
      call. = FALSE
    )
  }
  expected_tissue_counts <- c(MNG = 4408L, dCLN = 8534L)
  observed_tissue_counts <- table(metadata$TissueGroup)
  if (!identical(
    as.integer(observed_tissue_counts[names(expected_tissue_counts)]),
    as.integer(expected_tissue_counts)
  )) {
    stop(
      "Figure 2d1 tissue counts differ from the locked input. Observed: ",
      paste(names(observed_tissue_counts), observed_tissue_counts,
            collapse = ", "),
      call. = FALSE
    )
  }

  counts <- methods::as(Matrix::t(cell_by_gene), "dgCMatrix")
  rownames(counts) <- genes
  colnames(counts) <- cells
  rownames(metadata) <- metadata$cell

  cat("Reclustering validated B cells for Figure 2d1...\n")
  obj <- CreateSeuratObject(
    counts = counts,
    assay = "RNA",
    meta.data = metadata,
    project = "validated_B_cells_pc10_res02"
  )
  rm(cell_by_gene, counts)
  invisible(gc())
  DefaultAssay(obj) <- "RNA"

  set.seed(seed_use)
  obj <- NormalizeData(obj, assay = "RNA", verbose = FALSE)
  set.seed(seed_use)
  obj <- FindVariableFeatures(
    obj,
    assay = "RNA",
    selection.method = "vst",
    nfeatures = 2000,
    verbose = FALSE
  )
  set.seed(seed_use)
  obj <- ScaleData(
    obj,
    assay = "RNA",
    features = VariableFeatures(obj),
    verbose = FALSE
  )
  set.seed(seed_use)
  obj <- RunPCA(
    obj,
    assay = "RNA",
    features = VariableFeatures(obj),
    npcs = 10,
    seed.use = seed_use,
    verbose = FALSE
  )
  set.seed(seed_use)
  obj <- FindNeighbors(
    obj,
    reduction = "pca",
    dims = 1:10,
    k.param = 30,
    verbose = FALSE
  )
  set.seed(seed_use)
  obj <- FindClusters(
    obj,
    resolution = 0.2,
    random.seed = seed_use,
    verbose = FALSE
  )
  set.seed(seed_use)
  obj <- RunUMAP(
    obj,
    reduction = "pca",
    dims = 1:10,
    n.neighbors = 30,
    seed.use = seed_use,
    verbose = FALSE
  )

  cluster_column <- "bcell_cluster_pc10_k30_res02_seed1234"
  obj[[cluster_column]] <- as.character(Idents(obj))
  umap <- Embeddings(obj, reduction = "umap")
  plot_df <- data.frame(
    cell = rownames(umap),
    UMAP_1 = umap[, 1],
    UMAP_2 = umap[, 2],
    tissue = obj$TissueGroup,
    cluster = obj@meta.data[[cluster_column]],
    stringsAsFactors = FALSE
  )

  # Random plotting order avoids systematically drawing either tissue last.
  set.seed(seed_use)
  plot_df <- plot_df[sample.int(nrow(plot_df)), , drop = FALSE]
  tissue_colors <- c(MNG = "#FF6B6B", dCLN = "#5B9DF5")

  figure_2d1 <- ggplot(
    plot_df,
    aes(x = UMAP_1, y = UMAP_2, color = tissue)
  ) +
    geom_point(size = 2.00, alpha = 1, stroke = 0) +
    scale_color_manual(values = tissue_colors, drop = FALSE) +
    scale_x_continuous(expand = expansion(mult = 0.055)) +
    scale_y_continuous(expand = expansion(mult = 0.055)) +
    coord_cartesian() +
    labs(x = "UMAP 1", y = "UMAP 2") +
    theme_classic(base_size = 18) +
    theme(
      axis.title = element_text(size = 32, face = "plain"),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      panel.border = element_rect(
        color = "black", fill = NA, linewidth = 1.15
      ),
      aspect.ratio = 1,
      legend.position = "none",
      plot.margin = margin(12, 12, 12, 12)
    )

  figure_tiff <- file.path(
    figure_output_dir,
    "figure_2d1_validated_bcells_pc10_k30_res02_seed1234_tissue_umap.tiff"
  )
  object_rds <- file.path(
    output_dir, "validated_bcells_pc10_k30_res02_seed1234.rds"
  )
  metadata_csv <- file.path(
    output_dir, "validated_bcells_pc10_k30_res02_seed1234_metadata.csv"
  )
  coordinate_csv <- file.path(
    output_dir,
    "validated_bcells_pc10_k30_res02_seed1234_umap_coordinates.csv"
  )
  cluster_summary_csv <- file.path(
    output_dir,
    "validated_bcells_pc10_k30_res02_seed1234_cluster_summary.csv"
  )

  ggsave(
    filename = figure_tiff,
    plot = figure_2d1,
    device = "tiff",
    width = 8,
    height = 8,
    units = "in",
    dpi = 300,
    compression = "lzw",
    bg = "white"
  )
  saveRDS(obj, object_rds, compress = TRUE)
  write.csv(obj@meta.data, metadata_csv, row.names = TRUE)
  write.csv(plot_df, coordinate_csv, row.names = FALSE)

  cluster_summary <- obj@meta.data |>
    count(
      TissueGroup,
      .data[[cluster_column]],
      name = "n_cells"
    ) |>
    arrange(TissueGroup, .data[[cluster_column]])
  write.csv(cluster_summary, cluster_summary_csv, row.names = FALSE)

  cat("Saved Figure 2d1:\n", figure_tiff, "\n", sep = "")
  list(
    object = obj,
    metadata = obj@meta.data,
    umap_coordinates = plot_df,
    cluster_summary = cluster_summary,
    plot = figure_2d1,
    figure_tiff = figure_tiff
  )
})

################################################################################
# Figure 2d2: k=6 cNMF pairwise program similarity and tissue bias
################################################################################

# Pairwise similarity is the Jaccard index among each program's top 50 genes.
# Programs are ordered by average-linkage clustering of top-50 gene overlap.
# The annotation strip identifies the tissue with higher mean relative usage.

figure_2d2_cnmf_pairwise_tissue_bias <- local({

  suppressPackageStartupMessages({
    library(ComplexHeatmap)
    library(circlize)
    library(dplyr)
    library(grid)
  })

  ################################################################################
  # k=6 cNMF pairwise program-similarity matrix with tissue-bias annotation
  #
  # This follows the existing pairwise cNMF template:
  #   - top 50 genes per program
  #   - pairwise Jaccard similarity, displayed as percentage and clipped at 25%
  #   - average-linkage ordering using top50 minus overlap as distance
  #   - square matrix, Cartesian row orientation, and the original figure ratio
  #
  # The annotation strip is changed from DBSCAN bias to Tissue bias. A program is
  # labeled Meninges or dCLN according to the tissue with the larger mean relative
  # usage across all validated B cells. M1/M2 define Meninges and L1/L2 define
  # dCLN. The cNMF solution is the k=6, density-threshold=0.1 run after removing
  # structural ribosomal, TCR-locus, and all Ig-locus genes.
  ################################################################################

  base_dir <- "/Users/shingheimok/Desktop/Phd_Study/research/Wu_lab/b_cells"
  run_dir <- file.path(
    base_dir,
    "outs/output/cnmf_stable_bcells_no_ribo_tcr_ig_k6_k15_hvg5000",
    "stable_bcells_no_ribo_tcr_ig_k6_k15_hvg5000"
  )
  score_file <- file.path(
    run_dir,
    "stable_bcells_no_ribo_tcr_ig_k6_k15_hvg5000.gene_spectra_score.k_6.dt_0_1.txt"
  )
  usage_file <- file.path(
    run_dir,
    "k6_density_0_1_sample_panels",
    "k6_dt01_program_usage_by_cell_and_sample.csv"
  )
  output_dir <- file.path(run_dir, "k6_pairwise_similarity_tissue_bias")
  figure_dir <- file.path(base_dir, "figures", "figure_2")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

  required_files <- c(score_file, usage_file)
  missing_files <- required_files[!file.exists(required_files)]
  if (length(missing_files) > 0L) {
    stop(
      "Missing required k=6 cNMF input files:\n",
      paste(missing_files, collapse = "\n"),
      call. = FALSE
    )
  }

  figure_tiff <- file.path(
    figure_dir,
    "k6_cnmf_pairwise_similarity_tissue_bias.tiff"
  )
  pairwise_csv <- file.path(
    output_dir, "k6_cnmf_program_pairwise_similarity_top50.csv"
  )
  mapping_csv <- file.path(
    output_dir, "k6_cnmf_program_order_mapping.csv"
  )
  tissue_bias_csv <- file.path(
    output_dir, "k6_cnmf_program_tissue_bias.csv"
  )
  top_gene_csv <- file.path(
    output_dir, "k6_cnmf_program_top50_genes.csv"
  )

  top_n <- 50L
  jaccard_clip <- 25

  normalize_program_id <- function(x) {
    x <- as.character(x)
    x <- sub("^Program_", "", x)
    x <- sub("^GEP", "", x)
    as.character(as.integer(x))
  }

  cat("Reading k=6 cNMF gene scores...\n")
  score_raw <- read.table(
    score_file,
    header = TRUE,
    row.names = 1,
    sep = "\t",
    check.names = FALSE,
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE
  )
  score_raw[] <- lapply(
    score_raw,
    function(x) suppressWarnings(as.numeric(as.character(x)))
  )
  score_df <- as.data.frame(t(as.matrix(score_raw)), check.names = FALSE)
  score_df <- score_df[!duplicated(rownames(score_df)), , drop = FALSE]
  colnames(score_df) <- normalize_program_id(colnames(score_df))
  program_ids <- colnames(score_df)

  if (length(program_ids) != 6L || anyNA(program_ids)) {
    stop(
      "Expected exactly six identifiable cNMF programs; observed: ",
      paste(program_ids, collapse = ", "),
      call. = FALSE
    )
  }

  program_gene_list <- setNames(vector("list", length(program_ids)), program_ids)
  for (program_id in program_ids) {
    scores <- score_df[[program_id]]
    names(scores) <- rownames(score_df)
    scores <- scores[!is.na(scores)]
    ranked_genes <- names(sort(scores, decreasing = TRUE))
    program_gene_list[[program_id]] <- unique(ranked_genes)[seq_len(top_n)]
  }

  top_gene_df <- bind_rows(lapply(program_ids, function(program_id) {
    data.frame(
      program = paste0("GEP", program_id),
      rank = seq_along(program_gene_list[[program_id]]),
      gene = program_gene_list[[program_id]],
      stringsAsFactors = FALSE
    )
  }))
  write.csv(top_gene_df, top_gene_csv, row.names = FALSE)

  n_programs <- length(program_ids)
  overlap_mat <- matrix(
    top_n,
    nrow = n_programs,
    ncol = n_programs,
    dimnames = list(program_ids, program_ids)
  )
  jaccard_mat <- matrix(
    1,
    nrow = n_programs,
    ncol = n_programs,
    dimnames = list(program_ids, program_ids)
  )
  distance_mat <- matrix(
    0,
    nrow = n_programs,
    ncol = n_programs,
    dimnames = list(program_ids, program_ids)
  )

  pairwise_rows <- vector("list", choose(n_programs, 2L))
  pair_index <- 1L
  for (i in seq_len(n_programs - 1L)) {
    for (j in seq.int(i + 1L, n_programs)) {
      program_1 <- program_ids[[i]]
      program_2 <- program_ids[[j]]
      genes_1 <- program_gene_list[[program_1]]
      genes_2 <- program_gene_list[[program_2]]
      overlap_n <- length(intersect(genes_1, genes_2))
      union_n <- length(union(genes_1, genes_2))
      jaccard <- overlap_n / union_n

      overlap_mat[program_1, program_2] <- overlap_n
      overlap_mat[program_2, program_1] <- overlap_n
      jaccard_mat[program_1, program_2] <- jaccard
      jaccard_mat[program_2, program_1] <- jaccard
      distance_mat[program_1, program_2] <- top_n - overlap_n
      distance_mat[program_2, program_1] <- top_n - overlap_n

      pairwise_rows[[pair_index]] <- data.frame(
        program_1 = paste0("GEP", program_1),
        program_2 = paste0("GEP", program_2),
        overlap_n = overlap_n,
        union_n = union_n,
        jaccard = jaccard,
        jaccard_percent = 100 * jaccard,
        distance_top50_minus_overlap = top_n - overlap_n,
        stringsAsFactors = FALSE
      )
      pair_index <- pair_index + 1L
    }
  }
  pairwise_df <- bind_rows(pairwise_rows) |>
    arrange(distance_top50_minus_overlap, desc(overlap_n))
  write.csv(pairwise_df, pairwise_csv, row.names = FALSE)

  # Faithful template ordering: average-linkage clustering of top50-overlap
  # distance. Axis labels retain the original GEP numbers even though their
  # positions follow the clustered order.
  program_hclust <- hclust(as.dist(distance_mat), method = "average")
  ordered_ids <- program_hclust$labels[program_hclust$order]
  program_mapping <- data.frame(
    matrix_position = seq_along(ordered_ids),
    original_program = paste0("GEP", ordered_ids),
    stringsAsFactors = FALSE
  )
  write.csv(program_mapping, mapping_csv, row.names = FALSE)

  cat("Calculating tissue bias from relative cNMF usage...\n")
  usage_df <- read.csv(usage_file, stringsAsFactors = FALSE, check.names = FALSE)
  required_usage_columns <- c("cell", "sample", "program", "usage")
  missing_usage_columns <- setdiff(required_usage_columns, colnames(usage_df))
  if (length(missing_usage_columns) > 0L) {
    stop(
      "Missing usage columns: ",
      paste(missing_usage_columns, collapse = ", "),
      call. = FALSE
    )
  }
  if (!all(usage_df$sample %in% c("M1", "M2", "L1", "L2"))) {
    stop("Unexpected sample labels in the k=6 usage table.", call. = FALSE)
  }
  usage_df <- usage_df |>
    mutate(
      tissue = ifelse(sample %in% c("M1", "M2"), "Meninges", "dCLN"),
      program_id = normalize_program_id(program)
    )

  tissue_means_long <- usage_df |>
    group_by(program_id, tissue) |>
    summarise(
      n_cells = n(),
      mean_usage = mean(usage),
      median_usage = median(usage),
      .groups = "drop"
    )
  tissue_means_wide <- reshape(
    as.data.frame(tissue_means_long),
    idvar = "program_id",
    timevar = "tissue",
    direction = "wide"
  )
  names(tissue_means_wide) <- sub("mean_usage\\.", "mean_usage_", names(tissue_means_wide))
  names(tissue_means_wide) <- sub("median_usage\\.", "median_usage_", names(tissue_means_wide))
  names(tissue_means_wide) <- sub("n_cells\\.", "n_cells_", names(tissue_means_wide))
  tissue_means_wide$tissue_bias <- ifelse(
    tissue_means_wide$mean_usage_Meninges >= tissue_means_wide$mean_usage_dCLN,
    "Meninges",
    "dCLN"
  )
  tissue_means_wide$mean_difference_Meninges_minus_dCLN <-
    tissue_means_wide$mean_usage_Meninges - tissue_means_wide$mean_usage_dCLN
  tissue_means_wide$program <- paste0("GEP", tissue_means_wide$program_id)
  tissue_bias_df <- tissue_means_wide |>
    select(
      program,
      n_cells_Meninges,
      n_cells_dCLN,
      mean_usage_Meninges,
      mean_usage_dCLN,
      mean_difference_Meninges_minus_dCLN,
      median_usage_Meninges,
      median_usage_dCLN,
      tissue_bias
    ) |>
    arrange(as.integer(sub("GEP", "", program)))
  write.csv(tissue_bias_df, tissue_bias_csv, row.names = FALSE)

  tissue_bias_by_id <- setNames(
    tissue_means_wide$tissue_bias,
    tissue_means_wide$program_id
  )
  if (anyNA(tissue_bias_by_id[ordered_ids])) {
    stop("Tissue bias could not be assigned to every ordered program.", call. = FALSE)
  }

  # Cartesian orientation: columns increase left to right; rows are reversed so
  # displayed program numbers increase from bottom to top, as in the template.
  column_ids <- ordered_ids
  row_ids <- rev(ordered_ids)
  jaccard_percent_display <- pmin(
    100 * jaccard_mat[row_ids, column_ids, drop = FALSE],
    jaccard_clip
  )

  mng_color <- "#FF6B6B"
  dcln_color <- "#5B9DF5"
  tissue_bias_levels <- c("Meninges", "dCLN")
  top_annotation <- HeatmapAnnotation(
    `Tissue bias` = factor(
      tissue_bias_by_id[column_ids],
      levels = tissue_bias_levels
    ),
    col = list(
      `Tissue bias` = c("Meninges" = mng_color, "dCLN" = dcln_color)
    ),
    show_annotation_name = FALSE,
    show_legend = FALSE,
    gp = gpar(col = "grey70"),
    simple_anno_size = unit(2.5, "mm")
  )

  similarity_color_function <- circlize::colorRamp2(
    c(0, 5, 10, 15, 20, 25),
    c("white", "#FFF4A3", "#FDB366", "#F46D43", "#9C179E", "#313695")
  )

  pairwise_heatmap <- Heatmap(
    jaccard_percent_display,
    name = "Similarity\n(Jaccard index)",
    col = similarity_color_function,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    top_annotation = top_annotation,
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    show_heatmap_legend = FALSE,
    row_labels = row_ids,
    column_labels = column_ids,
    row_names_side = "left",
    column_names_side = "bottom",
    column_title_side = "bottom",
    row_names_gp = gpar(fontsize = 14),
    column_names_gp = gpar(fontsize = 14),
    column_names_rot = 0,
    rect_gp = gpar(col = "#D9D9D9", lwd = 0.5),
    border = TRUE,
    border_gp = gpar(col = "black", lwd = 2),
    width = unit(95, "mm"),
    height = unit(95, "mm"),
    row_title = "cNMF programs (GEP)",
    column_title = "cNMF programs (GEP)",
    row_title_gp = gpar(fontsize = 16),
    column_title_gp = gpar(fontsize = 16)
  )

  similarity_legend <- Legend(
    title = "Similarity\n(Jaccard index)",
    col_fun = similarity_color_function,
    at = c(5, 10, 15, 20, 25),
    labels = c("5", "10", "15", "20", "25"),
    direction = "horizontal",
    legend_width = unit(38, "mm"),
    title_gp = gpar(fontsize = 15, fontface = "bold"),
    labels_gp = gpar(fontsize = 13)
  )

  tissue_bias_legend <- Legend(
    at = tissue_bias_levels,
    labels = tissue_bias_levels,
    legend_gp = gpar(
      fill = c(mng_color, dcln_color),
      col = c(mng_color, dcln_color)
    ),
    title = "Tissue bias",
    ncol = 1,
    title_gp = gpar(fontsize = 15, fontface = "bold"),
    labels_gp = gpar(fontsize = 14),
    grid_width = unit(5.5, "mm"),
    grid_height = unit(5.5, "mm")
  )

  cat("Rendering the k=6 pairwise plot...\n")
  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_tiff(
      filename = figure_tiff,
      width = 7.7,
      height = 5.6,
      units = "in",
      res = 300,
      compression = "lzw",
      background = "white"
    )
  } else {
    grDevices::tiff(
      filename = figure_tiff,
      width = 7.7,
      height = 5.6,
      units = "in",
      res = 300,
      compression = "lzw",
      bg = "white"
    )
  }
  grid.newpage()
  pushViewport(
    viewport(
      layout = grid.layout(
        nrow = 1,
        ncol = 2,
        widths = unit.c(unit(1, "null"), unit(42, "mm"))
      )
    )
  )

  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
  draw(
    pairwise_heatmap,
    newpage = FALSE,
    padding = unit(c(4, 4, 4, 4), "mm")
  )
  popViewport()

  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2, clip = "off"))
  draw(
    tissue_bias_legend,
    x = unit(-18, "mm"),
    y = unit(0.82, "npc"),
    just = c("left", "top")
  )
  draw(
    similarity_legend,
    x = unit(-18, "mm"),
    y = unit(0.225, "npc"),
    just = c("left", "bottom")
  )
  popViewport(2)
  dev.off()

  cat("\nPairwise k=6 cNMF plot complete.\n")
  cat("Figure: ", figure_tiff, "\n", sep = "")
  cat("Program order: ", mapping_csv, "\n", sep = "")
  cat("Tissue bias: ", tissue_bias_csv, "\n", sep = "")
  print(tissue_bias_df, row.names = FALSE)

  list(
    pairwise_similarity = pairwise_df,
    program_order = program_mapping,
    tissue_bias = tissue_bias_df,
    top_genes = top_gene_df,
    figure_tiff = figure_tiff
  )
})

################################################################################
# Figure 2d3: k=6 cNMF per-cell program-usage UMAPs (GEP1-GEP5)
################################################################################

# Uses the Figure 2d1 UMAP and the density-threshold 0.1 k=6 cNMF usages.
# GEP6 is omitted because it is a contamination-associated nuisance component.

figure_2d3_cnmf_program_usage_umaps <- local({

  suppressPackageStartupMessages({
    library(cowplot)
    library(ggplot2)
    library(grid)
  })

  ################################################################################
  # k=6 cNMF usage UMAPs for GEP1-GEP5
  #
  # Uses the exact 12,942 validated B cells and PC10/k30/res0.2/seed1234 UMAP from
  # Figure 2d1. Per-cell usages come from the density-threshold 0.1 k=6 solution
  # generated after excluding structural ribosomal, TCR-locus, and all Ig-locus
  # genes. GEP6 is intentionally omitted.
  ################################################################################

  set.seed(1234)

  base_dir <- "/Users/shingheimok/Desktop/Phd_Study/research/Wu_lab/b_cells"
  coordinate_file <- file.path(
    base_dir,
    "outs/output/bcell_recluster_pc10_k30_res02_seed1234",
    "validated_bcells_pc10_k30_res02_seed1234_umap_coordinates.csv"
  )
  usage_file <- file.path(
    base_dir,
    "outs/output/cnmf_stable_bcells_no_ribo_tcr_ig_k6_k15_hvg5000",
    "stable_bcells_no_ribo_tcr_ig_k6_k15_hvg5000",
    "k6_density_0_1_sample_panels",
    "k6_dt01_program_usage_by_cell_and_sample.csv"
  )
  figure_dir <- file.path(base_dir, "figures", "figure_2")
  output_tiff <- file.path(
    figure_dir,
    "k6_cnmf_GEP1_GEP5_usage_umaps.tiff"
  )
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

  required_files <- c(coordinate_file, usage_file)
  missing_files <- required_files[!file.exists(required_files)]
  if (length(missing_files) > 0L) {
    stop(
      "Missing required input files:\n",
      paste(missing_files, collapse = "\n"),
      call. = FALSE
    )
  }

  coordinates <- read.csv(
    coordinate_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  usage_long <- read.csv(
    usage_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  required_coordinate_columns <- c("cell", "UMAP_1", "UMAP_2", "tissue")
  required_usage_columns <- c("cell", "sample", "program", "usage")
  missing_coordinate_columns <- setdiff(
    required_coordinate_columns, colnames(coordinates)
  )
  missing_usage_columns <- setdiff(required_usage_columns, colnames(usage_long))
  if (length(missing_coordinate_columns) > 0L) {
    stop(
      "Missing UMAP columns: ",
      paste(missing_coordinate_columns, collapse = ", "),
      call. = FALSE
    )
  }
  if (length(missing_usage_columns) > 0L) {
    stop(
      "Missing usage columns: ",
      paste(missing_usage_columns, collapse = ", "),
      call. = FALSE
    )
  }
  if (nrow(coordinates) != 12942L || anyDuplicated(coordinates$cell)) {
    stop(
      "Expected 12,942 unique cells in the Figure 2d1 UMAP coordinates.",
      call. = FALSE
    )
  }
  expected_tissue_counts <- c(MNG = 4408L, dCLN = 8534L)
  observed_tissue_counts <- table(coordinates$tissue)
  if (!identical(
    as.integer(observed_tissue_counts[names(expected_tissue_counts)]),
    as.integer(expected_tissue_counts)
  )) {
    stop("UMAP tissue counts do not match the locked B-cell input.", call. = FALSE)
  }

  program_order <- paste0("GEP", 1:6)
  usage_counts <- table(usage_long$program)
  if (
    !all(program_order %in% names(usage_counts)) ||
      any(usage_counts[program_order] != nrow(coordinates))
  ) {
    stop("Each of GEP1-GEP6 must have one usage per UMAP cell.", call. = FALSE)
  }
  if (
    length(setdiff(coordinates$cell, usage_long$cell)) > 0L ||
      length(setdiff(unique(usage_long$cell), coordinates$cell)) > 0L
  ) {
    stop("The UMAP and cNMF usage cell sets do not match.", call. = FALSE)
  }

  # Shared 0-1 display scale follows the supplied template. A small number of
  # cNMF relative-usage values exceed 1 because of numerical normalization; these
  # values are retained in the source table and squished to the top display color.
  usage_colors <- c(
    "#FFF5F0", "#FEE0D2", "#FCBBA1", "#FC9272",
    "#FB6A4A", "#DE2D26", "#99000D"
  )

  make_program_panel <- function(program_id, show_x_title = FALSE) {
    program_usage <- usage_long[
      usage_long$program == program_id,
      c("cell", "usage"),
      drop = FALSE
    ]
    program_usage <- program_usage[
      match(coordinates$cell, program_usage$cell),
      ,
      drop = FALSE
    ]
    if (
      anyNA(program_usage$cell) ||
        !identical(as.character(program_usage$cell), as.character(coordinates$cell))
    ) {
      stop("Usage alignment failed for ", program_id, ".", call. = FALSE)
    }

    plot_data <- cbind(
      coordinates[, c("cell", "UMAP_1", "UMAP_2"), drop = FALSE],
      usage = program_usage$usage
    )
    # Plot low-usage cells first and high-usage cells last so focal regions remain
    # visible; randomized tie-breaking avoids ordering artifacts among equal values.
    tie_break <- sample.int(nrow(plot_data))
    plot_data <- plot_data[order(plot_data$usage, tie_break), , drop = FALSE]

    ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2, color = usage)) +
      geom_point(size = 0.42, alpha = 1, stroke = 0) +
      scale_color_gradientn(
        colors = usage_colors,
        limits = c(0, 1),
        breaks = seq(0, 1, by = 0.2),
        labels = sprintf("%.1f", seq(0, 1, by = 0.2)),
        oob = scales::squish,
        name = "Program\nusage",
        guide = guide_colorbar(
          title.position = "top",
          title.hjust = 0,
          barwidth = unit(5.2, "mm"),
          barheight = unit(39, "mm"),
          frame.colour = "black",
          frame.linewidth = 0.9,
          ticks.colour = "black",
          ticks.linewidth = 0.65,
          label.position = "right"
        )
      ) +
      scale_x_continuous(expand = expansion(mult = 0.045)) +
      scale_y_continuous(expand = expansion(mult = 0.045)) +
      coord_fixed() +
      labs(
        title = sub("GEP", "Program ", program_id),
        # Reserve identical title space in both rows so every square panel has
        # exactly the same dimensions; the first-row title is transparent.
        x = "UMAP 1",
        y = NULL
      ) +
      theme_classic(base_size = 17) +
      theme(
        plot.title = element_text(
          size = 21,
          face = "plain",
          hjust = 0.5,
          margin = margin(b = 7)
        ),
        axis.title.x = element_text(
          size = 20,
          color = if (show_x_title) "black" else "transparent",
          margin = margin(t = 5)
        ),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.border = element_rect(
          color = "black",
          fill = NA,
          linewidth = 0.85
        ),
        aspect.ratio = 1,
        legend.title = element_text(size = 17, face = "plain"),
        legend.text = element_text(size = 13),
        legend.position = "right",
        plot.margin = margin(5, 7, 5, 7)
      )
  }

  panels <- list(
    GEP1 = make_program_panel("GEP1", show_x_title = FALSE),
    GEP2 = make_program_panel("GEP2", show_x_title = FALSE),
    GEP3 = make_program_panel("GEP3", show_x_title = FALSE),
    GEP4 = make_program_panel("GEP4", show_x_title = TRUE),
    GEP5 = make_program_panel("GEP5", show_x_title = TRUE)
  )

  legend_grob <- cowplot::get_legend(panels$GEP5)
  panels_without_legends <- lapply(
    panels,
    function(panel) panel + theme(legend.position = "none")
  )

  top_row <- cowplot::plot_grid(
    panels_without_legends$GEP1,
    panels_without_legends$GEP2,
    panels_without_legends$GEP3,
    nrow = 1,
    align = "hv",
    axis = "tblr",
    rel_widths = c(1, 1, 1)
  )
  bottom_row <- cowplot::plot_grid(
    panels_without_legends$GEP4,
    panels_without_legends$GEP5,
    cowplot::ggdraw(legend_grob),
    nrow = 1,
    align = "hv",
    axis = "tblr",
    rel_widths = c(1, 1, 1)
  )
  figure_program_usage_umaps <- cowplot::plot_grid(
    top_row,
    bottom_row,
    ncol = 1,
    rel_heights = c(1, 1)
  )

  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_tiff(
      filename = output_tiff,
      width = 10.2,
      height = 7.0,
      units = "in",
      res = 300,
      compression = "lzw",
      background = "white"
    )
    print(figure_program_usage_umaps)
    dev.off()
  } else {
    ggsave(
      filename = output_tiff,
      plot = figure_program_usage_umaps,
      device = "tiff",
      width = 10.2,
      height = 7.0,
      units = "in",
      dpi = 300,
      compression = "lzw",
      bg = "white"
    )
  }

  cat("Saved k=6 GEP1-GEP5 usage UMAPs:\n", output_tiff, "\n", sep = "")

  list(
    plot = figure_program_usage_umaps,
    panels = panels,
    coordinates = coordinates,
    usage = usage_long,
    figure_tiff = output_tiff
  )
})

################################################################################
# Disabled legacy panel: MNG-higher versus Kolz Th17 CNS-higher DEG Venn
################################################################################

# Both DEG sets use q_hurdle < 0.05 and log2FC > 0.2. The Kolz comparison is
# Th17 CNS versus LN (not Th1); spleen is excluded. Only representative genes
# are printed, while the region counts are calculated from all qualifying DEGs.

if (FALSE) {
figure_2e_deg_venn <- local({
  local_analysis_root <- Sys.getenv(
    "EAE_LOCAL_ANALYSIS_ROOT",
    unset = "/Users/shingheimok/Desktop/Phd_Study/research/Wu_lab/b_cells"
  )
  figure_output_dir <- file.path(local_analysis_root, "figures", "figure_2")
  output_dir <- file.path(
    local_analysis_root,
    "outs", "output", "DEG", "MNG_vs_Kolz_Th17_CNS_DEG_venn"
  )
  dir.create(figure_output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  our_deg_csv <- file.path(
    local_analysis_root,
    "outs", "output", "DEG", "validated_bcells_MNG_vs_dCLN_MAST_msigdbr",
    "validated_bcells_MNG_vs_dCLN_MAST_all_genes.csv"
  )
  kolz_deg_csv <- file.path(
    local_analysis_root,
    "outs", "output", "DEG", "Kolz_Th17_CNS_vs_LN_MAST_gprofiler2",
    "Kolz_Th17_CNS_vs_LN_MAST_all_genes.csv"
  )
  if (!file.exists(our_deg_csv)) stop("Our DEG table was not found: ", our_deg_csv)
  if (!file.exists(kolz_deg_csv)) stop("Kolz DEG table was not found: ", kolz_deg_csv)

  q_threshold <- 0.05
  minimum_log2fc <- 0.2
  our_deg <- read.csv(our_deg_csv, check.names = FALSE, stringsAsFactors = FALSE)
  kolz_deg <- read.csv(kolz_deg_csv, check.names = FALSE, stringsAsFactors = FALSE)
  our_set <- sort(unique(our_deg$gene[
    !is.na(our_deg$q_hurdle) & our_deg$q_hurdle < q_threshold &
      our_deg$avg_log2FC_MNG_vs_dCLN > minimum_log2fc
  ]))
  kolz_set <- sort(unique(kolz_deg$gene[
    !is.na(kolz_deg$q_hurdle) & kolz_deg$q_hurdle < q_threshold &
      kolz_deg$avg_log2FC_CNS_vs_LN > minimum_log2fc
  ]))
  our_specific <- setdiff(our_set, kolz_set)
  shared <- intersect(our_set, kolz_set)
  kolz_specific <- setdiff(kolz_set, our_set)

  selected_genes <- data.frame(
    gene = c(
      # Our core: OXPHOS and MHC-I antigen processing/presentation.
      "Ndufa6", "Ndufv3", "Psmb9", "Psmb10", "Tap1", "H2-K1", "H2-D1", "Psmb8",
      # Other representative programs in our MNG B cells.
      "Gpx4", "Aldoa", "Eno1",
      # Shared core: serine/one-carbon, arginine, and amino-acid anabolism.
      "Mthfd2", "Ass1", "Shmt2", "Psat1", "Phgdh", "Asl", "Slc7a5", "Sat1",
      # Additional shared metabolic genes.
      "Tpi1", "Pgk1",
      # Kolz core: IL-17A/NF-kappaB/chemokine response and glycolysis.
      "Traf6", "Nfkbia", "Ccl4", "Ccl5", "Cxcl10", "Ccr5",
      "Ldha", "Pgam1", "Gpi1", "Hk2",
      # Other Kolz CNS activation/secretory genes.
      "Il2ra", "Prdm1"
    ),
    region = c(
      rep("Our MNG specific", 11L),
      rep("Shared", 10L),
      rep("Kolz Th17 CNS specific", 12L)
    ),
    core_gene = c(
      rep(TRUE, 8L), rep(FALSE, 3L),
      rep(TRUE, 8L), rep(FALSE, 2L),
      rep(TRUE, 10L), rep(FALSE, 2L)
    ),
    biological_theme = c(
      "OXPHOS", "OXPHOS",
      "MHC-I immunoproteasome", "MHC-I immunoproteasome",
      "MHC-I peptide transport", "MHC-I presentation",
      "MHC-I presentation", "MHC-I immunoproteasome",
      "redox", "glycolysis", "glycolysis",
      "one-carbon metabolism", "arginine biosynthesis",
      "serine/one-carbon metabolism", "serine biosynthesis",
      "serine biosynthesis", "arginine biosynthesis",
      "amino-acid transport", "polyamine/amino-acid metabolism",
      "glycolysis", "glycolysis",
      "IL-17A/NF-kappaB signaling", "IL-17A/NF-kappaB signaling",
      "chemokine response", "chemokine response", "chemokine response",
      "chemokine response", "glycolysis", "glycolysis", "glycolysis",
      "glycolysis", "lymphocyte activation", "antibody secretion"
    ),
    stringsAsFactors = FALSE
  )

  expected_sets <- list(
    "Our MNG specific" = our_specific,
    "Shared" = shared,
    "Kolz Th17 CNS specific" = kolz_specific
  )
  invalid_selected <- unlist(lapply(names(expected_sets), function(region_name) {
    genes_now <- selected_genes$gene[selected_genes$region == region_name]
    setdiff(genes_now, expected_sets[[region_name]])
  }), use.names = FALSE)
  if (length(invalid_selected) > 0L) {
    stop(
      "Selected Venn genes are not in their required DEG region: ",
      paste(unique(invalid_selected), collapse = ", "),
      call. = FALSE
    )
  }

  our_stats <- our_deg[, c("gene", "avg_log2FC_MNG_vs_dCLN", "q_hurdle")]
  names(our_stats)[2:3] <- c("our_avg_log2FC_MNG_vs_dCLN", "our_q_hurdle")
  kolz_stats <- kolz_deg[, c("gene", "avg_log2FC_CNS_vs_LN", "q_hurdle")]
  names(kolz_stats)[2:3] <- c("kolz_avg_log2FC_CNS_vs_LN", "kolz_q_hurdle")
  selected_audit <- selected_genes |>
    dplyr::left_join(our_stats, by = "gene") |>
    dplyr::left_join(kolz_stats, by = "gene")
  write.csv(
    selected_audit,
    file.path(output_dir, "figure_2e_displayed_gene_audit.csv"),
    row.names = FALSE
  )

  count_audit <- data.frame(
    region = c("Our MNG specific", "Shared", "Kolz Th17 CNS specific"),
    n_genes = c(length(our_specific), length(shared), length(kolz_specific)),
    total_set_size = c(length(our_set), length(shared), length(kolz_set)),
    q_hurdle_threshold = q_threshold,
    minimum_positive_log2FC = minimum_log2fc,
    stringsAsFactors = FALSE
  )
  write.csv(
    count_audit,
    file.path(output_dir, "figure_2e_venn_region_counts.csv"),
    row.names = FALSE
  )

  our_color <- "#E64B35"
  kolz_color <- "#7E57C2"
  ellipse_grob <- function(cx, cy, rx, ry, color) {
    theta <- seq(0, 2 * pi, length.out = 500L)
    grid::polygonGrob(
      x = grid::unit(cx + rx * cos(theta), "npc"),
      y = grid::unit(cy + ry * sin(theta), "npc"),
      gp = grid::gpar(fill = NA, col = color, lwd = 4.2, linejoin = "round")
    )
  }
  draw_rows <- function(rows, x, y_top, step, core = TRUE) {
    for (i in seq_along(rows)) {
      grid::grid.text(
        rows[[i]],
        x = grid::unit(x, "npc"),
        y = grid::unit(y_top - (i - 1L) * step, "npc"),
        gp = grid::gpar(
          fontsize = if (core) 20.5 else 19.5,
          fontface = if (core) "bold.italic" else "italic",
          col = "black"
        )
      )
    }
  }

  figure_tiff <- file.path(
    figure_output_dir,
    "figure_2e_MNG_vs_Kolz_Th17_CNS_DEG_venn.tiff"
  )
  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_tiff(
      figure_tiff,
      width = 10,
      height = 8.8,
      units = "in",
      res = 300,
      compression = "lzw",
      background = "white"
    )
  } else {
    grDevices::tiff(
      figure_tiff,
      width = 10,
      height = 8.8,
      units = "in",
      res = 300,
      bg = "white"
    )
  }
  grid::grid.newpage()
  grid::grid.text(
    "Up-regulated genes",
    x = unit(0.50, "npc"), y = unit(0.975, "npc"),
    gp = grid::gpar(fontsize = 25, fontface = "plain")
  )
  grid::grid.lines(
    x = unit(c(0.285, 0.715), "npc"), y = unit(c(0.935, 0.935), "npc"),
    gp = grid::gpar(col = "black", lwd = 3.8, lineend = "butt")
  )
  grid::grid.text(
    "Th1 EAE",
    x = unit(0.355, "npc"), y = unit(0.865, "npc"),
    gp = grid::gpar(fontsize = 22)
  )
  grid::grid.text(
    "Th17 EAE\n(Kolz et al.)",
    x = unit(0.645, "npc"), y = unit(0.875, "npc"),
    gp = grid::gpar(fontsize = 21, lineheight = 0.95)
  )
  grid::grid.draw(ellipse_grob(0.36, 0.47, 0.315, 0.355, our_color))
  grid::grid.draw(ellipse_grob(0.64, 0.47, 0.315, 0.355, kolz_color))

  # Bold-italic core pathway genes are grouped above the secondary examples.
  draw_rows(
    c("Ndufa6  Ndufv3", "Psmb9  Psmb10", "Tap1  H2-K1", "H2-D1  Psmb8"),
    x = 0.220, y_top = 0.675, step = 0.055, core = TRUE
  )
  draw_rows(
    c("Gpx4  Aldoa", "Eno1"),
    x = 0.205, y_top = 0.375, step = 0.058, core = FALSE
  )
  draw_rows(
    c("Mthfd2  Ass1", "Shmt2  Psat1", "Phgdh  Asl", "Slc7a5  Sat1"),
    x = 0.500, y_top = 0.675, step = 0.055, core = TRUE
  )
  draw_rows(
    c("Tpi1  Pgk1"),
    x = 0.500, y_top = 0.375, step = 0.058, core = FALSE
  )
  draw_rows(
    c("Traf6  Nfkbia", "Ccl4  Ccl5", "Cxcl10  Ccr5", "Ldha  Pgam1", "Gpi1  Hk2"),
    x = 0.795, y_top = 0.675, step = 0.055, core = TRUE
  )
  draw_rows(
    c("Il2ra  Prdm1"),
    x = 0.795, y_top = 0.325, step = 0.058, core = FALSE
  )

  grDevices::dev.off()

  cat("Saved Figure 2e:\n", figure_tiff, "\n", sep = "")
  list(
    our_set = our_set,
    kolz_set = kolz_set,
    shared = shared,
    selected_gene_audit = selected_audit,
    region_counts = count_audit,
    figure_tiff = figure_tiff
  )
})
}

################################################################################
# Figure 2e: k=6 cNMF program usage and pseudobulk heatmaps by sample
################################################################################

# Uses the stabilized 12,942-cell validated-B-cell cNMF solution generated
# after excluding structural ribosomal, TCR-locus, and all Ig-locus genes.
# GEP6 is deliberately hidden in this main Figure 2e version. The five violin
# plotting areas are square and equal-width; each heatmap is aligned to the
# exact width of its violin. Heatmap x-axis labels are hidden but ticks remain.

figure_2e_cnmf_usage <- local({
  previous_hide_gep6 <- Sys.getenv("HIDE_GEP6", unset = NA_character_)
  on.exit({
    if (is.na(previous_hide_gep6)) {
      Sys.unsetenv("HIDE_GEP6")
    } else {
      Sys.setenv(HIDE_GEP6 = previous_hide_gep6)
    }
  }, add = TRUE)
  Sys.setenv(HIDE_GEP6 = "true")

  suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(cowplot)
    library(scales)
    library(grid)
  })

  # Adapted from the existing local violin + pseudobulk heatmap template.
  # Input: stabilized validated B-cell cNMF, k=6, density threshold=0.1,
  #        with structural ribosomal, TCR-locus, and all Ig-locus genes excluded.
  # Groups: M1, M2, L1, L2.

  base_dir <- "/Users/shingheimok/Desktop/Phd_Study/research/Wu_lab/b_cells"
  run_dir <- file.path(
    base_dir,
    "outs/output/cnmf_stable_bcells_no_ribo_tcr_ig_k6_k15_hvg5000",
    "stable_bcells_no_ribo_tcr_ig_k6_k15_hvg5000"
  )
  data_dir <- file.path(run_dir, "k6_density_0_1_sample_panels")
  figure_dir <- file.path(base_dir, "figures/figure_2")
  dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

  prep_script <- file.path(base_dir, "prepare_k6_program_usage_by_sample.py")
  usage_file <- file.path(data_dir, "k6_dt01_program_usage_by_cell_and_sample.csv")
  heat_file <- file.path(data_dir, "k6_dt01_pseudobulk_heatmap_by_sample.csv")
  sample_count_file <- file.path(data_dir, "k6_dt01_sample_cell_counts.csv")
  usage_summary_file <- file.path(data_dir, "k6_dt01_program_usage_summary_by_sample.csv")
  hide_gep6 <- tolower(Sys.getenv("HIDE_GEP6", unset = "false")) %in%
    c("1", "true", "yes")
  out_file <- file.path(
    figure_dir,
    if (hide_gep6) {
      "k6_program_usage_violin_pseudobulk_M1_M2_L1_L2_GEP6_hidden.tiff"
    } else {
      "k6_program_usage_violin_pseudobulk_M1_M2_L1_L2.tiff"
    }
  )

  if (!file.exists(prep_script)) {
    stop("Missing preparation script: ", prep_script, call. = FALSE)
  }

  required_tables <- c(usage_file, heat_file, sample_count_file, usage_summary_file)
  if (any(!file.exists(required_tables))) {
    status <- system2(Sys.which("python3"), prep_script)
    if (!identical(status, 0L)) {
      stop("Python preparation step failed with status ", status, call. = FALSE)
    }
  }

  usage_df <- read.csv(usage_file, stringsAsFactors = FALSE)
  heat_df <- read.csv(heat_file, stringsAsFactors = FALSE)
  sample_counts <- read.csv(sample_count_file, stringsAsFactors = FALSE)

  sample_levels <- c("M1", "M2", "L1", "L2")
  program_levels <- if (hide_gep6) paste0("GEP", 1:5) else paste0("GEP", 1:6)

  expected_counts <- c(M1 = 3492L, M2 = 916L, L1 = 3514L, L2 = 5020L)
  observed_counts <- setNames(sample_counts$n_cells, sample_counts$sample)
  if (!identical(as.integer(observed_counts[sample_levels]), as.integer(expected_counts))) {
    stop(
      "Sample cell counts differ from the expected stabilized B-cell input.",
      call. = FALSE
    )
  }

  program_titles <- c(
    GEP1 = "MAPK pathway/\ntrafficking",
    GEP2 = "BCR signaling/\nmitochondrial",
    GEP3 = "MHC-I antigen processing/\nglycolysis",
    GEP4 = "Amino acid\nbiosynthesis",
    GEP5 = "Interferon\nresponse",
    GEP6 = "Granulocytic/\ncell-cycle"
  )

  # Robust display limits retain the distributional shape without allowing a
  # handful of extreme cells to flatten the low-usage GEP5/GEP6 violins.
  program_y_limits <- list(
    GEP1 = c(0, 1.15),
    GEP2 = c(0, 1.00),
    GEP3 = c(0, 1.00),
    GEP4 = c(0, 1.10),
    GEP5 = c(0, 0.50),
    GEP6 = c(0, 0.25)
  )

  sample_cols <- c(
    M1 = "#E64B35FF",
    M2 = "#F39B7FFF",
    L1 = "#00A087FF",
    L2 = "#4DBBD5FF"
  )

  usage_df <- usage_df |>
    filter(program %in% program_levels) |>
    mutate(
      sample = factor(sample, levels = sample_levels),
      program = factor(program, levels = program_levels)
    )

  heat_df <- heat_df |>
    filter(program %in% program_levels) |>
    mutate(
      sample = factor(sample, levels = sample_levels),
      program = factor(program, levels = program_levels),
      gene_label = sub("^mt-", "mt_", gene)
    )

  if (anyNA(usage_df$sample) || anyNA(usage_df$program)) {
    stop("Unexpected sample or program label in usage table.", call. = FALSE)
  }

  plot_list <- vector("list", length(program_levels))
  top_plot_list <- vector("list", length(program_levels))
  bottom_plot_list <- vector("list", length(program_levels))
  names(plot_list) <- program_levels
  names(top_plot_list) <- program_levels
  names(bottom_plot_list) <- program_levels
  heat_legend_source <- NULL

  for (program_id in program_levels) {
    violin_program <- usage_df |>
      filter(program == program_id)

    heat_program <- heat_df |>
      filter(program == program_id) |>
      arrange(rank, sample) |>
      mutate(
        gene_label = factor(
          gene_label,
          levels = rev(unique(gene_label[order(rank)]))
        )
      )

    y_limits <- program_y_limits[[program_id]]
    y_breaks <- pretty(y_limits, n = 4)
    y_breaks <- y_breaks[y_breaks >= y_limits[1] & y_breaks <= y_limits[2]]

    p_top <- ggplot(
      violin_program,
      aes(x = sample, y = usage, fill = sample, group = sample)
    ) +
      geom_violin(
        width = 0.82,
        scale = "width",
        color = "#2f2f2f",
        linewidth = 0.25,
        trim = TRUE
      ) +
      scale_fill_manual(values = sample_cols, drop = FALSE, guide = "none") +
      scale_x_discrete(drop = FALSE, labels = rep("", length(sample_levels))) +
      scale_y_continuous(breaks = y_breaks, expand = expansion(mult = c(0, 0.02))) +
      coord_cartesian(ylim = y_limits) +
      labs(
        title = program_id,
        subtitle = unname(program_titles[program_id]),
        x = NULL,
        y = if (program_id == "GEP1") "Usage" else NULL
      ) +
      theme_bw(base_size = 12) +
      theme(
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1.6),
        axis.line = element_blank(),
        plot.title = element_text(
          size = 18, face = "bold", hjust = 0.5, margin = margin(b = 2)
        ),
        plot.subtitle = element_text(
          size = 13.3, face = "bold", hjust = 0.5,
          lineheight = 0.88, margin = margin(b = 5)
        ),
        axis.text.x = element_blank(),
        axis.ticks.x = element_line(color = "black", linewidth = 0.7),
        axis.ticks.y = element_line(color = "black", linewidth = 0.7),
        axis.ticks.length = unit(2.8, "pt"),
        axis.title.y = if (program_id == "GEP1") {
          element_text(size = 19, face = "bold", margin = margin(r = 2))
        } else {
          element_blank()
        },
        axis.text.y = element_text(size = 13.5, face = "bold"),
        aspect.ratio = 1,
        plot.margin = margin(0, 2, -2, 2)
      )

    p_bottom <- ggplot(heat_program, aes(x = sample, y = gene_label)) +
      geom_tile(
        aes(fill = z_expression),
        width = 1, height = 1, color = NA, linewidth = 0
      ) +
      geom_point(
        aes(size = percent_expressing),
        shape = 21, stroke = 0.45, color = "white", fill = NA
      ) +
      scale_fill_gradientn(
        colours = c(
          "#0b305c", "#3d6fa6", "#b9d7ee", "#ffffff",
          "#f3c3bc", "#bf3042", "#660321"
        ),
        values = scales::rescale(c(-1.5, -0.9, -0.3, 0, 0.3, 0.9, 1.5)),
        limits = c(-1.5, 1.5),
        oob = scales::squish,
        name = "Expression level\n(z-normal pseudobulk\nlog2(TP10K))"
      ) +
      scale_size_area(
        max_size = 6,
        breaks = c(100, 75, 50, 25, 0),
        limits = c(0, 100),
        name = "Expressing\ncells (%)"
      ) +
      scale_x_discrete(drop = FALSE, expand = c(0, 0)) +
      scale_y_discrete(expand = c(0, 0)) +
      labs(x = NULL, y = NULL) +
      guides(
        size = guide_legend(order = 1),
        fill = guide_colorbar(order = 2)
      ) +
      theme_bw(base_size = 11) +
      theme(
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1.6),
        axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_line(color = "black", linewidth = 0.7),
        axis.ticks.y = element_line(color = "black", linewidth = 0.7),
        axis.ticks.length = unit(2.8, "pt"),
        axis.text.y = element_text(size = 11.7, face = "bold"),
        plot.margin = margin(-4, 2, 2, 2)
      )

    top_plot_list[[program_id]] <- p_top + theme(legend.position = "none")
    bottom_plot_list[[program_id]] <- p_bottom + theme(legend.position = "none")

    if (program_id == "GEP1") {
      heat_legend_source <- p_bottom
    }
  }

  # Align all six violin panels and all six heatmaps in one operation.  Pairwise
  # alignment is insufficient here because differing y-axis labels and tick widths
  # otherwise leave each program with a different physical panel width.
  all_aligned <- cowplot::align_plots(
    plotlist = c(top_plot_list, bottom_plot_list),
    align = "v",
    axis = "lr"
  )
  aligned_top <- all_aligned[seq_along(program_levels)]
  aligned_bottom <- all_aligned[length(program_levels) + seq_along(program_levels)]

  for (i in seq_along(program_levels)) {
    plot_list[[i]] <- cowplot::plot_grid(
      aligned_top[[i]],
      cowplot::ggdraw(),
      aligned_bottom[[i]],
      ncol = 1,
      rel_heights = c(0.44, 0.045, 0.515),
      align = "v",
      axis = "lr"
    )
  }

  sample_legend_plot <- ggplot(
    data.frame(
      x = 1,
      y = seq_along(sample_levels),
      sample = factor(sample_levels, levels = sample_levels)
    ),
    aes(x = x, y = y, color = sample)
  ) +
    geom_point(size = 3.2) +
    scale_color_manual(
      values = sample_cols,
      breaks = sample_levels,
      labels = sample_levels,
      name = "Sample"
    ) +
    guides(
      color = guide_legend(
        override.aes = list(size = 6, shape = 16),
        order = 1
      )
    ) +
    theme_void(base_size = 10) +
    theme(
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 17),
      legend.text = element_text(size = 15, face = "bold"),
      legend.key.size = unit(0.42, "in"),
      legend.spacing.y = unit(0.08, "in"),
      legend.box.margin = margin(0, 0, 0, 0)
    )

  sample_legend <- cowplot::get_legend(sample_legend_plot)

  size_legend <- cowplot::get_legend(
    heat_legend_source +
      guides(
        fill = "none",
        size = guide_legend(
          order = 1,
          override.aes = list(
            shape = 16,
            color = "black",
            fill = "black",
            stroke = 0,
            size = c(8.4, 6.8, 5.2, 3.6, 1.8)
          )
        )
      ) +
      theme(
        legend.position = "right",
        legend.title = element_text(face = "bold", size = 17),
        legend.text = element_text(size = 15),
        legend.key.size = unit(0.46, "in"),
        legend.spacing.y = unit(0.08, "in"),
        legend.box.margin = margin(0, 0, 0, 0)
      )
  )

  heat_legend <- cowplot::get_legend(
    heat_legend_source +
      guides(
        size = "none",
        fill = guide_colorbar(
          order = 2,
          frame.colour = "black",
          frame.linewidth = 0.8,
          ticks.colour = "black",
          ticks.linewidth = 0.8,
          barwidth = unit(0.34, "in"),
          barheight = unit(1.55, "in")
        )
      ) +
      theme(
        legend.position = "right",
        legend.title = element_text(face = "bold", size = 17),
        legend.text = element_text(size = 15),
        legend.key.size = unit(0.45, "in"),
        legend.background = element_blank(),
        legend.box.margin = margin(2, 2, 2, 2)
      )
  )

  legend_col <- cowplot::plot_grid(
    sample_legend,
    size_legend,
    heat_legend,
    ncol = 1,
    rel_heights = c(0.31, 0.24, 0.45),
    align = "v"
  )

  panel_sequence <- list()
  for (i in seq_along(plot_list)) {
    panel_sequence[[length(panel_sequence) + 1L]] <- plot_list[[i]]
    if (i < length(plot_list)) {
      panel_sequence[[length(panel_sequence) + 1L]] <- cowplot::ggdraw()
    }
  }

  panel_rel_widths <- c(rep(c(1, 0.060), length(plot_list) - 1L), 1)
  p_panels <- cowplot::plot_grid(
    plotlist = panel_sequence,
    nrow = 1,
    rel_widths = panel_rel_widths,
    align = "h",
    axis = "tb"
  )

  p_final <- cowplot::plot_grid(
    p_panels,
    legend_col,
    nrow = 1,
    rel_widths = c(1, 0.135),
    align = "h"
  )

  ggsave(
    filename = out_file,
    plot = p_final,
    device = "tiff",
    width = if (hide_gep6) 24 else 28,
    height = 8.8,
    units = "in",
    dpi = 300,
    compression = "lzw",
    bg = "white"
  )

  cat("\nSaved 300-dpi TIFF:\n", out_file, "\n", sep = "")
  cat("Sample counts:\n")
  print(sample_counts, row.names = FALSE)

  list(
    plot = p_final,
    figure_tiff = out_file,
    sample_counts = sample_counts,
    displayed_programs = program_levels
  )
})

################################################################################
# Disabled legacy panel: independent no-immunoglobulin UMAPs of four datasets
################################################################################

# Disabled: retain this complete four-dataset UMAP workflow for reference, but
# do not execute or regenerate the panel when figure_2.R is run.
if (FALSE) {

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

figure_2c <- cowplot::plot_grid(
  p_mng, p_dcln,
  p_th1, p_th17,
  ncol = 2,
  align = "hv",
  axis = "tblr",
  rel_widths = c(1, 1),
  rel_heights = c(1, 1)
)

tiff_file <- file.path(output_dir, "figure_2c_stable_dbscan_umaps.tiff")
png_file <- file.path(output_dir, "figure_2c_stable_dbscan_umaps.png")

ggsave(
  tiff_file,
  figure_2c,
  width = 9.6,
  height = 9.4,
  dpi = 300,
  compression = "lzw",
  bg = "white"
)
ggsave(
  png_file,
  figure_2c,
  width = 9.6,
  height = 9.4,
  dpi = 300,
  bg = "white"
)

cat("Saved Figure 2c:\n", tiff_file, "\n", png_file, "\n", sep = "")
}
