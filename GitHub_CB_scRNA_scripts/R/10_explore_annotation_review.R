root <- Sys.getenv("CB_REVISION_ROOT")
if (!nzchar(root)) root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
Sys.setenv(CB_REVISION_ROOT = root)
project_lib <- file.path(root, "renv", "library", "R-4.4", "x86_64-w64-mingw32")
.libPaths(unique(c(project_lib, .Library, .Library.site)))
source(file.path(root, "R", "00_setup.R"), chdir = FALSE)

log_step("Annotation review: exploratory evidence tables")
base <- readRDS(file.path(root, "objects", "04_base_unintegrated.rds"))
rpca <- readRDS(file.path(root, "objects", "05_base_rpca_visualization_only.rds"))
base$cluster_rpca_visualization <- as.character(rpca$cluster_rpca_visualization[colnames(base)])

legacy_to_broad <- c(
  "FibroStromal" = "Mesenchymal",
  "LesionStromal" = "Mesenchymal",
  "Pericytes" = "Mesenchymal",
  "LEPR+ MSC-like" = "Mesenchymal",
  "VSMC" = "Mesenchymal",
  "ChondroFibroStromal" = "Mesenchymal",
  "Mature OB" = "Mesenchymal",
  "Act OB-Stromal" = "Mesenchymal",
  "CNCC-OCPs" = "Mesenchymal",
  "Resident Mac" = "Myeloid",
  "Inflam Mono" = "Myeloid",
  "TREM2+ Mac" = "Myeloid",
  "OC" = "Myeloid",
  "Neutrophils" = "Myeloid",
  "Mast cells" = "Myeloid",
  "T/NK" = "Lymphoid",
  "B cells" = "Lymphoid",
  "Endothelial" = "Endothelial",
  "Schwann" = "Neural"
)
base$legacy_broad <- unname(legacy_to_broad[as.character(base$legacy_celltype)])

lineage_markers <- list(
  Mesenchymal = c(
    "COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "COL6A1", "PDGFRA",
    "LEPR", "CXCL12", "RUNX2", "SP7", "ALPL", "BGLAP", "RGS5", "ACTA2"
  ),
  Myeloid = c(
    "PTPRC", "LST1", "TYROBP", "FCER1G", "AIF1", "CTSS", "LYZ", "CTSK",
    "ACP5", "MMP9", "CSF1R"
  ),
  Lymphoid = c("PTPRC", "CD3D", "CD3E", "TRBC1", "NKG7", "GNLY", "MS4A1", "CD79A", "CD74"),
  Endothelial = c("PECAM1", "VWF", "EMCN", "KDR", "CLDN5", "RAMP2", "ENG", "ESAM"),
  Neural = c("S100B", "SOX10", "PLP1", "MPZ", "SLC1A3", "NGFR", "GFAP")
)

data_mat <- SeuratObject::LayerData(base, assay = "RNA", layer = "data")
for (nm in names(lineage_markers)) {
  genes <- intersect(lineage_markers[[nm]], rownames(data_mat))
  raw_score <- Matrix::colMeans(data_mat[genes, , drop = FALSE])
  z <- as.numeric(scale(raw_score))
  z[!is.finite(z)] <- 0
  base[[paste0("score_", nm)]] <- z
}

score_cols <- paste0("score_", names(lineage_markers))
cluster_scores <- base@meta.data |>
  dplyr::group_by(cluster_rpca_visualization) |>
  dplyr::summarise(
    dplyr::across(dplyr::all_of(score_cols), mean),
    n_cells = dplyr::n(),
    .groups = "drop"
  )
score_matrix <- as.matrix(cluster_scores[, score_cols])
cluster_scores$marker_top_lineage <- names(lineage_markers)[max.col(score_matrix, ties.method = "first")]
sorted_scores <- t(apply(score_matrix, 1, sort, decreasing = TRUE))
cluster_scores$marker_score_margin <- sorted_scores[, 1] - sorted_scores[, 2]

legacy_counts <- base@meta.data |>
  dplyr::filter(!is.na(legacy_broad)) |>
  dplyr::count(cluster_rpca_visualization, legacy_broad, name = "n_legacy") |>
  dplyr::group_by(cluster_rpca_visualization) |>
  dplyr::mutate(legacy_fraction = n_legacy / sum(n_legacy)) |>
  dplyr::arrange(dplyr::desc(n_legacy), .by_group = TRUE) |>
  dplyr::slice_head(n = 1) |>
  dplyr::ungroup()
cluster_evidence <- dplyr::left_join(cluster_scores, legacy_counts, by = "cluster_rpca_visualization")
cluster_evidence$reviewed_broad <- ifelse(
  !is.na(cluster_evidence$legacy_fraction) & cluster_evidence$legacy_fraction >= 0.60,
  cluster_evidence$legacy_broad,
  cluster_evidence$marker_top_lineage
)
cluster_evidence$evidence_status <- ifelse(
  !is.na(cluster_evidence$legacy_broad) &
    cluster_evidence$legacy_broad == cluster_evidence$marker_top_lineage,
  "legacy_and_marker_agree",
  ifelse(
    !is.na(cluster_evidence$legacy_fraction) & cluster_evidence$legacy_fraction >= 0.85,
    "high_purity_legacy_marker_review_needed",
    "manual_review_required"
  )
)
readr::write_csv(
  cluster_evidence,
  file.path(root, "results", "qc", "broad_lineage_cluster_evidence.csv")
)

fine_purity <- base@meta.data |>
  dplyr::filter(!is.na(legacy_celltype)) |>
  dplyr::count(cluster_rpca_visualization, legacy_celltype, name = "n_cells") |>
  dplyr::group_by(cluster_rpca_visualization) |>
  dplyr::mutate(fraction = n_cells / sum(n_cells)) |>
  dplyr::arrange(dplyr::desc(n_cells), .by_group = TRUE) |>
  dplyr::ungroup()
readr::write_csv(
  fine_purity,
  file.path(root, "results", "qc", "legacy_fine_annotation_by_rpca_cluster.csv")
)

base$broad_lineage_reviewed <- cluster_evidence$reviewed_broad[
  match(base$cluster_rpca_visualization, cluster_evidence$cluster_rpca_visualization)
]

Seurat::Idents(base) <- factor(
  base$broad_lineage_reviewed,
  levels = c("Mesenchymal", "Myeloid", "Lymphoid", "Endothelial", "Neural")
)
broad_markers <- Seurat::FindAllMarkers(
  base,
  assay = "RNA",
  only.pos = TRUE,
  min.pct = 0.20,
  logfc.threshold = 0.25,
  test.use = "wilcox",
  verbose = FALSE
)
readr::write_csv(
  broad_markers,
  file.path(root, "results", "qc", "broad_lineage_markers_exploratory.csv")
)

stromal_states <- c("LesionStromal", "Act OB-Stromal")
stromal <- subset(base, subset = legacy_celltype %in% stromal_states & sample_id == "CB")
Seurat::Idents(stromal) <- stromal$legacy_celltype
de <- Seurat::FindMarkers(
  stromal,
  ident.1 = "Act OB-Stromal",
  ident.2 = "LesionStromal",
  assay = "RNA",
  test.use = "wilcox",
  min.pct = 0.10,
  logfc.threshold = 0,
  verbose = FALSE
)
de$gene <- rownames(de)
de$pct_difference <- de$pct.1 - de$pct.2
de$comparison <- "Act_OB_Stromal_vs_LesionStromal"
readr::write_csv(
  de,
  file.path(root, "results", "qc", "act_ob_vs_lesion_markers_corrected.csv")
)

known_genes <- intersect(
  c(
    "WNT5A", "FAP", "TNFSF11", "CSF1", "CXCL12", "IL6", "POSTN", "CTHRC1",
    "IBSP", "ALPL", "SPP1", "BGLAP", "RUNX2", "SP7", "COL1A1", "COL1A2",
    "COL3A1", "DCN", "LUM"
  ),
  rownames(data_mat)
)
candidate_genes <- unique(c(
  head(de$gene[order(de$avg_log2FC, decreasing = TRUE)], 75),
  head(de$gene[order(de$avg_log2FC, decreasing = FALSE)], 75),
  known_genes
))
candidate_genes <- intersect(candidate_genes, rownames(data_mat))
cells_act <- colnames(stromal)[stromal$legacy_celltype == "Act OB-Stromal"]
cells_lesion <- colnames(stromal)[stromal$legacy_celltype == "LesionStromal"]
n_balanced <- min(length(cells_act), length(cells_lesion))
set.seed(seed + 100)
boot <- replicate(100, {
  a <- sample(cells_act, n_balanced, replace = FALSE)
  l <- sample(cells_lesion, n_balanced, replace = FALSE)
  Matrix::rowMeans(data_mat[candidate_genes, a, drop = FALSE]) -
    Matrix::rowMeans(data_mat[candidate_genes, l, drop = FALSE])
})
bootstrap_summary <- data.frame(
  gene = candidate_genes,
  median_delta_log_normalized = apply(boot, 1, stats::median),
  q025 = apply(boot, 1, stats::quantile, probs = 0.025),
  q975 = apply(boot, 1, stats::quantile, probs = 0.975),
  sign_consistency = pmax(rowMeans(boot > 0), rowMeans(boot < 0)),
  n_per_group = n_balanced,
  iterations = 100
)
bootstrap_summary <- dplyr::left_join(
  bootstrap_summary,
  de[, c("gene", "avg_log2FC", "pct.1", "pct.2", "pct_difference", "p_val_adj")],
  by = "gene"
)
readr::write_csv(
  bootstrap_summary,
  file.path(root, "results", "qc", "act_ob_vs_lesion_balanced_subsampling.csv")
)

known_marker_table <- bootstrap_summary[bootstrap_summary$gene %in% known_genes, ]
readr::write_csv(
  known_marker_table,
  file.path(root, "results", "qc", "act_ob_vs_lesion_known_markers.csv")
)

log_step("Annotation exploratory tables complete")
