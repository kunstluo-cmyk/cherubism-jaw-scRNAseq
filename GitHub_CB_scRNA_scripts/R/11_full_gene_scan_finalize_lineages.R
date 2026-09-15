root <- Sys.getenv("CB_REVISION_ROOT")
if (!nzchar(root)) root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
Sys.setenv(CB_REVISION_ROOT = root)
project_lib <- file.path(root, "renv", "library", "R-4.4", "x86_64-w64-mingw32")
.libPaths(unique(c(project_lib, .Library, .Library.site)))
source(file.path(root, "R", "00_setup.R"), chdir = FALSE)

if (!requireNamespace("RANN", quietly = TRUE)) {
  stop("RANN is required for conservative nearest-neighbour transfer.")
}

log_step("Global annotation review: six-lineage adjudication and full-gene scans")
dir.create(file.path(root, "results", "annotation"), recursive = TRUE, showWarnings = FALSE)
base <- readRDS(file.path(root, "objects", "04_base_unintegrated.rds"))
rpca <- readRDS(file.path(root, "objects", "05_base_rpca_visualization_only.rds"))
base$condition <- as.character(base$group)
base$cluster_rpca_visualization <- as.character(rpca$cluster_rpca_visualization[colnames(base)])

legacy_to_broad6 <- c(
  "FibroStromal" = "Stromal/osteogenic",
  "LesionStromal" = "Stromal/osteogenic",
  "LEPR+ MSC-like" = "Stromal/osteogenic",
  "ChondroFibroStromal" = "Stromal/osteogenic",
  "Mature OB" = "Stromal/osteogenic",
  "Act OB-Stromal" = "Stromal/osteogenic",
  "CNCC-OCPs" = "Stromal/osteogenic",
  "Pericytes" = "Mural",
  "VSMC" = "Mural",
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
base$legacy_broad6 <- unname(legacy_to_broad6[as.character(base$legacy_celltype)])

lineage_markers6 <- list(
  `Stromal/osteogenic` = c(
    "COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "COL6A1", "PDGFRA",
    "LEPR", "CXCL12", "RUNX2", "SP7", "ALPL", "BGLAP", "IBSP"
  ),
  Mural = c(
    "RGS5", "CSPG4", "MCAM", "PDGFRB", "RBP1", "DES", "ACTA2", "TAGLN",
    "MYH11", "NOTCH3", "KCNJ8", "ABCC9"
  ),
  Myeloid = c(
    "PTPRC", "LST1", "TYROBP", "FCER1G", "AIF1", "CTSS", "LYZ", "CTSK",
    "ACP5", "MMP9", "CSF1R"
  ),
  Lymphoid = c(
    "PTPRC", "CD3D", "CD3E", "TRBC1", "NKG7", "GNLY", "MS4A1", "CD79A", "CD74"
  ),
  Endothelial = c("PECAM1", "VWF", "EMCN", "KDR", "CLDN5", "RAMP2", "ENG", "ESAM"),
  Neural = c("S100B", "SOX10", "PLP1", "MPZ", "SLC1A3", "NGFR", "GFAP")
)

data_mat <- SeuratObject::LayerData(base, assay = "RNA", layer = "data")
raw_mat <- SeuratObject::LayerData(base, assay = "RNA_raw", layer = "counts")

for (nm in names(lineage_markers6)) {
  genes <- intersect(lineage_markers6[[nm]], rownames(data_mat))
  if (!length(genes)) stop("No available markers for lineage: ", nm)
  raw_score <- Matrix::colMeans(data_mat[genes, , drop = FALSE])
  z <- as.numeric(scale(raw_score))
  z[!is.finite(z)] <- 0
  base[[paste0("score6_", make.names(nm))]] <- z
}

score_cols <- paste0("score6_", make.names(names(lineage_markers6)))
cluster_scores <- base@meta.data |>
  dplyr::group_by(cluster_rpca_visualization) |>
  dplyr::summarise(
    dplyr::across(dplyr::all_of(score_cols), mean),
    n_cells = dplyr::n(),
    .groups = "drop"
  )
score_matrix <- as.matrix(cluster_scores[, score_cols])
cluster_scores$marker_top_lineage <- names(lineage_markers6)[max.col(score_matrix, ties.method = "first")]
sorted_scores <- t(apply(score_matrix, 1, sort, decreasing = TRUE))
cluster_scores$marker_score_margin <- sorted_scores[, 1] - sorted_scores[, 2]

legacy_counts <- base@meta.data |>
  dplyr::filter(!is.na(legacy_broad6)) |>
  dplyr::count(cluster_rpca_visualization, legacy_broad6, name = "n_legacy") |>
  dplyr::group_by(cluster_rpca_visualization) |>
  dplyr::mutate(legacy_fraction = n_legacy / sum(n_legacy)) |>
  dplyr::arrange(dplyr::desc(n_legacy), .by_group = TRUE) |>
  dplyr::slice_head(n = 1) |>
  dplyr::ungroup()

cluster_evidence <- dplyr::left_join(cluster_scores, legacy_counts, by = "cluster_rpca_visualization")
cluster_evidence$reviewed_broad6 <- ifelse(
  !is.na(cluster_evidence$legacy_fraction) & cluster_evidence$legacy_fraction >= 0.60,
  cluster_evidence$legacy_broad6,
  cluster_evidence$marker_top_lineage
)
cluster_evidence$evidence_status <- dplyr::case_when(
  !is.na(cluster_evidence$legacy_broad6) &
    cluster_evidence$legacy_broad6 == cluster_evidence$marker_top_lineage ~ "legacy_and_marker_agree",
  !is.na(cluster_evidence$legacy_fraction) & cluster_evidence$legacy_fraction >= 0.85 ~
    "high_purity_legacy_marker_review_needed",
  TRUE ~ "manual_review_required"
)

# Full-transcriptome adjudication: cluster 18 is a small pDC population
# (GZMB/GPR183/LILRA4/IL3RA/HLA-II), not a B-cell cluster.
pdc_cluster <- cluster_evidence$cluster_rpca_visualization == "18"
cluster_evidence$reviewed_broad6[pdc_cluster] <- "Myeloid"
cluster_evidence$evidence_status[pdc_cluster] <- "manual_override_pDC_full_gene_scan"
readr::write_csv(
  cluster_evidence,
  file.path(root, "results", "annotation", "broad6_cluster_evidence.csv")
)

cluster_map <- stats::setNames(
  cluster_evidence$reviewed_broad6,
  cluster_evidence$cluster_rpca_visualization
)
base$broad_lineage_reviewed <- unname(cluster_map[base$cluster_rpca_visualization])

# Preserve only old fine labels concordant with the reviewed broad lineage.
# All missing or discordant cells are transferred within their reviewed lineage.
legacy_fine <- as.character(base$legacy_celltype)
legacy_fine_broad <- unname(legacy_to_broad6[legacy_fine])
legacy_concordant <- !is.na(legacy_fine) & legacy_fine_broad == base$broad_lineage_reviewed
fine <- ifelse(legacy_concordant, legacy_fine, NA_character_)
fine_confidence <- ifelse(legacy_concordant, 1, NA_real_)
fine_status <- ifelse(
  legacy_concordant,
  "legacy_label_broad_revalidated",
  ifelse(is.na(legacy_fine), "pending_knn_transfer", "legacy_discordant_pending_transfer")
)
names(fine) <- names(fine_confidence) <- names(fine_status) <- colnames(base)

emb <- Seurat::Embeddings(rpca, reduction = "integrated.rpca")
emb <- emb[colnames(base), seq_len(min(30, ncol(emb))), drop = FALSE]
for (broad in names(lineage_markers6)) {
  reference <- which(legacy_concordant & base$broad_lineage_reviewed == broad)
  query <- which(is.na(fine) & base$broad_lineage_reviewed == broad)
  if (!length(query)) next
  if (!length(reference)) stop("No concordant reference cells for lineage: ", broad)
  nn <- RANN::nn2(
    data = emb[reference, , drop = FALSE],
    query = emb[query, , drop = FALSE],
    k = min(31, length(reference))
  )$nn.idx
  neighbour_labels <- matrix(legacy_fine[reference][nn], nrow = nrow(nn), ncol = ncol(nn))
  transferred <- apply(neighbour_labels, 1, function(x) {
    tab <- sort(table(x), decreasing = TRUE)
    c(label = names(tab)[1], confidence = as.numeric(tab[1]) / length(x))
  })
  transferred <- t(transferred)
  fine[query] <- transferred[, "label"]
  fine_confidence[query] <- as.numeric(transferred[, "confidence"])
  low <- fine_confidence[query] < 0.50
  fine_status[query] <- ifelse(
    low,
    "knn_transfer_within_lineage_low_confidence",
    "knn_transfer_within_lineage_supported"
  )
}

# Two small populations are not represented correctly by the legacy 19-label set.
# Cluster 14 is antibody-secreting/plasma-like (JCHAIN/MZB1/DERL3/IGH), while
# cluster 18 is pDC-like (GZMB/LILRA4/IL3RA/GPR183/HLA-II).
plasma_cells <- base$cluster_rpca_visualization == "14"
pdc_cells <- base$cluster_rpca_visualization == "18"
fine[plasma_cells] <- "Plasma cells"
fine_confidence[plasma_cells] <- 1
fine_status[plasma_cells] <- "manual_full_gene_plasma_cell_annotation"
fine[pdc_cells] <- "pDC"
fine_confidence[pdc_cells] <- 1
fine_status[pdc_cells] <- "manual_full_gene_pDC_annotation"

base$celltype_reviewed_provisional <- fine[colnames(base)]
base$fine_annotation_confidence <- as.numeric(fine_confidence[colnames(base)])
base$fine_annotation_status <- fine_status[colnames(base)]

readr::write_csv(
  base@meta.data |>
    tibble::rownames_to_column("cell_id") |>
    dplyr::select(
      cell_id, sample_id, condition, cluster_rpca_visualization,
      broad_lineage_reviewed, legacy_celltype, celltype_reviewed_provisional,
      fine_annotation_confidence, fine_annotation_status
    ),
  file.path(root, "results", "annotation", "cell_annotation_audit.csv.gz")
)

composition <- base@meta.data |>
  dplyr::mutate(
    display_group = ifelse(
      broad_lineage_reviewed %in% c("Stromal/osteogenic", "Mural"),
      "Mesenchymal/mural", "Hematopoietic/endothelial/neural"
    )
  ) |>
  dplyr::count(sample_id, condition, display_group, broad_lineage_reviewed,
               celltype_reviewed_provisional, name = "n_cells") |>
  dplyr::group_by(sample_id) |>
  dplyr::mutate(sample_total = sum(n_cells), proportion_total = n_cells / sample_total) |>
  dplyr::group_by(sample_id, display_group) |>
  dplyr::mutate(group_total = sum(n_cells), proportion_within_display_group = n_cells / group_total) |>
  dplyr::ungroup()
readr::write_csv(
  composition,
  file.path(root, "results", "annotation", "composition_reviewed_by_sample.csv")
)

full_gene_scan <- function(labels, prefix, min_cells = 3L) {
  labels <- as.character(labels)
  keep <- !is.na(labels)
  labels <- labels[keep]
  x <- data_mat[, keep, drop = FALSE]
  counts <- raw_mat[, keep, drop = FALSE]
  genes <- rownames(x)
  total_sum <- Matrix::rowSums(x)
  total_detect <- Matrix::rowSums(counts > 0)
  n_total <- ncol(x)
  groups <- sort(unique(labels))
  out_path <- file.path(root, "results", "annotation", paste0(prefix, "_full_gene_scan.csv.gz"))
  top_path <- file.path(root, "results", "annotation", paste0(prefix, "_top50_markers.csv"))
  con <- gzfile(out_path, open = "wt")
  on.exit(close(con), add = TRUE)
  top_list <- list()
  wrote_header <- FALSE
  for (g in groups) {
    idx <- which(labels == g)
    n_in <- length(idx)
    n_out <- n_total - n_in
    if (n_in < min_cells || n_out < min_cells) next
    sum_in <- Matrix::rowSums(x[, idx, drop = FALSE])
    detect_in <- Matrix::rowSums(counts[, idx, drop = FALSE] > 0)
    mean_in <- sum_in / n_in
    mean_out <- (total_sum - sum_in) / n_out
    pct_in <- detect_in / n_in
    pct_out <- (total_detect - detect_in) / n_out
    tab <- data.frame(
      group = g,
      gene = genes,
      n_in = n_in,
      n_out = n_out,
      mean_log_normalized_in = mean_in,
      mean_log_normalized_out = mean_out,
      delta_mean_log_normalized = mean_in - mean_out,
      pct_in = pct_in,
      pct_out = pct_out,
      pct_difference = pct_in - pct_out,
      stringsAsFactors = FALSE
    )
    tab$ranking_score <- tab$delta_mean_log_normalized *
      sqrt(pmax(tab$pct_difference, 0)) * sqrt(pmax(tab$pct_in, 0))
    utils::write.table(
      tab, con, sep = ",", row.names = FALSE, col.names = !wrote_header,
      append = wrote_header, quote = TRUE, qmethod = "double"
    )
    wrote_header <- TRUE
    top_list[[g]] <- tab |>
      dplyr::filter(pct_in >= 0.10, delta_mean_log_normalized > 0) |>
      dplyr::arrange(dplyr::desc(ranking_score), dplyr::desc(delta_mean_log_normalized)) |>
      dplyr::slice_head(n = 50)
    log_step("Full-gene scan ", prefix, ": ", g, " (n=", n_in, ")")
  }
  close(con)
  on.exit(NULL, add = FALSE)
  readr::write_csv(dplyr::bind_rows(top_list), top_path)
  invisible(list(full = out_path, top = top_path))
}

# These three outputs explicitly retain all genes; top-50 files are convenience views only.
full_gene_scan(base$cluster_rpca_visualization, "rpca_cluster")
full_gene_scan(base$broad_lineage_reviewed, "reviewed_broad_lineage")
full_gene_scan(base$celltype_reviewed_provisional, "reviewed_fine_celltype")

rpca$broad_lineage_reviewed <- base$broad_lineage_reviewed[colnames(rpca)]
rpca$celltype_reviewed_provisional <- base$celltype_reviewed_provisional[colnames(rpca)]
rpca$fine_annotation_confidence <- base$fine_annotation_confidence[colnames(rpca)]
rpca$fine_annotation_status <- base$fine_annotation_status[colnames(rpca)]

saveRDS(base, file.path(root, "objects", "06_base_broad_lineages_reviewed.rds"), compress = "gzip")
saveRDS(rpca, file.path(root, "objects", "07_rpca_broad_lineages_reviewed.rds"), compress = "gzip")

summary_tab <- base@meta.data |>
  dplyr::count(broad_lineage_reviewed, celltype_reviewed_provisional,
               fine_annotation_status, name = "n_cells")
readr::write_csv(summary_tab, file.path(root, "results", "annotation", "annotation_review_summary.csv"))
log_step("Global annotation review full-gene scans complete")
