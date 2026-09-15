root <- Sys.getenv("CB_REVISION_ROOT")
if (!nzchar(root)) root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
Sys.setenv(CB_REVISION_ROOT = root)
project_lib <- file.path(root, "renv", "library", "R-4.4", "x86_64-w64-mingw32")
.libPaths(unique(c(project_lib, .Library, .Library.site)))
source(file.path(root, "R", "00_setup.R"), chdir = FALSE)

log_step("Reviewer-targeted expression audit")
out_dir <- file.path(root, "results", "annotation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
base <- readRDS(file.path(root, "objects", "06_base_broad_lineages_reviewed.rds"))
data_mat <- SeuratObject::LayerData(base, assay = "RNA", layer = "data")
corrected_counts <- SeuratObject::LayerData(base, assay = "RNA", layer = "counts")
raw_counts <- SeuratObject::LayerData(base, assay = "RNA_raw", layer = "counts")
meta <- base@meta.data
meta$condition <- as.character(meta$group)
meta$fine_state <- as.character(base$celltype_reviewed_provisional)

target_genes <- c(
  "SH3BP2", "TNF", "TNFRSF1A", "TNFRSF1B", "TLR2", "TLR4", "MYD88", "NFKB1", "RELA",
  "IL1B", "IL6", "CXCL8", "NLRP3", "CSF1", "TNFSF11", "TNFRSF11A", "TNFRSF11B",
  "CTSK", "ACP5", "MMP9", "SPP1", "TCIRG1", "DCSTAMP", "OCSTAMP", "TREM2",
  "ESR1", "ESR2", "CYP19A1"
)
target_genes <- intersect(target_genes, rownames(data_mat))

summarise_expression <- function(mat, genes, groups, layer_name) {
  groups <- as.character(groups)
  out <- lapply(sort(unique(groups)), function(g) {
    idx <- which(groups == g)
    data.frame(
      group = g,
      gene = genes,
      n_cells = length(idx),
      mean_expression = Matrix::rowMeans(mat[genes, idx, drop = FALSE]),
      percent_detected = 100 * Matrix::rowMeans(mat[genes, idx, drop = FALSE] > 0),
      layer = layer_name,
      stringsAsFactors = FALSE
    )
  }) |>
    dplyr::bind_rows()
  out
}

target_meta <- meta |>
  dplyr::transmute(
    cell_id = rownames(meta), condition, fine_state,
    broad_lineage = as.character(base$broad_lineage_reviewed)
  )
target_by_state <- lapply(sort(unique(target_meta$fine_state)), function(st) {
  lapply(sort(unique(target_meta$condition)), function(cond) {
    keep <- rownames(target_meta)[target_meta$fine_state == st & target_meta$condition == cond]
    if (!length(keep)) return(NULL)
    ans <- data.frame(
      fine_state = st, condition = cond, gene = target_genes,
      n_cells = length(keep),
      mean_log_normalized = Matrix::rowMeans(data_mat[target_genes, keep, drop = FALSE]),
      percent_detected = 100 * Matrix::rowMeans(corrected_counts[target_genes, keep, drop = FALSE] > 0),
      stringsAsFactors = FALSE
    )
    ans
  }) |>
    dplyr::bind_rows()
}) |>
  dplyr::bind_rows()
readr::write_csv(target_by_state, file.path(out_dir, "reviewer_target_genes_by_fine_state.csv"))

target_by_broad <- lapply(sort(unique(target_meta$broad_lineage)), function(st) {
  lapply(sort(unique(target_meta$condition)), function(cond) {
    keep <- rownames(target_meta)[target_meta$broad_lineage == st & target_meta$condition == cond]
    if (!length(keep)) return(NULL)
    data.frame(
      broad_lineage = st, condition = cond, gene = target_genes,
      n_cells = length(keep),
      mean_log_normalized = Matrix::rowMeans(data_mat[target_genes, keep, drop = FALSE]),
      percent_detected = 100 * Matrix::rowMeans(corrected_counts[target_genes, keep, drop = FALSE] > 0),
      stringsAsFactors = FALSE
    )
  }) |>
    dplyr::bind_rows()
}) |>
  dplyr::bind_rows()
readr::write_csv(target_by_broad, file.path(out_dir, "reviewer_target_genes_by_broad_lineage.csv"))

collagen_genes <- intersect(
  c("COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL5A2", "COL6A1", "COL6A2", "COL6A3", "COL11A1", "COL12A1"),
  rownames(data_mat)
)
collagen_states <- c(
  "FibroStromal", "LesionStromal", "Act OB-Stromal", "ChondroFibroStromal",
  "LEPR+ MSC-like", "Mature OB", "CNCC-OCPs", "Resident Mac", "Inflam Mono",
  "TREM2+ Mac", "OC", "Neutrophils", "Mast cells"
)
collagen_rows <- lapply(intersect(collagen_states, unique(meta$fine_state)), function(st) {
  keep <- rownames(meta)[meta$fine_state == st]
  lib_raw <- Matrix::colSums(raw_counts[, keep, drop = FALSE])
  lib_raw[lib_raw <= 0] <- 1
  raw_norm <- log1p(raw_counts[collagen_genes, keep, drop = FALSE] %*% Matrix::Diagonal(x = 10000 / lib_raw))
  data.frame(
    fine_state = st,
    gene = collagen_genes,
    n_cells = length(keep),
    corrected_mean_log_normalized = Matrix::rowMeans(data_mat[collagen_genes, keep, drop = FALSE]),
    corrected_percent_detected = 100 * Matrix::rowMeans(corrected_counts[collagen_genes, keep, drop = FALSE] > 0),
    raw_mean_log_normalized = Matrix::rowMeans(raw_norm),
    raw_percent_detected = 100 * Matrix::rowMeans(raw_counts[collagen_genes, keep, drop = FALSE] > 0),
    stringsAsFactors = FALSE
  )
}) |>
  dplyr::bind_rows()
readr::write_csv(collagen_rows, file.path(out_dir, "reviewer_collagen_expression_ambient_audit.csv"))

myeloid_states <- c("Resident Mac", "Inflam Mono", "TREM2+ Mac", "OC", "Neutrophils", "Mast cells", "pDC")
myeloid_composition <- meta |>
  dplyr::filter(fine_state %in% myeloid_states) |>
  dplyr::count(condition, fine_state, name = "n_cells") |>
  dplyr::group_by(condition) |>
  dplyr::mutate(myeloid_total = sum(n_cells), percent_within_myeloid = 100 * n_cells / myeloid_total) |>
  dplyr::ungroup()
readr::write_csv(myeloid_composition, file.path(out_dir, "reviewer_myeloid_composition_by_sample.csv"))

summary_lines <- data.frame(
  metric = c(
    "SH3BP2_available_in_matrix", "target_gene_count", "collagen_gene_count",
    "myeloid_states_in_comparison", "CB_cells", "CTRL_cells",
    "CB_myeloid_cells", "CTRL_myeloid_cells"
  ),
  value = c(
    "SH3BP2" %in% rownames(data_mat), length(target_genes), length(collagen_genes),
    length(myeloid_states), nrow(meta[meta$condition == "CB", ]), nrow(meta[meta$condition == "CTRL", ]),
    sum(meta$condition == "CB" & meta$fine_state %in% myeloid_states),
    sum(meta$condition == "CTRL" & meta$fine_state %in% myeloid_states)
  )
)
readr::write_csv(summary_lines, file.path(out_dir, "reviewer_targeted_audit_summary.csv"))
log_step("Reviewer-targeted expression audit complete")
