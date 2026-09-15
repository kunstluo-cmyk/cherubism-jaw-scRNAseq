root <- Sys.getenv("CB_REVISION_ROOT")
if (!nzchar(root)) root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
Sys.setenv(CB_REVISION_ROOT = root)
project_lib <- file.path(root, "renv", "library", "R-4.4", "x86_64-w64-mingw32")
.libPaths(unique(c(project_lib, .Library, .Library.site)))
source(file.path(root, "R", "00_setup.R"), chdir = FALSE)

log_step("Focused revalidation: LesionStromal versus Act OB-Stromal")
dir.create(file.path(root, "results", "annotation"), recursive = TRUE, showWarnings = FALSE)
base <- readRDS(file.path(root, "objects", "06_base_broad_lineages_reviewed.rds"))

target_states <- c("LesionStromal", "Act OB-Stromal")
cells <- colnames(base)[
  base$sample_id == "CB" & base$celltype_reviewed_provisional %in% target_states
]
corrected_counts <- SeuratObject::LayerData(base, assay = "RNA", layer = "counts")[, cells, drop = FALSE]
meta <- base@meta.data[cells, , drop = FALSE]

stromal <- Seurat::CreateSeuratObject(
  counts = corrected_counts,
  assay = "RNA",
  meta.data = meta,
  project = "CB_lesion_vs_act_ob"
)
stromal <- Seurat::NormalizeData(stromal, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
stromal <- Seurat::FindVariableFeatures(stromal, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
stromal <- Seurat::ScaleData(stromal, features = Seurat::VariableFeatures(stromal), verbose = FALSE)
stromal <- Seurat::RunPCA(stromal, features = Seurat::VariableFeatures(stromal), npcs = 30, verbose = FALSE)
stromal <- Seurat::FindNeighbors(stromal, reduction = "pca", dims = 1:20, verbose = FALSE)
stromal <- Seurat::FindClusters(stromal, resolution = 0.4, algorithm = 1, random.seed = seed, verbose = FALSE)
stromal <- Seurat::RunUMAP(
  stromal, reduction = "pca", dims = 1:20, seed.use = seed,
  reduction.name = "umap.stromal", verbose = FALSE
)

labels <- as.character(stromal$celltype_reviewed_provisional)
clusters <- as.character(stromal$seurat_clusters)
contingency <- as.data.frame.matrix(table(cluster = clusters, state = labels)) |>
  tibble::rownames_to_column("stromal_cluster") |>
  dplyr::mutate(cluster_total = rowSums(dplyr::across(dplyr::all_of(target_states))))
readr::write_csv(
  contingency,
  file.path(root, "results", "annotation", "stromal_recluster_contingency.csv")
)

nmi_score <- function(x, y) {
  tab <- table(x, y)
  joint <- tab / sum(tab)
  px <- rowSums(joint)
  py <- colSums(joint)
  expected <- outer(px, py)
  nz <- joint > 0
  mi <- sum(joint[nz] * log(joint[nz] / expected[nz]))
  hx <- -sum(px[px > 0] * log(px[px > 0]))
  hy <- -sum(py[py > 0] * log(py[py > 0]))
  if (hx == 0 || hy == 0) return(NA_real_)
  mi / sqrt(hx * hy)
}

data_mat <- SeuratObject::LayerData(stromal, assay = "RNA", layer = "data")
count_mat <- SeuratObject::LayerData(stromal, assay = "RNA", layer = "counts")
act <- which(labels == "Act OB-Stromal")
lesion <- which(labels == "LesionStromal")
mean_act <- Matrix::rowMeans(data_mat[, act, drop = FALSE])
mean_lesion <- Matrix::rowMeans(data_mat[, lesion, drop = FALSE])
pct_act <- Matrix::rowMeans(count_mat[, act, drop = FALSE] > 0)
pct_lesion <- Matrix::rowMeans(count_mat[, lesion, drop = FALSE] > 0)

direct <- data.frame(
  gene = rownames(data_mat),
  n_act_ob_stromal = length(act),
  n_lesion_stromal = length(lesion),
  mean_act_ob_stromal = mean_act,
  mean_lesion_stromal = mean_lesion,
  delta_act_minus_lesion = mean_act - mean_lesion,
  pct_act_ob_stromal = pct_act,
  pct_lesion_stromal = pct_lesion,
  pct_difference_act_minus_lesion = pct_act - pct_lesion,
  stringsAsFactors = FALSE
)
direct$ranking_score <- direct$delta_act_minus_lesion *
  sqrt(abs(direct$pct_difference_act_minus_lesion)) *
  sqrt(pmax(direct$pct_act_ob_stromal, direct$pct_lesion_stromal))
readr::write_csv(
  direct,
  file.path(root, "results", "annotation", "act_ob_vs_lesion_direct_full_gene_scan.csv.gz")
)

nontechnical <- !grepl("^(MT-|RPL|RPS)", direct$gene)
top_act <- direct |>
  dplyr::filter(
    nontechnical, pct_act_ob_stromal >= 0.10,
    delta_act_minus_lesion > 0, pct_difference_act_minus_lesion > 0
  ) |>
  dplyr::arrange(dplyr::desc(ranking_score)) |>
  dplyr::slice_head(n = 200)
top_lesion <- direct |>
  dplyr::filter(
    nontechnical, pct_lesion_stromal >= 0.10,
    delta_act_minus_lesion < 0, pct_difference_act_minus_lesion < 0
  ) |>
  dplyr::arrange(ranking_score) |>
  dplyr::slice_head(n = 200)

known <- c(
  "WNT5A", "FAP", "TNFSF11", "CSF1", "CXCL12", "IL6", "RUNX2", "SOX5",
  "PDZRN4", "BICC1", "GPC6", "PRKG1", "RORA", "VCAN", "COL11A1",
  "MYL9", "IGFBP7", "TAGLN", "POSTN", "CTHRC1", "SFRP4", "LUM", "DCN",
  "COL1A1", "COL1A2", "COL3A1", "ALPL", "SP7", "BGLAP", "IBSP", "SPP1"
)
candidates <- unique(c(top_act$gene, top_lesion$gene, intersect(known, rownames(data_mat))))
candidates <- intersect(candidates, rownames(data_mat))

set.seed(seed + 120)
n_balanced <- min(length(act), length(lesion))
iterations <- 200L
boot <- vapply(seq_len(iterations), function(i) {
  ai <- sample(act, n_balanced, replace = FALSE)
  li <- sample(lesion, n_balanced, replace = FALSE)
  Matrix::rowMeans(data_mat[candidates, ai, drop = FALSE]) -
    Matrix::rowMeans(data_mat[candidates, li, drop = FALSE])
}, numeric(length(candidates)))
rownames(boot) <- candidates

robust <- data.frame(
  gene = candidates,
  median_delta_act_minus_lesion = apply(boot, 1, stats::median),
  q025 = apply(boot, 1, stats::quantile, probs = 0.025),
  q975 = apply(boot, 1, stats::quantile, probs = 0.975),
  sign_consistency = pmax(rowMeans(boot > 0), rowMeans(boot < 0)),
  n_cells_per_state_per_iteration = n_balanced,
  iterations = iterations,
  stringsAsFactors = FALSE
) |>
  dplyr::left_join(direct, by = "gene") |>
  dplyr::arrange(dplyr::desc(abs(median_delta_act_minus_lesion)))
readr::write_csv(
  robust,
  file.path(root, "results", "annotation", "act_ob_vs_lesion_balanced_robustness.csv")
)

marker_evidence <- robust |>
  dplyr::filter(gene %in% known) |>
  dplyr::mutate(
    favored_state = dplyr::case_when(
      q025 > 0 ~ "Act OB-Stromal",
      q975 < 0 ~ "LesionStromal",
      TRUE ~ "Non-discriminating"
    )
  ) |>
  dplyr::arrange(favored_state, dplyr::desc(abs(median_delta_act_minus_lesion)))
readr::write_csv(
  marker_evidence,
  file.path(root, "results", "annotation", "act_ob_vs_lesion_marker_evidence.csv")
)

hvg <- intersect(Seurat::VariableFeatures(stromal), rownames(data_mat))
all_expressed <- rownames(data_mat)[(mean_act + mean_lesion) > 0]
profile_summary <- data.frame(
  metric = c(
    "all_expressed_gene_spearman_correlation",
    "variable_gene_spearman_correlation",
    "NMI_recluster_vs_state"
  ),
  value = c(
    stats::cor(mean_act[all_expressed], mean_lesion[all_expressed], method = "spearman"),
    stats::cor(mean_act[hvg], mean_lesion[hvg], method = "spearman"),
    nmi_score(labels, clusters)
  )
)
readr::write_csv(
  profile_summary,
  file.path(root, "results", "annotation", "stromal_state_separation_metrics.csv")
)

umap <- as.data.frame(Seurat::Embeddings(stromal, reduction = "umap.stromal")) |>
  tibble::rownames_to_column("cell_id") |>
  dplyr::mutate(
    state = labels[match(cell_id, colnames(stromal))],
    stromal_cluster = clusters[match(cell_id, colnames(stromal))]
  )
readr::write_csv(umap, file.path(root, "results", "annotation", "source_data_stromal_umap.csv"))

saveRDS(stromal, file.path(root, "objects", "08_cb_lesion_act_ob_reclustered.rds"), compress = "gzip")
log_step("Focused stromal-state revalidation complete")
