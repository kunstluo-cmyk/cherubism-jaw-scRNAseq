log_step("Stage 3: building primary unintegrated base object")

base <- readRDS(file.path(root, "objects", "03_decontx_corrected.rds"))
Seurat::DefaultAssay(base) <- "RNA"

legacy_path <- normalizePath(
  file.path(root, config$legacy_annotation_rdata), winslash = "/", mustWork = TRUE
)
legacy_env <- new.env(parent = emptyenv())
load(legacy_path, envir = legacy_env)
legacy_names <- ls(legacy_env)[vapply(ls(legacy_env), function(nm) {
  inherits(legacy_env[[nm]], "Seurat")
}, logical(1))]
if (length(legacy_names)) {
  legacy <- legacy_env[[legacy_names[1]]]
  if ("celltype_anno" %in% colnames(legacy@meta.data)) {
    annotation_map <- stats::setNames(as.character(legacy$celltype_anno), colnames(legacy))
    base$legacy_celltype <- unname(annotation_map[colnames(base)])
  } else {
    base$legacy_celltype <- NA_character_
  }
} else {
  base$legacy_celltype <- NA_character_
}
base$legacy_annotation_status <- ifelse(
  is.na(base$legacy_celltype), "not_in_legacy_final_object", "legacy_only_requires_revalidation"
)
rm(legacy_env)

set.seed(seed)
base <- Seurat::NormalizeData(
  base,
  normalization.method = config$normalization$method,
  scale.factor = as.numeric(config$normalization$scale_factor),
  verbose = FALSE
)
base <- Seurat::FindVariableFeatures(
  base,
  selection.method = "vst",
  nfeatures = as.integer(config$normalization$variable_features),
  verbose = FALSE
)
base <- Seurat::ScaleData(base, features = Seurat::VariableFeatures(base), verbose = FALSE)
base <- Seurat::RunPCA(
  base,
  features = Seurat::VariableFeatures(base),
  npcs = as.integer(config$dimension_reduction$npcs),
  seed.use = seed,
  verbose = FALSE
)
dims_use <- seq_len(as.integer(config$dimension_reduction$dimensions))
base <- Seurat::RunUMAP(
  base,
  reduction = "pca",
  dims = dims_use,
  reduction.name = "umap.unintegrated",
  reduction.key = "UMAPUNINT_",
  seed.use = seed,
  verbose = FALSE
)
base <- Seurat::FindNeighbors(
  base,
  reduction = "pca",
  dims = dims_use,
  graph.name = c("RNA_nn_unintegrated", "RNA_snn_unintegrated"),
  verbose = FALSE
)
base <- Seurat::FindClusters(
  base,
  graph.name = "RNA_snn_unintegrated",
  resolution = as.numeric(config$dimension_reduction$clustering_resolution),
  cluster.name = "cluster_unintegrated",
  random.seed = seed,
  verbose = FALSE
)
base@misc$analysis_scope <- list(
  biological_replicates = c(CB = 1, CTRL = 1),
  primary_expression_assay = "RNA (DecontX corrected)",
  raw_audit_assay = "RNA_raw",
  differential_expression = "Use unintegrated RNA; disease-level inference is not estimable with n=1 per group.",
  clinical_covariates = "Not used in the computational model; participant-level details are reported in the manuscript."
)
saveRDS(base, file.path(root, "objects", "04_base_unintegrated.rds"), compress = "gzip")

log_step("Building RPCA-aligned visualization-only object")
rpca <- readRDS(file.path(root, "objects", "03_decontx_corrected.rds"))
rpca@meta.data <- base@meta.data[colnames(rpca), , drop = FALSE]
Seurat::DefaultAssay(rpca) <- "RNA"
rpca[["RNA"]] <- split(rpca[["RNA"]], f = rpca$sample_id)
rpca <- Seurat::NormalizeData(
  rpca,
  normalization.method = config$normalization$method,
  scale.factor = as.numeric(config$normalization$scale_factor),
  verbose = FALSE
)
rpca <- Seurat::FindVariableFeatures(
  rpca,
  selection.method = "vst",
  nfeatures = as.integer(config$normalization$variable_features),
  verbose = FALSE
)
rpca <- Seurat::ScaleData(rpca, features = Seurat::VariableFeatures(rpca), verbose = FALSE)
rpca <- Seurat::RunPCA(
  rpca,
  features = Seurat::VariableFeatures(rpca),
  npcs = as.integer(config$dimension_reduction$npcs),
  seed.use = seed,
  verbose = FALSE
)
rpca <- Seurat::IntegrateLayers(
  object = rpca,
  method = Seurat::RPCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.rpca",
  dims = dims_use,
  verbose = FALSE
)
rpca <- Seurat::FindNeighbors(
  rpca,
  reduction = "integrated.rpca",
  dims = dims_use,
  graph.name = c("rpca_nn", "rpca_snn"),
  verbose = FALSE
)
rpca <- Seurat::FindClusters(
  rpca,
  graph.name = "rpca_snn",
  resolution = as.numeric(config$dimension_reduction$clustering_resolution),
  cluster.name = "cluster_rpca_visualization",
  random.seed = seed,
  verbose = FALSE
)
rpca <- Seurat::RunUMAP(
  rpca,
  reduction = "integrated.rpca",
  dims = dims_use,
  reduction.name = "umap.rpca",
  reduction.key = "UMAPRPCA_",
  seed.use = seed,
  verbose = FALSE
)
rpca[["RNA"]] <- SeuratObject::JoinLayers(rpca[["RNA"]])
rpca@misc$integration_notice <- paste(
  "RPCA is used only for visualization and annotation sensitivity.",
  "Because donor and disease group are perfectly confounded, it cannot identify or remove donor effects.",
  "Do not use integrated values for differential expression."
)
saveRDS(
  rpca,
  file.path(root, "objects", "05_base_rpca_visualization_only.rds"),
  compress = "gzip"
)

base_summary <- data.frame(
  object = c("04_base_unintegrated", "05_base_rpca_visualization_only"),
  n_cells = c(ncol(base), ncol(rpca)),
  n_genes = c(nrow(base), nrow(rpca)),
  reduction = c("pca + umap.unintegrated", "integrated.rpca + umap.rpca"),
  permitted_use = c(
    "primary expression analysis, markers, descriptive effect sizes",
    "visualization and annotation sensitivity only"
  )
)
readr::write_csv(base_summary, file.path(root, "results", "qc", "base_object_summary.csv"))

rm(base, rpca)
gc()
log_step("Stage 3 complete")
