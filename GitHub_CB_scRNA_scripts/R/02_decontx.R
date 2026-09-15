log_step("Stage 2: DecontX ambient-RNA estimation")

singlets <- readRDS(file.path(root, "objects", "02_qc_singlets_raw.rds"))
sample_ids <- c("CTRL", "CB")
corrected_list <- vector("list", length(sample_ids))
names(corrected_list) <- sample_ids
gene_audits <- vector("list", length(sample_ids))
contamination_summaries <- vector("list", length(sample_ids))

for (i in seq_along(sample_ids)) {
  sid <- sample_ids[i]
  log_step("DecontX sample ", sid)
  seu <- subset(singlets, subset = sample_id == sid)
  raw_counts <- SeuratObject::LayerData(seu, assay = "RNA", layer = "counts")
  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = raw_counts))

  sce <- celda::decontX(
    sce,
    maxIter = as.integer(config$ambient_rna$max_iterations),
    varGenes = as.integer(config$ambient_rna$variable_genes),
    seed = seed + i,
    logfile = file.path(root, "logs", paste0("decontX_", sid, ".log")),
    verbose = TRUE
  )

  corrected_counts <- celda::decontXcounts(sce)
  contamination <- SummarizedExperiment::colData(sce)$decontX_contamination
  decontx_cluster <- as.character(SummarizedExperiment::colData(sce)$decontX_clusters)

  meta <- seu@meta.data
  meta$nCount_raw <- seu$nCount_raw
  meta$nFeature_raw <- seu$nFeature_raw
  meta <- meta[, setdiff(colnames(meta), c("nCount_RNA", "nFeature_RNA")), drop = FALSE]
  meta$decontX_contamination <- contamination
  meta$decontX_cluster <- decontx_cluster
  meta$decontX_high_contamination <- contamination > as.numeric(config$ambient_rna$high_contamination_flag)
  meta$ambient_method <- "DecontX_filtered_cells_no_empty_droplets"

  corrected <- Seurat::CreateSeuratObject(
    counts = corrected_counts,
    assay = "RNA",
    project = paste0("CB_revision_", sid),
    meta.data = meta,
    min.cells = 0,
    min.features = 0
  )
  corrected[["RNA_raw"]] <- SeuratObject::CreateAssay5Object(counts = raw_counts)
  Seurat::DefaultAssay(corrected) <- "RNA"
  corrected_list[[sid]] <- corrected

  raw_sum <- Matrix::rowSums(raw_counts)
  corrected_sum <- Matrix::rowSums(corrected_counts)
  gene_audits[[sid]] <- data.frame(
    sample_id = sid,
    gene = rownames(raw_counts),
    raw_counts = as.numeric(raw_sum),
    corrected_counts = as.numeric(corrected_sum),
    removed_counts = as.numeric(raw_sum - corrected_sum),
    removed_fraction = ifelse(raw_sum > 0, as.numeric((raw_sum - corrected_sum) / raw_sum), NA_real_)
  )
  contamination_summaries[[sid]] <- data.frame(
    sample_id = sid,
    n_cells = length(contamination),
    median = stats::median(contamination),
    q25 = stats::quantile(contamination, 0.25, names = FALSE),
    q75 = stats::quantile(contamination, 0.75, names = FALSE),
    q95 = stats::quantile(contamination, 0.95, names = FALSE),
    percent_gt_0_25 = mean(contamination > 0.25) * 100
  )

  saveRDS(
    corrected,
    file.path(root, "objects", paste0("02a_", sid, "_decontx.rds")),
    compress = "gzip"
  )
  rm(seu, sce, raw_counts, corrected_counts, corrected)
  gc()
}

decontx <- merge(corrected_list[["CB"]], y = corrected_list[["CTRL"]], merge.data = FALSE)
decontx[["RNA"]] <- SeuratObject::JoinLayers(decontx[["RNA"]])
decontx[["RNA_raw"]] <- SeuratObject::JoinLayers(decontx[["RNA_raw"]])
decontx$sample_id <- factor(decontx$sample_id, levels = c("CTRL", "CB"))
decontx@misc$ambient_rna <- list(
  method = "DecontX",
  mode = "filtered cells without empty-droplet background",
  corrected_assay = "RNA",
  uncorrected_assay = "RNA_raw",
  warning = "Repeat with raw droplets if Cell Ranger raw matrix becomes available."
)
saveRDS(decontx, file.path(root, "objects", "03_decontx_corrected.rds"), compress = "gzip")

gene_audit <- dplyr::bind_rows(gene_audits)
top_removed <- gene_audit |>
  dplyr::group_by(sample_id) |>
  dplyr::slice_max(order_by = removed_counts, n = 100, with_ties = FALSE) |>
  dplyr::ungroup()
focus_genes <- c(
  "COL1A1", "COL1A2", "COL3A1", "COL6A1", "COL6A2", "COL6A3", "POSTN",
  "CTHRC1", "SPP1", "TNF", "SH3BP2", "TLR2", "TLR4", "MYD88"
)
focus_audit <- gene_audit[gene_audit$gene %in% focus_genes, ]
readr::write_csv(top_removed, file.path(root, "results", "qc", "decontx_top_removed_genes.csv"))
readr::write_csv(focus_audit, file.path(root, "results", "qc", "decontx_focus_genes.csv"))
readr::write_csv(
  dplyr::bind_rows(contamination_summaries),
  file.path(root, "results", "qc", "decontx_contamination_summary.csv")
)

rm(singlets, corrected_list, gene_audits, gene_audit, decontx)
gc()
log_step("Stage 2 complete")

