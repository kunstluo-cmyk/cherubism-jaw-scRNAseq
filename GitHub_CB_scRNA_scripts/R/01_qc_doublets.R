log_step("Stage 1: loading raw objects")

input_path <- normalizePath(file.path(root, config$input_rdata), winslash = "/", mustWork = TRUE)
e <- new.env(parent = emptyenv())
load(input_path, envir = e)
if (!all(c("CB", "CTRL") %in% ls(e))) stop("Input RData must contain Seurat objects named CB and CTRL.")

sample_meta <- readr::read_csv(
  file.path(root, "config", "sample_metadata.csv"), show_col_types = FALSE
)

prepare_sample <- function(obj, sample_id) {
  Seurat::DefaultAssay(obj) <- "RNA"
  if (!identical(SeuratObject::Layers(obj[["RNA"]]), "counts")) {
    obj[["RNA"]] <- SeuratObject::CreateAssay5Object(
      counts = SeuratObject::LayerData(obj, assay = "RNA", layer = "counts")
    )
  }
  row <- sample_meta[sample_meta$sample_id == sample_id, , drop = FALSE]
  if (nrow(row) != 1) stop("Missing or duplicate sample metadata for ", sample_id)

  obj$sample_id <- sample_id
  obj$cell_barcode_original <- colnames(obj)
  for (nm in setdiff(colnames(row), "sample_id")) obj[[nm]] <- row[[nm]][1]
  obj$nCount_raw <- obj$nCount_RNA
  obj$nFeature_raw <- obj$nFeature_RNA
  obj[["percent_mt_raw"]] <- Seurat::PercentageFeatureSet(obj, pattern = "^MT-")
  obj[["percent_ribo_raw"]] <- Seurat::PercentageFeatureSet(obj, pattern = "^RP[SL]")
  hb_features <- grep("^HB[ABDEGQZ]", rownames(obj), value = TRUE)
  obj[["percent_hb_raw"]] <- Seurat::PercentageFeatureSet(obj, features = hb_features)
  obj
}

cb <- prepare_sample(e$CB, "CB")
ctrl <- prepare_sample(e$CTRL, "CTRL")
if (length(intersect(colnames(cb), colnames(ctrl)))) {
  stop("CB and CTRL contain overlapping cell barcodes; explicit renaming is required.")
}
rm(e)

raw <- merge(cb, y = ctrl, merge.data = FALSE)
raw[["RNA"]] <- SeuratObject::JoinLayers(raw[["RNA"]])
raw$sample_id <- factor(raw$sample_id, levels = c("CTRL", "CB"))
rm(cb, ctrl)

qc_cfg <- config$qc
thresholds <- lapply(levels(raw$sample_id), function(sid) {
  d <- raw@meta.data[raw$sample_id == sid, , drop = FALSE]
  q <- as.numeric(qc_cfg$upper_quantile)
  upper_mad <- as.numeric(qc_cfg$upper_mad_multiplier)
  mt_mad <- as.numeric(qc_cfg$mitochondrial_mad_multiplier)
  data.frame(
    sample_id = sid,
    min_counts = as.numeric(qc_cfg$min_counts),
    min_features = as.numeric(qc_cfg$min_features),
    max_counts = min(
      stats::quantile(d$nCount_raw, q, na.rm = TRUE, names = FALSE),
      stats::median(d$nCount_raw) + upper_mad * safe_mad(d$nCount_raw)
    ),
    max_features = min(
      stats::quantile(d$nFeature_raw, q, na.rm = TRUE, names = FALSE),
      stats::median(d$nFeature_raw) + upper_mad * safe_mad(d$nFeature_raw)
    ),
    max_percent_mt = max(
      as.numeric(qc_cfg$mitochondrial_floor_percent),
      min(
        as.numeric(qc_cfg$mitochondrial_cap_percent),
        stats::median(d$percent_mt_raw) + mt_mad * safe_mad(d$percent_mt_raw)
      )
    )
  )
})
thresholds <- do.call(rbind, thresholds)
readr::write_csv(thresholds, file.path(root, "results", "qc", "qc_thresholds.csv"))

thr <- thresholds[match(as.character(raw$sample_id), thresholds$sample_id), ]
raw$qc_low_counts <- raw$nCount_raw < thr$min_counts
raw$qc_low_features <- raw$nFeature_raw < thr$min_features
raw$qc_high_counts <- raw$nCount_raw > thr$max_counts
raw$qc_high_features <- raw$nFeature_raw > thr$max_features
raw$qc_high_mito <- raw$percent_mt_raw > thr$max_percent_mt
flag_cols <- c(
  "qc_low_counts", "qc_low_features", "qc_high_counts", "qc_high_features", "qc_high_mito"
)
raw$qc_primary_pass <- rowSums(raw@meta.data[, flag_cols, drop = FALSE]) == 0
raw$qc_fail_reason <- apply(raw@meta.data[, flag_cols, drop = FALSE], 1, function(z) {
  hit <- sub("^qc_", "", flag_cols[as.logical(z)])
  if (length(hit)) paste(hit, collapse = ";") else "pass"
})

saveRDS(raw, file.path(root, "objects", "00_raw_qc_audited.rds"), compress = "gzip")
log_step("Primary QC retained ", sum(raw$qc_primary_pass), " of ", ncol(raw), " cells")

qc_input <- subset(raw, cells = colnames(raw)[raw$qc_primary_pass])
sce <- Seurat::as.SingleCellExperiment(qc_input, assay = "RNA")
set.seed(seed)
sce <- scDblFinder::scDblFinder(
  sce,
  samples = as.character(SummarizedExperiment::colData(sce)$sample_id),
  dbr = NULL,
  multiSampleMode = "split",
  BPPARAM = BiocParallel::SerialParam(progressbar = TRUE),
  verbose = TRUE
)

qc_input$scDblFinder_score <- SummarizedExperiment::colData(sce)$scDblFinder.score
qc_input$scDblFinder_class <- as.character(SummarizedExperiment::colData(sce)$scDblFinder.class)
saveRDS(qc_input, file.path(root, "objects", "01_qc_with_doublets.rds"), compress = "gzip")

singlets <- subset(qc_input, subset = scDblFinder_class == "singlet")
saveRDS(singlets, file.path(root, "objects", "02_qc_singlets_raw.rds"), compress = "gzip")

audit <- raw@meta.data
audit$cell_id <- rownames(audit)
audit$scDblFinder_score <- NA_real_
audit$scDblFinder_class <- "not_run_primary_qc_fail"
idx <- match(colnames(qc_input), audit$cell_id)
audit$scDblFinder_score[idx] <- qc_input$scDblFinder_score
audit$scDblFinder_class[idx] <- qc_input$scDblFinder_class
audit$final_qc_singlet <- audit$qc_primary_pass & audit$scDblFinder_class == "singlet"
readr::write_csv(audit, file.path(root, "results", "qc", "cell_qc_audit.csv"))

stage_counts <- rbind(
  data.frame(sample_id = levels(raw$sample_id), stage = "Raw", n_cells = as.integer(table(raw$sample_id))),
  data.frame(
    sample_id = levels(raw$sample_id), stage = "Primary QC",
    n_cells = as.integer(table(factor(raw$sample_id[raw$qc_primary_pass], levels = levels(raw$sample_id))))
  ),
  data.frame(
    sample_id = levels(raw$sample_id), stage = "Singlet",
    n_cells = as.integer(table(factor(singlets$sample_id, levels = levels(raw$sample_id))))
  )
)
raw_n <- stage_counts$n_cells[match(stage_counts$sample_id, stage_counts$sample_id[stage_counts$stage == "Raw"])]
stage_counts$retained_percent_of_raw <- 100 * stage_counts$n_cells / raw_n
readr::write_csv(stage_counts, file.path(root, "results", "qc", "qc_counts_by_stage.csv"))

doublet_summary <- qc_input@meta.data |>
  dplyr::count(sample_id, scDblFinder_class, name = "n_cells") |>
  dplyr::group_by(sample_id) |>
  dplyr::mutate(percent = 100 * n_cells / sum(n_cells)) |>
  dplyr::ungroup()
readr::write_csv(doublet_summary, file.path(root, "results", "qc", "doublet_summary.csv"))

rm(raw, qc_input, sce, singlets, audit)
gc()
log_step("Stage 1 complete")
