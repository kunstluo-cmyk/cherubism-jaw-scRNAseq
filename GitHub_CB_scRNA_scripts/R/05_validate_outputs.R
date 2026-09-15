log_step("Stage 5: validating objects and output bundle")

raw <- readRDS(file.path(root, "objects", "00_raw_qc_audited.rds"))
qc <- readRDS(file.path(root, "objects", "01_qc_with_doublets.rds"))
singlets <- readRDS(file.path(root, "objects", "02_qc_singlets_raw.rds"))
corrected <- readRDS(file.path(root, "objects", "03_decontx_corrected.rds"))
base <- readRDS(file.path(root, "objects", "04_base_unintegrated.rds"))
rpca <- readRDS(file.path(root, "objects", "05_base_rpca_visualization_only.rds"))

expected_final <- ncol(singlets)
checks <- data.frame(
  check = c(
    "raw_cell_count", "qc_subset_matches_flag", "singlets_match_calls",
    "corrected_cell_count", "base_cell_count", "rpca_cell_count",
    "corrected_and_raw_assay_cells_identical", "corrected_and_raw_assay_genes_identical",
    "base_has_unintegrated_umap", "rpca_has_integrated_reduction", "rpca_has_umap",
    "sample_metadata_complete", "contamination_finite", "rpca_warning_present"
  ),
  passed = c(
    ncol(raw) == 20655,
    ncol(qc) == sum(raw$qc_primary_pass),
    ncol(singlets) == sum(qc$scDblFinder_class == "singlet"),
    ncol(corrected) == expected_final,
    ncol(base) == expected_final,
    ncol(rpca) == expected_final,
    identical(colnames(corrected[["RNA"]]), colnames(corrected[["RNA_raw"]])),
    identical(rownames(corrected[["RNA"]]), rownames(corrected[["RNA_raw"]])),
    "umap.unintegrated" %in% names(base@reductions),
    "integrated.rpca" %in% names(rpca@reductions),
    "umap.rpca" %in% names(rpca@reductions),
    all(!is.na(corrected$sample_id)) && all(!is.na(corrected$group)),
    all(is.finite(corrected$decontX_contamination)),
    nzchar(rpca@misc$integration_notice)
  ),
  details = c(
    paste0("observed=", ncol(raw), "; expected=20655"),
    paste0("observed=", ncol(qc), "; flagged=", sum(raw$qc_primary_pass)),
    paste0("observed=", ncol(singlets), "; singlet_calls=", sum(qc$scDblFinder_class == "singlet")),
    paste0("observed=", ncol(corrected), "; expected=", expected_final),
    paste0("observed=", ncol(base), "; expected=", expected_final),
    paste0("observed=", ncol(rpca), "; expected=", expected_final),
    "cell order in RNA and RNA_raw",
    "gene order in RNA and RNA_raw",
    paste(names(base@reductions), collapse = ","),
    paste(names(rpca@reductions), collapse = ","),
    paste(names(rpca@reductions), collapse = ","),
    "non-identifying sample and group metadata present for every retained cell",
    "all DecontX estimates are finite",
    rpca@misc$integration_notice
  )
)

readr::write_csv(checks, file.path(root, "results", "qc", "validation_checks.csv"))

raw_counts <- SeuratObject::LayerData(base, assay = "RNA_raw", layer = "counts")
corrected_counts <- SeuratObject::LayerData(base, assay = "RNA", layer = "counts")
if (any(Matrix::colSums(corrected_counts) > Matrix::colSums(raw_counts) + 1e-6)) {
  stop("Corrected library size exceeds raw library size for at least one cell.")
}

focus_genes <- intersect(
  c(
    "COL1A1", "COL1A2", "COL3A1", "COL6A1", "COL6A2", "COL6A3", "POSTN",
    "CTHRC1", "SPP1", "SH3BP2", "TNF", "TLR2", "TLR4", "MYD88"
  ),
  rownames(base)
)
group <- interaction(
  base$sample_id,
  ifelse(is.na(base$legacy_celltype), "Unassigned_legacy", base$legacy_celltype),
  drop = TRUE,
  sep = "|"
)

ambient_by_group <- lapply(levels(group), function(g) {
  cells <- which(group == g)
  parts <- strsplit(g, "|", fixed = TRUE)[[1]]
  raw_lib <- Matrix::colSums(raw_counts[, cells, drop = FALSE])
  corrected_lib <- Matrix::colSums(corrected_counts[, cells, drop = FALSE])
  do.call(rbind, lapply(focus_genes, function(gene) {
    r <- as.numeric(raw_counts[gene, cells])
    c <- as.numeric(corrected_counts[gene, cells])
    data.frame(
      sample_id = parts[1],
      legacy_celltype = paste(parts[-1], collapse = "|"),
      n_cells = length(cells),
      gene = gene,
      raw_percent_detected = mean(r > 0) * 100,
      corrected_percent_detected = mean(c > 0) * 100,
      raw_mean_log_normalized = mean(log1p(10000 * r / pmax(raw_lib, 1))),
      corrected_mean_log_normalized = mean(log1p(10000 * c / pmax(corrected_lib, 1)))
    )
  }))
})
ambient_by_group <- dplyr::bind_rows(ambient_by_group)
readr::write_csv(
  ambient_by_group,
  file.path(root, "results", "qc", "ambient_focus_by_legacy_celltype.csv")
)

contamination_by_group <- base@meta.data |>
  dplyr::mutate(
    legacy_celltype = ifelse(is.na(legacy_celltype), "Unassigned_legacy", legacy_celltype)
  ) |>
  dplyr::group_by(sample_id, legacy_celltype) |>
  dplyr::summarise(
    n_cells = dplyr::n(),
    median_contamination = stats::median(decontX_contamination),
    q95_contamination = stats::quantile(decontX_contamination, 0.95),
    percent_gt_0_25 = mean(decontX_contamination > 0.25) * 100,
    .groups = "drop"
  )
readr::write_csv(
  contamination_by_group,
  file.path(root, "results", "qc", "contamination_by_legacy_celltype.csv")
)

bundle_files <- list.files(
  root,
  pattern = "[.](R|csv|yml|md|rds|svg|pdf|tiff|png|txt|lock)$",
  recursive = TRUE,
  full.names = TRUE
)
bundle_files <- bundle_files[!grepl("/renv/library/", gsub("\\\\", "/", bundle_files))]
bundle_files <- bundle_files[basename(bundle_files) != "output_manifest.csv"]
info <- file.info(bundle_files)
manifest <- data.frame(
  path = substring(normalizePath(bundle_files, winslash = "/"), nchar(root) + 2),
  size_bytes = info$size,
  md5 = unname(tools::md5sum(bundle_files))
)
readr::write_csv(manifest, file.path(root, "results", "qc", "output_manifest.csv"))

if (!all(checks$passed)) {
  stop("Validation failed: ", paste(checks$check[!checks$passed], collapse = ", "))
}

rm(raw, qc, singlets, corrected, base, rpca, raw_counts, corrected_counts)
gc()
log_step("Stage 5 complete; all ", nrow(checks), " validation checks passed")
