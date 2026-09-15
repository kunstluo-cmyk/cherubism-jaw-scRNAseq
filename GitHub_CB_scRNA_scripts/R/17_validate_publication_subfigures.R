root <- Sys.getenv("CB_REVISION_ROOT")
if (!nzchar(root)) root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
project_lib <- file.path(root, "renv", "library", "R-4.4", "x86_64-w64-mingw32")
.libPaths(unique(c(project_lib, .Library, .Library.site)))
out_dir <- file.path(root, "figures", "publication_subfigures")
report_dir <- file.path(root, "results", "annotation")

bases <- c(
  "UMAP_fine_cell_states", "UMAP_major_lineages",
  "Barplot_major_lineage_cell_counts",
  "UMAP_condition_CB_CTRL_split", "UMAP_cell_states_CB_CTRL_split",
  "Barplot_mesenchymal_mural", "Barplot_hematopoietic_endothelial_neural",
  "Barplot_stacked_cell_state_composition", "Bubbleplot_major_lineage_markers",
  "Bubbleplot_fine_cell_state_markers", "Bubbleplot_mesenchymal_mural_markers",
  "Bubbleplot_hematopoietic_endothelial_neural_markers",
  "Bubbleplot_mesenchymal_endothelial_other_markers",
  "Bubbleplot_immune_cell_markers",
  "Bubbleplot_stromal_state_discriminators",
  "FeaturePlot_10gene_2x5",
  "SH3BP2_expression_CB_CTRL_major_lineages",
  paste0("FeatureUMAP_", c("COL1A1","PDZRN4","WNT5A","CTHRC1","LUM","TREM2","ACP5","SH3BP2"))
)
exts <- c("svg","pdf","tiff","png")

rows <- lapply(bases, function(b) {
  paths <- file.path(out_dir, paste0(b, ".", exts))
  exists <- file.exists(paths)
  sizes <- rep(NA_real_, length(paths))
  sizes[exists] <- file.info(paths[exists])$size
  png_path <- file.path(out_dir, paste0(b, ".png"))
  wh <- c(NA_integer_, NA_integer_)
  if (file.exists(png_path) && requireNamespace("png", quietly=TRUE)) {
    d <- dim(png::readPNG(png_path, native=TRUE))
    wh <- c(d[2], d[1])
  }
  exact_nonimmune <- b %in% c(
    "Bubbleplot_mesenchymal_mural_markers",
    "Bubbleplot_mesenchymal_endothelial_other_markers"
  )
  exact_immune <- b %in% c(
    "Bubbleplot_hematopoietic_endothelial_neural_markers",
    "Bubbleplot_immune_cell_markers"
  )
  exact_pair <- b == "Bubbleplot_stromal_state_discriminators"
  exact_feature_grid <- b == "FeaturePlot_10gene_2x5"
  pixel_pass <- if (exact_nonimmune) {
    all(wh == c(8000,2900))
  } else if (exact_immune) {
    all(wh == c(6300,2900))
  } else if (exact_pair) {
    all(wh == c(6800,3600))
  } else if (exact_feature_grid) {
    all(wh == c(6107,2036))
  } else {
    all(wh >= c(650,450))
  }
  data.frame(
    figure=b,
    all_four_formats=all(exists),
    minimum_file_bytes=min(sizes, na.rm=TRUE),
    png_width_px=wh[1], png_height_px=wh[2],
    pass=all(exists) && all(sizes > 1000, na.rm=TRUE) && pixel_pass,
    stringsAsFactors=FALSE
  )
})
report <- do.call(rbind, rows)
readr::write_csv(report, file.path(report_dir, "publication_subfigure_validation.csv"))
if (!all(report$pass)) stop("Publication subfigure validation failed: ", paste(report$figure[!report$pass], collapse=", "))
message("Publication subfigure validation passed: ", nrow(report), " figures / ", nrow(report)*4, " files")
