root <- Sys.getenv("CB_REVISION_ROOT")
if (!nzchar(root)) root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
Sys.setenv(CB_REVISION_ROOT = root)
project_lib <- file.path(root, "renv", "library", "R-4.4", "x86_64-w64-mingw32")
.libPaths(unique(c(project_lib, .Library, .Library.site)))
source(file.path(root, "R", "00_setup.R"), chdir = FALSE)

log_step("Validation: annotation objects, full-gene evidence, and figure exports")
ann_dir <- file.path(root, "results", "annotation")
fig_dir <- file.path(root, "figures", "annotation")
audit <- readr::read_csv(file.path(ann_dir, "cell_annotation_audit.csv.gz"), show_col_types = FALSE)
top <- readr::read_csv(file.path(ann_dir, "reviewed_fine_celltype_top50_markers.csv"), show_col_types = FALSE)

fine_to_broad <- c(
  "FibroStromal" = "Stromal/osteogenic", "LesionStromal" = "Stromal/osteogenic",
  "Act OB-Stromal" = "Stromal/osteogenic", "ChondroFibroStromal" = "Stromal/osteogenic",
  "LEPR+ MSC-like" = "Stromal/osteogenic", "Mature OB" = "Stromal/osteogenic",
  "CNCC-OCPs" = "Stromal/osteogenic", "Pericytes" = "Mural", "VSMC" = "Mural",
  "Endothelial" = "Endothelial", "Resident Mac" = "Myeloid", "Inflam Mono" = "Myeloid",
  "TREM2+ Mac" = "Myeloid", "OC" = "Myeloid", "Neutrophils" = "Myeloid",
  "Mast cells" = "Myeloid", "pDC" = "Myeloid", "T/NK" = "Lymphoid",
  "B cells" = "Lymphoid", "Plasma cells" = "Lymphoid", "Schwann" = "Neural"
)

checks <- list()
add_check <- function(name, pass, observed, expected) {
  checks[[length(checks) + 1L]] <<- data.frame(
    check = name, pass = isTRUE(pass), observed = as.character(observed),
    expected = as.character(expected), stringsAsFactors = FALSE
  )
}

add_check("cell_count", nrow(audit) == 17967, nrow(audit), 17967)
add_check("cell_ids_unique", !anyDuplicated(audit$cell_id), nrow(unique(audit["cell_id"])), nrow(audit))
add_check("broad_no_missing", !anyNA(audit$broad_lineage_reviewed), sum(is.na(audit$broad_lineage_reviewed)), 0)
add_check("fine_no_missing", !anyNA(audit$celltype_reviewed_provisional), sum(is.na(audit$celltype_reviewed_provisional)), 0)
add_check("six_broad_lineages", dplyr::n_distinct(audit$broad_lineage_reviewed) == 6,
          dplyr::n_distinct(audit$broad_lineage_reviewed), 6)
add_check("twenty_one_fine_states", dplyr::n_distinct(audit$celltype_reviewed_provisional) == 21,
          dplyr::n_distinct(audit$celltype_reviewed_provisional), 21)
expected_broad <- unname(fine_to_broad[audit$celltype_reviewed_provisional])
discordant <- sum(expected_broad != audit$broad_lineage_reviewed, na.rm = TRUE)
add_check("fine_broad_concordance", discordant == 0, discordant, 0)
add_check("cluster14_plasma", all(audit$celltype_reviewed_provisional[audit$cluster_rpca_visualization == 14] == "Plasma cells"),
          paste(unique(audit$celltype_reviewed_provisional[audit$cluster_rpca_visualization == 14]), collapse = ";"), "Plasma cells")
add_check("cluster18_pDC", all(audit$celltype_reviewed_provisional[audit$cluster_rpca_visualization == 18] == "pDC"),
          paste(unique(audit$celltype_reviewed_provisional[audit$cluster_rpca_visualization == 18]), collapse = ";"), "pDC")

expected_markers <- list(
  Endothelial = c("PECAM1", "VWF", "EMCN"),
  `Mature OB` = c("BGLAP", "IBSP"),
  OC = c("ACP5", "MMP9"),
  Schwann = c("S100B", "MPZ"),
  `Plasma cells` = c("JCHAIN", "MZB1"),
  pDC = c("GZMB", "LILRA4")
)
for (state in names(expected_markers)) {
  observed <- top$gene[top$group == state]
  hits <- intersect(expected_markers[[state]], observed)
  add_check(
    paste0("top50_marker_support_", make.names(state)),
    length(hits) >= 1, paste(hits, collapse = ";"),
    paste(expected_markers[[state]], collapse = " or ")
  )
}

full_files <- c(
  "rpca_cluster_full_gene_scan.csv.gz",
  "reviewed_broad_lineage_full_gene_scan.csv.gz",
  "reviewed_fine_celltype_full_gene_scan.csv.gz",
  "act_ob_vs_lesion_direct_full_gene_scan.csv.gz"
)
for (nm in full_files) {
  info <- file.info(file.path(ann_dir, nm))
  add_check(paste0("full_scan_exists_", nm), isTRUE(info$size > 1000000), info$size, ">1000000 bytes")
}

figure_bases <- c(
  "Figure1_revised_global_annotation",
  "FigureS_global_fine_marker_audit",
  "FigureS_stromal_state_revalidation"
)
for (b in figure_bases) {
  for (ext in c("svg", "pdf", "tiff", "png")) {
    p <- file.path(fig_dir, paste0(b, ".", ext))
    info <- file.info(p)
    add_check(paste0("figure_exists_", b, "_", ext), isTRUE(info$size > 10000), info$size, ">10000 bytes")
  }
}

if (requireNamespace("magick", quietly = TRUE)) {
  for (b in figure_bases) {
    info <- magick::image_info(magick::image_read(file.path(fig_dir, paste0(b, ".png"))))
    add_check(paste0("preview_width_", b), info$width >= 1400, info$width, ">=1400 px")
    add_check(paste0("preview_height_", b), info$height >= 1000, info$height, ">=1000 px")
  }
}

checks <- dplyr::bind_rows(checks)
readr::write_csv(checks, file.path(ann_dir, "annotation_figure_validation_checks.csv"))
if (any(!checks$pass)) {
  print(checks[!checks$pass, ])
  stop("Annotation/figure validation failed: ", sum(!checks$pass), " checks")
}
log_step("Annotation/figure validation passed: ", nrow(checks), " checks")
