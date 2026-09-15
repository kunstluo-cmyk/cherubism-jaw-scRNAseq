root <- Sys.getenv("CB_REVISION_ROOT")
if (!nzchar(root)) root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

suppressPackageStartupMessages({
  library(magick)
  library(dplyr)
  library(tibble)
})

fig_dir <- file.path(root, "figures", "publication_subfigures", "cellchat_figure6_rebuild")
res_dir <- file.path(root, "results", "cellchat_figure6_rebuild")

stems <- c(
  A = "A_CTRL_common_lineage_CellChat_interaction_strength",
  B = "B_CB_common_lineage_CellChat_interaction_strength",
  C = "C_common_source_target_CB_minus_CTRL_heatmap",
  D = "D_CB_stromal_to_TREM2_OC_pathway_strength_heatmap",
  E = "E_CB_key_ligand_receptor_bubbleplot",
  F = "F_CB_sender_ligand_receiver_receptor_expression",
  G = "G_NicheNet_v2_candidate_ligand_activity_ranking",
  H = "H_NicheNet_v2_ligand_target_heatmap"
)
expected <- tribble(
  ~panel, ~width_px, ~height_px,
  "A", 3300L, 2850L,
  "B", 3300L, 2850L,
  "C", 3900L, 3150L,
  "D", 4500L, 3300L,
  "E", 5400L, 3900L,
  "F", 6250L, 1700L,
  "G", 3900L, 3300L,
  "H", 5400L, 3600L
)

formats <- c("svg", "pdf", "tiff", "png")
file_audit <- bind_rows(lapply(names(stems), function(panel_i) {
  bind_rows(lapply(formats, function(ext_i) {
    f <- file.path(fig_dir, paste0(stems[[panel_i]], ".", ext_i))
    tibble(
      panel = panel_i,
      format = ext_i,
      file = normalizePath(f, winslash = "/", mustWork = FALSE),
      exists = file.exists(f),
      bytes = if (file.exists(f)) file.info(f)$size else NA_real_
    )
  }))
}))

raster_audit <- bind_rows(lapply(names(stems), function(panel_i) {
  bind_rows(lapply(c("png", "tiff"), function(ext_i) {
    f <- file.path(fig_dir, paste0(stems[[panel_i]], ".", ext_i))
    info <- image_info(image_read(f))[1, ]
    tibble(
      panel = panel_i,
      format = ext_i,
      width_px = info$width,
      height_px = info$height,
      colorspace = info$colorspace
    )
  }))
})) |>
  left_join(expected, by = "panel", suffix = c("_observed", "_expected")) |>
  mutate(
    dimension_pass = width_px_observed == width_px_expected &
      height_px_observed == height_px_expected
  )

write.csv(file_audit, file.path(res_dir, "Figure6_file_format_audit.csv"), row.names = FALSE)
write.csv(raster_audit, file.path(res_dir, "Figure6_raster_dimension_audit.csv"), row.names = FALSE)

if (!all(file_audit$exists) || any(file_audit$bytes <= 0, na.rm = TRUE)) {
  stop("Figure 6 file-format audit failed.")
}
if (!all(raster_audit$dimension_pass)) stop("Figure 6 raster-dimension audit failed.")

prepare_panel <- function(panel_i, width = 2400L) {
  f <- file.path(fig_dir, paste0(stems[[panel_i]], ".png"))
  img <- image_read(f) |>
    image_scale(paste0(width, "x")) |>
    image_border("white", "55x55")
  image_annotate(
    img, panel_i, gravity = "northwest", location = "+18+12",
    size = 68, font = "Arial", weight = 700, color = "black"
  )
}

pad_pair <- function(img1, img2) {
  i1 <- image_info(img1)[1, ]
  i2 <- image_info(img2)[1, ]
  h <- max(i1$height, i2$height)
  img1 <- image_extent(img1, geometry = paste0(i1$width, "x", h), gravity = "north", color = "white")
  img2 <- image_extent(img2, geometry = paste0(i2$width, "x", h), gravity = "north", color = "white")
  image_append(c(img1, img2), stack = FALSE)
}

rows <- list(
  pad_pair(prepare_panel("A"), prepare_panel("B")),
  pad_pair(prepare_panel("C"), prepare_panel("D")),
  pad_pair(prepare_panel("E"), prepare_panel("F")),
  pad_pair(prepare_panel("G"), prepare_panel("H"))
)
row_width <- max(vapply(rows, function(x) image_info(x)$width[[1]], numeric(1)))
rows <- lapply(rows, function(x) image_extent(x, paste0(row_width, "x", image_info(x)$height[[1]]),
                                               gravity = "center", color = "white"))
preview <- image_append(do.call(c, rows), stack = TRUE)
image_write(preview, file.path(fig_dir, "Figure6_A-H_composite_preview.png"), format = "png", density = "300x300")

source_tables <- c(
  "A_broad_strength_matrix_CTRL.csv",
  "B_broad_strength_matrix_CB.csv",
  "C_common_lineage_descriptive_log2_CB_CTRL.csv",
  "D_CB_pathway_strength_source_data.csv",
  "E_CB_key_ligand_receptor_source_data.csv",
  "F_ligand_receptor_expression_source_data.csv",
  "G_NicheNet_v2_candidate_ligand_activities.csv",
  "H_NicheNet_v2_top_ligand_target_links.csv"
)
table_audit <- tibble(
  file = source_tables,
  exists = file.exists(file.path(res_dir, source_tables)),
  rows = vapply(source_tables, function(x) {
    f <- file.path(res_dir, x)
    if (!file.exists(f)) return(NA_integer_)
    nrow(read.csv(f, check.names = FALSE))
  }, integer(1))
)
write.csv(table_audit, file.path(res_dir, "Figure6_source_table_audit.csv"), row.names = FALSE)
if (!all(table_audit$exists) || any(table_audit$rows <= 0, na.rm = TRUE)) {
  stop("Figure 6 source-data audit failed.")
}

writeLines(
  c(
    "Recommended layout: A-B / C-D / E-F / G-H.",
    "A and B use the same broad lineages, node order, and global edge-width maximum.",
    "C is a descriptive log2 strength ratio and is not a donor-level statistical comparison.",
    "D-F are CB-focused, DecontX-corrected sender-receiver evidence.",
    "G-H use the official human NicheNet v2 ligand-target prior and pooled descriptive receiver-state contrasts.",
    "All A-H individual panels passed SVG/PDF/TIFF/PNG presence and raster-dimension checks."
  ),
  file.path(res_dir, "Figure6_layout_and_QA_notes.txt")
)
