root <- Sys.getenv("CB_REVISION_ROOT")
if (!nzchar(root)) root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
config <- yaml::read_yaml(file.path(root, "config", "analysis_config.yml"))
seed <- as.integer(config$seed)
set.seed(seed)

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggrastr)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

out_dir <- file.path(root, "figures", "publication_subfigures", "monocle2_joint_stromal")
src_dir <- file.path(root, "results", "trajectory", "monocle2_joint_stromal")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(src_dir, recursive = TRUE, showWarnings = FALSE)

state_order <- c(
  "LEPR+ MSC-like", "Mature OB", "CNCC-OCPs", "FibroStromal",
  "ECM-remodeling stromal state", "ChondroFibroStromal", "LesionStromal"
)
state_colors <- c(
  "FibroStromal" = "#8DD3C7", "LesionStromal" = "#BEBADA",
  "ECM-remodeling stromal state" = "#BC80BD", "ChondroFibroStromal" = "#FCCDE5",
  "LEPR+ MSC-like" = "#FB8072", "Mature OB" = "#FDB462", "CNCC-OCPs" = "#FFED6F"
)
condition_colors <- c("CTRL" = "#4DBBD5", "CB" = "#E64B35")

# OC support is intentionally restricted to positive stromal support factors:
# TNFSF11/RANKL drives osteoclast differentiation; CSF1 and IL34 provide CSF1R
# survival/differentiation signals; CCL2 supports monocyte/osteoclast recruitment.
# TNFRSF11B/OPG is excluded because it antagonizes RANKL.
programs <- list(
  "MSC niche" = c("LEPR", "CXCL12", "CFD"),
  "Osteogenic program" = c("DLX5", "SP7", "ALPL"),
  "Mineralization" = c("BGLAP", "DMP1", "IFITM5", "IBSP"),
  "Fibro-chondrogenic" = c("ASPN", "COL11A1", "WNT5A", "ACAN", "COL10A1", "COMP"),
  "Lesion / ECM remodeling" = c("FAP", "POSTN", "CTHRC1", "PDPN", "MMP2"),
  "OC support" = c("TNFSF11", "CSF1", "IL34", "CCL2")
)
program_order <- names(programs)
program_colors <- c(
  "MSC niche" = "#FB8072", "Osteogenic program" = "#E6AB02",
  "Mineralization" = "#FDB462", "Fibro-chondrogenic" = "#8DD3C7",
  "Lesion / ECM remodeling" = "#BC80BD", "OC support" = "#6A3D9A"
)
gene_program <- tibble(
  program = rep(program_order, lengths(programs)),
  gene = unlist(programs, use.names = FALSE)
)

theme_pub <- theme_classic(base_size = 11.5, base_family = "Arial") +
  theme(
    text = element_text(family = "Arial", face = "bold", colour = "black"),
    axis.title = element_text(size = 11.5, face = "bold"),
    axis.text = element_text(size = 9.5, face = "bold", colour = "black"),
    axis.line = element_line(linewidth = 0.5, colour = "black"),
    axis.ticks = element_line(linewidth = 0.45, colour = "black"),
    strip.background = element_blank(),
    strip.text = element_text(size = 10.5, face = "bold.italic", colour = "black"),
    legend.title = element_text(size = 10.5, face = "bold"),
    legend.text = element_text(size = 9.2, face = "bold"),
    plot.title = element_text(size = 13, face = "bold"),
    panel.grid = element_blank(),
    plot.margin = margin(7, 9, 7, 7)
  )

save_plot <- function(plot, filename, width_mm, height_mm, dpi = 600) {
  path <- file.path(out_dir, filename)
  wi <- width_mm / 25.4; hi <- height_mm / 25.4
  svglite::svglite(paste0(path, ".svg"), width = wi, height = hi); print(plot); dev.off()
  cairo_pdf(paste0(path, ".pdf"), width = wi, height = hi, family = "Arial"); print(plot); dev.off()
  ragg::agg_tiff(paste0(path, ".tiff"), width = wi, height = hi, units = "in",
                 res = dpi, compression = "lzw"); print(plot); dev.off()
  ragg::agg_png(paste0(path, ".png"), width = wi, height = hi, units = "in", res = 300);
  print(plot); dev.off()
}

save_plot_exact_pixels <- function(
  plot, filename, width_px, height_px, dpi = 600
) {
  path <- file.path(out_dir, filename)
  wi <- width_px / dpi
  hi <- height_px / dpi

  svglite::svglite(paste0(path, ".svg"), width = wi, height = hi)
  print(plot)
  dev.off()

  cairo_pdf(
    paste0(path, ".pdf"), width = wi, height = hi, family = "Arial"
  )
  print(plot)
  dev.off()

  ragg::agg_tiff(
    paste0(path, ".tiff"), width = width_px, height = height_px,
    units = "px", res = dpi, compression = "lzw"
  )
  print(plot)
  dev.off()

  ragg::agg_png(
    paste0(path, ".png"), width = width_px, height = height_px,
    units = "px", res = dpi
  )
  print(plot)
  dev.off()
}

save_heatmap <- function(ht, filename, width_mm, height_mm, dpi = 600,
                         bottom_pad_mm = 5, left_pad_mm = 14) {
  path <- file.path(out_dir, filename)
  wi <- width_mm / 25.4; hi <- height_mm / 25.4
  draw_ht <- function() ComplexHeatmap::draw(
    ht, heatmap_legend_side = "right", annotation_legend_side = "right",
    merge_legends = FALSE,
    padding = unit(c(5, 8, bottom_pad_mm, left_pad_mm), "mm")
  )
  svglite::svglite(paste0(path, ".svg"), width = wi, height = hi); draw_ht(); dev.off()
  cairo_pdf(paste0(path, ".pdf"), width = wi, height = hi, family = "Arial"); draw_ht(); dev.off()
  ragg::agg_tiff(paste0(path, ".tiff"), width = wi, height = hi, units = "in",
                 res = dpi, compression = "lzw"); draw_ht(); dev.off()
  ragg::agg_png(paste0(path, ".png"), width = wi, height = hi, units = "in", res = 300);
  draw_ht(); dev.off()
}

mode_value <- function(x) {
  x <- x[!is.na(x)]
  if (!length(x)) return(NA_character_)
  names(which.max(table(x)))[1]
}

zscore_rows <- function(m, limit = 2.5) {
  z <- t(scale(t(m)))
  z[!is.finite(z)] <- 0
  pmax(pmin(z, limit), -limit)
}

wrap_state <- function(x) {
  x <- as.character(x)
  x[x == "ECM-remodeling stromal state"] <- "ECM-remodeling\nstromal state"
  x[x == "ChondroFibroStromal"] <- "ChondroFibro\nStromal"
  x[x == "LEPR+ MSC-like"] <- "LEPR+\nMSC-like"
  x[x == "CNCC-OCPs"] <- "CNCC-\nOCPs"
  x
}

message("Loading corrected expression and existing Monocle2 coordinates...")
coord <- read.csv(file.path(src_dir, "source_data_B_Monocle2_DDRTree_pseudotime.csv"),
                  check.names = FALSE, stringsAsFactors = FALSE)
obj <- readRDS(file.path(root, "objects", "06_base_broad_lineages_reviewed.rds"))
cells <- coord$cell_id[coord$cell_id %in% colnames(obj)]
obj <- subset(obj, cells = cells)
obj$cell_state_display <- as.character(obj$celltype_reviewed_provisional)
obj$cell_state_display[obj$cell_state_display == "Act OB-Stromal"] <-
  "ECM-remodeling stromal state"
obj$cell_state_display <- factor(obj$cell_state_display, levels = state_order)
obj$condition_display <- factor(as.character(obj$group), levels = c("CTRL", "CB"))

coord <- coord[match(colnames(obj), coord$cell_id), , drop = FALSE]
stopifnot(identical(coord$cell_id, colnames(obj)))
log_data <- SeuratObject::LayerData(obj, assay = "RNA", layer = "data")

genes_present <- intersect(gene_program$gene, rownames(log_data))
genes_missing <- setdiff(gene_program$gene, genes_present)
gene_program_present <- gene_program |> filter(gene %in% genes_present)
write.csv(
  tibble(
    program = program_order,
    requested_genes = vapply(programs, paste, collapse = ";", character(1)),
    present_genes = vapply(programs, function(g) paste(intersect(g, genes_present), collapse = ";"), character(1)),
    n_present = vapply(programs, function(g) length(intersect(g, genes_present)), integer(1))
  ),
  file.path(src_dir, "functional_program_gene_definitions_with_OC_support.csv"), row.names = FALSE
)
write.csv(data.frame(gene = genes_missing),
          file.path(src_dir, "functional_program_genes_missing_after_OC_update.csv"), row.names = FALSE)

message("Recomputing the six-program state heatmap...")
state_cells <- split(colnames(obj), as.character(obj$cell_state_display))[state_order]
state_avg <- vapply(state_cells, function(cc) {
  Matrix::rowMeans(log_data[genes_present, cc, drop = FALSE])
}, numeric(length(genes_present)))
rownames(state_avg) <- genes_present
colnames(state_avg) <- state_order
state_z_all <- zscore_rows(state_avg, 2.5)
state_z <- state_z_all[gene_program_present$gene, state_order, drop = FALSE]
row_program <- factor(gene_program_present$program, levels = program_order)

write.csv(as.data.frame(state_avg) |> rownames_to_column("gene"),
          file.path(src_dir, "source_data_D_state_average_log_expression.csv"), row.names = FALSE)
write.csv(as.data.frame(state_z) |> rownames_to_column("gene"),
          file.path(src_dir, "source_data_D_state_gene_zscores.csv"), row.names = FALSE)

heat_col <- circlize::colorRamp2(c(-2.5, 0, 2.5), c("#3C8DBC", "#F7F7F7", "#B30000"))
legend_gp <- list(
  title_gp = gpar(fontfamily = "Arial", fontface = "bold", fontsize = 10),
  labels_gp = gpar(fontfamily = "Arial", fontface = "bold", fontsize = 8.5)
)
top_state <- HeatmapAnnotation(
  `Cell state` = factor(state_order, levels = names(state_colors)),
  col = list(`Cell state` = state_colors), show_annotation_name = FALSE,
  annotation_legend_param = list(`Cell state` = legend_gp)
)
left_program <- rowAnnotation(
  Program = row_program, col = list(Program = program_colors), show_annotation_name = FALSE,
  annotation_legend_param = list(Program = legend_gp), width = unit(4, "mm")
)
ht_state <- Heatmap(
  state_z, name = "Gene-wise\nZ-score", col = heat_col,
  cluster_rows = FALSE, cluster_columns = FALSE, row_split = row_program,
  row_title = NULL, row_gap = unit(1.8, "mm"), top_annotation = top_state,
  left_annotation = left_program, column_names_rot = 35,
  column_labels = wrap_state(state_order),
  column_names_gp = gpar(fontfamily = "Arial", fontface = "bold", fontsize = 9),
  row_names_side = "right",
  row_names_gp = gpar(fontfamily = "Arial", fontface = "bold.italic", fontsize = 9),
  column_title = "Functional programs by cell state",
  column_title_gp = gpar(fontfamily = "Arial", fontface = "bold", fontsize = 12),
  rect_gp = gpar(col = "white", lwd = 0.7), heatmap_legend_param = legend_gp
)
save_heatmap(
  ht_state, "D_cell_state_by_functional_program_gene_heatmap", 215, 174,
  bottom_pad_mm = 26
)

message("Recomputing pseudotime heatmap with OC-support genes...")
pt_df <- data.frame(
  cell_id = coord$cell_id,
  pseudotime = coord$Pseudotime,
  pseudotime_scaled = coord$pseudotime_scaled,
  cell_state = as.character(coord$cell_state),
  condition = as.character(coord$condition), stringsAsFactors = FALSE
) |>
  arrange(pseudotime)
pt_df$bin <- pmin(90L, ceiling(seq_len(nrow(pt_df)) / nrow(pt_df) * 90L))
bin_info <- pt_df |>
  group_by(bin) |>
  summarise(
    pseudotime = mean(pseudotime), pseudotime_scaled = mean(pseudotime_scaled),
    cell_state = mode_value(cell_state), condition = mode_value(condition),
    n_cells = n(), .groups = "drop"
  )
bin_levels <- sort(unique(pt_df$bin))
expr_pt <- log_data[genes_present, pt_df$cell_id, drop = FALSE]
bin_avg <- vapply(bin_levels, function(b) {
  cc <- pt_df$cell_id[pt_df$bin == b]
  Matrix::rowMeans(expr_pt[, cc, drop = FALSE])
}, numeric(length(genes_present)))
rownames(bin_avg) <- genes_present
colnames(bin_avg) <- paste0("B", sprintf("%02d", bin_levels))
smoothed <- t(apply(bin_avg, 1, function(y) {
  x <- bin_info$pseudotime_scaled
  if (sd(y, na.rm = TRUE) == 0) return(y)
  fit <- smooth.spline(x = x, y = y, spar = 0.48)
  predict(fit, x = x)$y
}))
rownames(smoothed) <- rownames(bin_avg)
colnames(smoothed) <- colnames(bin_avg)
pt_z_all <- zscore_rows(smoothed, 2.5)
pt_z <- pt_z_all[gene_program_present$gene, , drop = FALSE]

write.csv(bin_info, file.path(src_dir, "source_data_E_pseudotime_bin_annotations.csv"), row.names = FALSE)
write.csv(as.data.frame(bin_avg) |> rownames_to_column("gene"),
          file.path(src_dir, "source_data_E_pseudotime_binned_log_expression.csv"), row.names = FALSE)
write.csv(as.data.frame(pt_z) |> rownames_to_column("gene"),
          file.path(src_dir, "source_data_E_pseudotime_smoothed_gene_zscores.csv"), row.names = FALSE)

pt_annotation <- HeatmapAnnotation(
  Condition = factor(bin_info$condition, levels = c("CTRL", "CB")),
  `Cell state` = factor(bin_info$cell_state, levels = names(state_colors)),
  Pseudotime = bin_info$pseudotime_scaled,
  col = list(
    Condition = condition_colors, `Cell state` = state_colors,
    Pseudotime = circlize::colorRamp2(c(0, 0.5, 1), c("#3C8DBC", "#FFF7BC", "#B30000"))
  ),
  show_annotation_name = TRUE,
  annotation_name_gp = gpar(fontfamily = "Arial", fontface = "bold", fontsize = 9),
  annotation_legend_param = list(Condition = legend_gp, `Cell state` = legend_gp,
                                 Pseudotime = legend_gp)
)
left_program_pt <- rowAnnotation(
  Program = row_program, col = list(Program = program_colors), show_annotation_name = FALSE,
  annotation_legend_param = list(Program = legend_gp), width = unit(4.5, "mm")
)
ht_pt <- Heatmap(
  pt_z, name = "Smoothed\nZ-score", col = heat_col,
  cluster_rows = FALSE, cluster_columns = FALSE, row_split = row_program,
  row_title = NULL, row_gap = unit(1.8, "mm"), top_annotation = pt_annotation,
  left_annotation = left_program_pt, show_column_names = FALSE,
  row_names_side = "right",
  row_names_gp = gpar(fontfamily = "Arial", fontface = "bold.italic", fontsize = 9.3),
  column_title = "Functional programs over exploratory pseudotime",
  column_title_gp = gpar(fontfamily = "Arial", fontface = "bold", fontsize = 11),
  rect_gp = gpar(col = NA), heatmap_legend_param = legend_gp
)
save_heatmap(
  ht_pt, "E_Monocle2_pseudotime_functional_program_gradient_heatmap", 215, 162,
  bottom_pad_mm = 5
)

message("Building descriptive program-score bar plot...")
# Program scores are the mean gene-wise standardized state-average expression.
# They are descriptive state scores, not independent replicate-level estimates.
score_df <- lapply(program_order, function(pr) {
  gg <- intersect(programs[[pr]], rownames(state_z_all))
  tibble(
    program = pr, cell_state = state_order,
    score = colMeans(state_z_all[gg, state_order, drop = FALSE]),
    n_genes = length(gg)
  )
}) |>
  bind_rows()
score_df$program <- factor(score_df$program, levels = program_order)
score_df$cell_state <- factor(score_df$cell_state, levels = state_order)
program_labeller <- ggplot2::as_labeller(c(
  "MSC niche" = "MSC niche",
  "Osteogenic program" = "Osteogenic program",
  "Mineralization" = "Mineralization",
  "Fibro-chondrogenic" = "Fibro-chondrogenic",
  "Lesion / ECM remodeling" = "Lesion / ECM\nremodeling",
  "OC support" = "OC support"
))
write.csv(score_df, file.path(src_dir, "source_data_F_functional_program_scores_by_cell_state.csv"),
          row.names = FALSE)

p_score <- ggplot(score_df, aes(cell_state, score, fill = cell_state)) +
  geom_hline(yintercept = 0, colour = "#555555", linewidth = 0.4) +
  geom_col(width = 0.72) +
  facet_wrap(
    ~program, ncol = 2, nrow = 3, scales = "free_y",
    labeller = program_labeller
  ) +
  scale_fill_manual(values = state_colors, labels = wrap_state, drop = FALSE) +
  scale_x_discrete(labels = wrap_state) +
  guides(
    fill = guide_legend(
      title = "Cell state", nrow = 3, byrow = TRUE,
      title.position = "top", title.hjust = 0.5,
      override.aes = list(colour = NA)
    )
  ) +
  theme_pub +
  theme(
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    axis.title.x = element_blank(), legend.position = "bottom",
    legend.direction = "horizontal", legend.box = "vertical",
    legend.key.width = unit(4.8, "mm"), legend.key.height = unit(4.8, "mm"),
    legend.spacing.x = unit(2.3, "mm"),
    panel.spacing.x = unit(5, "mm"), panel.spacing.y = unit(7, "mm")
  ) +
  labs(
    y = "Mean program score (gene-wise Z-score)", title = NULL
  )
save_plot_exact_pixels(
  p_score,
  "F_functional_program_score_barplots_by_cell_state",
  width_px = 2700,
  height_px = 2980,
  dpi = 600
)

message("Building key-gene pseudotime curves...")
curve_genes <- c(
  "LEPR", "CXCL12", "DLX5", "ALPL", "BGLAP", "DMP1",
  "ACAN", "COL10A1", "POSTN", "CTHRC1", "TNFSF11", "CSF1"
)
curve_genes <- intersect(curve_genes, rownames(log_data))
curve_levels <- curve_genes

# Plot a deterministic, state-stratified subset of individual cells for visual
# density and overlay a trend estimated from all 90 equal-cell pseudotime bins.
cell_plot_ids <- pt_df |>
  group_by(cell_state) |>
  group_modify(~slice_sample(.x, n = min(900L, nrow(.x)), replace = FALSE)) |>
  ungroup() |>
  pull(cell_id)
cell_plot_meta <- pt_df[match(cell_plot_ids, pt_df$cell_id), ]
point_mat <- as.matrix(log_data[curve_genes, cell_plot_ids, drop = FALSE])
point_df <- as.data.frame(t(point_mat)) |>
  rownames_to_column("cell_id") |>
  left_join(cell_plot_meta[, c("cell_id", "pseudotime_scaled", "cell_state")], by = "cell_id") |>
  pivot_longer(cols = all_of(curve_genes), names_to = "gene", values_to = "expression")

curve_bin_df <- as.data.frame(t(bin_avg[curve_genes, , drop = FALSE])) |>
  mutate(pseudotime_scaled = bin_info$pseudotime_scaled) |>
  pivot_longer(cols = all_of(curve_genes), names_to = "gene", values_to = "binned_expression")
curve_line_df <- curve_bin_df |>
  group_by(gene) |>
  group_modify(~{
    fit <- mgcv::gam(
      binned_expression ~ s(pseudotime_scaled, k = 6, bs = "cs"),
      data = .x, method = "REML"
    )
    xx <- seq(0, 1, length.out = 160)
    tibble(
      pseudotime_scaled = xx,
      expression = pmax(0, as.numeric(predict(fit, newdata = data.frame(pseudotime_scaled = xx))))
    )
  }) |>
  ungroup()

point_df$gene <- factor(point_df$gene, levels = curve_levels)
curve_line_df$gene <- factor(curve_line_df$gene, levels = curve_levels)
write.csv(curve_bin_df, file.path(src_dir, "source_data_G_key_gene_binned_expression_curves.csv"),
          row.names = FALSE)
write.csv(curve_line_df, file.path(src_dir, "source_data_G_key_gene_smoothed_expression_curves.csv"),
          row.names = FALSE)
write.csv(
  point_df |> group_by(gene, cell_state) |> summarise(n_plotted_cells = n(), .groups = "drop"),
  file.path(src_dir, "source_data_G_curve_point_sampling_summary.csv"), row.names = FALSE
)

p_curve <- ggplot(point_df, aes(pseudotime_scaled, expression, colour = cell_state)) +
  ggrastr::geom_point_rast(size = 0.25, alpha = 0.16, raster.dpi = 600) +
  geom_line(
    data = curve_line_df,
    aes(pseudotime_scaled, expression, group = gene),
    inherit.aes = FALSE, colour = "black", linewidth = 0.8, lineend = "round"
  ) +
  facet_wrap(~gene, nrow = 2, scales = "free_y") +
  scale_colour_manual(values = state_colors, drop = FALSE) +
  guides(
    colour = guide_legend(
      override.aes = list(size = 4.5, alpha = 1),
      keyheight = unit(5.5, "mm"),
      keywidth = unit(5.5, "mm")
    )
  ) +
  scale_x_continuous(breaks = c(0, 0.5, 1), expand = expansion(mult = c(0.02, 0.03))) +
  theme_pub +
  theme(
    legend.position = "right", axis.text = element_text(size = 8.3),
    strip.text = element_text(size = 10.5, face = "bold.italic"),
    panel.spacing.x = unit(4.3, "mm"), panel.spacing.y = unit(5.2, "mm")
  ) +
  labs(
    x = "Exploratory Monocle2 pseudotime (scaled)",
    y = "Log-normalized expression", colour = "Cell state",
    title = "Key-gene expression patterns over exploratory pseudotime"
  )
save_plot_exact_pixels(
  p_curve,
  "G_key_gene_expression_curves_over_Monocle2_pseudotime",
  width_px = 6300,
  height_px = 2750,
  dpi = 600
)

qa <- tibble(
  check = c(
    "n_cells", "n_programs", "n_requested_genes", "n_present_genes",
    "OC_support_definition", "score_interpretation", "curve_interpretation"
  ),
  value = c(
    ncol(obj), length(programs), nrow(gene_program), length(genes_present),
    "TNFSF11;CSF1;IL34;CCL2 (positive support factors; TNFRSF11B excluded)",
    "Descriptive state-average gene-wise Z-score; no donor-level inference",
    "Black lines are smooth-spline trends over 90 equal-cell bins; points are a deterministic state-stratified display subset"
  )
)
write.csv(qa, file.path(src_dir, "OC_support_program_and_curve_QA_manifest.csv"), row.names = FALSE)

message("OC-support update complete: ", out_dir)
