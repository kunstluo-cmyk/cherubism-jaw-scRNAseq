root <- Sys.getenv("CB_REVISION_ROOT")
if (!nzchar(root)) root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
config <- yaml::read_yaml(file.path(root, "config", "analysis_config.yml"))
seed <- as.integer(config$seed)

suppressPackageStartupMessages({
  library(monocle)
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(ggrepel)
  library(ggrastr)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

out_dir <- file.path(root, "figures", "publication_subfigures", "monocle2_joint_stromal")
src_dir <- file.path(root, "results", "trajectory", "monocle2_joint_stromal")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

state_display <- c(
  "FibroStromal", "LesionStromal", "ECM-remodeling stromal state",
  "ChondroFibroStromal", "LEPR+ MSC-like", "Mature OB", "CNCC-OCPs"
)
state_colors <- c(
  "FibroStromal" = "#8DD3C7", "LesionStromal" = "#BEBADA",
  "ECM-remodeling stromal state" = "#BC80BD", "ChondroFibroStromal" = "#FCCDE5",
  "LEPR+ MSC-like" = "#FB8072", "Mature OB" = "#FDB462", "CNCC-OCPs" = "#FFED6F"
)
condition_colors <- c("CTRL" = "#4DBBD5", "CB" = "#E64B35")
program_colors <- c(
  "MSC niche" = "#FB8072", "Osteogenic program" = "#E6AB02",
  "Mineralization" = "#FDB462", "Fibro-chondrogenic" = "#8DD3C7",
  "Lesion / ECM remodeling" = "#BC80BD", "OC support" = "#6A3D9A"
)
programs <- list(
  "MSC niche" = c("LEPR", "CXCL12", "CFD"),
  "Osteogenic program" = c("DLX5", "SP7", "ALPL"),
  "Mineralization" = c("BGLAP", "DMP1", "IFITM5", "IBSP"),
  "Fibro-chondrogenic" = c("ASPN", "COL11A1", "WNT5A", "ACAN", "COL10A1", "COMP"),
  "Lesion / ECM remodeling" = c("FAP", "POSTN", "CTHRC1", "PDPN", "MMP2"),
  "OC support" = c("TNFSF11", "CSF1", "IL34", "CCL2")
)
program_order <- names(programs)
gene_program_map <- setNames(rep(program_order, lengths(programs)), unlist(programs, use.names = FALSE))

wrap_state <- function(x) {
  x <- as.character(x)
  x[x == "ECM-remodeling stromal state"] <- "ECM-remodeling\nstromal state"
  x[x == "ChondroFibroStromal"] <- "ChondroFibro\nStromal"
  x[x == "LEPR+ MSC-like"] <- "LEPR+\nMSC-like"
  x[x == "CNCC-OCPs"] <- "CNCC-\nOCPs"
  x
}

theme_pub <- theme_classic(base_size = 12, base_family = "Arial") +
  theme(
    text = element_text(family = "Arial", face = "bold", colour = "black"),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10.5, face = "bold", colour = "black"),
    axis.line = element_line(linewidth = 0.55, colour = "black"),
    axis.ticks = element_line(linewidth = 0.5, colour = "black"),
    legend.title = element_text(size = 10.5, face = "bold"),
    legend.text = element_text(size = 9.5, face = "bold"),
    plot.title = element_text(size = 13, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 9.5, face = "bold", colour = "#333333"),
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

save_heatmap <- function(ht, filename, width_mm, height_mm, dpi = 600,
                         bottom_pad_mm = 5, left_pad_mm = 8) {
  path <- file.path(out_dir, filename)
  wi <- width_mm / 25.4; hi <- height_mm / 25.4
  draw_ht <- function() ComplexHeatmap::draw(
    ht, heatmap_legend_side = "right", annotation_legend_side = "right",
    merge_legends = FALSE, padding = unit(c(5, 8, bottom_pad_mm, left_pad_mm), "mm")
  )
  svglite::svglite(paste0(path, ".svg"), width = wi, height = hi); draw_ht(); dev.off()
  cairo_pdf(paste0(path, ".pdf"), width = wi, height = hi, family = "Arial"); draw_ht(); dev.off()
  ragg::agg_tiff(paste0(path, ".tiff"), width = wi, height = hi, units = "in",
                 res = dpi, compression = "lzw"); draw_ht(); dev.off()
  ragg::agg_png(paste0(path, ".png"), width = wi, height = hi, units = "in", res = 300);
  draw_ht(); dev.off()
}

joint_umap <- read.csv(file.path(src_dir, "source_data_A_joint_RPCA_UMAP.csv"), check.names = FALSE)
joint_umap$cell_state <- factor(joint_umap$cell_state, levels = state_display)
label_umap <- joint_umap |>
  group_by(cell_state) |>
  summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2), .groups = "drop") |>
  mutate(label = wrap_state(cell_state))

p_umap <- ggplot(joint_umap, aes(UMAP_1, UMAP_2, colour = cell_state)) +
  ggrastr::geom_point_rast(size = 0.32, alpha = 0.86, raster.dpi = 600) +
  scale_colour_manual(values = state_colors, drop = FALSE) +
  ggrepel::geom_text_repel(
    data = label_umap, aes(UMAP_1, UMAP_2, label = label), inherit.aes = FALSE,
    family = "Arial", fontface = "bold", size = 3.35, colour = "black",
    box.padding = 0.26, point.padding = 0.07, seed = seed, max.overlaps = Inf,
    min.segment.length = 0, segment.colour = "#555555", segment.size = 0.35
  ) +
  coord_equal() + theme_pub +
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
        axis.line = element_blank(), legend.position = "none") +
  annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf,
           arrow = arrow(length = unit(2.4, "mm")), linewidth = 0.55) +
  annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf,
           arrow = arrow(length = unit(2.4, "mm")), linewidth = 0.55) +
  annotate("text", x = Inf, y = -Inf, label = "UMAP1", hjust = 1.1, vjust = 1.7,
           family = "Arial", fontface = "bold", size = 3.8) +
  annotate("text", x = -Inf, y = Inf, label = "UMAP2", angle = 90, hjust = 1.1,
           vjust = -0.5, family = "Arial", fontface = "bold", size = 3.8) +
  labs(title = "Joint stromal/osteogenic landscape")
save_plot(p_umap, "A_joint_stromal_RPCA_UMAP_cell_state_no_trajectory", 165, 108)

cds <- readRDS(file.path(root, "objects", "09_joint_stromal_monocle2_exploratory.rds"))
coord <- read.csv(file.path(src_dir, "source_data_B_Monocle2_DDRTree_pseudotime.csv"), check.names = FALSE)
label_traj <- coord |>
  group_by(cell_state) |>
  summarise(Component_1 = median(Component_1), Component_2 = median(Component_2), .groups = "drop") |>
  mutate(label = wrap_state(cell_state))
root_state <- read.csv(file.path(src_dir, "analysis_manifest.csv"), check.names = FALSE) |>
  filter(item == "root_state") |>
  pull(value)

p_traj <- monocle::plot_cell_trajectory(
  cds, color_by = "Pseudotime", show_branch_points = FALSE,
  show_state_number = FALSE, cell_size = 0.68
) +
  scale_colour_gradientn(
    colours = c("#3C8DBC", "#8DD3C7", "#FFF7BC", "#F46D43", "#B30000"),
    limits = range(pData(cds)$Pseudotime, na.rm = TRUE), name = "Pseudotime"
  ) +
  ggrepel::geom_text_repel(
    data = label_traj, aes(Component_1, Component_2, label = label), inherit.aes = FALSE,
    family = "Arial", fontface = "bold", size = 3.15, colour = "black",
    box.padding = 0.28, point.padding = 0.06, seed = seed, max.overlaps = Inf,
    min.segment.length = 0, segment.colour = "#555555", segment.size = 0.35
  ) +
  theme_pub +
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
        axis.line = element_blank(), legend.position = "right") +
  labs(title = "Exploratory Monocle2 ordering",
       subtitle = paste0("Rooted in DDRTree state ", root_state, " (LEPR+ MSC-like enriched)"))
save_plot(p_traj, "B_Monocle2_DDRTree_pseudotime_with_cell_state_labels", 175, 112)

count_df <- read.csv(file.path(src_dir, "source_data_C_stromal_cell_counts.csv"), check.names = FALSE)
count_df$condition <- factor(count_df$condition, levels = c("CTRL", "CB"))
count_df$cell_state <- factor(count_df$cell_state, levels = state_display)
p_count <- ggplot(count_df, aes(n_cells, factor(cell_state, levels = rev(state_display)), fill = condition)) +
  geom_col(position = position_dodge(width = 0.78), width = 0.64) +
  geom_text(aes(label = scales::comma(n_cells)), position = position_dodge(width = 0.78),
            hjust = -0.16, family = "Arial", fontface = "bold", size = 3.4) +
  scale_fill_manual(values = condition_colors, breaks = c("CTRL", "CB"), drop = FALSE) +
  scale_x_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.18))) +
  theme_pub + theme(axis.title.y = element_blank(), legend.position = "top") +
  labs(x = "Number of cells", fill = "Condition", title = "Stromal/osteogenic cell-state counts")
save_plot(p_count, "C_stromal_cell_state_counts_CB_CTRL", 165, 108)

stage_df <- read.csv(file.path(src_dir, "source_data_C2_pseudotime_stage_composition.csv"), check.names = FALSE)
stage_df$condition <- factor(stage_df$condition, levels = c("CTRL", "CB"))
stage_df$pseudotime_stage <- factor(stage_df$pseudotime_stage, levels = c("Early", "Middle", "Late"))
stage_cols <- c("Early" = "#3C8DBC", "Middle" = "#FFF7BC", "Late" = "#B30000")
p_stage <- ggplot(stage_df, aes(condition, percent, fill = pseudotime_stage)) +
  geom_col(width = 0.66) +
  geom_text(aes(label = ifelse(percent >= 2, sprintf("%.1f%%", percent), "")),
            position = position_stack(vjust = 0.5), family = "Arial", fontface = "bold", size = 3.6) +
  scale_fill_manual(values = stage_cols, breaks = c("Early", "Middle", "Late"), drop = FALSE) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  theme_pub + theme(axis.title.x = element_blank(), legend.position = "top") +
  labs(y = "Percentage of stromal/osteogenic cells", fill = "Exploratory stage",
       title = "Monocle2 pseudotime-stage composition")
save_plot(p_stage, "C2_Monocle2_pseudotime_stage_composition", 122, 105)

heat_col <- circlize::colorRamp2(c(-2.5, 0, 2.5), c("#3C8DBC", "#F7F7F7", "#B30000"))
legend_gp <- list(
  title_gp = gpar(fontfamily = "Arial", fontface = "bold", fontsize = 10),
  labels_gp = gpar(fontfamily = "Arial", fontface = "bold", fontsize = 8.5)
)

state_z_df <- read.csv(file.path(src_dir, "source_data_D_state_gene_zscores.csv"), check.names = FALSE)
state_z <- as.matrix(state_z_df[, setdiff(colnames(state_z_df), "gene"), drop = FALSE])
rownames(state_z) <- state_z_df$gene
gene_program <- factor(unname(gene_program_map[rownames(state_z)]), levels = names(program_colors))
state_order_pt <- colnames(state_z)

top_state <- HeatmapAnnotation(
  `Cell state` = factor(state_order_pt, levels = names(state_colors)),
  col = list(`Cell state` = state_colors), show_annotation_name = FALSE,
  annotation_legend_param = list(`Cell state` = legend_gp)
)
left_program <- rowAnnotation(
  Program = gene_program, col = list(Program = program_colors), show_annotation_name = FALSE,
  annotation_legend_param = list(Program = legend_gp), width = unit(4, "mm")
)
ht_state <- Heatmap(
  state_z, name = "Gene-wise\nZ-score", col = heat_col,
  cluster_rows = FALSE, cluster_columns = FALSE, row_split = gene_program,
  row_title = NULL, row_gap = unit(2.1, "mm"), top_annotation = top_state,
  left_annotation = left_program, column_names_rot = 35,
  column_labels = wrap_state(state_order_pt),
  column_names_gp = gpar(fontfamily = "Arial", fontface = "bold", fontsize = 9),
  row_names_side = "right",
  row_names_gp = gpar(fontfamily = "Arial", fontface = "bold.italic", fontsize = 9),
  column_title = "Functional programs by cell state",
  column_title_gp = gpar(fontfamily = "Arial", fontface = "bold", fontsize = 12),
  rect_gp = gpar(col = "white", lwd = 0.7),
  heatmap_legend_param = legend_gp
)
save_heatmap(
  ht_state, "D_cell_state_by_functional_program_gene_heatmap", 215, 154,
  bottom_pad_mm = 26, left_pad_mm = 14
)

pt_z_df <- read.csv(file.path(src_dir, "source_data_E_pseudotime_smoothed_gene_zscores.csv"), check.names = FALSE)
pt_z <- as.matrix(pt_z_df[, setdiff(colnames(pt_z_df), "gene"), drop = FALSE])
rownames(pt_z) <- pt_z_df$gene
bin_info <- read.csv(file.path(src_dir, "source_data_E_pseudotime_bin_annotations.csv"), check.names = FALSE)
gene_program_pt <- factor(unname(gene_program_map[rownames(pt_z)]), levels = names(program_colors))

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
  annotation_legend_param = list(Condition = legend_gp, `Cell state` = legend_gp, Pseudotime = legend_gp)
)
left_program_pt <- rowAnnotation(
  Program = gene_program_pt, col = list(Program = program_colors), show_annotation_name = FALSE,
  annotation_legend_param = list(Program = legend_gp), width = unit(4.5, "mm")
)
ht_pt <- Heatmap(
  pt_z, name = "Smoothed\nZ-score", col = heat_col,
  cluster_rows = FALSE, cluster_columns = FALSE, row_split = gene_program_pt,
  row_title = NULL, row_gap = unit(2.2, "mm"), top_annotation = pt_annotation,
  left_annotation = left_program_pt, show_column_names = FALSE,
  row_names_side = "right",
  row_names_gp = gpar(fontfamily = "Arial", fontface = "bold.italic", fontsize = 9.5),
  column_title = "Functional programs over exploratory pseudotime",
  column_title_gp = gpar(fontfamily = "Arial", fontface = "bold", fontsize = 11),
  rect_gp = gpar(col = NA), heatmap_legend_param = legend_gp
)
save_heatmap(
  ht_pt, "E_Monocle2_pseudotime_functional_program_gradient_heatmap", 215, 145,
  bottom_pad_mm = 5, left_pad_mm = 14
)

message("Render-only refresh complete: ", out_dir)
