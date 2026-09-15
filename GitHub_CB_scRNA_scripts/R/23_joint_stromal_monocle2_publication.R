root <- Sys.getenv("CB_REVISION_ROOT")
if (!nzchar(root)) root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
Sys.setenv(CB_REVISION_ROOT = root)
config <- yaml::read_yaml(file.path(root, "config", "analysis_config.yml"))
seed <- as.integer(config$seed)
log_file <- file.path(root, "logs", "pipeline.log")
log_step <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

suppressPackageStartupMessages({
  library(Seurat)
  library(monocle)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(Matrix)
  library(ggrepel)
  library(ggrastr)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

set.seed(seed)

out_dir <- file.path(root, "figures", "publication_subfigures", "monocle2_joint_stromal")
src_dir <- file.path(root, "results", "trajectory", "monocle2_joint_stromal")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(src_dir, recursive = TRUE, showWarnings = FALSE)

state_raw <- c(
  "FibroStromal", "LesionStromal", "Act OB-Stromal", "ChondroFibroStromal",
  "LEPR+ MSC-like", "Mature OB", "CNCC-OCPs"
)
display_name <- c("Act OB-Stromal" = "ECM-remodeling stromal state")
rename_state <- function(x) {
  x <- as.character(x)
  x[x %in% names(display_name)] <- unname(display_name[x[x %in% names(display_name)]])
  x
}
state_display <- rename_state(state_raw)

state_colors <- c(
  "FibroStromal" = "#8DD3C7",
  "LesionStromal" = "#BEBADA",
  "ECM-remodeling stromal state" = "#BC80BD",
  "ChondroFibroStromal" = "#FCCDE5",
  "LEPR+ MSC-like" = "#FB8072",
  "Mature OB" = "#FDB462",
  "CNCC-OCPs" = "#FFED6F"
)
condition_colors <- c("CTRL" = "#4DBBD5", "CB" = "#E64B35")
program_colors <- c(
  "MSC niche" = "#FB8072",
  "Osteogenic program" = "#E6AB02",
  "Mineralization" = "#FDB462",
  "Fibro-chondrogenic" = "#8DD3C7",
  "Lesion / ECM remodeling" = "#BC80BD",
  "OC support" = "#6A3D9A"
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
program_gene_table <- tibble(
  program = rep(program_order, lengths(programs)),
  gene = unlist(programs, use.names = FALSE)
)

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
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4

  svglite::svglite(paste0(path, ".svg"), width = width_in, height = height_in)
  print(plot)
  grDevices::dev.off()

  grDevices::cairo_pdf(paste0(path, ".pdf"), width = width_in, height = height_in,
                       family = "Arial")
  print(plot)
  grDevices::dev.off()

  ragg::agg_tiff(paste0(path, ".tiff"), width = width_in, height = height_in,
                 units = "in", res = dpi, compression = "lzw")
  print(plot)
  grDevices::dev.off()

  ragg::agg_png(paste0(path, ".png"), width = width_in, height = height_in,
                units = "in", res = 300)
  print(plot)
  grDevices::dev.off()
}

save_heatmap <- function(ht, filename, width_mm, height_mm, dpi = 600) {
  path <- file.path(out_dir, filename)
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4
  draw_ht <- function() ComplexHeatmap::draw(
    ht, heatmap_legend_side = "right", annotation_legend_side = "right",
    merge_legends = FALSE, padding = unit(c(4, 8, 4, 4), "mm")
  )

  svglite::svglite(paste0(path, ".svg"), width = width_in, height = height_in)
  draw_ht()
  grDevices::dev.off()

  grDevices::cairo_pdf(paste0(path, ".pdf"), width = width_in, height = height_in,
                       family = "Arial")
  draw_ht()
  grDevices::dev.off()

  ragg::agg_tiff(paste0(path, ".tiff"), width = width_in, height = height_in,
                 units = "in", res = dpi, compression = "lzw")
  draw_ht()
  grDevices::dev.off()

  ragg::agg_png(paste0(path, ".png"), width = width_in, height = height_in,
                units = "in", res = 300)
  draw_ht()
  grDevices::dev.off()
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

log_step("Monocle2 publication rebuild: loading RPCA visualization object")
rpca <- readRDS(file.path(root, "objects", "07_rpca_broad_lineages_reviewed.rds"))
keep_rpca <- which(as.character(rpca$broad_lineage_reviewed) == "Stromal/osteogenic")
joint_umap <- as.data.frame(Seurat::Embeddings(rpca, reduction = "umap.rpca")[keep_rpca, , drop = FALSE]) |>
  rownames_to_column("cell_id")
colnames(joint_umap)[2:3] <- c("UMAP_1", "UMAP_2")
joint_umap$cell_state <- factor(
  rename_state(rpca$celltype_reviewed_provisional[keep_rpca]),
  levels = state_display
)
joint_umap$condition <- factor(as.character(rpca$group[keep_rpca]), levels = c("CTRL", "CB"))
joint_umap <- joint_umap |> filter(!is.na(cell_state))
write.csv(joint_umap, file.path(src_dir, "source_data_A_joint_RPCA_UMAP.csv"), row.names = FALSE)

label_umap <- joint_umap |>
  group_by(cell_state) |>
  summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2), .groups = "drop")

p_umap <- ggplot(joint_umap, aes(UMAP_1, UMAP_2, colour = cell_state)) +
  ggrastr::geom_point_rast(size = 0.36, alpha = 0.88, raster.dpi = 600) +
  scale_colour_manual(values = state_colors, drop = FALSE) +
  ggrepel::geom_label_repel(
    data = label_umap,
    aes(UMAP_1, UMAP_2, label = cell_state),
    inherit.aes = FALSE, family = "Arial", fontface = "bold", size = 3.15,
    colour = "black", fill = scales::alpha("white", 0.78), label.size = 0,
    box.padding = 0.28, point.padding = 0.08, seed = seed, max.overlaps = Inf,
    min.segment.length = 0
  ) +
  coord_equal() +
  theme_pub +
  theme(
    axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
    axis.line = element_blank(), legend.position = "none"
  ) +
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

rm(rpca)
invisible(gc())

log_step("Monocle2 publication rebuild: loading corrected unintegrated object")
obj <- readRDS(file.path(root, "objects", "06_base_broad_lineages_reviewed.rds"))
cells <- joint_umap$cell_id
cells <- cells[cells %in% colnames(obj)]
obj <- subset(obj, cells = cells)
obj$cell_state_display <- factor(rename_state(obj$celltype_reviewed_provisional), levels = state_display)
obj$condition_display <- factor(as.character(obj$group), levels = c("CTRL", "CB"))

counts <- SeuratObject::LayerData(obj, assay = "RNA", layer = "counts")
counts <- counts[, colnames(obj), drop = FALSE]
library_size <- Matrix::colSums(counts)
library_size_factor <- library_size / stats::median(library_size[library_size > 0])

# Select ordering genes within the stromal/osteogenic subset without using the
# condition label. All cells are retained. Library-size factors are calculated
# from the complete corrected count matrix before reducing the gene dimension.
obj <- Seurat::FindVariableFeatures(
  obj, assay = "RNA", selection.method = "vst", nfeatures = 1800, verbose = FALSE
)
ordering_genes <- unique(c(
  Seurat::VariableFeatures(obj), intersect(unlist(programs), rownames(counts))
))
gene_detection <- Matrix::rowSums(counts[ordering_genes, , drop = FALSE] > 0)
ordering_genes <- ordering_genes[gene_detection >= 20]
counts_traj <- counts[ordering_genes, , drop = FALSE]

cell_meta <- obj@meta.data[, c("condition_display", "cell_state_display", "sample_id"), drop = FALSE]
colnames(cell_meta) <- c("condition", "cell_state", "sample_id")
cell_meta$cell_state <- factor(as.character(cell_meta$cell_state), levels = state_display)
cell_meta$condition <- factor(as.character(cell_meta$condition), levels = c("CTRL", "CB"))

pd <- new("AnnotatedDataFrame", data = cell_meta)
fd <- new("AnnotatedDataFrame", data = data.frame(
  gene_short_name = rownames(counts_traj), row.names = rownames(counts_traj), stringsAsFactors = FALSE
))
cds <- monocle::newCellDataSet(
  counts_traj, phenoData = pd, featureData = fd,
  lowerDetectionLimit = 0.5,
  expressionFamily = VGAM::negbinomial.size()
)
sizeFactors(cds) <- library_size_factor[colnames(cds)]
cds <- monocle::detectGenes(cds, min_expr = 0.1)
if (length(ordering_genes) < 200) stop("Too few Monocle2 ordering genes after expression filtering.")
cds <- monocle::setOrderingFilter(cds, ordering_genes)

log_step("Monocle2 publication rebuild: running DDRTree on ", ncol(cds), " stromal/osteogenic cells")
cds <- monocle::reduceDimension(
  cds, max_components = 2, method = "DDRTree", num_dim = 20,
  norm_method = "log", pseudo_expr = 1, verbose = TRUE
)
cds <- monocle::orderCells(cds)

state_enrichment <- as.data.frame.matrix(table(pData(cds)$State, pData(cds)$cell_state))
if (!"LEPR+ MSC-like" %in% colnames(state_enrichment)) stop("LEPR+ MSC-like state missing.")
root_state <- rownames(state_enrichment)[which.max(state_enrichment[, "LEPR+ MSC-like"])]
cds <- monocle::orderCells(cds, root_state = as.numeric(as.character(root_state)))

pData(cds)$pseudotime_scaled <- pData(cds)$Pseudotime / max(pData(cds)$Pseudotime, na.rm = TRUE)
pt_breaks <- stats::quantile(pData(cds)$pseudotime_scaled, probs = c(0, 1/3, 2/3, 1),
                            na.rm = TRUE, names = FALSE)
pt_breaks <- unique(pt_breaks)
if (length(pt_breaks) == 4) {
  pData(cds)$pseudotime_stage <- cut(
    pData(cds)$pseudotime_scaled, breaks = pt_breaks, include.lowest = TRUE,
    labels = c("Early", "Middle", "Late")
  )
} else {
  pData(cds)$pseudotime_stage <- cut(
    rank(pData(cds)$pseudotime_scaled, ties.method = "first") / ncol(cds),
    breaks = c(0, 1/3, 2/3, 1), include.lowest = TRUE,
    labels = c("Early", "Middle", "Late")
  )
}

coord <- as.data.frame(t(monocle::reducedDimS(cds))) |>
  rownames_to_column("cell_id")
colnames(coord)[2:3] <- c("Component_1", "Component_2")
coord <- bind_cols(coord, as.data.frame(pData(cds)[coord$cell_id,
  c("condition", "cell_state", "State", "Pseudotime", "pseudotime_scaled", "pseudotime_stage"),
  drop = FALSE]))
write.csv(coord, file.path(src_dir, "source_data_B_Monocle2_DDRTree_pseudotime.csv"), row.names = FALSE)
write.csv(data.frame(
  ordering_gene = ordering_genes,
  detected_cells = as.integer(Matrix::rowSums(counts_traj[ordering_genes, , drop = FALSE] > 0))
),
          file.path(src_dir, "monocle2_ordering_genes.csv"), row.names = FALSE)
write.csv(state_enrichment, file.path(src_dir, "monocle2_state_by_cell_state_counts.csv"))

label_traj <- coord |>
  group_by(cell_state) |>
  summarise(Component_1 = median(Component_1), Component_2 = median(Component_2), .groups = "drop")

p_traj <- monocle::plot_cell_trajectory(
  cds, color_by = "Pseudotime", show_branch_points = FALSE,
  show_state_number = FALSE, cell_size = 0.72
) +
  scale_colour_gradientn(
    colours = c("#3C8DBC", "#8DD3C7", "#FFF7BC", "#F46D43", "#B30000"),
    limits = range(pData(cds)$Pseudotime, na.rm = TRUE), name = "Pseudotime"
  ) +
  ggrepel::geom_label_repel(
    data = label_traj,
    aes(Component_1, Component_2, label = cell_state),
    inherit.aes = FALSE, family = "Arial", fontface = "bold", size = 3.0,
    colour = "black", fill = scales::alpha("white", 0.82), label.size = 0,
    box.padding = 0.28, point.padding = 0.08, seed = seed,
    max.overlaps = Inf, min.segment.length = 0
  ) +
  theme_pub +
  theme(
    axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
    axis.line = element_blank(), legend.position = "right"
  ) +
  labs(
    title = "Exploratory Monocle2 ordering",
    subtitle = paste0("Rooted in DDRTree state ", root_state, " (LEPR+ MSC-like enriched)")
  )
save_plot(p_traj, "B_Monocle2_DDRTree_pseudotime_with_cell_state_labels", 175, 112)

count_df <- as.data.frame(table(
  condition = factor(obj$condition_display, levels = c("CTRL", "CB")),
  cell_state = factor(obj$cell_state_display, levels = state_display)
), stringsAsFactors = FALSE)
colnames(count_df)[3] <- "n_cells"
count_df <- count_df |>
  group_by(condition) |>
  mutate(total_stromal = sum(n_cells), percent_stromal = 100 * n_cells / total_stromal) |>
  ungroup()
write.csv(count_df, file.path(src_dir, "source_data_C_stromal_cell_counts.csv"), row.names = FALSE)

p_count <- ggplot(count_df, aes(n_cells, factor(cell_state, levels = rev(state_display)), fill = condition)) +
  geom_col(position = position_dodge(width = 0.78), width = 0.64) +
  geom_text(
    aes(label = scales::comma(n_cells)), position = position_dodge(width = 0.78),
    hjust = -0.16, family = "Arial", fontface = "bold", size = 3.4
  ) +
  scale_fill_manual(values = condition_colors, drop = FALSE) +
  scale_x_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.18))) +
  theme_pub +
  theme(axis.title.y = element_blank(), legend.position = "top") +
  labs(x = "Number of cells", fill = "Condition", title = "Stromal/osteogenic cell-state counts")
save_plot(p_count, "C_stromal_cell_state_counts_CB_CTRL", 165, 108)

stage_df <- as.data.frame(table(
  condition = pData(cds)$condition,
  pseudotime_stage = pData(cds)$pseudotime_stage
), stringsAsFactors = FALSE) |>
  group_by(condition) |>
  mutate(percent = 100 * Freq / sum(Freq)) |>
  ungroup()
write.csv(stage_df, file.path(src_dir, "source_data_C2_pseudotime_stage_composition.csv"), row.names = FALSE)

stage_cols <- c("Early" = "#3C8DBC", "Middle" = "#FFF7BC", "Late" = "#B30000")
p_stage <- ggplot(stage_df, aes(condition, percent, fill = pseudotime_stage)) +
  geom_col(width = 0.66) +
  geom_text(aes(label = sprintf("%.1f%%", percent)), position = position_stack(vjust = 0.5),
            family = "Arial", fontface = "bold", size = 3.6) +
  scale_fill_manual(values = stage_cols, drop = FALSE) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  theme_pub +
  theme(axis.title.x = element_blank(), legend.position = "top") +
  labs(y = "Percentage of stromal/osteogenic cells", fill = "Exploratory stage",
       title = "Monocle2 pseudotime-stage composition")
save_plot(p_stage, "C2_Monocle2_pseudotime_stage_composition", 122, 105)

log_data <- SeuratObject::LayerData(obj, assay = "RNA", layer = "data")
genes_present <- intersect(program_gene_table$gene, rownames(log_data))
genes_missing <- setdiff(program_gene_table$gene, genes_present)
write.csv(data.frame(gene = genes_missing), file.path(src_dir, "requested_program_genes_missing.csv"),
          row.names = FALSE)
if (length(genes_present) < 15) stop("Too few requested functional-program genes are present.")

program_gene_table_present <- program_gene_table |> filter(gene %in% genes_present)
state_cells <- split(colnames(obj), as.character(obj$cell_state_display))
state_cells <- state_cells[state_display]
state_avg <- vapply(state_cells, function(cc) {
  if (!length(cc)) return(rep(NA_real_, length(genes_present)))
  Matrix::rowMeans(log_data[genes_present, cc, drop = FALSE])
}, numeric(length(genes_present)))
rownames(state_avg) <- genes_present
colnames(state_avg) <- state_display
state_z <- zscore_rows(state_avg, 2.5)

state_pt <- coord |>
  group_by(cell_state) |>
  summarise(median_pseudotime = median(Pseudotime), .groups = "drop") |>
  arrange(median_pseudotime)
state_order_pt <- as.character(state_pt$cell_state)
state_order_pt <- state_order_pt[state_order_pt %in% colnames(state_z)]
state_z <- state_z[program_gene_table_present$gene, state_order_pt, drop = FALSE]
gene_program <- factor(program_gene_table_present$program, levels = program_order)

write.csv(
  as.data.frame(state_avg) |> rownames_to_column("gene"),
  file.path(src_dir, "source_data_D_state_average_log_expression.csv"), row.names = FALSE
)
write.csv(
  as.data.frame(state_z) |> rownames_to_column("gene"),
  file.path(src_dir, "source_data_D_state_gene_zscores.csv"), row.names = FALSE
)
write.csv(state_pt, file.path(src_dir, "source_data_state_median_pseudotime.csv"), row.names = FALSE)

heat_col <- circlize::colorRamp2(c(-2.5, 0, 2.5), c("#3C8DBC", "#F7F7F7", "#B30000"))
top_state <- HeatmapAnnotation(
  `Cell state` = factor(state_order_pt, levels = names(state_colors)),
  col = list(`Cell state` = state_colors),
  show_annotation_name = FALSE,
  annotation_legend_param = list(
    `Cell state` = list(title_gp = gpar(fontfamily = "Arial", fontface = "bold", fontsize = 10),
                        labels_gp = gpar(fontfamily = "Arial", fontface = "bold", fontsize = 8.5))
  )
)
left_program <- rowAnnotation(
  Program = gene_program,
  col = list(Program = program_colors),
  show_annotation_name = FALSE,
  annotation_legend_param = list(
    Program = list(title_gp = gpar(fontfamily = "Arial", fontface = "bold", fontsize = 10),
                   labels_gp = gpar(fontfamily = "Arial", fontface = "bold", fontsize = 8.5))
  ),
  width = unit(4, "mm")
)

ht_state <- Heatmap(
  state_z, name = "Gene-wise\nZ-score", col = heat_col,
  cluster_rows = FALSE, cluster_columns = FALSE,
  row_split = gene_program, row_gap = unit(2.1, "mm"),
  top_annotation = top_state, left_annotation = left_program,
  column_names_rot = 35, column_names_gp = gpar(fontfamily = "Arial", fontface = "bold", fontsize = 9),
  row_names_side = "right", row_names_gp = gpar(fontfamily = "Arial", fontface = "bold.italic", fontsize = 9),
  row_title_gp = gpar(fontfamily = "Arial", fontface = "bold", fontsize = 10),
  column_title = "Cell states ordered by median exploratory pseudotime",
  column_title_gp = gpar(fontfamily = "Arial", fontface = "bold", fontsize = 12),
  rect_gp = gpar(col = "white", lwd = 0.7),
  heatmap_legend_param = list(
    title_gp = gpar(fontfamily = "Arial", fontface = "bold", fontsize = 10),
    labels_gp = gpar(fontfamily = "Arial", fontface = "bold", fontsize = 9)
  )
)
save_heatmap(ht_state, "D_cell_state_by_functional_program_gene_heatmap", 190, 145)

log_step("Monocle2 publication rebuild: constructing pseudotime-gradient heatmap")
pt_df <- data.frame(
  cell_id = rownames(pData(cds)),
  pseudotime = pData(cds)$Pseudotime,
  pseudotime_scaled = pData(cds)$pseudotime_scaled,
  cell_state = as.character(pData(cds)$cell_state),
  condition = as.character(pData(cds)$condition),
  stringsAsFactors = FALSE
) |>
  arrange(pseudotime)
pt_df$bin <- pmin(90L, ceiling(seq_len(nrow(pt_df)) / nrow(pt_df) * 90L))

bin_info <- pt_df |>
  group_by(bin) |>
  summarise(
    pseudotime = mean(pseudotime),
    pseudotime_scaled = mean(pseudotime_scaled),
    cell_state = mode_value(cell_state),
    condition = mode_value(condition),
    n_cells = n(), .groups = "drop"
  )

expr_pt <- log_data[genes_present, pt_df$cell_id, drop = FALSE]
bin_levels <- sort(unique(pt_df$bin))
bin_avg <- vapply(bin_levels, function(b) {
  cc <- pt_df$cell_id[pt_df$bin == b]
  Matrix::rowMeans(expr_pt[, cc, drop = FALSE])
}, numeric(length(genes_present)))
rownames(bin_avg) <- genes_present
colnames(bin_avg) <- paste0("B", sprintf("%02d", bin_levels))

smoothed <- t(apply(bin_avg, 1, function(y) {
  x <- bin_info$pseudotime_scaled
  if (sum(is.finite(y)) < 5 || stats::sd(y, na.rm = TRUE) == 0) return(y)
  fit <- stats::smooth.spline(x = x, y = y, spar = 0.48)
  stats::predict(fit, x = x)$y
}))
rownames(smoothed) <- rownames(bin_avg)
colnames(smoothed) <- colnames(bin_avg)
pt_z <- zscore_rows(smoothed, 2.5)
pt_z <- pt_z[program_gene_table_present$gene, , drop = FALSE]

write.csv(bin_info, file.path(src_dir, "source_data_E_pseudotime_bin_annotations.csv"), row.names = FALSE)
write.csv(as.data.frame(bin_avg) |> rownames_to_column("gene"),
          file.path(src_dir, "source_data_E_pseudotime_binned_log_expression.csv"), row.names = FALSE)
write.csv(as.data.frame(pt_z) |> rownames_to_column("gene"),
          file.path(src_dir, "source_data_E_pseudotime_smoothed_gene_zscores.csv"), row.names = FALSE)

pt_annotation <- HeatmapAnnotation(
  `Cell state` = factor(bin_info$cell_state, levels = names(state_colors)),
  Pseudotime = bin_info$pseudotime_scaled,
  col = list(
    `Cell state` = state_colors,
    Pseudotime = circlize::colorRamp2(c(0, 0.5, 1), c("#3C8DBC", "#FFF7BC", "#B30000"))
  ),
  show_annotation_name = TRUE,
  annotation_name_gp = gpar(fontfamily = "Arial", fontface = "bold", fontsize = 9),
  annotation_legend_param = list(
    `Cell state` = list(title_gp = gpar(fontfamily = "Arial", fontface = "bold", fontsize = 10),
                        labels_gp = gpar(fontfamily = "Arial", fontface = "bold", fontsize = 8.5)),
    Pseudotime = list(title_gp = gpar(fontfamily = "Arial", fontface = "bold", fontsize = 10),
                      labels_gp = gpar(fontfamily = "Arial", fontface = "bold", fontsize = 9))
  )
)

left_program_pt <- rowAnnotation(
  Program = gene_program,
  col = list(Program = program_colors),
  show_annotation_name = FALSE,
  annotation_legend_param = list(
    Program = list(title_gp = gpar(fontfamily = "Arial", fontface = "bold", fontsize = 10),
                   labels_gp = gpar(fontfamily = "Arial", fontface = "bold", fontsize = 8.5))
  ), width = unit(4.5, "mm")
)

ht_pt <- Heatmap(
  pt_z, name = "Smoothed\nZ-score", col = heat_col,
  cluster_rows = FALSE, cluster_columns = FALSE,
  row_split = gene_program, row_gap = unit(2.2, "mm"),
  top_annotation = pt_annotation, left_annotation = left_program_pt,
  show_column_names = FALSE,
  row_names_side = "right", row_names_gp = gpar(fontfamily = "Arial", fontface = "bold.italic", fontsize = 9.5),
  row_title_gp = gpar(fontfamily = "Arial", fontface = "bold", fontsize = 10),
  column_title = "Exploratory Monocle2 pseudotime (early to late)",
  column_title_gp = gpar(fontfamily = "Arial", fontface = "bold", fontsize = 12),
  rect_gp = gpar(col = NA),
  heatmap_legend_param = list(
    title_gp = gpar(fontfamily = "Arial", fontface = "bold", fontsize = 10),
    labels_gp = gpar(fontfamily = "Arial", fontface = "bold", fontsize = 9)
  )
)
save_heatmap(ht_pt, "E_Monocle2_pseudotime_functional_program_gradient_heatmap", 205, 145)

state_summary <- coord |>
  group_by(condition, cell_state) |>
  summarise(
    n_cells = n(), median_pseudotime = median(Pseudotime),
    q25_pseudotime = quantile(Pseudotime, 0.25), q75_pseudotime = quantile(Pseudotime, 0.75),
    .groups = "drop"
  )
write.csv(state_summary, file.path(src_dir, "descriptive_state_pseudotime_summary.csv"), row.names = FALSE)

saveRDS(cds, file.path(root, "objects", "09_joint_stromal_monocle2_exploratory.rds"), compress = FALSE)

manifest <- tibble(
  item = c(
    "analysis_space_A", "analysis_space_B_E", "trajectory_method", "root_state",
    "ordering_gene_selection", "n_cells", "n_ordering_genes", "n_program_genes_present", "missing_program_genes",
    "interpretation_boundary"
  ),
  value = c(
    "RPCA integrated UMAP (visualization only)",
    "DecontX-corrected unintegrated RNA expression",
    "Monocle2 DDRTree exploratory ordering",
    root_state,
    "Top 1800 variable genes within the stromal/osteogenic subset plus requested program genes; condition label not used; genes detected in at least 20 cells",
    as.character(ncol(cds)), as.character(length(ordering_genes)),
    as.character(length(genes_present)), paste(genes_missing, collapse = ";"),
    "One donor per condition; descriptive only; inferred pseudotime does not establish lineage direction"
  )
)
write.csv(manifest, file.path(src_dir, "analysis_manifest.csv"), row.names = FALSE)

log_step("Monocle2 publication rebuild complete. Root state=", root_state,
         "; cells=", ncol(cds), "; output=", out_dir)
