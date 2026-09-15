root <- Sys.getenv("CB_REVISION_ROOT")
if (!nzchar(root)) root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
config <- yaml::read_yaml(file.path(root, "config", "analysis_config.yml"))
seed <- as.integer(config$seed)
set.seed(seed)

suppressPackageStartupMessages({
  library(Seurat)
  library(monocle)
  library(Matrix)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggrepel)
  library(ggrastr)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

out_dir <- file.path(
  root, "figures", "publication_subfigures", "myeloid_TREM2_OC_axis_revision"
)
src_dir <- file.path(
  root, "results", "trajectory", "monocle2_joint_myeloid_trem2_oc"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(src_dir, recursive = TRUE, showWarnings = FALSE)

# This analysis is deliberately restricted to the monocyte/macrophage/OC
# compartment. Neutrophils, mast cells and pDCs are not forced onto an OC axis.
state_order <- c("Inflam Mono", "Resident Mac", "TREM2+ Mac", "OC")
state_colors <- c(
  "Inflam Mono" = "#E68472",
  "Resident Mac" = "#79B8A9",
  "TREM2+ Mac" = "#9A78B4",
  "OC" = "#527BA8"
)
condition_colors <- c("CTRL" = "#4DBBD5", "CB" = "#E64B35")
stage_colors <- c("Early" = "#4DBBD5", "Middle" = "#F3D89A", "Late" = "#6A3D9A")
pseudotime_colors <- c("#4DBBD5", "#A8D8CE", "#F3D89A", "#E68472", "#6A3D9A")

# Transcript programs are descriptive. In particular, the TNF/TLR/MYD88 set
# is not labelled as pathway activation because RNA abundance alone is
# insufficient to establish signaling activity.
programs <- list(
  "Inflammatory monocyte" = c("FCN1", "S100A8", "S100A9", "IL1B", "CXCL8"),
  "Resident macrophage" = c("C1QA", "C1QB", "C1QC", "CD163", "MRC1", "LYVE1"),
  "TREM2 / lipid remodeling" = c(
    "TREM2", "APOE", "GPNMB", "SPP1", "HMOX1", "LIPA", "CTSD"
  ),
  "OC commitment / fusion" = c(
    "TNFRSF11A", "NFATC1", "DCSTAMP", "OCSTAMP", "ATP6V0D2"
  ),
  "OC resorption" = c("ACP5", "CTSK", "MMP9", "CA2", "TCIRG1", "CLCN7"),
  "TNF / TLR / MYD88 transcripts" = c("TNF", "TLR2", "TLR4", "MYD88", "NFKBIA")
)
program_order <- names(programs)
program_colors <- c(
  "Inflammatory monocyte" = "#E68472",
  "Resident macrophage" = "#79B8A9",
  "TREM2 / lipid remodeling" = "#9A78B4",
  "OC commitment / fusion" = "#DDA85C",
  "OC resorption" = "#527BA8",
  "TNF / TLR / MYD88 transcripts" = "#B55D78"
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
    plot.title = element_text(size = 13, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 9.3, face = "bold", colour = "#3A3A3A"),
    panel.grid = element_blank(),
    plot.margin = margin(7, 9, 7, 7)
  )

save_plot <- function(plot, filename, width_mm, height_mm, dpi = 600) {
  path <- file.path(out_dir, filename)
  wi <- width_mm / 25.4
  hi <- height_mm / 25.4

  svglite::svglite(paste0(path, ".svg"), width = wi, height = hi)
  print(plot)
  grDevices::dev.off()

  grDevices::cairo_pdf(
    paste0(path, ".pdf"), width = wi, height = hi, family = "Arial"
  )
  print(plot)
  grDevices::dev.off()

  ragg::agg_tiff(
    paste0(path, ".tiff"), width = wi, height = hi, units = "in",
    res = dpi, compression = "lzw"
  )
  print(plot)
  grDevices::dev.off()

  ragg::agg_png(
    paste0(path, ".png"), width = wi, height = hi, units = "in", res = 300
  )
  print(plot)
  grDevices::dev.off()
}

save_plot_exact_pixels <- function(
  plot, filename, width_px, height_px, dpi = 600
) {
  path <- file.path(out_dir, filename)
  wi <- width_px / dpi
  hi <- height_px / dpi

  svglite::svglite(paste0(path, ".svg"), width = wi, height = hi)
  print(plot)
  grDevices::dev.off()

  grDevices::cairo_pdf(
    paste0(path, ".pdf"), width = wi, height = hi, family = "Arial"
  )
  print(plot)
  grDevices::dev.off()

  ragg::agg_tiff(
    paste0(path, ".tiff"), width = width_px, height = height_px,
    units = "px", res = dpi, compression = "lzw"
  )
  print(plot)
  grDevices::dev.off()

  ragg::agg_png(
    paste0(path, ".png"), width = width_px, height = height_px,
    units = "px", res = dpi
  )
  print(plot)
  grDevices::dev.off()
}

save_heatmap <- function(
  ht, filename, width_mm, height_mm, dpi = 600,
  bottom_pad_mm = 5, left_pad_mm = 12
) {
  path <- file.path(out_dir, filename)
  wi <- width_mm / 25.4
  hi <- height_mm / 25.4
  draw_ht <- function() ComplexHeatmap::draw(
    ht,
    heatmap_legend_side = "right",
    annotation_legend_side = "right",
    merge_legends = FALSE,
    padding = unit(c(5, 8, bottom_pad_mm, left_pad_mm), "mm")
  )

  svglite::svglite(paste0(path, ".svg"), width = wi, height = hi)
  draw_ht()
  grDevices::dev.off()

  grDevices::cairo_pdf(
    paste0(path, ".pdf"), width = wi, height = hi, family = "Arial"
  )
  draw_ht()
  grDevices::dev.off()

  ragg::agg_tiff(
    paste0(path, ".tiff"), width = wi, height = hi, units = "in",
    res = dpi, compression = "lzw"
  )
  draw_ht()
  grDevices::dev.off()

  ragg::agg_png(
    paste0(path, ".png"), width = wi, height = hi, units = "in", res = 300
  )
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

wrap_state <- function(x) {
  x <- as.character(x)
  x[x == "Inflam Mono"] <- "Inflam\nMono"
  x[x == "Resident Mac"] <- "Resident\nMac"
  x[x == "TREM2+ Mac"] <- "TREM2+\nMac"
  x
}

message("Loading the reviewed RPCA visualization object...")
rpca <- readRDS(file.path(root, "objects", "07_rpca_broad_lineages_reviewed.rds"))
rpca_state <- as.character(rpca$celltype_reviewed_provisional)
axis_cells <- colnames(rpca)[rpca_state %in% state_order]
if (length(axis_cells) < 200) stop("Too few cells in the prespecified myeloid TREM2-OC subset.")

joint_umap <- as.data.frame(
  Seurat::Embeddings(rpca, reduction = "umap.rpca")[axis_cells, , drop = FALSE]
) |>
  rownames_to_column("cell_id")
colnames(joint_umap)[2:3] <- c("UMAP_1", "UMAP_2")
joint_umap$cell_state <- factor(
  as.character(rpca$celltype_reviewed_provisional[axis_cells]),
  levels = state_order
)
joint_umap$condition <- factor(as.character(rpca$group[axis_cells]), levels = c("CTRL", "CB"))
joint_umap <- joint_umap |> filter(!is.na(cell_state), !is.na(condition))
write.csv(
  joint_umap, file.path(src_dir, "source_data_A_joint_myeloid_RPCA_UMAP.csv"),
  row.names = FALSE
)

label_umap <- joint_umap |>
  group_by(cell_state) |>
  summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2), .groups = "drop")

p_umap <- ggplot(joint_umap, aes(UMAP_1, UMAP_2, colour = cell_state)) +
  ggrastr::geom_point_rast(size = 0.55, alpha = 0.90, raster.dpi = 600) +
  scale_colour_manual(values = state_colors, drop = FALSE) +
  ggrepel::geom_label_repel(
    data = label_umap,
    aes(UMAP_1, UMAP_2, label = cell_state),
    inherit.aes = FALSE,
    family = "Arial", fontface = "bold", size = 3.25,
    colour = "black", fill = scales::alpha("white", 0.82), label.size = 0,
    box.padding = 0.32, point.padding = 0.09, seed = seed,
    max.overlaps = Inf, min.segment.length = 0
  ) +
  coord_equal() +
  theme_pub +
  theme(
    axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
    axis.line = element_blank(), legend.position = "none"
  ) +
  annotate(
    "segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf,
    arrow = arrow(length = unit(2.4, "mm")), linewidth = 0.55
  ) +
  annotate(
    "segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf,
    arrow = arrow(length = unit(2.4, "mm")), linewidth = 0.55
  ) +
  annotate(
    "text", x = Inf, y = -Inf, label = "UMAP1", hjust = 1.1, vjust = 1.7,
    family = "Arial", fontface = "bold", size = 3.8
  ) +
  annotate(
    "text", x = -Inf, y = Inf, label = "UMAP2", angle = 90,
    hjust = 1.1, vjust = -0.5, family = "Arial", fontface = "bold", size = 3.8
  ) +
  labs(title = "Joint monocyte-macrophage-osteoclast landscape")
save_plot(
  p_umap, "A_joint_myeloid_RPCA_UMAP_cell_state_no_trajectory", 170, 112
)

message("Computing the integrated-RPCA distance to the pooled OC centroid...")
dims_use <- seq_len(min(
  as.integer(config$dimension_reduction$dimensions),
  ncol(Seurat::Embeddings(rpca, reduction = "integrated.rpca"))
))
axis_emb <- Seurat::Embeddings(rpca, reduction = "integrated.rpca")[
  joint_umap$cell_id, dims_use, drop = FALSE
]
oc_rows <- as.character(joint_umap$cell_state) == "OC"
if (sum(oc_rows) < 20) stop("Too few OC cells to define a stable pooled centroid.")
oc_centroid <- colMeans(axis_emb[oc_rows, , drop = FALSE])
distance_to_oc <- sqrt(rowSums(sweep(axis_emb, 2, oc_centroid, FUN = "-")^2))
distance_df <- joint_umap[, c("cell_id", "condition", "cell_state")]
distance_df$distance_to_OC_centroid <- as.numeric(distance_to_oc[distance_df$cell_id])
distance_df$cell_state <- factor(as.character(distance_df$cell_state), levels = state_order)
distance_df$condition <- factor(as.character(distance_df$condition), levels = c("CTRL", "CB"))

distance_summary <- distance_df |>
  group_by(condition, cell_state) |>
  summarise(
    n_cells = n(),
    median_distance = median(distance_to_OC_centroid),
    q25_distance = quantile(distance_to_OC_centroid, 0.25),
    q75_distance = quantile(distance_to_OC_centroid, 0.75),
    .groups = "drop"
  )
write.csv(
  distance_df, file.path(src_dir, "source_data_H_distance_to_OC_centroid_per_cell.csv"),
  row.names = FALSE
)
write.csv(
  distance_summary,
  file.path(src_dir, "source_data_H_distance_to_OC_centroid_summary.csv"),
  row.names = FALSE
)
write.csv(
  data.frame(dimension = colnames(axis_emb), OC_centroid = as.numeric(oc_centroid)),
  file.path(src_dir, "source_data_H_pooled_OC_centroid_integrated_RPCA.csv"),
  row.names = FALSE
)

p_distance <- ggplot(
  distance_df,
  aes(
    cell_state, distance_to_OC_centroid, fill = condition,
    group = interaction(cell_state, condition)
  )
) +
  geom_violin(
    position = position_dodge(width = 0.82), width = 0.78,
    scale = "width", trim = TRUE, alpha = 0.88,
    colour = "#333333", linewidth = 0.38
  ) +
  geom_boxplot(
    position = position_dodge(width = 0.82), width = 0.14,
    outlier.shape = NA, fill = "white", alpha = 0.92,
    colour = "#333333", linewidth = 0.42
  ) +
  scale_fill_manual(
    values = condition_colors, breaks = c("CTRL", "CB"), drop = FALSE
  ) +
  scale_x_discrete(labels = wrap_state) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
  theme_pub +
  theme(legend.position = "top", axis.title.x = element_blank()) +
  labs(
    y = "Distance to pooled OC centroid\n(integrated-RPCA space)",
    fill = "Condition",
    title = "Myeloid distance to the OC centroid",
    subtitle = paste0(length(dims_use), " RPCA dimensions; descriptive per-cell distributions")
  )
save_plot(p_distance, "H_Myeloid_distance_to_OC_centroid_violin", 178, 112)

rm(rpca, axis_emb)
invisible(gc())

message("Loading corrected unintegrated expression for Monocle2...")
obj <- readRDS(file.path(root, "objects", "06_base_broad_lineages_reviewed.rds"))
cells <- joint_umap$cell_id[joint_umap$cell_id %in% colnames(obj)]
obj <- subset(obj, cells = cells)
obj$cell_state_display <- factor(
  as.character(obj$celltype_reviewed_provisional), levels = state_order
)
obj$condition_display <- factor(as.character(obj$group), levels = c("CTRL", "CB"))

counts <- SeuratObject::LayerData(obj, assay = "RNA", layer = "counts")
counts <- counts[, colnames(obj), drop = FALSE]
library_size <- Matrix::colSums(counts)
library_size_factor <- library_size / stats::median(library_size[library_size > 0])

obj <- Seurat::FindVariableFeatures(
  obj, assay = "RNA", selection.method = "vst", nfeatures = 1800, verbose = FALSE
)
ordering_genes <- unique(c(
  Seurat::VariableFeatures(obj), intersect(unlist(programs), rownames(counts))
))
gene_detection <- Matrix::rowSums(counts[ordering_genes, , drop = FALSE] > 0)
ordering_genes <- ordering_genes[gene_detection >= 10]
if (length(ordering_genes) < 200) stop("Too few Monocle2 ordering genes after filtering.")
counts_traj <- counts[ordering_genes, , drop = FALSE]

cell_meta <- obj@meta.data[, c(
  "condition_display", "cell_state_display", "sample_id"
), drop = FALSE]
colnames(cell_meta) <- c("condition", "cell_state", "sample_id")
cell_meta$cell_state <- factor(as.character(cell_meta$cell_state), levels = state_order)
cell_meta$condition <- factor(as.character(cell_meta$condition), levels = c("CTRL", "CB"))

pd <- new("AnnotatedDataFrame", data = cell_meta)
fd <- new("AnnotatedDataFrame", data = data.frame(
  gene_short_name = rownames(counts_traj),
  row.names = rownames(counts_traj), stringsAsFactors = FALSE
))
cds <- monocle::newCellDataSet(
  counts_traj,
  phenoData = pd,
  featureData = fd,
  lowerDetectionLimit = 0.5,
  expressionFamily = VGAM::negbinomial.size()
)
sizeFactors(cds) <- library_size_factor[colnames(cds)]
cds <- monocle::detectGenes(cds, min_expr = 0.1)
cds <- monocle::setOrderingFilter(cds, ordering_genes)

message("Running Monocle2 DDRTree on ", ncol(cds), " cells...")
cds <- monocle::reduceDimension(
  cds, max_components = 2, method = "DDRTree", num_dim = 20,
  norm_method = "log", pseudo_expr = 1, verbose = TRUE
)
cds <- monocle::orderCells(cds)

state_tab <- as.matrix(table(pData(cds)$State, pData(cds)$cell_state))
if (!"Inflam Mono" %in% colnames(state_tab)) stop("Inflam Mono is absent after DDRTree.")
root_audit <- tibble(
  State = rownames(state_tab),
  n_total = as.numeric(rowSums(state_tab)),
  inflam_mono_n = as.numeric(state_tab[, "Inflam Mono"]),
  inflam_mono_fraction = inflam_mono_n / n_total
) |>
  arrange(desc(inflam_mono_fraction), desc(inflam_mono_n), desc(n_total))
eligible_roots <- root_audit |> filter(n_total >= 20)
if (!nrow(eligible_roots)) eligible_roots <- root_audit
root_state <- eligible_roots$State[1]
cds <- monocle::orderCells(cds, root_state = as.numeric(as.character(root_state)))

pData(cds)$pseudotime_scaled <- pData(cds)$Pseudotime /
  max(pData(cds)$Pseudotime, na.rm = TRUE)
rank_fraction <- rank(pData(cds)$pseudotime_scaled, ties.method = "first") / ncol(cds)
pData(cds)$pseudotime_stage <- cut(
  rank_fraction,
  breaks = c(0, 1 / 3, 2 / 3, 1), include.lowest = TRUE,
  labels = c("Early", "Middle", "Late")
)

coord <- as.data.frame(t(monocle::reducedDimS(cds))) |>
  rownames_to_column("cell_id")
colnames(coord)[2:3] <- c("Component_1", "Component_2")
coord <- bind_cols(
  coord,
  as.data.frame(pData(cds)[
    coord$cell_id,
    c(
      "condition", "cell_state", "State", "Pseudotime",
      "pseudotime_scaled", "pseudotime_stage"
    ),
    drop = FALSE
  ])
)
write.csv(
  coord, file.path(src_dir, "source_data_B_Monocle2_DDRTree_pseudotime.csv"),
  row.names = FALSE
)
write.csv(root_audit, file.path(src_dir, "monocle2_root_state_audit.csv"), row.names = FALSE)
write.csv(
  data.frame(
    ordering_gene = ordering_genes,
    detected_cells = as.integer(Matrix::rowSums(counts_traj > 0))
  ),
  file.path(src_dir, "monocle2_ordering_genes.csv"), row.names = FALSE
)
write.csv(
  as.data.frame.matrix(state_tab) |> rownames_to_column("State"),
  file.path(src_dir, "monocle2_DDRTree_state_by_cell_state_counts.csv"),
  row.names = FALSE
)

label_traj <- coord |>
  group_by(cell_state) |>
  summarise(
    Component_1 = median(Component_1), Component_2 = median(Component_2),
    .groups = "drop"
  )

p_traj <- monocle::plot_cell_trajectory(
  cds,
  color_by = "Pseudotime",
  show_branch_points = FALSE,
  show_state_number = FALSE,
  cell_size = 0.95
) +
  scale_colour_gradientn(
    colours = pseudotime_colors,
    limits = range(pData(cds)$Pseudotime, na.rm = TRUE),
    name = "Pseudotime"
  ) +
  ggrepel::geom_label_repel(
    data = label_traj,
    aes(Component_1, Component_2, label = cell_state),
    inherit.aes = FALSE,
    family = "Arial", fontface = "bold", size = 3.2,
    colour = "black", fill = scales::alpha("white", 0.84), label.size = 0,
    box.padding = 0.32, point.padding = 0.09, seed = seed,
    max.overlaps = Inf, min.segment.length = 0
  ) +
  theme_pub +
  theme(
    axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
    axis.line = element_blank(), legend.position = "right"
  ) +
  labs(
    title = "Exploratory Monocle2 ordering of the myeloid TREM2-OC compartment",
    subtitle = paste0(
      "Rooted in DDRTree state ", root_state,
      " (highest eligible Inflam Mono fraction)"
    )
  )
save_plot(
  p_traj, "B_Monocle2_DDRTree_pseudotime_with_cell_state_labels", 182, 112
)

count_df <- as.data.frame(table(
  condition = factor(obj$condition_display, levels = c("CTRL", "CB")),
  cell_state = factor(obj$cell_state_display, levels = state_order)
), stringsAsFactors = FALSE)
colnames(count_df)[3] <- "n_cells"
count_df <- count_df |>
  group_by(condition) |>
  mutate(total_axis_cells = sum(n_cells), percent_axis_cells = 100 * n_cells / total_axis_cells) |>
  ungroup()
count_df$condition <- factor(as.character(count_df$condition), levels = c("CTRL", "CB"))
count_df$cell_state <- factor(as.character(count_df$cell_state), levels = state_order)
write.csv(
  count_df, file.path(src_dir, "source_data_C_myeloid_cell_state_counts.csv"),
  row.names = FALSE
)

p_count <- ggplot(
  count_df,
  aes(n_cells, factor(cell_state, levels = rev(state_order)), fill = condition)
) +
  geom_col(position = position_dodge(width = 0.78), width = 0.64) +
  geom_text(
    aes(label = scales::comma(n_cells)),
    position = position_dodge(width = 0.78), hjust = -0.16,
    family = "Arial", fontface = "bold", size = 3.5
  ) +
  scale_fill_manual(
    values = condition_colors, breaks = c("CTRL", "CB"), drop = FALSE
  ) +
  scale_x_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.20))) +
  scale_y_discrete(labels = wrap_state) +
  theme_pub +
  theme(axis.title.y = element_blank(), legend.position = "top") +
  labs(
    x = "Number of cells", fill = "Condition",
    title = "Myeloid TREM2-OC cell-state counts"
  )
save_plot(p_count, "C_myeloid_cell_state_counts_CB_CTRL", 165, 105)

stage_df <- as.data.frame(table(
  condition = pData(cds)$condition,
  pseudotime_stage = pData(cds)$pseudotime_stage
), stringsAsFactors = FALSE) |>
  group_by(condition) |>
  mutate(percent = 100 * Freq / sum(Freq)) |>
  ungroup()
stage_df$condition <- factor(as.character(stage_df$condition), levels = c("CTRL", "CB"))
stage_df$pseudotime_stage <- factor(
  as.character(stage_df$pseudotime_stage), levels = c("Early", "Middle", "Late")
)
write.csv(
  stage_df, file.path(src_dir, "source_data_C2_pseudotime_stage_composition.csv"),
  row.names = FALSE
)

p_stage <- ggplot(stage_df, aes(condition, percent, fill = pseudotime_stage)) +
  geom_col(width = 0.66) +
  geom_text(
    aes(label = sprintf("%.1f%%", percent)),
    position = position_stack(vjust = 0.5),
    family = "Arial", fontface = "bold", size = 3.6
  ) +
  scale_fill_manual(
    values = stage_colors, breaks = c("Early", "Middle", "Late"), drop = FALSE
  ) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  theme_pub +
  theme(axis.title.x = element_blank(), legend.position = "top") +
  labs(
    y = "Percentage of TREM2-OC axis cells",
    fill = "Exploratory stage",
    title = "Monocle2 pseudotime-stage composition"
  )
save_plot(p_stage, "C2_Monocle2_pseudotime_stage_composition", 128, 105)

message("Building state and pseudotime functional-program heatmaps...")
log_data <- SeuratObject::LayerData(obj, assay = "RNA", layer = "data")
genes_present <- intersect(gene_program$gene, rownames(log_data))
genes_missing <- setdiff(gene_program$gene, genes_present)
gene_program_present <- gene_program |> filter(gene %in% genes_present)
write.csv(
  tibble(
    program = program_order,
    requested_genes = vapply(programs, paste, collapse = ";", character(1)),
    present_genes = vapply(
      programs,
      function(g) paste(intersect(g, genes_present), collapse = ";"),
      character(1)
    ),
    n_present = vapply(
      programs, function(g) length(intersect(g, genes_present)), integer(1)
    )
  ),
  file.path(src_dir, "functional_program_gene_definitions.csv"), row.names = FALSE
)
write.csv(
  data.frame(gene = genes_missing),
  file.path(src_dir, "requested_program_genes_missing.csv"), row.names = FALSE
)
if (length(genes_present) < 20) stop("Too few requested myeloid program genes are present.")

state_cells <- split(colnames(obj), as.character(obj$cell_state_display))[state_order]
state_avg <- vapply(state_cells, function(cc) {
  Matrix::rowMeans(log_data[genes_present, cc, drop = FALSE])
}, numeric(length(genes_present)))
rownames(state_avg) <- genes_present
colnames(state_avg) <- state_order
state_z_all <- zscore_rows(state_avg, 2.5)
state_z <- state_z_all[gene_program_present$gene, state_order, drop = FALSE]
row_program <- factor(gene_program_present$program, levels = program_order)

write.csv(
  as.data.frame(state_avg) |> rownames_to_column("gene"),
  file.path(src_dir, "source_data_D_state_average_log_expression.csv"),
  row.names = FALSE
)
write.csv(
  as.data.frame(state_z) |> rownames_to_column("gene"),
  file.path(src_dir, "source_data_D_state_gene_zscores.csv"),
  row.names = FALSE
)

heat_col_state <- circlize::colorRamp2(
  c(-1.5, 0, 1.5), c("#527BA8", "#F7F7F7", "#B55D78")
)
heat_col_pt <- circlize::colorRamp2(
  c(-2.5, 0, 2.5), c("#527BA8", "#F7F7F7", "#B55D78")
)
legend_gp <- list(
  title_gp = gpar(fontfamily = "Arial", fontface = "bold", fontsize = 10),
  labels_gp = gpar(fontfamily = "Arial", fontface = "bold", fontsize = 8.5)
)
top_state <- HeatmapAnnotation(
  `Cell state` = factor(state_order, levels = names(state_colors)),
  col = list(`Cell state` = state_colors),
  show_annotation_name = FALSE,
  annotation_legend_param = list(`Cell state` = legend_gp)
)
left_program <- rowAnnotation(
  Program = row_program,
  col = list(Program = program_colors),
  show_annotation_name = FALSE,
  annotation_legend_param = list(Program = legend_gp),
  width = unit(4.5, "mm")
)
ht_state <- Heatmap(
  state_z,
  name = "Gene-wise\nZ-score",
  col = heat_col_state,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_split = row_program,
  row_title = NULL,
  row_gap = unit(1.8, "mm"),
  top_annotation = top_state,
  left_annotation = left_program,
  column_labels = wrap_state(state_order),
  column_names_rot = 0,
  column_names_gp = gpar(
    fontfamily = "Arial", fontface = "bold", fontsize = 10
  ),
  row_names_side = "right",
  row_names_gp = gpar(
    fontfamily = "Arial", fontface = "bold.italic", fontsize = 8.7
  ),
  column_title = "Myeloid transcript programs by annotated state",
  column_title_gp = gpar(
    fontfamily = "Arial", fontface = "bold", fontsize = 11.5
  ),
  rect_gp = gpar(col = "white", lwd = 0.65),
  heatmap_legend_param = c(legend_gp, list(at = c(-1.5, 0, 1.5)))
)
save_heatmap(
  ht_state, "D_cell_state_by_TREM2_OC_program_gene_heatmap", 200, 170,
  bottom_pad_mm = 14
)

pt_df <- data.frame(
  cell_id = coord$cell_id,
  pseudotime = coord$Pseudotime,
  pseudotime_scaled = coord$pseudotime_scaled,
  cell_state = as.character(coord$cell_state),
  condition = as.character(coord$condition),
  stringsAsFactors = FALSE
) |>
  arrange(pseudotime)
pt_df$bin <- pmin(75L, ceiling(seq_len(nrow(pt_df)) / nrow(pt_df) * 75L))
bin_info <- pt_df |>
  group_by(bin) |>
  summarise(
    pseudotime = mean(pseudotime),
    pseudotime_scaled = mean(pseudotime_scaled),
    cell_state = mode_value(cell_state),
    condition = mode_value(condition),
    n_cells = n(),
    .groups = "drop"
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
  if (stats::sd(y, na.rm = TRUE) == 0) return(y)
  fit <- smooth.spline(x = x, y = y, spar = 0.50)
  predict(fit, x = x)$y
}))
rownames(smoothed) <- rownames(bin_avg)
colnames(smoothed) <- colnames(bin_avg)
pt_z_all <- zscore_rows(smoothed, 2.5)
pt_z <- pt_z_all[gene_program_present$gene, , drop = FALSE]

write.csv(
  bin_info, file.path(src_dir, "source_data_E_pseudotime_bin_annotations.csv"),
  row.names = FALSE
)
write.csv(
  as.data.frame(bin_avg) |> rownames_to_column("gene"),
  file.path(src_dir, "source_data_E_pseudotime_binned_log_expression.csv"),
  row.names = FALSE
)
write.csv(
  as.data.frame(pt_z) |> rownames_to_column("gene"),
  file.path(src_dir, "source_data_E_pseudotime_smoothed_gene_zscores.csv"),
  row.names = FALSE
)

pt_annotation <- HeatmapAnnotation(
  Condition = factor(bin_info$condition, levels = c("CTRL", "CB")),
  `Cell state` = factor(bin_info$cell_state, levels = state_order),
  Pseudotime = bin_info$pseudotime_scaled,
  col = list(
    Condition = condition_colors,
    `Cell state` = state_colors,
    Pseudotime = circlize::colorRamp2(
      c(0, 0.5, 1), c("#4DBBD5", "#F3D89A", "#6A3D9A")
    )
  ),
  show_annotation_name = TRUE,
  annotation_name_gp = gpar(
    fontfamily = "Arial", fontface = "bold", fontsize = 9
  ),
  annotation_legend_param = list(
    Condition = legend_gp, `Cell state` = legend_gp, Pseudotime = legend_gp
  )
)
left_program_pt <- rowAnnotation(
  Program = row_program,
  col = list(Program = program_colors),
  show_annotation_name = FALSE,
  annotation_legend_param = list(Program = legend_gp),
  width = unit(4.5, "mm")
)
ht_pt <- Heatmap(
  pt_z,
  name = "Smoothed\nZ-score",
  col = heat_col_pt,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_split = row_program,
  row_title = NULL,
  row_gap = unit(1.8, "mm"),
  top_annotation = pt_annotation,
  left_annotation = left_program_pt,
  show_column_names = FALSE,
  row_names_side = "right",
  row_names_gp = gpar(
    fontfamily = "Arial", fontface = "bold.italic", fontsize = 8.8
  ),
  column_title = "Myeloid transcript programs over exploratory pseudotime",
  column_title_gp = gpar(
    fontfamily = "Arial", fontface = "bold", fontsize = 11.5
  ),
  rect_gp = gpar(col = NA),
  heatmap_legend_param = c(legend_gp, list(at = c(-2.5, 0, 2.5)))
)
save_heatmap(
  ht_pt, "E_Monocle2_pseudotime_TREM2_OC_program_gradient_heatmap",
  215, 172
)

message("Building descriptive program-score bar plots...")
score_df <- lapply(program_order, function(pr) {
  gg <- intersect(programs[[pr]], rownames(state_z_all))
  tibble(
    program = pr,
    cell_state = state_order,
    score = colMeans(state_z_all[gg, state_order, drop = FALSE]),
    n_genes = length(gg)
  )
}) |>
  bind_rows()
score_df$program <- factor(score_df$program, levels = program_order)
score_df$cell_state <- factor(score_df$cell_state, levels = state_order)
program_labeller <- ggplot2::as_labeller(c(
  "Inflammatory monocyte" = "Inflammatory\nmonocyte",
  "Resident macrophage" = "Resident\nmacrophage",
  "TREM2 / lipid remodeling" = "TREM2 / lipid\nremodeling",
  "OC commitment / fusion" = "OC commitment /\nfusion",
  "OC resorption" = "OC resorption",
  "TNF / TLR / MYD88 transcripts" = "TNF / TLR / MYD88\ntranscripts"
))
write.csv(
  score_df, file.path(src_dir, "source_data_F_program_scores_by_cell_state.csv"),
  row.names = FALSE
)

p_score <- ggplot(score_df, aes(cell_state, score, fill = cell_state)) +
  geom_hline(yintercept = 0, colour = "#555555", linewidth = 0.4) +
  geom_col(width = 0.72) +
  facet_wrap(
    ~program, ncol = 2, nrow = 3, scales = "free_y",
    labeller = program_labeller
  ) +
  scale_fill_manual(values = state_colors, breaks = state_order, drop = FALSE) +
  scale_x_discrete(labels = wrap_state) +
  guides(
    fill = guide_legend(
      title = "Cell state", nrow = 2, byrow = TRUE,
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
  "F_TREM2_OC_program_score_barplots_by_cell_state",
  width_px = 2700,
  height_px = 2980,
  dpi = 600
)

message("Building key-gene pseudotime curves...")
curve_genes <- c(
  "FCN1", "IL1B", "C1QA", "MRC1", "TREM2", "GPNMB",
  "SPP1", "TNFRSF11A", "NFATC1", "ACP5", "CTSK", "MMP9"
)
curve_genes <- intersect(curve_genes, rownames(log_data))
curve_levels <- curve_genes

cell_plot_ids <- pt_df |>
  group_by(cell_state) |>
  group_modify(~slice_sample(.x, n = min(600L, nrow(.x)), replace = FALSE)) |>
  ungroup() |>
  pull(cell_id)
cell_plot_meta <- pt_df[match(cell_plot_ids, pt_df$cell_id), ]
point_mat <- as.matrix(log_data[curve_genes, cell_plot_ids, drop = FALSE])
point_df <- as.data.frame(t(point_mat)) |>
  rownames_to_column("cell_id") |>
  left_join(
    cell_plot_meta[, c("cell_id", "pseudotime_scaled", "cell_state")],
    by = "cell_id"
  ) |>
  pivot_longer(
    cols = all_of(curve_genes), names_to = "gene", values_to = "expression"
  )

curve_bin_df <- as.data.frame(t(bin_avg[curve_genes, , drop = FALSE])) |>
  mutate(pseudotime_scaled = bin_info$pseudotime_scaled) |>
  pivot_longer(
    cols = all_of(curve_genes),
    names_to = "gene", values_to = "binned_expression"
  )
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
      expression = pmax(
        0,
        as.numeric(predict(fit, newdata = data.frame(pseudotime_scaled = xx)))
      )
    )
  }) |>
  ungroup()

point_df$gene <- factor(point_df$gene, levels = curve_levels)
point_df$cell_state <- factor(as.character(point_df$cell_state), levels = state_order)
curve_line_df$gene <- factor(curve_line_df$gene, levels = curve_levels)
write.csv(
  curve_bin_df, file.path(src_dir, "source_data_G_key_gene_binned_expression.csv"),
  row.names = FALSE
)
write.csv(
  curve_line_df,
  file.path(src_dir, "source_data_G_key_gene_smoothed_expression_curves.csv"),
  row.names = FALSE
)
write.csv(
  point_df |>
    group_by(gene, cell_state) |>
    summarise(n_plotted_cells = n(), .groups = "drop"),
  file.path(src_dir, "source_data_G_curve_point_sampling_summary.csv"),
  row.names = FALSE
)

p_curve <- ggplot(
  point_df, aes(pseudotime_scaled, expression, colour = cell_state)
) +
  ggrastr::geom_point_rast(size = 0.34, alpha = 0.20, raster.dpi = 600) +
  geom_line(
    data = curve_line_df,
    aes(pseudotime_scaled, expression, group = gene),
    inherit.aes = FALSE,
    colour = "black", linewidth = 0.82, lineend = "round"
  ) +
  facet_wrap(~gene, nrow = 2, scales = "free_y") +
  scale_colour_manual(values = state_colors, breaks = state_order, drop = FALSE) +
  guides(
    colour = guide_legend(
      override.aes = list(size = 4.5, alpha = 1),
      keyheight = unit(5.5, "mm"), keywidth = unit(5.5, "mm")
    )
  ) +
  scale_x_continuous(
    breaks = c(0, 0.5, 1), expand = expansion(mult = c(0.02, 0.03))
  ) +
  theme_pub +
  theme(
    legend.position = "right",
    axis.text = element_text(size = 8.3),
    strip.text = element_text(size = 10.5, face = "bold.italic"),
    panel.spacing.x = unit(4.3, "mm"),
    panel.spacing.y = unit(5.2, "mm")
  ) +
  labs(
    x = "Exploratory Monocle2 pseudotime (scaled)",
    y = "Log-normalized expression",
    colour = "Cell state",
    title = "Key-gene expression patterns over exploratory myeloid pseudotime"
  )
save_plot_exact_pixels(
  p_curve,
  "G_key_gene_expression_curves_over_Monocle2_pseudotime",
  width_px = 6300,
  height_px = 2750,
  dpi = 600
)

state_summary <- coord |>
  group_by(condition, cell_state) |>
  summarise(
    n_cells = n(),
    median_pseudotime = median(Pseudotime),
    q25_pseudotime = quantile(Pseudotime, 0.25),
    q75_pseudotime = quantile(Pseudotime, 0.75),
    median_scaled_pseudotime = median(pseudotime_scaled),
    .groups = "drop"
  )
write.csv(
  state_summary, file.path(src_dir, "descriptive_state_pseudotime_summary.csv"),
  row.names = FALSE
)
write.csv(
  coord |>
    group_by(cell_state) |>
    summarise(
      n_cells = n(), median_pseudotime = median(Pseudotime),
      median_scaled_pseudotime = median(pseudotime_scaled), .groups = "drop"
    ),
  file.path(src_dir, "source_data_state_median_pseudotime.csv"),
  row.names = FALSE
)

saveRDS(
  cds,
  file.path(root, "objects", "10_joint_myeloid_trem2_oc_monocle2_exploratory.rds"),
  compress = FALSE
)

manifest <- tibble(
  item = c(
    "analysis_space_A_H", "analysis_space_B_G", "included_states",
    "excluded_myeloid_states", "trajectory_method", "root_state",
    "root_rule", "ordering_gene_selection", "n_cells", "n_ordering_genes",
    "n_program_genes_present", "missing_program_genes",
    "OC_centroid_definition", "condition_comparison",
    "interpretation_boundary"
  ),
  value = c(
    "RPCA integrated space (visualization and geometric distance only)",
    "DecontX-corrected unintegrated RNA expression",
    paste(state_order, collapse = ";"),
    "Neutrophils;Mast cells;pDC",
    "Monocle2 DDRTree exploratory ordering",
    as.character(root_state),
    "Eligible DDRTree state (>=20 cells) with the highest Inflam Mono fraction",
    "Top 1800 variable genes in the four-state subset plus requested program genes; condition not used; detected in >=10 cells",
    as.character(ncol(cds)),
    as.character(length(ordering_genes)),
    as.character(length(genes_present)),
    paste(genes_missing, collapse = ";"),
    paste0("Pooled OC centroid in integrated-RPCA dimensions 1-", length(dims_use)),
    "Descriptive only; no donor-level inferential test",
    "One donor per condition; pseudotime and centroid distance do not establish lineage direction, causality or a universal disease mechanism"
  )
)
write.csv(manifest, file.path(src_dir, "analysis_manifest.csv"), row.names = FALSE)

figure_stems <- c(
  "A_joint_myeloid_RPCA_UMAP_cell_state_no_trajectory",
  "B_Monocle2_DDRTree_pseudotime_with_cell_state_labels",
  "C_myeloid_cell_state_counts_CB_CTRL",
  "C2_Monocle2_pseudotime_stage_composition",
  "D_cell_state_by_TREM2_OC_program_gene_heatmap",
  "E_Monocle2_pseudotime_TREM2_OC_program_gradient_heatmap",
  "F_TREM2_OC_program_score_barplots_by_cell_state",
  "G_key_gene_expression_curves_over_Monocle2_pseudotime",
  "H_Myeloid_distance_to_OC_centroid_violin"
)
qa_files <- tidyr::crossing(stem = figure_stems, extension = c("svg", "pdf", "tiff", "png")) |>
  mutate(
    path = file.path(out_dir, paste0(stem, ".", extension)),
    exists = file.exists(path),
    bytes = ifelse(exists, file.info(path)$size, NA_real_)
  )
write.csv(qa_files, file.path(src_dir, "figure_export_QA_manifest.csv"), row.names = FALSE)
if (!all(qa_files$exists) || any(qa_files$bytes <= 0, na.rm = TRUE)) {
  stop("At least one requested publication export is missing or empty.")
}

writeLines(capture.output(sessionInfo()), file.path(src_dir, "sessionInfo.txt"))
message(
  "Joint myeloid TREM2-OC Monocle2 analysis complete. Root state=", root_state,
  "; cells=", ncol(cds), "; output=", out_dir
)
