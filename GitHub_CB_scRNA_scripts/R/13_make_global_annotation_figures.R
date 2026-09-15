root <- Sys.getenv("CB_REVISION_ROOT")
if (!nzchar(root)) root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
Sys.setenv(CB_REVISION_ROOT = root)
project_lib <- file.path(root, "renv", "library", "R-4.4", "x86_64-w64-mingw32")
.libPaths(unique(c(project_lib, .Library, .Library.site)))
source(file.path(root, "R", "00_setup.R"), chdir = FALSE)

if (!requireNamespace("ggrastr", quietly = TRUE)) stop("ggrastr is required for rasterized scatter layers.")
if (!requireNamespace("ggrepel", quietly = TRUE)) stop("ggrepel is required for UMAP labels.")

log_step("Figure generation: revised global annotation and stromal-state evidence")
fig_dir <- file.path(root, "figures", "annotation")
src_dir <- file.path(root, "results", "annotation")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

rpca <- readRDS(file.path(root, "objects", "07_rpca_broad_lineages_reviewed.rds"))
stromal <- readRDS(file.path(root, "objects", "08_cb_lesion_act_ob_reclustered.rds"))

condition_colors <- read_condition_colors()
pal <- readr::read_csv(file.path(root, "config", "palette.csv"), show_col_types = FALSE)
fine_colors <- stats::setNames(pal$color[pal$category == "celltype"], pal$label[pal$category == "celltype"])
broad_colors <- c(
  "Stromal/osteogenic" = "#8DD3C7",
  "Mural" = "#7393B3",
  "Myeloid" = "#33A02C",
  "Lymphoid" = "#CCEBC5",
  "Endothelial" = "#43AA8B",
  "Neural" = "#FB9A99"
)

display_name <- c(
  "Act OB-Stromal" = "ECM-remodeling stromal state",
  "TREM2+ Mac" = "TREM2+ Mac",
  "LEPR+ MSC-like" = "LEPR+ MSC-like",
  "CNCC-OCPs" = "CNCC-OCPs",
  "ChondroFibroStromal" = "ChondroFibroStromal"
)
rename_display <- function(x) {
  x <- as.character(x)
  hit <- x %in% names(display_name)
  x[hit] <- unname(display_name[x[hit]])
  x
}
fine_colors_display <- fine_colors
names(fine_colors_display) <- rename_display(names(fine_colors_display))

fine_order <- c(
  "FibroStromal", "LesionStromal", "Act OB-Stromal", "ChondroFibroStromal",
  "LEPR+ MSC-like", "Mature OB", "CNCC-OCPs", "Pericytes", "VSMC",
  "Endothelial", "Resident Mac", "Inflam Mono", "TREM2+ Mac", "OC",
  "Neutrophils", "Mast cells", "pDC", "T/NK", "B cells", "Plasma cells", "Schwann"
)
fine_order_display <- rename_display(fine_order)
broad_order <- c("Stromal/osteogenic", "Mural", "Endothelial", "Myeloid", "Lymphoid", "Neural")

umap <- as.data.frame(Seurat::Embeddings(rpca, reduction = "umap.rpca")) |>
  tibble::rownames_to_column("cell_id")
colnames(umap)[2:3] <- c("UMAP_1", "UMAP_2")
umap$condition <- factor(as.character(rpca$group), levels = c("CTRL", "CB"))
umap$broad_lineage <- factor(as.character(rpca$broad_lineage_reviewed), levels = broad_order)
umap$celltype <- factor(rename_display(rpca$celltype_reviewed_provisional), levels = fine_order_display)
readr::write_csv(umap, file.path(src_dir, "source_data_global_annotation_umap.csv.gz"))

label_position <- function(df, group_col) {
  df |>
    dplyr::filter(!is.na(.data[[group_col]])) |>
    dplyr::group_by(.data[[group_col]]) |>
    dplyr::summarise(UMAP_1 = stats::median(UMAP_1), UMAP_2 = stats::median(UMAP_2), .groups = "drop")
}

theme_umap <- theme_cb(7) +
  ggplot2::theme(
    axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
    axis.line = ggplot2::element_line(linewidth = 0.35),
    legend.key.height = grid::unit(2.8, "mm"), legend.key.width = grid::unit(2.8, "mm")
  )

p_broad <- ggplot2::ggplot(umap, ggplot2::aes(UMAP_1, UMAP_2, colour = broad_lineage)) +
  ggrastr::geom_point_rast(size = 0.22, alpha = 0.85, raster.dpi = 500) +
  ggrepel::geom_text_repel(
    data = label_position(umap, "broad_lineage"),
    ggplot2::aes(x = UMAP_1, y = UMAP_2, label = broad_lineage), inherit.aes = FALSE,
    size = 2.4, fontface = "bold", colour = "black", segment.color = NA,
    box.padding = 0.2, max.overlaps = Inf
  ) +
  ggplot2::scale_colour_manual(values = broad_colors, drop = FALSE) +
  ggplot2::labs(x = "UMAP1", y = "UMAP2", colour = NULL, title = "Reviewed major lineages") +
  theme_umap + ggplot2::theme(legend.position = "none")

p_fine <- ggplot2::ggplot(umap, ggplot2::aes(UMAP_1, UMAP_2, colour = celltype)) +
  ggrastr::geom_point_rast(size = 0.17, alpha = 0.82, raster.dpi = 500) +
  ggplot2::scale_colour_manual(values = fine_colors_display, drop = FALSE, na.value = "grey80") +
  ggplot2::labs(x = "UMAP1", y = "UMAP2", colour = NULL, title = "Transcriptome-reviewed cell states") +
  theme_umap +
  ggplot2::theme(legend.position = "right", legend.text = ggplot2::element_text(size = 5.2)) +
  ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 1.8, alpha = 1), ncol = 1))

p_condition <- ggplot2::ggplot(umap, ggplot2::aes(UMAP_1, UMAP_2)) +
  ggrastr::geom_point_rast(colour = "grey90", size = 0.12, raster.dpi = 500) +
  ggrastr::geom_point_rast(
    ggplot2::aes(colour = condition), size = 0.18, alpha = 0.85, raster.dpi = 500
  ) +
  ggplot2::facet_wrap(~condition, nrow = 1) +
  ggplot2::scale_colour_manual(values = condition_colors, drop = FALSE) +
  ggplot2::labs(x = "UMAP1", y = "UMAP2", colour = NULL, title = "Condition-resolved UMAP") +
  theme_umap +
  ggplot2::theme(legend.position = "none", strip.background = ggplot2::element_blank())

composition <- umap |>
  dplyr::filter(!is.na(celltype)) |>
  dplyr::mutate(
    display_group = ifelse(
      broad_lineage %in% c("Stromal/osteogenic", "Mural"),
      "Mesenchymal / mural", "Hematopoietic / endothelial / neural"
    )
  ) |>
  dplyr::count(condition, display_group, celltype, name = "n_cells") |>
  dplyr::group_by(condition) |>
  dplyr::mutate(sample_total = sum(n_cells), percent_of_all_cells = 100 * n_cells / sample_total) |>
  dplyr::ungroup()
readr::write_csv(composition, file.path(src_dir, "source_data_composition_barplots.csv"))

p_composition <- ggplot2::ggplot(
  composition,
  ggplot2::aes(x = percent_of_all_cells, y = celltype, fill = condition)
) +
  ggplot2::geom_col(position = ggplot2::position_dodge2(width = 0.78, preserve = "single"), width = 0.70) +
  ggplot2::facet_grid(display_group ~ ., scales = "free_y", space = "free_y") +
  ggplot2::scale_fill_manual(values = condition_colors) +
  ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.04))) +
  ggplot2::labs(x = "Cells (% of all cells in each sample)", y = NULL, fill = NULL,
                title = "Separated cell-state composition") +
  theme_cb(7) +
  ggplot2::theme(
    legend.position = "top", strip.background = ggplot2::element_rect(fill = "grey95", colour = NA),
    strip.text.y.right = ggplot2::element_text(angle = 90, size = 5.6),
    axis.text.y = ggplot2::element_text(size = 5.6)
  )

summarise_dot <- function(object, grouping, genes) {
  mat <- SeuratObject::LayerData(object, assay = "RNA", layer = "data")
  counts <- SeuratObject::LayerData(object, assay = "RNA", layer = "counts")
  genes <- intersect(genes, rownames(mat))
  grouping <- as.character(grouping)
  groups <- unique(grouping)
  out <- lapply(groups, function(g) {
    idx <- which(grouping == g)
    data.frame(
      group = g,
      gene = genes,
      average_expression = Matrix::rowMeans(mat[genes, idx, drop = FALSE]),
      percent_expressed = 100 * Matrix::rowMeans(counts[genes, idx, drop = FALSE] > 0),
      stringsAsFactors = FALSE
    )
  }) |>
    dplyr::bind_rows() |>
    dplyr::group_by(gene) |>
    dplyr::mutate(
      scaled_average = {
        s <- stats::sd(average_expression)
        if (is.finite(s) && s > 0) as.numeric(scale(average_expression)) else rep(0, dplyr::n())
      },
      scaled_average = pmax(-1, pmin(2, scaled_average))
    ) |>
    dplyr::ungroup()
  out
}

broad_marker_groups <- list(
  `Stromal/osteogenic` = c("COL1A1", "DCN", "PDGFRA", "RUNX2"),
  Mural = c("RGS5", "PDGFRB", "ACTA2", "MYH11"),
  Endothelial = c("PECAM1", "EMCN", "VWF", "CLDN5"),
  Myeloid = c("LST1", "TYROBP", "FCER1G", "LYZ"),
  Lymphoid = c("PTPRC", "CD3D", "NKG7", "CD79A"),
  Neural = c("S100B", "PLP1", "MPZ", "CDH19")
)
broad_genes <- unlist(broad_marker_groups, use.names = FALSE)
broad_dot <- summarise_dot(rpca, rpca$broad_lineage_reviewed, broad_genes)
broad_dot$group <- factor(broad_dot$group, levels = rev(broad_order))
broad_dot$gene <- factor(broad_dot$gene, levels = broad_genes)
readr::write_csv(broad_dot, file.path(src_dir, "source_data_broad_lineage_dotplot.csv"))

p_dot <- ggplot2::ggplot(broad_dot, ggplot2::aes(gene, group)) +
  ggplot2::geom_point(ggplot2::aes(size = percent_expressed, colour = scaled_average)) +
  ggplot2::scale_colour_gradient(low = "#FDE0DD", high = "#D7301F", limits = c(-1, 2), oob = scales::squish) +
  ggplot2::scale_size(range = c(0.2, 4.3), limits = c(0, 100)) +
  ggplot2::labs(x = NULL, y = NULL, colour = "Scaled average\nexpression", size = "% expressed",
                title = "Canonical lineage markers") +
  theme_cb(7) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 60, hjust = 1, vjust = 1, size = 5.8),
    axis.text.y = ggplot2::element_text(size = 6),
    legend.position = "right"
  )

figure1 <- ((p_broad | p_fine) / (p_condition | p_composition) / p_dot) +
  patchwork::plot_layout(heights = c(1.0, 1.12, 0.72), guides = "keep") +
  patchwork::plot_annotation(tag_levels = "A")
save_pub_r(figure1, file.path(fig_dir, "Figure1_revised_global_annotation"), 183, 235, dpi = 600)

# Supplementary 21-state marker audit.
fine_marker_groups <- list(
  FibroStromal = c("COL11A1", "ASPN"),
  LesionStromal = c("PDZRN4", "BICC1", "GPC6", "RUNX2"),
  `ECM-remodeling stromal state` = c("MYL9", "IGFBP7", "LUM", "CTHRC1"),
  ChondroFibroStromal = c("ITGBL1", "THBS2"),
  `LEPR+ MSC-like` = c("LEPR", "CXCL12"),
  `Mature OB` = c("BGLAP", "IBSP"),
  `CNCC-OCPs` = c("SATB2", "CXCL14"),
  Pericytes = c("RGS5", "PDGFRB"),
  VSMC = c("ACTA2", "MYH11"),
  Endothelial = c("PECAM1", "EMCN"),
  `Resident Mac` = c("C1QA", "MRC1"),
  `Inflam Mono` = c("IL1B", "S100A9"),
  `TREM2+ Mac` = c("TREM2", "HMOX1"),
  OC = c("ACP5", "MMP9"),
  Neutrophils = c("FCGR3B", "CSF3R"),
  `Mast cells` = c("TPSAB1", "CPA3"),
  pDC = c("GZMB", "LILRA4"),
  `T/NK` = c("CD3D", "NKG7"),
  `B cells` = c("MS4A1", "CD79A"),
  `Plasma cells` = c("JCHAIN", "MZB1"),
  Schwann = c("S100B", "MPZ")
)
fine_genes <- unique(unlist(fine_marker_groups, use.names = FALSE))
fine_dot <- summarise_dot(rpca, rename_display(rpca$celltype_reviewed_provisional), fine_genes)
fine_dot$group <- factor(fine_dot$group, levels = rev(fine_order_display))
fine_dot$gene <- factor(fine_dot$gene, levels = fine_genes)
readr::write_csv(fine_dot, file.path(src_dir, "source_data_fine_annotation_dotplot.csv"))

p_fine_dot <- ggplot2::ggplot(fine_dot, ggplot2::aes(gene, group)) +
  ggplot2::geom_point(ggplot2::aes(size = percent_expressed, colour = scaled_average)) +
  ggplot2::scale_colour_gradient(low = "#FDE0DD", high = "#D7301F", limits = c(-1, 2), oob = scales::squish) +
  ggplot2::scale_size(range = c(0.15, 4), limits = c(0, 100)) +
  ggplot2::labs(x = NULL, y = NULL, colour = "Scaled average\nexpression", size = "% expressed") +
  theme_cb(7) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 60, hjust = 1, size = 5.4),
    axis.text.y = ggplot2::element_text(size = 5.5)
  )
save_pub_r(p_fine_dot, file.path(fig_dir, "FigureS_global_fine_marker_audit"), 183, 135, dpi = 600)

# Focused stromal-state figure.
stromal_umap <- as.data.frame(Seurat::Embeddings(stromal, reduction = "umap.stromal")) |>
  tibble::rownames_to_column("cell_id")
colnames(stromal_umap)[2:3] <- c("UMAP_1", "UMAP_2")
stromal_umap$state_raw <- as.character(stromal$celltype_reviewed_provisional)
stromal_umap$state <- factor(
  rename_display(stromal_umap$state_raw),
  levels = c("LesionStromal", "ECM-remodeling stromal state")
)
stromal_umap$cluster <- factor(as.character(stromal$seurat_clusters))

p_state <- ggplot2::ggplot(stromal_umap, ggplot2::aes(UMAP_1, UMAP_2, colour = state)) +
  ggrastr::geom_point_rast(size = 0.28, alpha = 0.85, raster.dpi = 500) +
  ggplot2::scale_colour_manual(values = fine_colors_display[c("LesionStromal", "ECM-remodeling stromal state")]) +
  ggplot2::labs(x = "UMAP1", y = "UMAP2", colour = NULL, title = "CB stromal states") +
  theme_umap + ggplot2::theme(legend.position = "bottom", legend.text = ggplot2::element_text(size = 5.6))

cluster_colors <- stats::setNames(c("#8DD3C7", "#BEBADA", "#80B1D3", "#B3DE69", "#FDB462", "#CAB2D6"), 0:5)
p_cluster <- ggplot2::ggplot(stromal_umap, ggplot2::aes(UMAP_1, UMAP_2, colour = cluster)) +
  ggrastr::geom_point_rast(size = 0.28, alpha = 0.85, raster.dpi = 500) +
  ggplot2::scale_colour_manual(values = cluster_colors) +
  ggplot2::labs(x = "UMAP1", y = "UMAP2", colour = "Recluster", title = "Independent reclustering") +
  theme_umap + ggplot2::theme(legend.position = "bottom")

stromal_markers <- c(
  "PDZRN4", "BICC1", "GPC6", "RUNX2", "FAP", "WNT5A", "VCAN", "SOX5",
  "MYL9", "IGFBP7", "CTHRC1", "LUM", "DCN", "POSTN", "COL1A1", "COL3A1",
  "ALPL", "BGLAP", "TNFSF11", "IBSP", "SP7"
)
stromal_dot <- summarise_dot(stromal, rename_display(stromal$celltype_reviewed_provisional), stromal_markers)
stromal_dot$group <- factor(stromal_dot$group, levels = rev(c("LesionStromal", "ECM-remodeling stromal state")))
stromal_dot$gene <- factor(stromal_dot$gene, levels = stromal_markers)
readr::write_csv(stromal_dot, file.path(src_dir, "source_data_stromal_marker_dotplot.csv"))

p_stromal_dot <- ggplot2::ggplot(stromal_dot, ggplot2::aes(gene, group)) +
  ggplot2::geom_point(ggplot2::aes(size = percent_expressed, colour = average_expression)) +
  ggplot2::scale_colour_gradient(low = "#FDE0DD", high = "#D7301F") +
  ggplot2::scale_size(range = c(0.4, 5), limits = c(0, 100)) +
  ggplot2::labs(x = NULL, y = NULL, colour = "Mean log-normalized\nexpression", size = "% expressed",
                title = "Markers that separate the related stromal states") +
  theme_cb(7) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 55, hjust = 1, size = 5.8))

effect <- readr::read_csv(file.path(src_dir, "act_ob_vs_lesion_marker_evidence.csv"), show_col_types = FALSE) |>
  dplyr::filter(gene %in% stromal_markers) |>
  dplyr::mutate(
    direction = dplyr::case_when(
      q025 > 0 ~ "ECM-remodeling stromal state",
      q975 < 0 ~ "LesionStromal",
      TRUE ~ "Non-discriminating"
    ),
    gene = factor(gene, levels = rev(stromal_markers))
  )
readr::write_csv(effect, file.path(src_dir, "source_data_stromal_effect_plot.csv"))
effect_colors <- c(
  "ECM-remodeling stromal state" = fine_colors_display[["ECM-remodeling stromal state"]],
  "LesionStromal" = fine_colors_display[["LesionStromal"]],
  "Non-discriminating" = "grey60"
)
p_effect <- ggplot2::ggplot(effect, ggplot2::aes(median_delta_act_minus_lesion, gene, colour = direction)) +
  ggplot2::geom_vline(xintercept = 0, linetype = 2, linewidth = 0.35, colour = "grey45") +
  ggplot2::geom_errorbar(
    ggplot2::aes(xmin = q025, xmax = q975),
    orientation = "y", width = 0, linewidth = 0.45
  ) +
  ggplot2::geom_point(size = 1.7) +
  ggplot2::scale_colour_manual(values = effect_colors) +
  ggplot2::labs(
    x = "Balanced mean-expression difference\n(positive = ECM-remodeling stromal state)",
    y = NULL, colour = NULL, title = "Cell-balanced robustness (200 iterations)"
  ) +
  theme_cb(7) +
  ggplot2::theme(legend.position = "bottom", legend.text = ggplot2::element_text(size = 5.5))

cluster_comp <- stromal_umap |>
  dplyr::count(cluster, state, name = "n_cells") |>
  dplyr::group_by(cluster) |>
  dplyr::mutate(percent = 100 * n_cells / sum(n_cells)) |>
  dplyr::ungroup()
readr::write_csv(cluster_comp, file.path(src_dir, "source_data_stromal_recluster_composition.csv"))
p_cluster_comp <- ggplot2::ggplot(cluster_comp, ggplot2::aes(cluster, percent, fill = state)) +
  ggplot2::geom_col(width = 0.72, colour = "white", linewidth = 0.15) +
  ggplot2::scale_fill_manual(values = fine_colors_display[c("LesionStromal", "ECM-remodeling stromal state")]) +
  ggplot2::scale_y_continuous(labels = function(x) paste0(x, "%"), expand = ggplot2::expansion(mult = c(0, 0.03))) +
  ggplot2::labs(x = "Independent stromal cluster", y = "State composition", fill = NULL,
                title = "State distribution across reclusters") +
  theme_cb(7) + ggplot2::theme(legend.position = "bottom", legend.text = ggplot2::element_text(size = 5.5))

figure_stromal <- ((p_state | p_cluster) / p_stromal_dot / (p_effect | p_cluster_comp)) +
  patchwork::plot_layout(heights = c(1, 0.70, 1.1)) +
  patchwork::plot_annotation(tag_levels = "A")
save_pub_r(figure_stromal, file.path(fig_dir, "FigureS_stromal_state_revalidation"), 183, 215, dpi = 600)

log_step("Revised global annotation and stromal-state figures complete")
