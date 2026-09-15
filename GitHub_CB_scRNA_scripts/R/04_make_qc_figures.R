log_step("Stage 4: making QC figures with the R backend")

condition_colors <- read_condition_colors()
stage_levels <- c("Raw", "Primary QC", "Singlet")

counts <- readr::read_csv(
  file.path(root, "results", "qc", "qc_counts_by_stage.csv"), show_col_types = FALSE
)
counts$stage <- factor(counts$stage, levels = stage_levels)
counts$sample_id <- factor(counts$sample_id, levels = c("CTRL", "CB"))

raw <- readRDS(file.path(root, "objects", "00_raw_qc_audited.rds"))
with_doublets <- readRDS(file.path(root, "objects", "01_qc_with_doublets.rds"))
corrected <- readRDS(file.path(root, "objects", "03_decontx_corrected.rds"))

p_a <- ggplot2::ggplot(
  counts,
  ggplot2::aes(x = stage, y = n_cells, group = sample_id, colour = sample_id)
) +
  ggplot2::geom_line(linewidth = 0.65) +
  ggplot2::geom_point(size = 2) +
  ggplot2::geom_text(
    ggplot2::aes(label = scales::comma(n_cells)), vjust = -0.7,
    size = 2.2, show.legend = FALSE
  ) +
  ggplot2::scale_colour_manual(values = condition_colors, name = NULL) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.14))) +
  ggplot2::labs(x = NULL, y = "Cells retained", colour = NULL, title = "Cell retention") +
  theme_cb() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1))

thr <- readr::read_csv(
  file.path(root, "results", "qc", "qc_thresholds.csv"), show_col_types = FALSE
)
thr$sample_id <- factor(thr$sample_id, levels = c("CTRL", "CB"))
scatter_df <- raw@meta.data
scatter_df$sample_id <- factor(scatter_df$sample_id, levels = c("CTRL", "CB"))
p_b <- ggplot2::ggplot(
  scatter_df,
  ggplot2::aes(x = nCount_raw, y = nFeature_raw, colour = qc_primary_pass)
) +
  ggplot2::geom_point(size = 0.18, alpha = 0.35) +
  ggplot2::geom_vline(
    data = thr, ggplot2::aes(xintercept = max_counts), linewidth = 0.35,
    linetype = 2
  ) +
  ggplot2::geom_hline(
    data = thr, ggplot2::aes(yintercept = max_features), linewidth = 0.35,
    linetype = 2
  ) +
  ggplot2::facet_wrap(~sample_id, scales = "free") +
  ggplot2::scale_colour_manual(
    values = c(`FALSE` = "#C9C9C9", `TRUE` = "#3B3B3B"), guide = "none"
  ) +
  ggplot2::scale_x_log10(labels = scales::label_number()) +
  ggplot2::labs(
    x = "UMI count (log scale)", y = "Detected genes",
    colour = "Primary QC", title = "Sample-aware complexity thresholds"
  ) +
  theme_cb()

mito_before <- data.frame(
  sample_id = raw$sample_id,
  percent_mt = raw$percent_mt_raw,
  stage = "Raw"
)
mito_after <- raw@meta.data[raw$qc_primary_pass, , drop = FALSE]
mito_after <- data.frame(
  sample_id = mito_after$sample_id,
  percent_mt = mito_after$percent_mt_raw,
  stage = "Primary QC"
)
mito_df <- rbind(mito_before, mito_after)
mito_df$stage <- factor(mito_df$stage, levels = c("Raw", "Primary QC"))
p_c <- ggplot2::ggplot(
  mito_df,
  ggplot2::aes(x = stage, y = percent_mt, fill = sample_id)
) +
  ggplot2::geom_violin(scale = "width", trim = TRUE, linewidth = 0.25) +
  ggplot2::geom_boxplot(width = 0.13, outlier.shape = NA, linewidth = 0.3, fill = "white") +
  ggplot2::scale_fill_manual(values = condition_colors, guide = "none") +
  ggplot2::coord_cartesian(ylim = c(0, 25)) +
  ggplot2::labs(x = NULL, y = "Mitochondrial reads (%)", fill = NULL, title = "Mitochondrial fraction") +
  theme_cb() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1))

dbl <- with_doublets@meta.data
dbl$sample_id <- factor(dbl$sample_id, levels = c("CTRL", "CB"))
p_d <- ggplot2::ggplot(
  dbl,
  ggplot2::aes(x = scDblFinder_class, y = scDblFinder_score, fill = sample_id)
) +
  ggplot2::geom_violin(scale = "width", trim = TRUE, linewidth = 0.25) +
  ggplot2::geom_boxplot(width = 0.13, outlier.shape = NA, linewidth = 0.3, fill = "white") +
  ggplot2::scale_fill_manual(values = condition_colors, guide = "none") +
  ggplot2::scale_y_continuous(trans = "log1p") +
  ggplot2::labs(x = NULL, y = "scDblFinder score (log1p)", fill = NULL, title = "Doublet detection") +
  theme_cb()

cont <- corrected@meta.data
cont$sample_id <- factor(cont$sample_id, levels = c("CTRL", "CB"))
p_e <- ggplot2::ggplot(
  cont,
  ggplot2::aes(x = sample_id, y = decontX_contamination, fill = sample_id)
) +
  ggplot2::geom_violin(scale = "width", trim = TRUE, linewidth = 0.25) +
  ggplot2::geom_boxplot(width = 0.13, outlier.shape = NA, linewidth = 0.3, fill = "white") +
  ggplot2::geom_hline(yintercept = 0.25, linetype = 2, linewidth = 0.35) +
  ggplot2::scale_fill_manual(values = condition_colors, guide = "none") +
  ggplot2::labs(
    x = NULL, y = "Estimated contamination fraction", fill = NULL,
    title = "Ambient-RNA estimate"
  ) +
  theme_cb()

qc_figure <- (p_a | p_b) / (p_c | p_d | p_e) +
  patchwork::plot_layout(heights = c(1, 1)) +
  patchwork::plot_annotation(
    tag_levels = "a",
    subtitle = "One donor per group; all cell-level QC summaries are descriptive."
  ) &
  ggplot2::theme(legend.position = "bottom")

save_pub_r(
  qc_figure,
  file.path(root, "figures", "qc", "FigureS_QC_rebuild"),
  width_mm = 183,
  height_mm = 138
)

base <- readRDS(file.path(root, "objects", "04_base_unintegrated.rds"))
rpca <- readRDS(file.path(root, "objects", "05_base_rpca_visualization_only.rds"))

umap_un <- as.data.frame(Seurat::Embeddings(base, "umap.unintegrated"))
umap_un$sample_id <- factor(base$sample_id, levels = c("CTRL", "CB"))
colnames(umap_un)[1:2] <- c("UMAP_1", "UMAP_2")
umap_rpca <- as.data.frame(Seurat::Embeddings(rpca, "umap.rpca"))
umap_rpca$sample_id <- factor(rpca$sample_id, levels = c("CTRL", "CB"))
colnames(umap_rpca)[1:2] <- c("UMAP_1", "UMAP_2")

umap_theme <- theme_cb() +
  ggplot2::theme(
    axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
    axis.line = ggplot2::element_blank()
  )
p_un <- ggplot2::ggplot(
  umap_un, ggplot2::aes(UMAP_1, UMAP_2, colour = sample_id)
) +
  ggplot2::geom_point(size = 0.18, alpha = 0.6) +
  ggplot2::scale_colour_manual(values = condition_colors) +
  ggplot2::labs(
    x = "UMAP 1", y = "UMAP 2", colour = NULL,
    title = "Unintegrated (primary representation)"
  ) + umap_theme
p_rpca <- ggplot2::ggplot(
  umap_rpca, ggplot2::aes(UMAP_1, UMAP_2, colour = sample_id)
) +
  ggplot2::geom_point(size = 0.18, alpha = 0.6) +
  ggplot2::scale_colour_manual(values = condition_colors) +
  ggplot2::labs(
    x = "UMAP 1", y = "UMAP 2", colour = NULL,
    title = "RPCA alignment (visualization only)"
  ) + umap_theme

embedding_figure <- (p_un | p_rpca) +
  patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(
    tag_levels = "a",
    subtitle = paste0(
      "RPCA is for visualization and annotation sensitivity only.\n",
      "Donor and disease group remain inseparable (one donor per group)."
    )
  ) &
  ggplot2::theme(legend.position = "bottom")
save_pub_r(
  embedding_figure,
  file.path(root, "figures", "qc", "FigureS_base_embeddings"),
  width_mm = 183,
  height_mm = 82
)

readr::write_csv(umap_un, file.path(root, "results", "qc", "source_data_umap_unintegrated.csv"))
readr::write_csv(umap_rpca, file.path(root, "results", "qc", "source_data_umap_rpca.csv"))
readr::write_csv(mito_df, file.path(root, "results", "qc", "source_data_mitochondrial_qc.csv"))
readr::write_csv(
  dbl[, c("sample_id", "scDblFinder_score", "scDblFinder_class")],
  file.path(root, "results", "qc", "source_data_doublets.csv")
)
readr::write_csv(
  cont[, c("sample_id", "decontX_contamination", "decontX_high_contamination")],
  file.path(root, "results", "qc", "source_data_decontx.csv")
)

rm(raw, with_doublets, corrected, base, rpca)
gc()
log_step("Stage 4 complete")
