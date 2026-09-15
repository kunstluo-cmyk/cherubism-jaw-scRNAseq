root <- Sys.getenv("CB_REVISION_ROOT")
if (!nzchar(root)) root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
Sys.setenv(CB_REVISION_ROOT = root)

project_lib <- file.path(root, "renv", "library", "R-4.4", "x86_64-w64-mingw32")
.libPaths(unique(c(project_lib, .Library, .Library.site)))
source(file.path(root, "R", "00_setup.R"), chdir = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(readr)
})

# Figure contract
# Core conclusion: TNF is descriptively more frequently detected in CB myeloid
# cells, whereas TLR2, TLR4 and MYD88 do not show a coordinated increase.
# The figure is descriptive because there is one biological donor per condition.
# No cell-level P values or significance symbols are shown.

out_dir <- file.path(root, "figures", "publication_subfigures")
result_dir <- file.path(root, "results", "annotation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

input_file <- file.path(result_dir, "reviewer_target_genes_by_broad_lineage.csv")
if (!file.exists(input_file)) stop("Missing source table: ", input_file)

genes <- c("TNF", "TLR2", "TLR4", "MYD88")
gene_labels <- c(
  TNF = "TNF",
  TLR2 = "TLR2",
  TLR4 = "TLR4",
  MYD88 = "MYD88"
)
broad_order <- c(
  "Stromal/osteogenic", "Mural", "Endothelial",
  "Myeloid", "Lymphoid", "Neural"
)
condition_order <- c("CTRL", "CB")
condition_colors <- read_condition_colors()[condition_order]

axis_df <- readr::read_csv(input_file, show_col_types = FALSE) |>
  filter(
    gene %in% genes,
    broad_lineage %in% broad_order,
    condition %in% condition_order
  ) |>
  mutate(
    gene = factor(gene, levels = genes),
    broad_lineage = factor(broad_lineage, levels = broad_order),
    condition = factor(condition, levels = condition_order),
    mean_log_normalized = as.numeric(mean_log_normalized),
    percent_detected = as.numeric(percent_detected),
    n_cells = as.integer(n_cells)
  )

if (nrow(axis_df) == 0L) stop("No target-gene data were found")
if (anyDuplicated(axis_df[c("broad_lineage", "condition", "gene")])) {
  stop("Target-gene source table contains duplicated lineage-condition-gene rows")
}

# Explicitly preserve combinations with no cells, particularly CB Neural.
axis_df <- axis_df |>
  tidyr::complete(
    broad_lineage = factor(broad_order, levels = broad_order),
    condition = factor(condition_order, levels = condition_order),
    gene = factor(genes, levels = genes)
  ) |>
  mutate(
    n_cells = tidyr::replace_na(n_cells, 0L),
    has_cells = n_cells > 0L
  )

# Two condition blocks, with deliberate separation matching the established
# manuscript dot-plot style.
gene_position <- tibble::tribble(
  ~condition, ~gene, ~x,
  "CTRL", "TNF", 1.0,
  "CTRL", "TLR2", 2.0,
  "CTRL", "TLR4", 3.0,
  "CTRL", "MYD88", 4.0,
  "CB", "TNF", 5.6,
  "CB", "TLR2", 6.6,
  "CB", "TLR4", 7.6,
  "CB", "MYD88", 8.6
) |>
  mutate(
    condition = factor(condition, levels = condition_order),
    gene = factor(gene, levels = genes)
  )

plot_df <- axis_df |>
  left_join(gene_position, by = c("condition", "gene")) |>
  mutate(y = match(as.character(broad_lineage), rev(broad_order)))

readr::write_csv(
  plot_df |>
    transmute(
      broad_lineage = as.character(broad_lineage),
      condition = as.character(condition),
      gene = as.character(gene),
      n_cells,
      mean_log_normalized,
      percent_detected
    ),
  file.path(result_dir, "source_data_TNF_TLR_MYD88_major_lineages.csv")
)

base_theme <- theme_classic(base_size = 11.5, base_family = "Arial") +
  theme(
    text = element_text(family = "Arial", face = "bold", colour = "black"),
    axis.title = element_text(size = 12.5, face = "bold"),
    axis.text = element_text(size = 11.2, face = "bold", colour = "black"),
    axis.line = element_line(linewidth = 0.65, colour = "black"),
    axis.ticks = element_line(linewidth = 0.55, colour = "black"),
    legend.title = element_text(size = 11.2, face = "bold"),
    legend.text = element_text(size = 10.4, face = "bold"),
    panel.grid = element_blank(),
    plot.tag = element_text(size = 15, face = "bold"),
    plot.margin = margin(7, 9, 7, 7)
  )

header_df <- tibble::tibble(
  condition = factor(condition_order, levels = condition_order),
  xmin = c(0.5, 5.1),
  xmax = c(4.5, 9.1),
  xmid = c(2.5, 7.1),
  label = condition_order
)

p_header <- ggplot(header_df) +
  geom_rect(
    aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = 1, fill = condition),
    colour = NA
  ) +
  geom_text(
    aes(x = xmid, y = 0.5, label = label),
    family = "Arial", fontface = "bold", size = 4.5, colour = "white"
  ) +
  scale_fill_manual(values = condition_colors, drop = FALSE) +
  coord_cartesian(xlim = c(0.5, 9.1), ylim = c(0, 1), expand = FALSE) +
  theme_void(base_family = "Arial") +
  theme(
    legend.position = "none",
    plot.margin = margin(0, 9, 0, 7),
    plot.tag = element_text(family = "Arial", size = 15, face = "bold", colour = "black"),
    plot.tag.position = c(0.005, 0.98)
  ) +
  labs(tag = "A")

no_cell_label <- tibble::tibble(x = 7.1, y = match("Neural", rev(broad_order)))

p_dot_core <- ggplot(filter(plot_df, has_cells), aes(x = x, y = y)) +
  geom_point(aes(size = percent_detected, colour = mean_log_normalized)) +
  geom_text(
    data = no_cell_label,
    aes(x = x, y = y, label = "No cells"),
    inherit.aes = FALSE,
    family = "Arial", fontface = "bold", size = 3.5, colour = "#5A5A5A"
  ) +
  scale_x_continuous(
    breaks = gene_position$x,
    labels = gene_labels[as.character(gene_position$gene)],
    limits = c(0.5, 9.1),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = seq_along(broad_order),
    labels = rev(broad_order),
    limits = c(0.55, length(broad_order) + 0.45),
    expand = c(0, 0)
  ) +
  scale_size(
    range = c(2.5, 10.5),
    limits = c(0, 80),
    breaks = c(0, 20, 40, 60, 80),
    name = "Percent expressed"
  ) +
  scale_colour_gradientn(
    colours = c("#FCE4D6", "#F46D43", "#B30000"),
    limits = c(0, max(plot_df$mean_log_normalized, na.rm = TRUE)),
    oob = scales::squish,
    name = "Average expression\n(log-normalized)"
  ) +
  base_theme +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10.8),
    axis.text.y = element_text(size = 11.5, margin = margin(r = 5)),
    legend.position = "right",
    plot.margin = margin(1, 9, 7, 7)
  )

p_dot <- p_header / p_dot_core +
  plot_layout(heights = c(0.105, 1))

myeloid_df <- axis_df |>
  filter(broad_lineage == "Myeloid", has_cells) |>
  mutate(
    gene = factor(gene, levels = rev(genes)),
    condition = factor(condition, levels = condition_order)
  )

readr::write_csv(
  myeloid_df |>
    transmute(
      condition = as.character(condition),
      gene = as.character(gene),
      n_cells,
      mean_log_normalized,
      percent_detected
    ),
  file.path(result_dir, "source_data_TNF_TLR_MYD88_myeloid_focus.csv")
)

connector_df <- myeloid_df |>
  select(gene, condition, mean_log_normalized, percent_detected) |>
  pivot_wider(
    names_from = condition,
    values_from = c(mean_log_normalized, percent_detected)
  )

p_percent <- ggplot(myeloid_df, aes(x = percent_detected, y = gene)) +
  geom_segment(
    data = connector_df,
    aes(
      x = percent_detected_CTRL,
      xend = percent_detected_CB,
      y = gene,
      yend = gene
    ),
    inherit.aes = FALSE,
    linewidth = 1.0,
    colour = "#C9C9C9"
  ) +
  geom_point(aes(fill = condition), shape = 21, size = 5.1, stroke = 0.7, colour = "white") +
  geom_text(
    aes(
      label = sprintf("%.1f", percent_detected), colour = condition,
      vjust = ifelse(condition == "CTRL", -1.15, 1.85)
    ),
    family = "Arial", fontface = "bold", size = 3.5, show.legend = FALSE
  ) +
  scale_fill_manual(values = condition_colors, drop = FALSE) +
  scale_colour_manual(values = condition_colors, drop = FALSE) +
  scale_x_continuous(
    limits = c(0, 60), breaks = seq(0, 60, 20),
    expand = expansion(mult = c(0, 0.03))
  ) +
  base_theme +
  theme(
    legend.position = "top",
    legend.title = element_text(size = 11.5),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 11.8, margin = margin(r = 5)),
    plot.title = element_text(size = 13, hjust = 0.5, face = "bold")
  ) +
  labs(
    x = "Positive myeloid cells (%)",
    fill = "Condition",
    title = "Detection frequency",
    tag = "B"
  ) +
  guides(colour = "none")

p_mean <- ggplot(myeloid_df, aes(x = mean_log_normalized, y = gene)) +
  geom_segment(
    data = connector_df,
    aes(
      x = mean_log_normalized_CTRL,
      xend = mean_log_normalized_CB,
      y = gene,
      yend = gene
    ),
    inherit.aes = FALSE,
    linewidth = 1.0,
    colour = "#C9C9C9"
  ) +
  geom_point(aes(fill = condition), shape = 21, size = 5.1, stroke = 0.7, colour = "white") +
  geom_text(
    aes(
      label = sprintf("%.2f", mean_log_normalized), colour = condition,
      vjust = ifelse(condition == "CTRL", -1.15, 1.85)
    ),
    family = "Arial", fontface = "bold", size = 3.5, show.legend = FALSE
  ) +
  scale_fill_manual(values = condition_colors, drop = FALSE) +
  scale_colour_manual(values = condition_colors, drop = FALSE) +
  scale_x_continuous(
    limits = c(0, 1.36), breaks = c(0, 0.4, 0.8, 1.2),
    expand = expansion(mult = c(0, 0.03))
  ) +
  base_theme +
  theme(
    legend.position = "top",
    legend.title = element_text(size = 11.5),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 11.8, margin = margin(r = 5)),
    plot.title = element_text(size = 13, hjust = 0.5, face = "bold")
  ) +
  labs(
    x = "Mean myeloid expression\n(log-normalized)",
    fill = "Condition",
    title = "Average abundance",
    tag = "C"
  ) +
  guides(colour = "none")

combined <- (p_dot | p_percent | p_mean) +
  plot_layout(widths = c(1.65, 1, 1), guides = "keep") &
  theme(plot.background = element_rect(fill = "white", colour = NA))

save_exact <- function(plot, path, width_px, height_px, dpi = 600) {
  width_in <- width_px / dpi
  height_in <- height_px / dpi

  svglite::svglite(paste0(path, ".svg"), width = width_in, height = height_in)
  print(plot)
  grDevices::dev.off()

  grDevices::cairo_pdf(
    paste0(path, ".pdf"), width = width_in, height = height_in,
    family = "Arial"
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

save_exact(
  p_dot,
  file.path(out_dir, "Bubbleplot_TNF_TLR_MYD88_major_lineages"),
  5000L, 3600L
)
save_exact(
  p_percent,
  file.path(out_dir, "Lollipop_TNF_TLR_MYD88_myeloid_detection"),
  3600L, 3000L
)
save_exact(
  p_mean,
  file.path(out_dir, "Lollipop_TNF_TLR_MYD88_myeloid_mean"),
  3600L, 3000L
)
save_exact(
  combined,
  file.path(out_dir, "FigureS_TNF_TLR_MYD88_expression"),
  8000L, 4200L
)

# Machine-readable QA and the exact values used in the reviewer response.
qa <- tibble::tibble(
  check = c(
    "one_biological_donor_per_condition",
    "inferential_statistics_shown",
    "target_gene_count",
    "major_lineage_count",
    "combined_width_px",
    "combined_height_px",
    "source_table"
  ),
  value = c(
    "TRUE", "FALSE", length(genes), length(broad_order),
    8000L, 4200L, basename(input_file)
  )
)
readr::write_csv(
  qa,
  file.path(result_dir, "validation_TNF_TLR_MYD88_figure.csv")
)

message("TNF/TLR2/TLR4/MYD88 reviewer-response figure exported")
