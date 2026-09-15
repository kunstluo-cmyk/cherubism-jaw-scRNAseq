root <- Sys.getenv("CB_REVISION_ROOT")
if (!nzchar(root)) root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
config <- yaml::read_yaml(file.path(root, "config", "analysis_config.yml"))
seed <- as.integer(config$seed)
set.seed(seed)

suppressPackageStartupMessages({
  library(Seurat)
  library(monocle)
  library(Matrix)
  library(UCell)
  library(BiocParallel)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggrepel)
  library(patchwork)
})

out_dir <- file.path(
  root, "figures", "publication_subfigures", "myeloid_TREM2_OC_axis_revision"
)
src_dir <- file.path(
  root, "results", "trajectory",
  "trem2_mac_oc_differentiation_propensity"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(src_dir, recursive = TRUE, showWarnings = FALSE)

# Evidence contract ---------------------------------------------------------
# The central question is whether TREM2+ macrophages occupy a more
# osteoclast-primed transcriptional state than the other non-OC myeloid
# states. TREM2 is deliberately excluded from every score below so that the
# conclusion cannot be produced simply by reusing the annotation marker.
# Because the study contains one participant per condition, all displays are
# descriptive and no donor-level inferential P values are reported.

state_order <- c("Inflam Mono", "Resident Mac", "TREM2+ Mac", "OC")
state_colors <- c(
  "Inflam Mono" = "#E68472",
  "Resident Mac" = "#79B8A9",
  "TREM2+ Mac" = "#9A78B4",
  "OC" = "#527BA8"
)
condition_colors <- c("CTRL" = "#4DBBD5", "CB" = "#E64B35")
score_colors <- c(
  "GO multinuclear OC differentiation" = "#A883B7",
  "OC commitment / fusion" = "#DDA85C",
  "OC resorption" = "#527BA8"
)

# GO set fixed to the MSigDB v2025.1.Hs GO:BP definition available when the
# revision was run. BBLN is retained in the definition table even if it is not
# detected in this dataset. The two curated sets represent distinct stages of
# OC differentiation and function; neither contains TREM2.
signatures_requested <- list(
  "GO multinuclear OC differentiation" = c(
    "BBLN", "CD109", "CD81", "DCSTAMP", "OCSTAMP", "SBNO2",
    "SH3PXD2A", "TCTA", "TNFRSF11A"
  ),
  "OC commitment / fusion" = c(
    "TNFRSF11A", "NFATC1", "DCSTAMP", "OCSTAMP", "ATP6V0D2",
    "OSCAR", "SIGLEC15", "TMEM64"
  ),
  "OC resorption" = c(
    "ACP5", "CTSK", "MMP9", "CA2", "TCIRG1", "CLCN7",
    "ATP6V1C1", "OSTM1", "SLC4A2", "SNX10"
  )
)
signature_order <- names(signatures_requested)
signature_labeller <- ggplot2::as_labeller(c(
  "GO multinuclear OC differentiation" = "GO multinuclear OC\ndifferentiation",
  "OC commitment / fusion" = "OC commitment / fusion",
  "OC resorption" = "OC resorption"
))

marker_programs <- list(
  "Commitment" = c("TNFRSF11A", "NFATC1", "OSCAR", "SIGLEC15"),
  "Fusion" = c("DCSTAMP", "OCSTAMP", "ATP6V0D2", "TMEM64"),
  "Resorption" = c("ACP5", "CTSK", "MMP9", "CA2", "TCIRG1", "CLCN7")
)
marker_order <- unlist(marker_programs, use.names = FALSE)
marker_program_colors <- c(
  "Commitment" = "#DDA85C",
  "Fusion" = "#9A78B4",
  "Resorption" = "#527BA8"
)

theme_pub <- theme_classic(base_size = 11.5, base_family = "Arial") +
  theme(
    text = element_text(family = "Arial", face = "bold", colour = "black"),
    axis.title = element_text(size = 11.5, face = "bold"),
    axis.text = element_text(size = 9.5, face = "bold", colour = "black"),
    axis.line = element_line(linewidth = 0.5, colour = "black"),
    axis.ticks = element_line(linewidth = 0.45, colour = "black"),
    strip.background = element_blank(),
    strip.text = element_text(size = 10.5, face = "bold", colour = "black"),
    legend.title = element_text(size = 10.5, face = "bold"),
    legend.text = element_text(size = 9.2, face = "bold"),
    plot.title = element_text(size = 13, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 9.1, face = "bold", colour = "#3A3A3A"),
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

median_iqr <- function(x) {
  tibble(
    median = median(x, na.rm = TRUE),
    q25 = quantile(x, 0.25, na.rm = TRUE, names = FALSE),
    q75 = quantile(x, 0.75, na.rm = TRUE, names = FALSE),
    mean = mean(x, na.rm = TRUE),
    pct_positive = 100 * mean(x > 0, na.rm = TRUE),
    n_cells = sum(is.finite(x))
  )
}

message("Loading reviewed corrected expression object...")
obj <- readRDS(file.path(root, "objects", "06_base_broad_lineages_reviewed.rds"))
axis_cells <- colnames(obj)[as.character(obj$celltype_reviewed_provisional) %in% state_order]
obj <- subset(obj, cells = axis_cells)
obj$cell_state <- factor(
  as.character(obj$celltype_reviewed_provisional), levels = state_order
)
obj$condition <- factor(as.character(obj$group), levels = c("CTRL", "CB"))
DefaultAssay(obj) <- "RNA"

counts <- SeuratObject::LayerData(obj, assay = "RNA", layer = "counts")
counts <- counts[, colnames(obj), drop = FALSE]

signature_definition <- tibble(
  signature = rep(signature_order, lengths(signatures_requested)),
  gene = unlist(signatures_requested, use.names = FALSE)
) |>
  mutate(
    contains_TREM2 = gene == "TREM2",
    present_in_dataset = gene %in% rownames(counts)
  )
if (any(signature_definition$contains_TREM2)) {
  stop("TREM2 must be excluded from every differentiation-propensity signature.")
}
write.csv(
  signature_definition,
  file.path(src_dir, "signature_gene_definitions_and_availability.csv"),
  row.names = FALSE
)

signatures_present <- lapply(
  signatures_requested, function(x) intersect(x, rownames(counts))
)
if (any(lengths(signatures_present) < 3)) {
  stop("At least one OC propensity signature has fewer than three detected genes.")
}

message("Calculating TREM2-independent UCell scores...")
score_matrix <- UCell::ScoreSignatures_UCell(
  matrix = counts,
  features = signatures_present,
  assay = "counts",
  BPPARAM = BiocParallel::SerialParam(),
  ncores = 1
)
score_matrix <- as.data.frame(score_matrix) |>
  rownames_to_column("cell_id")
score_cols <- setdiff(colnames(score_matrix), "cell_id")
if (length(score_cols) != length(signature_order)) {
  stop("Unexpected number of UCell score columns.")
}
names(score_cols) <- signature_order

cell_meta <- obj@meta.data |>
  as.data.frame() |>
  rownames_to_column("cell_id") |>
  transmute(
    cell_id,
    condition = factor(as.character(condition), levels = c("CTRL", "CB")),
    cell_state = factor(as.character(cell_state), levels = state_order)
  )

score_long <- cell_meta |>
  left_join(score_matrix, by = "cell_id") |>
  pivot_longer(
    cols = all_of(unname(score_cols)),
    names_to = "score_column", values_to = "score"
  ) |>
  mutate(
    signature = names(score_cols)[match(score_column, unname(score_cols))],
    signature = factor(signature, levels = signature_order)
  )
write.csv(
  score_long,
  file.path(src_dir, "source_data_A_B_UCell_scores_per_cell.csv"),
  row.names = FALSE
)

score_summary <- score_long |>
  group_by(condition, cell_state, signature) |>
  group_modify(~median_iqr(.x$score)) |>
  ungroup()
write.csv(
  score_summary,
  file.path(src_dir, "source_data_B_condition_state_UCell_score_summary.csv"),
  row.names = FALSE
)

pooled_summary <- score_long |>
  group_by(cell_state, signature) |>
  group_modify(~median_iqr(.x$score)) |>
  ungroup()
write.csv(
  pooled_summary,
  file.path(src_dir, "source_data_A_pooled_state_UCell_score_summary.csv"),
  row.names = FALSE
)

# A. Three independent TREM2-free score distributions ----------------------
state_n <- cell_meta |>
  count(cell_state, name = "n_cells") |>
  complete(cell_state = factor(state_order, levels = state_order), fill = list(n_cells = 0))
state_axis_labels <- setNames(
  paste0(as.character(state_n$cell_state), "\n(n=", state_n$n_cells, ")"),
  as.character(state_n$cell_state)
)

p_score_violin <- ggplot(
  score_long,
  aes(cell_state, score, fill = cell_state, colour = cell_state)
) +
  geom_violin(
    scale = "width", trim = TRUE, alpha = 0.78,
    linewidth = 0.42, width = 0.84
  ) +
  geom_boxplot(
    width = 0.13, outlier.shape = NA, fill = "white", colour = "black",
    linewidth = 0.42, alpha = 0.92
  ) +
  facet_wrap(
    ~signature, nrow = 1, scales = "free_y", labeller = signature_labeller
  ) +
  scale_fill_manual(values = state_colors, drop = FALSE) +
  scale_colour_manual(values = state_colors, drop = FALSE) +
  scale_x_discrete(labels = state_axis_labels, drop = FALSE) +
  theme_pub +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 26, hjust = 1, vjust = 1),
    panel.spacing.x = unit(11, "pt")
  ) +
  labs(
    x = NULL, y = "UCell score",
    title = "TREM2+ macrophages show an osteoclast-primed transcript pattern",
    subtitle = paste0(
      "TREM2 is excluded from all three signatures; pooled cells are shown descriptively ",
      "and no donor-level test is performed"
    )
  )
save_plot(
  p_score_violin,
  "A_TREM2_independent_OC_UCell_scores_violin",
  220, 118
)

# B. Condition-aware score summaries ---------------------------------------
score_summary$cell_state <- factor(score_summary$cell_state, levels = state_order)
score_summary$signature <- factor(score_summary$signature, levels = signature_order)
dodge <- position_dodge(width = 0.58)

p_condition <- ggplot(
  score_summary,
  aes(cell_state, median, colour = condition, group = condition)
) +
  geom_linerange(
    aes(ymin = q25, ymax = q75), position = dodge,
    linewidth = 0.85, alpha = 0.92
  ) +
  geom_point(
    aes(size = n_cells), position = dodge, stroke = 0.5, alpha = 0.96
  ) +
  facet_wrap(
    ~signature, nrow = 1, scales = "free_y", labeller = signature_labeller
  ) +
  scale_colour_manual(values = condition_colors, drop = FALSE) +
  scale_size_area(
    name = "Cells", max_size = 6.4,
    breaks = c(3, 25, 100, 250), limits = c(1, max(score_summary$n_cells))
  ) +
  scale_x_discrete(drop = FALSE) +
  theme_pub +
  theme(
    axis.text.x = element_text(angle = 28, hjust = 1, vjust = 1),
    panel.spacing.x = unit(11, "pt"),
    legend.position = "right"
  ) +
  labs(
    x = NULL, y = "Median UCell score (IQR)", colour = "Condition",
    title = "Condition-aware osteoclast-propensity score summaries",
    subtitle = paste0(
      "Point size reports the number of cells; CTRL TREM2+ Mac contains only three cells, ",
      "so all contrasts remain descriptive"
    )
  )
save_plot(
  p_condition,
  "B_condition_aware_OC_UCell_score_summary",
  224, 123
)

# C. Marker-level commitment, fusion and resorption evidence ---------------
marker_present <- marker_order[marker_order %in% rownames(counts)]
marker_missing <- setdiff(marker_order, marker_present)
write.csv(
  data.frame(gene = marker_missing),
  file.path(src_dir, "requested_marker_genes_missing.csv"),
  row.names = FALSE
)
if (length(marker_present) < 8) stop("Too few marker genes are present for the dot plot.")

lib_size <- Matrix::colSums(counts)
norm_mat <- Matrix::t(Matrix::t(counts[marker_present, , drop = FALSE]) /
  pmax(lib_size, 1)) * 1e4
norm_mat@x <- log1p(norm_mat@x)

state_cells <- split(colnames(obj), as.character(obj$cell_state))[state_order]
dot_df <- bind_rows(lapply(state_order, function(st) {
  cells <- state_cells[[st]]
  x <- norm_mat[, cells, drop = FALSE]
  tibble(
    cell_state = st,
    gene = marker_present,
    average_expression = Matrix::rowMeans(x),
    percent_expressed = 100 * Matrix::rowMeans(counts[marker_present, cells, drop = FALSE] > 0)
  )
})) |>
  group_by(gene) |>
  mutate(
    scaled_average_expression = if (sd(average_expression) > 0) {
      as.numeric(scale(average_expression))
    } else {
      0
    },
    scaled_average_expression = pmax(pmin(scaled_average_expression, 1.5), -1.5)
  ) |>
  ungroup() |>
  mutate(
    cell_state = factor(cell_state, levels = rev(state_order)),
    gene = factor(gene, levels = marker_present),
    x = as.numeric(gene)
  )
write.csv(
  dot_df,
  file.path(src_dir, "source_data_C_OC_marker_dotplot.csv"),
  row.names = FALSE
)

header_df <- tibble(
  program = rep(names(marker_programs), lengths(marker_programs)),
  gene = unlist(marker_programs, use.names = FALSE)
) |>
  filter(gene %in% marker_present) |>
  mutate(x = match(gene, marker_present)) |>
  group_by(program) |>
  summarise(
    xmin = min(x) - 0.5, xmax = max(x) + 0.5,
    x = mean(c(xmin, xmax)), .groups = "drop"
  )

p_header <- ggplot(header_df) +
  geom_rect(
    aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = 1, fill = program),
    colour = NA
  ) +
  geom_text(
    aes(x, 0.5, label = program), family = "Arial", fontface = "bold",
    colour = "white", size = 4.1
  ) +
  scale_fill_manual(values = marker_program_colors) +
  scale_x_continuous(limits = c(0.5, length(marker_present) + 0.5), expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme_void(base_family = "Arial") +
  theme(legend.position = "none", plot.margin = margin(0, 8, 0, 6))

p_dot <- ggplot(
  dot_df,
  aes(x, cell_state, size = percent_expressed, colour = scaled_average_expression)
) +
  geom_point(alpha = 0.98) +
  scale_x_continuous(
    breaks = seq_along(marker_present), labels = marker_present,
    limits = c(0.5, length(marker_present) + 0.5), expand = c(0, 0)
  ) +
  scale_y_discrete(drop = FALSE) +
  scale_size_area(
    name = "Percent expressed", max_size = 8.0,
    breaks = c(0, 25, 50, 75, 100), limits = c(0, 100)
  ) +
  scale_colour_gradient2(
    low = "#527BA8", mid = "#F7F7F7", high = "#B33A2B", midpoint = 0,
    limits = c(-1.5, 1.5), oob = scales::squish,
    name = "Average expression\n(scaled by gene)"
  ) +
  theme_pub +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 48, hjust = 1, vjust = 1, size = 9.1),
    legend.position = "right",
    plot.margin = margin(0, 8, 7, 6)
  )

p_marker <- p_header / p_dot +
  plot_layout(heights = c(0.14, 1)) +
  plot_annotation(
    title = "Osteoclast commitment, fusion and resorption markers",
    subtitle = paste0(
      "Dot size indicates detected-cell fraction; colour indicates gene-wise scaled mean ",
      "log-normalized expression"
    ),
    theme = theme(
      text = element_text(family = "Arial", face = "bold", colour = "black"),
      plot.title = element_text(size = 13, face = "bold"),
      plot.subtitle = element_text(size = 9.1, face = "bold", colour = "#3A3A3A"),
      plot.margin = margin(7, 7, 3, 7)
    )
  )
save_plot(
  p_marker,
  "C_OC_commitment_fusion_resorption_marker_dotplot",
  213, 116
)

# D. Occupancy of the OC-enriched DDRTree state ----------------------------
traj_dir <- file.path(root, "results", "trajectory", "monocle2_joint_myeloid_trem2_oc")
coord <- read.csv(
  file.path(traj_dir, "source_data_B_Monocle2_DDRTree_pseudotime.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
coord$cell_state <- factor(coord$cell_state, levels = state_order)
state_composition <- coord |>
  count(State, cell_state, name = "n_cells") |>
  group_by(State) |>
  mutate(
    state_total = sum(n_cells),
    fraction_within_DDRTree_state = n_cells / state_total
  ) |>
  ungroup()
oc_enriched_state <- state_composition |>
  filter(cell_state == "OC") |>
  arrange(desc(fraction_within_DDRTree_state), desc(n_cells)) |>
  slice(1) |>
  pull(State)
oc_state_total <- state_composition |>
  filter(State == oc_enriched_state) |>
  summarise(n = sum(n_cells)) |>
  pull(n)
oc_state_oc_n <- state_composition |>
  filter(State == oc_enriched_state, cell_state == "OC") |>
  pull(n_cells)

state_occupancy <- coord |>
  count(cell_state, State, name = "n_cells") |>
  group_by(cell_state) |>
  mutate(
    state_total_for_cell_state = sum(n_cells),
    percent_in_state = 100 * n_cells / state_total_for_cell_state
  ) |>
  ungroup() |>
  filter(State == oc_enriched_state) |>
  complete(
    cell_state = factor(state_order, levels = state_order),
    fill = list(n_cells = 0, percent_in_state = 0)
  ) |>
  mutate(
    cell_state = factor(cell_state, levels = state_order),
    label = sprintf("%.1f%%", percent_in_state)
  )
write.csv(
  state_composition,
  file.path(src_dir, "source_data_D_DDRTree_state_composition.csv"),
  row.names = FALSE
)
write.csv(
  state_occupancy,
  file.path(src_dir, "source_data_D_OC_enriched_state_occupancy.csv"),
  row.names = FALSE
)

p_occupancy <- ggplot(
  state_occupancy,
  aes(cell_state, percent_in_state, fill = cell_state)
) +
  geom_col(width = 0.68, colour = "black", linewidth = 0.38) +
  geom_text(
    aes(label = label), vjust = -0.45, family = "Arial", fontface = "bold",
    size = 3.7, colour = "black"
  ) +
  scale_fill_manual(values = state_colors, drop = FALSE) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(
    limits = c(0, max(state_occupancy$percent_in_state) * 1.18 + 1),
    expand = expansion(mult = c(0, 0))
  ) +
  theme_pub +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 24, hjust = 1, vjust = 1)
  ) +
  labs(
    x = NULL,
    y = paste0("Cells in OC-enriched DDRTree state ", oc_enriched_state, " (%)"),
    title = "TREM2+ macrophages preferentially occupy the OC-enriched branch state",
    subtitle = paste0(
      "DDRTree state ", oc_enriched_state, " contains ", oc_state_oc_n, "/",
      oc_state_total, " OC cells (",
      sprintf("%.1f", 100 * oc_state_oc_n / oc_state_total), "%); topology is exploratory"
    )
  )
save_plot(
  p_occupancy,
  "D_OC_enriched_DDRTree_state_occupancy",
  177, 108
)

# E. Map the independent GO score onto the existing Monocle2 ordering -------
cds <- readRDS(file.path(
  root, "objects", "10_joint_myeloid_trem2_oc_monocle2_exploratory.rds"
))
go_score_col <- score_cols[["GO multinuclear OC differentiation"]]
go_score_by_cell <- score_matrix[[go_score_col]]
names(go_score_by_cell) <- score_matrix$cell_id
if (!all(colnames(cds) %in% names(go_score_by_cell))) {
  stop("Some Monocle2 cells are missing from the UCell score table.")
}
pData(cds)$GO_multinuclear_OC_UCell <- go_score_by_cell[colnames(cds)]

coord_cds <- as.data.frame(t(monocle::reducedDimS(cds))) |>
  rownames_to_column("cell_id")
colnames(coord_cds)[2:3] <- c("Component_1", "Component_2")
coord_cds <- bind_cols(
  coord_cds,
  as.data.frame(pData(cds)[
    coord_cds$cell_id,
    c("cell_state", "State", "Pseudotime", "GO_multinuclear_OC_UCell"),
    drop = FALSE
  ])
)
label_state <- coord_cds |>
  group_by(cell_state) |>
  summarise(
    Component_1 = median(Component_1),
    Component_2 = median(Component_2),
    .groups = "drop"
  )
write.csv(
  coord_cds,
  file.path(src_dir, "source_data_E_GO_OC_score_on_Monocle2_trajectory.csv"),
  row.names = FALSE
)

p_score_trajectory <- monocle::plot_cell_trajectory(
  cds,
  color_by = "GO_multinuclear_OC_UCell",
  show_branch_points = FALSE,
  show_state_number = FALSE,
  cell_size = 1.0
) +
  scale_colour_gradientn(
    colours = c("#D8E6E2", "#79B8A9", "#DDA85C", "#B55D78", "#713D68"),
    limits = range(pData(cds)$GO_multinuclear_OC_UCell, na.rm = TRUE),
    name = "GO multinuclear OC\nUCell score"
  ) +
  ggrepel::geom_label_repel(
    data = label_state,
    aes(Component_1, Component_2, label = cell_state),
    inherit.aes = FALSE,
    family = "Arial", fontface = "bold", size = 3.25,
    colour = "black", fill = scales::alpha("white", 0.86),
    label.size = 0, box.padding = 0.34, point.padding = 0.1,
    seed = seed, max.overlaps = Inf, min.segment.length = 0
  ) +
  theme_pub +
  theme(
    axis.title = element_blank(), axis.text = element_blank(),
    axis.ticks = element_blank(), axis.line = element_blank(),
    legend.position = "right"
  ) +
  labs(
    title = "TREM2-independent OC differentiation score follows the OC arm",
    subtitle = paste0(
      "GO multinuclear osteoclast differentiation score mapped onto the existing ",
      "exploratory Monocle2 DDRTree ordering"
    )
  )
save_plot(
  p_score_trajectory,
  "E_GO_multinuclear_OC_score_on_Monocle2_trajectory",
  184, 114
)

# Reproducibility and export QA ---------------------------------------------
analysis_manifest <- tibble(
  field = c(
    "analysis_scope", "comparison", "primary_signature_source",
    "annotation_marker_excluded", "inference_level", "object",
    "monocle2_object", "seed", "R_version", "UCell_version"
  ),
  value = c(
    "Inflam Mono, Resident Mac, TREM2+ Mac and OC only",
    "TREM2+ macrophage OC differentiation propensity versus other non-OC myeloid states",
    "MSigDB v2025.1.Hs GOBP_MULTINUCLEAR_OSTEOCLAST_DIFFERENTIATION",
    "TREM2",
    "Descriptive per-cell analysis; one participant per condition; no donor-level P values",
    "06_base_broad_lineages_reviewed.rds",
    "10_joint_myeloid_trem2_oc_monocle2_exploratory.rds",
    as.character(seed),
    R.version.string,
    as.character(packageVersion("UCell"))
  )
)
write.csv(
  analysis_manifest,
  file.path(src_dir, "analysis_manifest.csv"),
  row.names = FALSE
)

figure_stems <- c(
  "A_TREM2_independent_OC_UCell_scores_violin",
  "B_condition_aware_OC_UCell_score_summary",
  "C_OC_commitment_fusion_resorption_marker_dotplot",
  "D_OC_enriched_DDRTree_state_occupancy",
  "E_GO_multinuclear_OC_score_on_Monocle2_trajectory"
)
exts <- c("svg", "pdf", "tiff", "png")
qa <- tidyr::crossing(figure = figure_stems, extension = exts) |>
  mutate(
    path = file.path(out_dir, paste0(figure, ".", extension)),
    exists = file.exists(path),
    bytes = ifelse(exists, file.info(path)$size, NA_real_),
    nonempty = exists & bytes > 0
  )
write.csv(
  qa,
  file.path(src_dir, "figure_export_QA_manifest.csv"),
  row.names = FALSE
)
if (!all(qa$nonempty)) stop("At least one requested publication export is missing or empty.")

capture.output(sessionInfo(), file = file.path(src_dir, "sessionInfo.txt"))
message("Completed TREM2+ macrophage OC differentiation-propensity figures: ", out_dir)
