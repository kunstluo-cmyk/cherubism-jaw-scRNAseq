root <- Sys.getenv("CB_REVISION_ROOT")
if (!nzchar(root)) root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(svglite)
  library(ragg)
})

set.seed(20260712)
future::plan("sequential")

out_dir <- file.path(root, "figures", "publication_subfigures", "cellchat_figure6_rebuild")
res_dir <- file.path(root, "results", "cellchat_figure6_rebuild")
resource_dir <- file.path(root, "resources", "nichenet_v2")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)

fine_colors <- c(
  "FibroStromal" = "#8DD3C7",
  "ChondroFibroStromal" = "#FCCDE5",
  "LesionStromal" = "#BEBADA",
  "ECM-remodeling stromal state" = "#BC80BD",
  "TREM2+ Mac" = "#33A02C",
  "OC" = "#1F78B4",
  "Resident Mac" = "#80B1D3"
)

theme_pub <- theme_classic(base_size = 11.5, base_family = "Arial") +
  theme(
    text = element_text(family = "Arial", face = "bold", colour = "black"),
    axis.title = element_text(size = 11.5, face = "bold"),
    axis.text = element_text(size = 9.5, face = "bold", colour = "black"),
    axis.line = element_line(linewidth = 0.5, colour = "black"),
    axis.ticks = element_line(linewidth = 0.45, colour = "black"),
    strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.5),
    strip.text = element_text(size = 10.5, face = "bold", colour = "black"),
    legend.title = element_text(size = 10.5, face = "bold"),
    legend.text = element_text(size = 9.2, face = "bold"),
    plot.title = element_text(size = 13, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 9.2, face = "bold", colour = "#3A3A3A"),
    panel.grid = element_blank(),
    plot.margin = margin(8, 10, 8, 8)
  )

save_plot_pixels <- function(plot, stem, width_px, height_px, dpi = 600) {
  path <- file.path(out_dir, stem)
  wi <- width_px / dpi
  hi <- height_px / dpi
  svglite::svglite(paste0(path, ".svg"), width = wi, height = hi)
  print(plot)
  dev.off()
  grDevices::cairo_pdf(paste0(path, ".pdf"), width = wi, height = hi, family = "Arial")
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

get_rna_data <- function(seu) {
  DefaultAssay(seu) <- "RNA"
  LayerData(seu, assay = "RNA", layer = "data")
}

fraction_expressing <- function(mat, cells) {
  if (!length(cells)) return(setNames(rep(0, nrow(mat)), rownames(mat)))
  Matrix::rowMeans(mat[, cells, drop = FALSE] > 0)
}

# Average precision and AUROC are computed directly from the official NicheNet
# v2 ligand-target prior. This avoids a package-version dependency while retaining
# the standard NicheNet single-ligand ranking quantities.
average_precision <- function(score, response) {
  keep <- is.finite(score) & !is.na(response)
  score <- score[keep]
  response <- as.logical(response[keep])
  if (!length(response) || !any(response) || all(response)) return(NA_real_)
  ord <- order(score, decreasing = TRUE, na.last = NA)
  y <- response[ord]
  precision <- cumsum(y) / seq_along(y)
  mean(precision[y])
}

rank_auroc <- function(score, response) {
  keep <- is.finite(score) & !is.na(response)
  score <- score[keep]
  response <- as.logical(response[keep])
  n_pos <- sum(response)
  n_neg <- sum(!response)
  if (n_pos == 0L || n_neg == 0L) return(NA_real_)
  ranks <- rank(score, ties.method = "average")
  (sum(ranks[response]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

score_ligand_activities <- function(ltm, geneset, background, ligands, context) {
  background <- intersect(background, rownames(ltm))
  geneset <- intersect(geneset, background)
  ligands <- intersect(ligands, colnames(ltm))
  response <- background %in% geneset
  baseline <- mean(response)
  bind_rows(lapply(ligands, function(ligand_i) {
    score <- as.numeric(ltm[background, ligand_i])
    tibble(
      receiver_context = context,
      ligand = ligand_i,
      auroc = rank_auroc(score, response),
      aupr = average_precision(score, response),
      aupr_corrected = average_precision(score, response) - baseline,
      pearson = suppressWarnings(cor(score, as.numeric(response), method = "pearson")),
      geneset_size = length(geneset),
      background_size = length(background),
      positive_fraction = baseline
    )
  }))
}

seu <- readRDS(file.path(root, "objects", "06_base_broad_lineages_reviewed.rds"))
DefaultAssay(seu) <- "RNA"
seu$receiver_state <- as.character(seu$celltype_reviewed_provisional)
seu$receiver_state[seu$receiver_state == "Resident Mac"] <- "Resident Mac"

state_order <- c("Resident Mac", "TREM2+ Mac", "OC")
cells_receiver <- colnames(seu)[seu$receiver_state %in% state_order]
seu_receiver <- subset(seu, cells = cells_receiver)
seu_receiver$receiver_state <- factor(seu_receiver$receiver_state, levels = state_order)
Idents(seu_receiver) <- "receiver_state"

receiver_counts <- as.data.frame(table(seu_receiver$condition, seu_receiver$receiver_state))
colnames(receiver_counts) <- c("condition", "receiver_state", "n_cells")
write.csv(receiver_counts, file.path(res_dir, "GH_receiver_cell_count_audit.csv"), row.names = FALSE)

# Pooled, descriptive cell-state contrasts. Cell-level adjusted P values are used
# only for ranking targets and are not interpreted as donor-level inference.
marker_cache <- file.path(res_dir, "GH_pooled_receiver_positive_markers.rds")
if (file.exists(marker_cache)) {
  marker_list <- readRDS(marker_cache)
} else {
  marker_list <- list(
    "TREM2+ Mac program" = FindMarkers(
      seu_receiver, ident.1 = "TREM2+ Mac", ident.2 = "Resident Mac",
      assay = "RNA", slot = "data", test.use = "wilcox",
      only.pos = TRUE, min.pct = 0.10, logfc.threshold = 0.25, verbose = FALSE
    ) |> rownames_to_column("gene"),
    "OC program" = FindMarkers(
      seu_receiver, ident.1 = "OC", ident.2 = "Resident Mac",
      assay = "RNA", slot = "data", test.use = "wilcox",
      only.pos = TRUE, min.pct = 0.10, logfc.threshold = 0.25, verbose = FALSE
    ) |> rownames_to_column("gene")
  )
  saveRDS(marker_list, marker_cache)
}

marker_long <- bind_rows(lapply(names(marker_list), function(nm) {
  marker_list[[nm]] |> mutate(receiver_context = nm)
}))
write.csv(marker_long, file.path(res_dir, "GH_pooled_receiver_positive_markers.csv"), row.names = FALSE)

mat_receiver <- get_rna_data(seu_receiver)
cells_by_state <- split(colnames(seu_receiver), as.character(seu_receiver$receiver_state))
expr_frac <- bind_rows(lapply(names(cells_by_state), function(state_i) {
  tibble(
    gene = rownames(mat_receiver),
    receiver_state = state_i,
    pct_expressed = fraction_expressing(mat_receiver, cells_by_state[[state_i]])
  )
}))
write.csv(expr_frac, file.path(res_dir, "GH_receiver_gene_detection_fraction.csv"), row.names = FALSE)

matrix_file <- file.path(resource_dir, "ligand_target_matrix_nsga2r_final.rds")
if (!file.exists(matrix_file)) {
  stop("Official NicheNet v2 ligand-target matrix is missing: ", matrix_file)
}
ligand_target_matrix <- readRDS(matrix_file)
stopifnot(is.matrix(ligand_target_matrix) || inherits(ligand_target_matrix, "Matrix"))

candidate_ligands <- scan(
  file.path(res_dir, "candidate_ligands_for_nichenet.txt"),
  what = character(), quiet = TRUE
)
candidate_ligands <- unique(intersect(candidate_ligands, colnames(ligand_target_matrix)))

# Require expression in at least one CB stromal sender and CellChat support.
sender_states <- c(
  "FibroStromal", "ChondroFibroStromal", "LesionStromal",
  "ECM-remodeling stromal state"
)
seu$fine_nichenet <- as.character(seu$celltype_reviewed_provisional)
seu$fine_nichenet[seu$fine_nichenet == "Act OB-Stromal"] <- "ECM-remodeling stromal state"
sender_cells <- colnames(seu)[seu$condition == "CB" & seu$fine_nichenet %in% sender_states]
sender_mat <- get_rna_data(seu)[, sender_cells, drop = FALSE]
sender_fraction <- fraction_expressing(sender_mat, sender_cells)
candidate_ligands <- candidate_ligands[
  !is.na(sender_fraction[candidate_ligands]) & sender_fraction[candidate_ligands] >= 0.05
]

activity_list <- list()
geneset_list <- list()
background_list <- list()
for (context_i in names(marker_list)) {
  focal_state <- if (context_i == "TREM2+ Mac program") "TREM2+ Mac" else "OC"
  markers_i <- marker_list[[context_i]] |>
    filter(p_val_adj < 0.05, avg_log2FC > 0.25, pct.1 >= 0.10) |>
    arrange(desc(avg_log2FC), p_val_adj)
  geneset_i <- head(markers_i$gene, 600)
  background_i <- expr_frac |>
    filter(receiver_state == focal_state, pct_expressed >= 0.10) |>
    pull(gene)
  geneset_i <- intersect(geneset_i, background_i)
  geneset_list[[context_i]] <- geneset_i
  background_list[[context_i]] <- background_i
  activity_list[[context_i]] <- score_ligand_activities(
    ligand_target_matrix, geneset_i, background_i, candidate_ligands, context_i
  )
}

activities <- bind_rows(activity_list) |>
  mutate(
    receiver_context = factor(
      receiver_context,
      levels = c("TREM2+ Mac program", "OC program")
    ),
    sender_detection = sender_fraction[ligand],
    aupr_corrected = replace_na(aupr_corrected, -Inf),
    pearson = replace_na(pearson, -Inf)
  ) |>
  group_by(receiver_context) |>
  arrange(desc(aupr_corrected), desc(pearson), .by_group = TRUE) |>
  mutate(activity_rank = row_number()) |>
  ungroup()
write.csv(activities, file.path(res_dir, "G_NicheNet_v2_candidate_ligand_activities.csv"), row.names = FALSE)

top_activities <- activities |>
  group_by(receiver_context) |>
  slice_head(n = 10) |>
  ungroup() |>
  mutate(
    ligand_facet = paste(receiver_context, ligand, sep = "___"),
    ligand_facet = factor(ligand_facet, levels = rev(unique(ligand_facet)))
  )

p_g <- ggplot(top_activities, aes(x = aupr_corrected, y = ligand_facet, fill = receiver_context)) +
  geom_col(width = 0.72, colour = "black", linewidth = 0.25) +
  geom_text(
    aes(label = sprintf("%.3f", aupr_corrected)),
    hjust = ifelse(top_activities$aupr_corrected >= 0, -0.10, 1.10),
    size = 2.8, family = "Arial", fontface = "bold"
  ) +
  facet_grid(receiver_context ~ ., scales = "free_y", space = "free_y") +
  scale_y_discrete(labels = function(x) sub(".*___", "", x)) +
  scale_fill_manual(values = c("TREM2+ Mac program" = fine_colors[["TREM2+ Mac"]],
                               "OC program" = fine_colors[["OC"]])) +
  scale_x_continuous(expand = expansion(mult = c(0.08, 0.20))) +
  labs(
    title = "NicheNet v2 candidate-ligand activity",
    subtitle = "CB stromal ligands ranked against pooled descriptive receiver-state gene sets",
    x = "Corrected area under the precision-recall curve",
    y = NULL
  ) +
  theme_pub +
  theme(legend.position = "none")

save_plot_pixels(p_g, "G_NicheNet_v2_candidate_ligand_activity_ranking", 3900, 3300)

# Extract top ligand-target links. Regulatory potential is displayed after
# within-ligand scaling because absolute NicheNet weights vary substantially.
links_all <- list()
for (context_i in names(geneset_list)) {
  ligands_i <- activities |>
    filter(receiver_context == context_i) |>
    slice_head(n = 6) |>
    pull(ligand)
  geneset_i <- geneset_list[[context_i]]
  for (ligand_i in ligands_i) {
    weights_i <- ligand_target_matrix[, ligand_i]
    names(weights_i) <- rownames(ligand_target_matrix)
    top250 <- names(sort(weights_i, decreasing = TRUE))[seq_len(min(250, length(weights_i)))]
    targets_i <- intersect(top250, geneset_i)
    if (!length(targets_i)) next
    targets_i <- head(targets_i[order(weights_i[targets_i], decreasing = TRUE)], 8)
    links_all[[paste(context_i, ligand_i, sep = "__")]] <- tibble(
      receiver_context = context_i,
      ligand = ligand_i,
      target = targets_i,
      regulatory_potential = as.numeric(weights_i[targets_i])
    )
  }
}
links <- bind_rows(links_all) |>
  group_by(receiver_context, ligand) |>
  mutate(regulatory_potential_scaled = regulatory_potential / max(regulatory_potential, na.rm = TRUE)) |>
  ungroup()
write.csv(links, file.path(res_dir, "H_NicheNet_v2_top_ligand_target_links.csv"), row.names = FALSE)

# Limit each receiver panel to the most informative nonredundant target genes.
target_order <- links |>
  group_by(receiver_context, target) |>
  summarise(max_weight = max(regulatory_potential, na.rm = TRUE), .groups = "drop") |>
  group_by(receiver_context) |>
  slice_max(max_weight, n = 22, with_ties = FALSE) |>
  arrange(receiver_context, max_weight) |>
  ungroup()

heat_df <- links |>
  semi_join(target_order, by = c("receiver_context", "target")) |>
  group_by(receiver_context, ligand, target) |>
  summarise(regulatory_potential_scaled = max(regulatory_potential_scaled), .groups = "drop") |>
  group_by(receiver_context) |>
  complete(ligand, target, fill = list(regulatory_potential_scaled = 0)) |>
  ungroup() |>
  left_join(target_order, by = c("receiver_context", "target")) |>
  mutate(
    target_facet = paste(receiver_context, target, sep = "___"),
    ligand_facet = paste(receiver_context, ligand, sep = "___")
  )

target_levels <- target_order |>
  mutate(target_facet = paste(receiver_context, target, sep = "___")) |>
  pull(target_facet)
ligand_levels <- activities |>
  group_by(receiver_context) |>
  slice_head(n = 6) |>
  ungroup() |>
  mutate(ligand_facet = paste(receiver_context, ligand, sep = "___")) |>
  pull(ligand_facet)
heat_df$target_facet <- factor(heat_df$target_facet, levels = target_levels)
heat_df$ligand_facet <- factor(heat_df$ligand_facet, levels = rev(ligand_levels))

make_target_heatmap <- function(context_i, show_x_title = TRUE) {
  df_i <- heat_df |> filter(receiver_context == context_i)
  ggplot(df_i, aes(x = target_facet, y = ligand_facet, fill = regulatory_potential_scaled)) +
    geom_tile(colour = "white", linewidth = 0.35) +
    scale_x_discrete(labels = function(x) sub(".*___", "", x), expand = c(0, 0), drop = TRUE) +
    scale_y_discrete(labels = function(x) sub(".*___", "", x), expand = c(0, 0), drop = TRUE) +
    scale_fill_gradientn(
      colours = c("#F7FBFF", "#C6DBEF", "#6BAED6", "#2171B5", "#08306B"),
      limits = c(0, 1), name = "Regulatory\npotential\n(within ligand)"
    ) +
    labs(
      title = context_i,
      x = if (show_x_title) "Predicted receiver target gene" else NULL,
      y = "Candidate stromal ligand"
    ) +
    theme_pub +
    theme(
      plot.title = element_text(
        hjust = 0.5, size = 11.5,
        colour = if (context_i == "TREM2+ Mac program") fine_colors[["TREM2+ Mac"]] else fine_colors[["OC"]]
      ),
      axis.text.x = element_text(angle = 50, hjust = 1, vjust = 1, size = 8.4),
      legend.position = "right"
    )
}

p_h <- make_target_heatmap("TREM2+ Mac program", show_x_title = FALSE) /
  make_target_heatmap("OC program", show_x_title = TRUE) +
  plot_layout(guides = "collect", heights = c(1, 1)) &
  theme(legend.position = "right")
p_h <- p_h + plot_annotation(
  title = "Top NicheNet ligand-target links",
  subtitle = "Targets are positively enriched in pooled descriptive receiver-state contrasts",
  theme = theme(
    plot.title = element_text(family = "Arial", face = "bold", size = 13, colour = "black"),
    plot.subtitle = element_text(family = "Arial", face = "bold", size = 9.2, colour = "#3A3A3A")
  )
)

save_plot_pixels(p_h, "H_NicheNet_v2_ligand_target_heatmap", 5400, 3600)

writeLines(capture.output(sessionInfo()), file.path(res_dir, "sessionInfo_31b.txt"))

summary_lines <- c(
  paste0("Official NicheNet v2 matrix: ", normalizePath(matrix_file, winslash = "/")),
  paste0("Candidate ligands after CellChat/sender-expression filtering: ", length(candidate_ligands)),
  paste0("TREM2+ Mac target-gene set in model/background: ",
         activity_list[["TREM2+ Mac program"]]$geneset_size[[1]]),
  paste0("OC target-gene set in model/background: ",
         activity_list[["OC program"]]$geneset_size[[1]]),
  "Important: pooled cell-state contrasts are descriptive and are not donor-level tests."
)
writeLines(summary_lines, file.path(res_dir, "GH_analysis_notes.txt"))
