root <- Sys.getenv("CB_REVISION_ROOT")
if (!nzchar(root)) root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
  library(scales)
  library(svglite)
  library(ragg)
})

set.seed(20260712)
future::plan("sequential")

out_dir <- file.path(root, "figures", "publication_subfigures", "cellchat_figure6_rebuild")
res_dir <- file.path(root, "results", "cellchat_figure6_rebuild")
obj_dir <- file.path(root, "objects", "cellchat_figure6_rebuild")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(obj_dir, recursive = TRUE, showWarnings = FALSE)

condition_colors <- c("CTRL" = "#4DBBD5", "CB" = "#E64B35")
broad_order <- c("Stromal/osteogenic", "Mural", "Endothelial", "Myeloid", "Lymphoid")
broad_colors <- c(
  "Stromal/osteogenic" = "#8DD3C7",
  "Mural" = "#7393B3",
  "Endothelial" = "#43AA8B",
  "Myeloid" = "#33A02C",
  "Lymphoid" = "#CCEBC5"
)
fine_order_cb <- c(
  "FibroStromal", "ChondroFibroStromal", "LesionStromal",
  "ECM-remodeling stromal state", "TREM2+ Mac", "OC"
)
fine_colors <- c(
  "FibroStromal" = "#8DD3C7",
  "ChondroFibroStromal" = "#FCCDE5",
  "LesionStromal" = "#BEBADA",
  "ECM-remodeling stromal state" = "#BC80BD",
  "TREM2+ Mac" = "#33A02C",
  "OC" = "#1F78B4"
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

save_base_pixels <- function(stem, width_px, height_px, draw_fun, dpi = 600) {
  path <- file.path(out_dir, stem)
  wi <- width_px / dpi
  hi <- height_px / dpi
  svglite::svglite(paste0(path, ".svg"), width = wi, height = hi)
  draw_fun()
  dev.off()
  grDevices::cairo_pdf(paste0(path, ".pdf"), width = wi, height = hi, family = "Arial")
  draw_fun()
  dev.off()
  ragg::agg_tiff(
    paste0(path, ".tiff"), width = width_px, height = height_px,
    units = "px", res = dpi, compression = "lzw"
  )
  draw_fun()
  dev.off()
  ragg::agg_png(
    paste0(path, ".png"), width = width_px, height = height_px,
    units = "px", res = dpi
  )
  draw_fun()
  dev.off()
}

get_rna_data <- function(seu) {
  DefaultAssay(seu) <- "RNA"
  LayerData(seu, assay = "RNA", layer = "data")
}

build_cellchat <- function(seu, labels, min_cells = 30L, nboot = 200L) {
  stopifnot(length(labels) == ncol(seu))
  data_input <- get_rna_data(seu)
  meta_input <- data.frame(labels = droplevels(labels), row.names = colnames(seu))
  cc <- createCellChat(object = data_input, meta = meta_input, group.by = "labels")
  db_secreted <- subsetDB(CellChatDB.human, search = "Secreted Signaling", key = "annotation")
  cc@DB <- db_secreted
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  cc <- computeCommunProb(
    cc,
    type = "triMean",
    raw.use = TRUE,
    population.size = FALSE,
    distance.use = FALSE,
    nboot = nboot,
    seed.use = 20260712
  )
  cc <- filterCommunication(cc, min.cells = min_cells)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  cc
}

extract_communication_fdr <- function(cc) {
  x <- subsetCommunication(cc)
  if (is.null(x) || nrow(x) == 0L) return(tibble())
  as_tibble(x) |>
    mutate(
      qval = p.adjust(pval, method = "BH"),
      source = as.character(source),
      target = as.character(target)
    )
}

seu <- readRDS(file.path(root, "objects", "06_base_broad_lineages_reviewed.rds"))
stopifnot(all(c("condition", "broad_lineage_reviewed", "celltype_reviewed_provisional") %in%
                colnames(seu@meta.data)))

seu$fine_cellchat <- as.character(seu$celltype_reviewed_provisional)
seu$fine_cellchat[seu$fine_cellchat == "Act OB-Stromal"] <- "ECM-remodeling stromal state"

cell_count_audit <- seu@meta.data |>
  as_tibble(rownames = "cell") |>
  count(condition, broad_lineage_reviewed, fine_cellchat, name = "n_cells")
write.csv(cell_count_audit, file.path(res_dir, "cell_count_audit.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 1. Comparable broad-lineage CellChat analysis
# -----------------------------------------------------------------------------
broad_cells <- colnames(seu)[seu$condition %in% c("CTRL", "CB") &
                              seu$broad_lineage_reviewed %in% broad_order]
seu_broad <- subset(seu, cells = broad_cells)
seu_broad$broad_cellchat <- factor(seu_broad$broad_lineage_reviewed, levels = broad_order)

broad_counts <- as.data.frame(with(seu_broad@meta.data, table(condition, broad_cellchat)))
colnames(broad_counts) <- c("condition", "cell_state", "n_cells")
write.csv(broad_counts, file.path(res_dir, "broad_cell_counts_by_condition.csv"), row.names = FALSE)

cc_broad <- list()
for (cond in c("CTRL", "CB")) {
  cache_file <- file.path(obj_dir, paste0("cellchat_broad_", cond, ".rds"))
  if (file.exists(cache_file)) {
    message("Loading cached broad CellChat: ", cond)
    cc_broad[[cond]] <- readRDS(cache_file)
  } else {
    message("Running broad CellChat: ", cond)
    cells_i <- colnames(seu_broad)[seu_broad$condition == cond]
    obj_i <- subset(seu_broad, cells = cells_i)
    obj_i$broad_cellchat <- factor(obj_i$broad_cellchat, levels = broad_order)
    cc_broad[[cond]] <- build_cellchat(obj_i, obj_i$broad_cellchat, min_cells = 30, nboot = 200)
    saveRDS(cc_broad[[cond]], cache_file)
  }
  df_i <- extract_communication_fdr(cc_broad[[cond]])
  write.csv(df_i, file.path(res_dir, paste0("broad_interactions_", cond, ".csv")), row.names = FALSE)
}

get_weight_matrix <- function(cc, order_use) {
  m <- cc@net$weight
  out <- matrix(0, nrow = length(order_use), ncol = length(order_use),
                dimnames = list(order_use, order_use))
  rr <- intersect(rownames(m), order_use)
  cc2 <- intersect(colnames(m), order_use)
  out[rr, cc2] <- m[rr, cc2, drop = FALSE]
  out
}

mat_ctrl <- get_weight_matrix(cc_broad$CTRL, broad_order)
mat_cb <- get_weight_matrix(cc_broad$CB, broad_order)
write.csv(as.data.frame(mat_ctrl) |> rownames_to_column("source"),
          file.path(res_dir, "A_broad_strength_matrix_CTRL.csv"), row.names = FALSE)
write.csv(as.data.frame(mat_cb) |> rownames_to_column("source"),
          file.path(res_dir, "B_broad_strength_matrix_CB.csv"), row.names = FALSE)

global_edge_max <- max(c(mat_ctrl, mat_cb), na.rm = TRUE)
draw_circle <- function(mat, title_text) {
  display_names <- c("Stromal/\nosteogenic", "Mural", "Endothelial", "Myeloid", "Lymphoid")
  mat_display <- mat
  rownames(mat_display) <- display_names
  colnames(mat_display) <- display_names
  par(family = "Arial", font = 2, mar = c(1.0, 1.0, 3.8, 1.0), xpd = TRUE)
  netVisual_circle(
    mat_display,
    color.use = broad_colors[broad_order],
    title.name = title_text,
    weight.scale = TRUE,
    vertex.weight = rep(1, length(broad_order)),
    vertex.weight.max = 1,
    vertex.size.max = 20,
    vertex.label.cex = 1.18,
    vertex.label.color = "black",
    edge.weight.max = global_edge_max,
    edge.width.max = 8,
    alpha.edge = 0.62,
    label.edge = FALSE,
    arrow.size = 0.16,
    margin = 0.23
  )
}

save_base_pixels(
  "A_CTRL_common_lineage_CellChat_interaction_strength", 3300, 2850,
  function() draw_circle(mat_ctrl, "CTRL: inferred interaction strength")
)
save_base_pixels(
  "B_CB_common_lineage_CellChat_interaction_strength", 3300, 2850,
  function() draw_circle(mat_cb, "CB: inferred interaction strength")
)

# Descriptive effect size; pseudocount is tied to positive values across both matrices.
positive_values <- c(mat_ctrl[mat_ctrl > 0], mat_cb[mat_cb > 0])
pseudo <- if (length(positive_values)) 0.05 * median(positive_values) else 1e-6
log2_ratio <- log2((mat_cb + pseudo) / (mat_ctrl + pseudo))
delta_df <- as.data.frame(as.table(log2_ratio), stringsAsFactors = FALSE) |>
  rename(source = Var1, target = Var2, log2_CB_CTRL = Freq) |>
  mutate(
    source = factor(source, levels = rev(broad_order)),
    target = factor(target, levels = broad_order)
  )
write.csv(delta_df, file.path(res_dir, "C_common_lineage_descriptive_log2_CB_CTRL.csv"), row.names = FALSE)

lim_delta <- max(abs(delta_df$log2_CB_CTRL), na.rm = TRUE)
p_c <- ggplot(delta_df, aes(target, source, fill = log2_CB_CTRL)) +
  geom_tile(colour = "white", linewidth = 0.9) +
  geom_text(aes(label = sprintf("%.1f", log2_CB_CTRL)), size = 3.25, family = "Arial", fontface = "bold") +
  scale_fill_gradient2(
    low = condition_colors[["CTRL"]], mid = "white", high = condition_colors[["CB"]],
    midpoint = 0, limits = c(-lim_delta, lim_delta), oob = squish,
    name = expression(log[2]~"(CB / CTRL)")
  ) +
  coord_equal() +
  labs(
    title = "Common-lineage communication: CB versus CTRL",
    subtitle = "Positive values indicate higher inferred strength in CB; no donor-level inference",
    x = "Receiver", y = "Sender"
  ) +
  theme_pub +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "right")
save_plot_pixels(p_c, "C_common_source_target_CB_minus_CTRL_heatmap", 3900, 3150)

# -----------------------------------------------------------------------------
# 2. CB lesion-stromal to TREM2+ macrophage / OC CellChat analysis
# -----------------------------------------------------------------------------
cb_cells <- colnames(seu)[seu$condition == "CB" & seu$fine_cellchat %in% fine_order_cb]
seu_cb <- subset(seu, cells = cb_cells)
seu_cb$fine_cellchat <- factor(seu_cb$fine_cellchat, levels = fine_order_cb)

fine_counts <- as.data.frame(table(seu_cb$fine_cellchat))
colnames(fine_counts) <- c("cell_state", "n_cells")
write.csv(fine_counts, file.path(res_dir, "CB_fine_cell_counts.csv"), row.names = FALSE)

fine_cache <- file.path(obj_dir, "cellchat_CB_lesion_stromal_TREM2_OC.rds")
if (file.exists(fine_cache)) {
  message("Loading cached CB fine-state CellChat")
  cc_cb_fine <- readRDS(fine_cache)
} else {
  message("Running CB fine-state CellChat")
  cc_cb_fine <- build_cellchat(seu_cb, seu_cb$fine_cellchat, min_cells = 30, nboot = 300)
  saveRDS(cc_cb_fine, fine_cache)
}
df_fine <- extract_communication_fdr(cc_cb_fine)
write.csv(df_fine, file.path(res_dir, "CB_fine_all_interactions.csv"), row.names = FALSE)

source_states <- fine_order_cb[1:4]
target_states <- fine_order_cb[5:6]
focus_df <- df_fine |>
  filter(source %in% source_states, target %in% target_states) |>
  arrange(qval, desc(prob))
write.csv(focus_df, file.path(res_dir, "CB_stromal_to_TREM2_OC_interactions.csv"), row.names = FALSE)

# Pathway-level matrix from CellChat netP.
path_array <- cc_cb_fine@netP$prob
pathways_all <- dimnames(path_array)[[3]]
path_long <- bind_rows(lapply(pathways_all, function(pw) {
  m <- path_array[, , pw, drop = TRUE]
  as.data.frame(as.table(m), stringsAsFactors = FALSE) |>
    rename(source = Var1, target = Var2, strength = Freq) |>
    mutate(pathway = pw)
})) |>
  filter(source %in% source_states, target %in% target_states, strength > 0)

pathway_rank <- path_long |>
  group_by(pathway) |>
  summarise(total_strength = sum(strength), .groups = "drop") |>
  arrange(desc(total_strength))
top_data_pathways <- head(pathway_rank$pathway[pathway_rank$total_strength > 1e-4], 5)
mechanistic_pathways <- intersect(c("CXCL", "TGFb", "CSF", "RANKL"), pathway_rank$pathway)
selected_pathways <- head(unique(c(top_data_pathways, mechanistic_pathways)), 8)
write.csv(pathway_rank, file.path(res_dir, "CB_focus_pathway_ranking.csv"), row.names = FALSE)
writeLines(selected_pathways, file.path(res_dir, "CB_selected_pathways.txt"))

path_plot <- path_long |>
  filter(pathway %in% selected_pathways) |>
  mutate(
    source_target = paste(source, target, sep = "  ->  "),
    source_target = factor(source_target,
                           levels = as.vector(outer(source_states, target_states,
                                                    paste, sep = "  ->  "))),
    pathway = factor(pathway, levels = rev(selected_pathways))
  )
write.csv(path_plot, file.path(res_dir, "D_CB_pathway_strength_source_data.csv"), row.names = FALSE)

p_d <- ggplot(path_plot, aes(source_target, pathway, fill = sqrt(strength))) +
  geom_tile(colour = "white", linewidth = 0.7) +
  scale_fill_gradientn(
    colours = c("#FFF3EC", "#F7B89C", "#E96B4C", "#B7191C"),
    name = "Inferred strength\n(square-root scale)"
  ) +
  labs(
    title = "CB stromal-to-myeloid signaling pathways",
    subtitle = "Lesion stromal sources to TREM2+ macrophages and osteoclasts",
    x = NULL, y = NULL
  ) +
  theme_pub +
  theme(
    axis.text.x = element_text(angle = 42, hjust = 1, size = 8.4),
    axis.text.y = element_text(size = 9.2),
    legend.position = "right"
  )
save_plot_pixels(p_d, "D_CB_stromal_to_TREM2_OC_pathway_strength_heatmap", 4500, 3300)

# Ligand-receptor bubble plot: retain the strongest non-redundant interactions.
lr_rank <- focus_df |>
  filter(pathway_name %in% selected_pathways, qval < 0.05) |>
  group_by(interaction_name_2, pathway_name) |>
  summarise(max_prob = max(prob), min_q = min(qval), .groups = "drop") |>
  arrange(desc(max_prob), min_q)
selected_lr <- head(lr_rank$interaction_name_2, 18)
if (length(selected_lr) < 8) {
  selected_lr <- head(unique(focus_df$interaction_name_2[order(focus_df$prob, decreasing = TRUE)]), 18)
}
lr_plot <- focus_df |>
  filter(interaction_name_2 %in% selected_lr) |>
  mutate(
    interaction = factor(interaction_name_2, levels = rev(selected_lr)),
    source_target = paste(source, target, sep = "  ->  "),
    source_target = factor(source_target,
                           levels = as.vector(outer(source_states, target_states,
                                                    paste, sep = "  ->  "))),
    # CellChat permutation P values have finite resolution (nboot = 300 here).
    # Bound zero-valued estimates at 1/(nboot+1) to avoid false precision.
    neglog10_q = -log10(pmax(qval, 1 / 301))
  )
write.csv(lr_plot, file.path(res_dir, "E_CB_key_ligand_receptor_source_data.csv"), row.names = FALSE)

p_e <- ggplot(lr_plot, aes(source_target, interaction)) +
  geom_point(aes(size = prob, colour = neglog10_q), alpha = 0.92) +
  scale_size_continuous(range = c(1.5, 9), name = "Interaction\nstrength") +
  scale_colour_gradientn(
    colours = c("#FDE0D2", "#F68A67", "#D7301F", "#8E0000"),
    name = expression(-log[10]~"FDR"~"(bounded)")
  ) +
  labs(
    title = "Prioritized ligand-receptor interactions in the CB lesion",
    subtitle = "CellChat-predicted stromal signals received by TREM2+ macrophages or osteoclasts",
    x = NULL, y = NULL
  ) +
  theme_pub +
  theme(axis.text.x = element_text(angle = 42, hjust = 1, size = 8.2),
        axis.text.y = element_text(size = 8.8), legend.position = "right")
save_plot_pixels(p_e, "E_CB_key_ligand_receptor_bubbleplot", 5400, 3900)

# Expression coverage of prespecified biologically interpretable ligands/receptors.
ligand_genes <- c("MIF", "PTN", "MDK", "POSTN", "CXCL12", "CSF1", "TNFSF11",
                  "TGFB1", "TGFB3", "SPP1")
receptor_genes <- c("CD74", "CD44", "CXCR4", "NCL", "ITGAV", "ITGB3", "ITGA4", "ITGB1",
                    "CSF1R", "TNFRSF11A", "TGFBR1", "TGFBR2")
genes_plot <- intersect(c(ligand_genes, receptor_genes), rownames(seu_cb))
expr_mat <- get_rna_data(seu_cb)[genes_plot, , drop = FALSE]
groups <- as.character(seu_cb$fine_cellchat)

expr_summary <- bind_rows(lapply(levels(seu_cb$fine_cellchat), function(g) {
  idx <- which(groups == g)
  if (!length(idx)) return(NULL)
  tibble(
    cell_state = g,
    gene = genes_plot,
    average_expression = Matrix::rowMeans(expr_mat[, idx, drop = FALSE]),
    percent_expressed = Matrix::rowMeans(expr_mat[, idx, drop = FALSE] > 0) * 100
  )
})) |>
  group_by(gene) |>
  mutate(scaled_average = as.numeric(scale(average_expression))) |>
  ungroup() |>
  mutate(
    scaled_average = ifelse(is.finite(scaled_average), scaled_average, 0),
    cell_state = factor(cell_state, levels = rev(fine_order_cb)),
    gene = factor(gene, levels = genes_plot),
    feature_class = ifelse(as.character(gene) %in% ligand_genes, "Ligands", "Receptors")
  )
write.csv(expr_summary, file.path(res_dir, "F_ligand_receptor_expression_source_data.csv"), row.names = FALSE)

p_f <- ggplot(expr_summary, aes(gene, cell_state)) +
  geom_point(aes(size = percent_expressed, colour = scaled_average)) +
  facet_grid(. ~ feature_class, scales = "free_x", space = "free_x") +
  scale_size_continuous(
    range = c(0.45, 6.2), limits = c(0, 100), breaks = c(0, 50, 100),
    name = "Percent expressed"
  ) +
  scale_colour_gradientn(
    colours = c("#FDE2D5", "#F98E6A", "#D7301F", "#8E0000"),
    values = rescale(c(-1, 0, 1, 2)),
    name = "Average expression\n(scaled)"
  ) +
  labs(
    title = "Candidate ligand and receptor expression across CB cell states",
    subtitle = "DecontX-corrected RNA expression in CB stromal sources and myeloid receivers",
    x = NULL, y = NULL
  ) +
  theme_pub +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1,
                               face = "bold.italic", size = 8.8),
    axis.text.y = element_text(size = 9.2),
    strip.text = element_text(size = 10.2, face = "bold"),
    panel.spacing.x = unit(8, "mm"),
    legend.position = "right",
    legend.key.height = unit(6.2, "mm"),
    legend.spacing.y = unit(1.2, "mm"),
    plot.margin = margin(5, 8, 4, 7)
  ) +
  guides(
    colour = guide_colourbar(barheight = unit(15, "mm"), barwidth = unit(4.5, "mm")),
    size = guide_legend(override.aes = list(colour = "black"))
  )
save_plot_pixels(p_f, "F_CB_sender_ligand_receiver_receptor_expression", 6250, 1700)

# Inputs for NicheNet / receiver-program validation.
candidate_ligands <- unique(focus_df$ligand[focus_df$qval < 0.05])
candidate_ligands <- intersect(candidate_ligands, rownames(seu_cb))
writeLines(candidate_ligands, file.path(res_dir, "candidate_ligands_for_nichenet.txt"))
saveRDS(seu_cb, file.path(obj_dir, "seurat_CB_sender_receiver_for_nichenet.rds"))

capture.output(sessionInfo(), file = file.path(res_dir, "sessionInfo_31a.txt"))
message("Stage 31a complete. Outputs: ", out_dir)
