#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)
set.seed(20260911)

suppressPackageStartupMessages({
  library(Seurat)
  library(monocle3)
  library(BiocNeighbors)
  library(igraph)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

root_dir <- Sys.getenv("CB_REVISION_ROOT")
if (!nzchar(root_dir)) {
  root_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
object_file <- file.path(root_dir, "objects", "06_base_broad_lineages_reviewed.rds")
result_dir <- file.path(root_dir, "results", "trajectory")
figure_dir <- file.path(root_dir, "figures", "trajectory_audit")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

state_col <- "celltype_reviewed_provisional"
condition_col <- "condition"

state_levels <- c(
  "LEPR+ MSC-like",
  "CNCC-OCPs",
  "Mature OB",
  "FibroStromal",
  "ChondroFibroStromal",
  "LesionStromal",
  "ECM-remodeling stromal state"
)

state_colors <- c(
  "LEPR+ MSC-like" = "#FB8072",
  "CNCC-OCPs" = "#FFED6F",
  "Mature OB" = "#FDB462",
  "FibroStromal" = "#8DD3C7",
  "ChondroFibroStromal" = "#FCCDE5",
  "LesionStromal" = "#BEBADA",
  "ECM-remodeling stromal state" = "#BC80BD"
)

condition_colors <- c("CTRL" = "#4DBBD5", "CB" = "#E64B35")

theme_audit <- function(base_size = 12) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.6, colour = "black"),
      axis.ticks = element_line(linewidth = 0.5, colour = "black"),
      axis.text = element_text(face = "bold", colour = "black"),
      axis.title = element_text(face = "bold", colour = "black"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", colour = "black"),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0),
      panel.grid = element_blank()
    )
}

save_audit <- function(plot, stem, width = 14, height = 10, dpi = 450) {
  ggsave(file.path(figure_dir, paste0(stem, ".pdf")), plot, width = width, height = height,
         device = cairo_pdf, bg = "white")
  ggsave(file.path(figure_dir, paste0(stem, ".png")), plot, width = width, height = height,
         dpi = dpi, bg = "white")
}

message("Loading reviewed base object")
obj <- readRDS(object_file)
stopifnot(state_col %in% colnames(obj[[]]), condition_col %in% colnames(obj[[]]))

keep_cells <- rownames(obj[[]])[obj[["broad_lineage_reviewed", drop = TRUE]] == "Stromal/osteogenic"]
stromal <- subset(obj, cells = keep_cells)
stromal$stromal_state <- as.character(stromal[[state_col, drop = TRUE]])
stromal$stromal_state[stromal$stromal_state == "Act OB-Stromal"] <-
  "ECM-remodeling stromal state"
stromal$stromal_state <- factor(stromal$stromal_state, levels = state_levels)
stromal$condition <- factor(stromal[[condition_col, drop = TRUE]], levels = c("CTRL", "CB"))

stopifnot(!anyNA(stromal$stromal_state), !anyNA(stromal$condition))

message("Recomputing a joint, non-integrated stromal state space")
DefaultAssay(stromal) <- "RNA"
stromal <- NormalizeData(stromal, normalization.method = "LogNormalize", scale.factor = 1e4,
                         verbose = FALSE)
stromal <- FindVariableFeatures(stromal, selection.method = "vst", nfeatures = 3000,
                                verbose = FALSE)
stromal <- ScaleData(stromal, features = VariableFeatures(stromal), verbose = FALSE)
stromal <- RunPCA(stromal, features = VariableFeatures(stromal), npcs = 30,
                  seed.use = 20260911, verbose = FALSE)
stromal <- RunUMAP(stromal, reduction = "pca", dims = 1:20, n.neighbors = 30,
                   min.dist = 0.30, metric = "cosine", seed.use = 20260911,
                   reduction.name = "umap.stromal_joint", verbose = FALSE)

saveRDS(stromal, file.path(root_dir, "objects", "08_joint_stromal_topology_audit.rds"),
        compress = FALSE)

meta <- stromal[[]] %>%
  tibble::rownames_to_column("cell") %>%
  select(cell, condition, stromal_state)

counts_condition_state <- meta %>%
  count(condition, stromal_state, name = "n_cells", .drop = FALSE) %>%
  group_by(condition) %>%
  mutate(fraction_within_condition = n_cells / sum(n_cells)) %>%
  ungroup()
write.csv(counts_condition_state,
          file.path(result_dir, "joint_stromal_state_counts_by_condition.csv"), row.names = FALSE)

umap <- Embeddings(stromal, "umap.stromal_joint")
pca <- Embeddings(stromal, "pca")[, 1:20, drop = FALSE]
umap_df <- as.data.frame(umap) %>%
  setNames(c("UMAP1", "UMAP2")) %>%
  tibble::rownames_to_column("cell") %>%
  left_join(meta, by = "cell")
write.csv(umap_df, file.path(result_dir, "source_data_joint_stromal_umap.csv"), row.names = FALSE)

label_df <- umap_df %>%
  group_by(stromal_state) %>%
  summarise(UMAP1 = median(UMAP1), UMAP2 = median(UMAP2), .groups = "drop")

p_state <- ggplot(umap_df, aes(UMAP1, UMAP2, colour = stromal_state)) +
  geom_point(size = 0.35, alpha = 0.82, stroke = 0) +
  ggrepel::geom_text_repel(
    data = label_df,
    aes(label = stromal_state),
    colour = "black", fontface = "bold", size = 3.5,
    box.padding = 0.35, point.padding = 0.2, segment.color = NA,
    show.legend = FALSE
  ) +
  scale_colour_manual(values = state_colors, drop = FALSE) +
  coord_equal() +
  labs(title = "Joint stromal/osteogenic state space", x = "UMAP1", y = "UMAP2",
       colour = "Cell state") +
  theme_audit() +
  theme(legend.position = "none")

p_condition <- ggplot(umap_df, aes(UMAP1, UMAP2, colour = condition)) +
  geom_point(size = 0.35, alpha = 0.82, stroke = 0) +
  facet_wrap(~condition, nrow = 1) +
  scale_colour_manual(values = condition_colors, drop = FALSE) +
  coord_equal() +
  labs(title = "Identical coordinates shown separately by specimen",
       x = "UMAP1", y = "UMAP2", colour = "Condition") +
  theme_audit() +
  theme(legend.position = "none")

message("Computing local-neighbour mixing and state adjacency")
k_use <- 30L
knn <- BiocNeighbors::findKNN(pca, k = k_use, BNPARAM = BiocNeighbors::KmknnParam())
knn_idx <- knn$index
condition_vec <- as.character(stromal$condition)
state_vec <- as.character(stromal$stromal_state)

same_condition <- rowMeans(matrix(condition_vec[knn_idx], nrow = nrow(knn_idx)) == condition_vec)
same_state <- rowMeans(matrix(state_vec[knn_idx], nrow = nrow(knn_idx)) == state_vec)

mixing_cell <- data.frame(
  cell = colnames(stromal),
  condition = condition_vec,
  stromal_state = state_vec,
  same_condition_neighbor_fraction = same_condition,
  same_state_neighbor_fraction = same_state
)
write.csv(mixing_cell, file.path(result_dir, "joint_stromal_knn_mixing_per_cell.csv"),
          row.names = FALSE)

mixing_summary <- mixing_cell %>%
  group_by(condition, stromal_state) %>%
  summarise(
    n_cells = n(),
    median_same_condition = median(same_condition_neighbor_fraction),
    q25_same_condition = quantile(same_condition_neighbor_fraction, 0.25),
    q75_same_condition = quantile(same_condition_neighbor_fraction, 0.75),
    median_same_state = median(same_state_neighbor_fraction),
    .groups = "drop"
  )
write.csv(mixing_summary, file.path(result_dir, "joint_stromal_knn_mixing_summary.csv"),
          row.names = FALSE)

source_state <- rep(state_vec, times = k_use)
target_state <- state_vec[as.vector(knn_idx)]
transition <- as.data.frame(table(source_state, target_state), stringsAsFactors = FALSE) %>%
  rename(source_state = source_state, target_state = target_state, n_edges = Freq) %>%
  group_by(source_state) %>%
  mutate(edge_fraction = n_edges / sum(n_edges)) %>%
  ungroup()
write.csv(transition, file.path(result_dir, "joint_stromal_state_knn_transition.csv"),
          row.names = FALSE)

p_transition <- ggplot(transition,
                       aes(x = target_state, y = source_state, fill = edge_fraction)) +
  geom_tile(colour = "white", linewidth = 0.2) +
  scale_fill_gradient(low = "#F7FBFF", high = "#2166AC",
                      labels = scales::percent_format(accuracy = 1)) +
  labs(title = "Root-free local adjacency between annotated states",
       x = "Neighbour state", y = "Index state", fill = "kNN edge\nfraction") +
  theme_audit(10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line = element_blank(), axis.ticks = element_blank())

transition_reciprocal <- transition %>%
  select(source_state, target_state, edge_fraction) %>%
  inner_join(
    transition %>%
      select(source_state, target_state, edge_fraction) %>%
      rename(source_state_rev = source_state,
             target_state_rev = target_state,
             edge_fraction_rev = edge_fraction),
    by = c("source_state" = "target_state_rev",
           "target_state" = "source_state_rev")
  ) %>%
  filter(source_state < target_state) %>%
  mutate(
    reciprocal_edge_fraction = pmin(edge_fraction, edge_fraction_rev),
    mean_edge_fraction = (edge_fraction + edge_fraction_rev) / 2
  ) %>%
  arrange(desc(reciprocal_edge_fraction))
write.csv(transition_reciprocal,
          file.path(result_dir, "joint_stromal_state_reciprocal_knn_edges.csv"),
          row.names = FALSE)

message("Estimating root-free centroid-MST edge stability")
calc_centroids <- function(indices_by_state, max_per_state = 300L) {
  do.call(rbind, lapply(names(indices_by_state), function(s) {
    idx <- indices_by_state[[s]]
    idx <- sample(idx, min(length(idx), max_per_state), replace = FALSE)
    colMeans(pca[idx, , drop = FALSE])
  }))
}

indices_by_state <- split(seq_along(state_vec), state_vec)
indices_by_state <- indices_by_state[state_levels]

mst_edges <- function(centroid_matrix) {
  d <- as.matrix(dist(centroid_matrix))
  g <- graph_from_adjacency_matrix(d, mode = "undirected", weighted = TRUE, diag = FALSE)
  m <- mst(g, weights = E(g)$weight)
  ed <- as_edgelist(m, names = TRUE)
  t(apply(ed, 1, sort))
}

n_boot <- 200L
edge_keys <- character()
for (b in seq_len(n_boot)) {
  cent <- calc_centroids(indices_by_state)
  rownames(cent) <- names(indices_by_state)
  ed <- mst_edges(cent)
  edge_keys <- c(edge_keys, apply(ed, 1, paste, collapse = " || "))
}
edge_support <- as.data.frame(table(edge_keys), stringsAsFactors = FALSE) %>%
  rename(edge = edge_keys, bootstrap_count = Freq) %>%
  mutate(bootstrap_support = bootstrap_count / n_boot) %>%
  separate(edge, into = c("state1", "state2"), sep = " \\|\\| ") %>%
  arrange(desc(bootstrap_support))
write.csv(edge_support, file.path(result_dir, "joint_stromal_centroid_mst_edge_stability.csv"),
          row.names = FALSE)

centroid_umap <- umap_df %>%
  group_by(stromal_state) %>%
  summarise(x = median(UMAP1), y = median(UMAP2), .groups = "drop") %>%
  mutate(stromal_state = as.character(stromal_state))

edge_plot <- edge_support %>%
  filter(bootstrap_support >= 0.50) %>%
  left_join(centroid_umap, by = c("state1" = "stromal_state")) %>%
  rename(x1 = x, y1 = y) %>%
  left_join(centroid_umap, by = c("state2" = "stromal_state")) %>%
  rename(x2 = x, y2 = y)

p_mst <- ggplot() +
  geom_point(data = umap_df, aes(UMAP1, UMAP2), colour = "grey82", size = 0.22,
             alpha = 0.35) +
  geom_segment(data = edge_plot,
               aes(x = x1, y = y1, xend = x2, yend = y2,
                   linewidth = bootstrap_support), colour = "black", alpha = 0.8) +
  geom_point(data = centroid_umap,
             aes(x, y, fill = stromal_state), shape = 21, size = 5.2,
             colour = "black", stroke = 0.6) +
  ggrepel::geom_text_repel(data = centroid_umap,
                           aes(x, y, label = stromal_state),
                           fontface = "bold", size = 3.3,
                           box.padding = 0.35, point.padding = 0.35,
                           segment.color = "grey50") +
  scale_fill_manual(values = state_colors, drop = FALSE) +
  scale_linewidth(range = c(0.5, 2.2), limits = c(0.5, 1),
                  labels = scales::percent_format(accuracy = 1)) +
  coord_equal() +
  labs(title = "Undirected state adjacency with bootstrap support",
       x = "UMAP1", y = "UMAP2", linewidth = "Edge support") +
  theme_audit() +
  theme(legend.position = "right")

# The reciprocal kNN network is the interpretable root-free display. Unlike an
# MST, it is allowed to remain disconnected and therefore does not manufacture
# a bridge between transcriptionally separated state systems.
network_edges <- transition_reciprocal %>%
  filter(reciprocal_edge_fraction >= 0.008) %>%
  left_join(centroid_umap, by = c("source_state" = "stromal_state")) %>%
  rename(x1 = x, y1 = y) %>%
  left_join(centroid_umap, by = c("target_state" = "stromal_state")) %>%
  rename(x2 = x, y2 = y)

node_counts <- meta %>%
  count(stromal_state, name = "n_cells") %>%
  mutate(stromal_state = as.character(stromal_state))
network_nodes <- centroid_umap %>% left_join(node_counts, by = "stromal_state")

p_network <- ggplot() +
  geom_point(data = umap_df, aes(UMAP1, UMAP2), colour = "grey86", size = 0.20,
             alpha = 0.30) +
  geom_segment(data = network_edges,
               aes(x = x1, y = y1, xend = x2, yend = y2,
                   linewidth = reciprocal_edge_fraction),
               colour = "black", alpha = 0.80, lineend = "round") +
  geom_point(data = network_nodes,
             aes(x, y, fill = stromal_state, size = n_cells),
             shape = 21, colour = "black", stroke = 0.6) +
  ggrepel::geom_text_repel(data = network_nodes,
                           aes(x, y, label = stromal_state),
                           fontface = "bold", size = 3.3,
                           box.padding = 0.35, point.padding = 0.35,
                           segment.color = "grey50") +
  scale_fill_manual(values = state_colors, drop = FALSE) +
  scale_linewidth(range = c(0.6, 2.4), guide = "none") +
  scale_size_continuous(range = c(3.5, 7.5), breaks = c(500, 1500, 3000)) +
  coord_equal() +
  labs(title = "Root-free reciprocal kNN adjacency (no forced bridge)",
       x = "UMAP1", y = "UMAP2", size = "Cells") +
  theme_audit() +
  theme(legend.position = "right")

message("Learning a partition-aware Monocle3 graph without forcing disconnected regions together")
counts_mat <- GetAssayData(stromal, assay = "RNA", layer = "counts")
gene_metadata <- data.frame(gene_short_name = rownames(counts_mat), row.names = rownames(counts_mat))
cell_metadata <- stromal[[]]
cds <- new_cell_data_set(counts_mat, cell_metadata = cell_metadata,
                         gene_metadata = gene_metadata)
cds <- estimate_size_factors(cds)
reducedDims(cds)$UMAP <- umap
cds <- cluster_cells(cds, reduction_method = "UMAP", cluster_method = "leiden",
                     k = 30, random_seed = 20260911)
cds_partitioned <- learn_graph(cds, use_partition = TRUE, close_loop = FALSE)

partition_vec <- as.character(partitions(cds_partitioned))
names(partition_vec) <- colnames(cds_partitioned)
partition_df <- data.frame(
  cell = colnames(cds_partitioned),
  partition = partition_vec,
  condition = condition_vec[colnames(cds_partitioned) %>% match(colnames(stromal))],
  stromal_state = state_vec[colnames(cds_partitioned) %>% match(colnames(stromal))]
)
write.csv(partition_df, file.path(result_dir, "joint_stromal_monocle_partitions_per_cell.csv"),
          row.names = FALSE)
write.csv(as.data.frame.matrix(table(partition_df$partition, partition_df$condition)),
          file.path(result_dir, "joint_stromal_monocle_partition_by_condition.csv"))
write.csv(as.data.frame.matrix(table(partition_df$partition, partition_df$stromal_state)),
          file.path(result_dir, "joint_stromal_monocle_partition_by_state.csv"))

message("Auditing root sensitivity on one deliberately connected graph")
cds_connected <- learn_graph(cds, use_partition = FALSE, close_loop = FALSE)

candidate_roots <- c(
  "LEPR+ MSC-like",
  "FibroStromal",
  "ChondroFibroStromal",
  "LesionStromal",
  "ECM-remodeling stromal state"
)

root_medoids <- vapply(candidate_roots, function(s) {
  idx <- which(state_vec == s)
  center <- colMeans(pca[idx, , drop = FALSE])
  idx[which.min(rowSums((pca[idx, , drop = FALSE] -
                          matrix(center, nrow = length(idx), ncol = ncol(pca), byrow = TRUE))^2))]
}, integer(1))
root_cells <- colnames(stromal)[root_medoids]
names(root_cells) <- candidate_roots

pt_list <- lapply(candidate_roots, function(s) {
  ordered <- order_cells(cds_connected, root_cells = root_cells[[s]])
  pt <- pseudotime(ordered)
  pt[!is.finite(pt)] <- NA_real_
  rng <- range(pt, na.rm = TRUE)
  if (diff(rng) > 0) pt <- (pt - rng[1]) / diff(rng)
  pt
})
names(pt_list) <- candidate_roots

old_multi_root_cells <- colnames(stromal)[state_vec %in% c("FibroStromal", "ChondroFibroStromal")]
cds_old_root <- order_cells(cds_connected, root_cells = old_multi_root_cells)
pt_old <- pseudotime(cds_old_root)
pt_old[!is.finite(pt_old)] <- NA_real_
old_rng <- range(pt_old, na.rm = TRUE)
if (diff(old_rng) > 0) pt_old <- (pt_old - old_rng[1]) / diff(old_rng)
pt_list[["Old multi-root definition"]] <- pt_old

pt_matrix <- do.call(cbind, pt_list)
write.csv(data.frame(cell = rownames(pt_matrix), pt_matrix, check.names = FALSE),
          file.path(result_dir, "joint_stromal_root_sensitivity_pseudotime.csv"), row.names = FALSE)

root_cor <- cor(pt_matrix, method = "spearman", use = "pairwise.complete.obs")
write.csv(root_cor, file.path(result_dir, "joint_stromal_root_sensitivity_spearman.csv"))

root_cor_df <- as.data.frame(as.table(root_cor), stringsAsFactors = FALSE) %>%
  rename(root_1 = Var1, root_2 = Var2, spearman_rho = Freq)

p_root_cor <- ggplot(root_cor_df, aes(root_2, root_1, fill = spearman_rho)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", spearman_rho)), fontface = "bold", size = 3) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-1, 1)) +
  coord_equal() +
  labs(title = "Pseudotime ordering changes with root choice",
       x = NULL, y = NULL, fill = "Spearman rho") +
  theme_audit(9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line = element_blank(), axis.ticks = element_blank())

pt_long <- data.frame(cell = rownames(pt_matrix), pt_matrix, check.names = FALSE) %>%
  pivot_longer(-cell, names_to = "root_definition", values_to = "scaled_pseudotime") %>%
  left_join(meta, by = "cell")

pt_summary <- pt_long %>%
  group_by(root_definition, condition, stromal_state) %>%
  summarise(
    n_cells = sum(!is.na(scaled_pseudotime)),
    median = median(scaled_pseudotime, na.rm = TRUE),
    q25 = quantile(scaled_pseudotime, 0.25, na.rm = TRUE),
    q75 = quantile(scaled_pseudotime, 0.75, na.rm = TRUE),
    fraction_low_0_1 = mean(scaled_pseudotime <= 0.1, na.rm = TRUE),
    .groups = "drop"
  )
write.csv(pt_summary, file.path(result_dir, "joint_stromal_root_sensitivity_summary.csv"),
          row.names = FALSE)

p_pt_state <- ggplot(pt_long,
                     aes(x = stromal_state, y = scaled_pseudotime, fill = stromal_state)) +
  geom_violin(scale = "width", trim = TRUE, linewidth = 0.25) +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", linewidth = 0.3) +
  facet_wrap(~root_definition, ncol = 2) +
  scale_fill_manual(values = state_colors, drop = FALSE) +
  labs(title = "State ordering is not invariant to the selected root",
       x = NULL, y = "Scaled pseudotime") +
  theme_audit(9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

saveRDS(cds_partitioned,
        file.path(root_dir, "objects", "08_joint_stromal_monocle_partitioned.rds"),
        compress = FALSE)
saveRDS(cds_connected,
        file.path(root_dir, "objects", "08_joint_stromal_monocle_connected_audit.rds"),
        compress = FALSE)

figure_overview <- (p_state | p_condition) / (p_network | p_transition) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 18))
save_audit(figure_overview, "FigureS_stromal_joint_topology_audit", width = 16, height = 11)

figure_root <- p_root_cor / p_pt_state +
  plot_layout(heights = c(0.9, 1.5)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 18))
save_audit(figure_root, "FigureS_stromal_root_sensitivity_audit", width = 14, height = 12)

summary_metrics <- data.frame(
  metric = c(
    "n_stromal_cells",
    "n_states",
    "n_monocle_partitions",
    "median_same_condition_neighbor_fraction",
    "median_same_state_neighbor_fraction",
    "old_multi_root_fraction_pseudotime_le_0_1"
  ),
  value = c(
    ncol(stromal),
    nlevels(stromal$stromal_state),
    length(unique(partition_vec)),
    median(same_condition),
    median(same_state),
    mean(pt_old <= 0.1, na.rm = TRUE)
  )
)
write.csv(summary_metrics, file.path(result_dir, "joint_stromal_topology_audit_summary.csv"),
          row.names = FALSE)

capture.output(sessionInfo(), file = file.path(result_dir, "sessionInfo_trajectory_audit.txt"))
message("Trajectory/topology audit completed")
