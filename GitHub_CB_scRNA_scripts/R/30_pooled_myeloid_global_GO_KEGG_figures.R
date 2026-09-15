root <- Sys.getenv("CB_REVISION_ROOT")
if (!nzchar(root)) root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
config <- yaml::read_yaml(file.path(root, "config", "analysis_config.yml"))
seed <- as.integer(config$seed)
set.seed(seed)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(stringr)
  library(purrr)
  library(ggplot2)
  library(ggrepel)
  library(scales)
  library(svglite)
  library(ragg)
  library(circlize)
})

result_dir <- file.path(root, "results", "enrichment", "pooled_myeloid_pairwise")
figure_dir <- file.path(
  root, "figures", "publication_subfigures",
  "pooled_myeloid_pairwise_global_GO_KEGG"
)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

de_all <- read_csv(file.path(result_dir, "DE_all_pairwise.csv.gz"), show_col_types = FALSE)
gsea_go <- read_csv(file.path(result_dir, "GSEA_GO_all_pairwise.csv.gz"), show_col_types = FALSE)
gsea_kegg <- read_csv(file.path(result_dir, "GSEA_KEGG_all_pairwise.csv.gz"), show_col_types = FALSE)
comparison_table <- read_csv(file.path(result_dir, "comparison_definitions.csv"), show_col_types = FALSE)

state_order <- c("Resident Mac", "TREM2+ Mac", "OC")
state_colors <- c(
  "Resident Mac" = "#79B8A9",
  "TREM2+ Mac" = "#9A78B4",
  "OC" = "#527BA8"
)
family_colors <- c(
  "Osteoclast / acidification" = "#527BA8",
  "Mitochondrial energy" = "#79B8A9",
  "ECM / adhesion / migration" = "#E3A85B",
  "Innate inflammation" = "#9A78B4",
  "Interferon / antiviral" = "#D96C75",
  "Translation" = "#74A6C6",
  "Metabolic remodeling" = "#C58BBE",
  "Other bone-related signaling" = "#B8A6D1"
)

comparison_axis_labels <- c(
  "TREM2_vs_Resident" = "Resident Mac enriched  <----  |  ---->  TREM2+ Mac enriched",
  "OC_vs_Resident" = "Resident Mac enriched  <----  |  ---->  OC enriched",
  "OC_vs_TREM2" = "TREM2+ Mac enriched  <----  |  ---->  OC enriched"
)
comparison_facet_labels <- c(
  "TREM2_vs_Resident" = "TREM2+ Mac vs\nResident Mac",
  "OC_vs_Resident" = "OC vs\nResident Mac",
  "OC_vs_TREM2" = "OC vs\nTREM2+ Mac"
)
comparison_order <- names(comparison_axis_labels)

theme_pub <- theme_classic(base_size = 11.5, base_family = "Arial") +
  theme(
    text = element_text(family = "Arial", face = "bold", colour = "black"),
    axis.title = element_text(size = 11.5, face = "bold"),
    axis.text = element_text(size = 9.4, face = "bold", colour = "black"),
    axis.line = element_line(linewidth = 0.5, colour = "black"),
    axis.ticks = element_line(linewidth = 0.45, colour = "black"),
    strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.55),
    strip.text = element_text(size = 10.4, face = "bold", colour = "black"),
    plot.title = element_text(size = 13.2, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 9.5, face = "bold", colour = "grey25"),
    plot.caption = element_text(size = 8.2, face = "bold", hjust = 0, colour = "grey25"),
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 8.8, face = "bold"),
    panel.grid = element_blank(),
    plot.margin = margin(7, 10, 7, 7)
  )

save_plot <- function(plot, stem, width_mm, height_mm) {
  path <- file.path(figure_dir, stem)
  wi <- width_mm / 25.4
  hi <- height_mm / 25.4
  svglite::svglite(paste0(path, ".svg"), width = wi, height = hi)
  print(plot)
  dev.off()
  grDevices::cairo_pdf(paste0(path, ".pdf"), width = wi, height = hi, family = "Arial")
  print(plot)
  dev.off()
  ragg::agg_tiff(
    paste0(path, ".tiff"), width = wi, height = hi, units = "in",
    res = 600, compression = "lzw", background = "white"
  )
  print(plot)
  dev.off()
  ragg::agg_png(
    paste0(path, ".png"), width = wi, height = hi, units = "in",
    res = 600, background = "white"
  )
  print(plot)
  dev.off()
}

save_base <- function(draw_fun, stem, width_mm, height_mm) {
  path <- file.path(figure_dir, stem)
  wi <- width_mm / 25.4
  hi <- height_mm / 25.4
  svglite::svglite(paste0(path, ".svg"), width = wi, height = hi)
  draw_fun()
  dev.off()
  grDevices::cairo_pdf(paste0(path, ".pdf"), width = wi, height = hi, family = "Arial")
  draw_fun()
  dev.off()
  ragg::agg_tiff(
    paste0(path, ".tiff"), width = wi, height = hi, units = "in",
    res = 600, compression = "lzw", background = "white"
  )
  draw_fun()
  dev.off()
  ragg::agg_png(
    paste0(path, ".png"), width = wi, height = hi, units = "in",
    res = 600, background = "white"
  )
  draw_fun()
  dev.off()
}

parse_genes <- function(x) {
  if (length(x) == 0 || is.na(x) || !nzchar(x)) return(character())
  unique(strsplit(x, "/", fixed = TRUE)[[1]])
}

top_leading <- function(leading_genes, comparison_id, n = 4L) {
  genes <- parse_genes(leading_genes)
  if (!length(genes)) return(character())
  de_all |>
    filter(comparison_id == !!comparison_id, gene %in% genes) |>
    arrange(desc(abs(avg_log2FC)), p_val_adj, gene) |>
    distinct(gene, .keep_all = TRUE) |>
    slice_head(n = n) |>
    pull(gene)
}

# -----------------------------------------------------------------------------
# A. Pairwise differential-expression landscape across all pooled cells.
# -----------------------------------------------------------------------------
manual_markers <- list(
  "TREM2_vs_Resident" = c("TREM2", "MSR1", "APOE", "LCOR", "VSIG4", "C1QA"),
  "OC_vs_Resident" = c("ACP5", "ATP6V0D2", "TCIRG1", "SIGLEC15", "SPP1", "C1QA"),
  "OC_vs_TREM2" = c("ACP5", "ATP6V0D2", "TCIRG1", "SIGLEC15", "TREM2", "TNF", "NFKBIA")
)

de_plot <- de_all |>
  mutate(
    comparison = factor(comparison, levels = comparison_table$comparison_label),
    neglog10_padj = pmin(-log10(pmax(p_val_adj, 1e-300)), 80),
    status = if_else(candidate_state_DEG, higher_in, "Not selected")
  )

label_de <- de_plot |>
  group_by(comparison_id, higher_in) |>
  arrange(p_val_adj, desc(abs(avg_log2FC))) |>
  mutate(
    rank_within_direction = row_number(),
    manual = map2_lgl(comparison_id, gene, ~ .y %in% manual_markers[[.x]])
  ) |>
  filter((candidate_state_DEG & rank_within_direction <= 4) | manual) |>
  ungroup() |>
  distinct(comparison_id, gene, .keep_all = TRUE)

p_de <- ggplot(de_plot, aes(avg_log2FC, neglog10_padj)) +
  geom_point(aes(colour = status), size = 0.85, alpha = 0.72) +
  geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", linewidth = 0.42, colour = "grey35") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.42, colour = "grey35") +
  ggrepel::geom_text_repel(
    data = label_de, aes(label = gene), family = "Arial", fontface = "bold",
    size = 3.05, colour = "black", segment.colour = "grey35", segment.size = 0.35,
    max.overlaps = Inf, box.padding = 0.35, point.padding = 0.15,
    min.segment.length = 0, seed = seed, show.legend = FALSE
  ) +
  facet_wrap(~comparison, nrow = 1, scales = "free") +
  scale_colour_manual(values = c(state_colors, "Not selected" = "#D2D2D2"), drop = FALSE) +
  labs(
    title = "Pairwise transcript differences among pooled myeloid cell states",
    subtitle = "All retained OC, TREM2+ Mac and Resident Mac cells pooled across CB and CTRL",
    x = "Average log2 fold change (first state vs second state)",
    y = expression(-log[10](cell-level~adjusted~P)), colour = "Higher in",
    caption = "Candidate state-associated genes: BH-adjusted P < 0.05, |avg log2FC| >= 0.25 and >=10% detection.\nCell-level P values are descriptive because state and specimen composition are partly confounded."
  ) +
  theme_pub +
  theme(legend.position = "bottom", legend.direction = "horizontal")

save_plot(p_de, "A_global_pairwise_DE_volcano", 245, 112)

for (cid in comparison_order) {
  comp_label <- comparison_table |>
    filter(comparison_id == cid) |>
    pull(comparison_label)
  p_one <- p_de %+% filter(de_plot, comparison_id == cid) +
    facet_null() +
    labs(
      title = comp_label,
      subtitle = "Pooled cell-state differential-expression landscape",
      caption = "Candidate genes: BH FDR < 0.05; |avg log2FC| >= 0.25; detection >= 10%.\nCell-level results are descriptive; state and specimen are partly confounded."
    ) +
    theme(legend.position = "bottom")
  p_one$layers[[4]]$data <- filter(label_de, comparison_id == cid)
  save_plot(p_one, paste0("A_", cid, "_DE_volcano"), 145, 118)
}

# -----------------------------------------------------------------------------
# B. Whole-GO scan: statistically supported, non-redundant representatives.
# These representatives were selected after reviewing the complete GSEA result,
# rather than restricting the analysis to the original eight mechanism terms.
# -----------------------------------------------------------------------------
go_representatives <- tribble(
  ~comparison_id, ~enriched_in, ~ID, ~family, ~display_order,
  "TREM2_vs_Resident", "Resident Mac", "GO:0051607", "Interferon / antiviral", 1,
  "TREM2_vs_Resident", "Resident Mac", "GO:0034340", "Interferon / antiviral", 2,
  "TREM2_vs_Resident", "TREM2+ Mac", "GO:0032543", "Translation", 3,
  "TREM2_vs_Resident", "TREM2+ Mac", "GO:0001755", "ECM / adhesion / migration", 4,
  "TREM2_vs_Resident", "TREM2+ Mac", "GO:0043648", "Metabolic remodeling", 5,
  "TREM2_vs_Resident", "TREM2+ Mac", "GO:0085029", "ECM / adhesion / migration", 6,
  "OC_vs_Resident", "Resident Mac", "GO:0071222", "Innate inflammation", 7,
  "OC_vs_Resident", "Resident Mac", "GO:0001819", "Innate inflammation", 8,
  "OC_vs_Resident", "Resident Mac", "GO:0034340", "Interferon / antiviral", 9,
  "OC_vs_Resident", "Resident Mac", "GO:0007249", "Innate inflammation", 10,
  "OC_vs_Resident", "OC", "GO:0006119", "Mitochondrial energy", 11,
  "OC_vs_Resident", "OC", "GO:0032543", "Translation", 12,
  "OC_vs_Resident", "OC", "GO:0030198", "ECM / adhesion / migration", 13,
  "OC_vs_Resident", "OC", "GO:0036035", "Osteoclast / acidification", 14,
  "OC_vs_Resident", "OC", "GO:0007042", "Osteoclast / acidification", 15,
  "OC_vs_TREM2", "TREM2+ Mac", "GO:0001819", "Innate inflammation", 16,
  "OC_vs_TREM2", "TREM2+ Mac", "GO:0071222", "Innate inflammation", 17,
  "OC_vs_TREM2", "TREM2+ Mac", "GO:0007249", "Innate inflammation", 18,
  "OC_vs_TREM2", "OC", "GO:0006119", "Mitochondrial energy", 19,
  "OC_vs_TREM2", "OC", "GO:0032543", "Translation", 20,
  "OC_vs_TREM2", "OC", "GO:0030198", "ECM / adhesion / migration", 21
)

go_plot_data <- gsea_go |>
  inner_join(go_representatives, by = c("comparison_id", "enriched_in", "ID")) |>
  mutate(
    comparison_axis = factor(comparison_axis_labels[comparison_id], levels = comparison_axis_labels[comparison_order]),
    comparison_facet = factor(comparison_facet_labels[comparison_id], levels = comparison_facet_labels[comparison_order]),
    signed_GeneRatio = if_else(enriched_in == ident_1, GeneRatio_numeric, -GeneRatio_numeric),
    neglog10_FDR = -log10(pmax(p.adjust, 1e-300)),
    top_genes = map2(leading_genes, comparison_id, top_leading, n = 3L),
    top_gene_label = map_chr(top_genes, ~ paste(.x, collapse = ", ")),
    display_label = paste0(str_wrap(str_to_sentence(Description), 34), "\n", top_gene_label),
    label_key = paste(comparison_id, display_order, display_label, sep = "|||"),
    label_key = factor(label_key, levels = rev(unique(label_key[order(display_order)]))),
    family = factor(family, levels = names(family_colors))
  ) |>
  arrange(match(comparison_id, comparison_order), display_order)

write_csv(go_plot_data, file.path(result_dir, "global_GO_representatives_figure_source.csv"))

xmax_go <- max(abs(go_plot_data$signed_GeneRatio), na.rm = TRUE) * 1.12
p_go_global <- ggplot(go_plot_data, aes(signed_GeneRatio, label_key)) +
  geom_vline(xintercept = 0, linewidth = 0.48, colour = "black") +
  geom_point(
    aes(size = Count, fill = neglog10_FDR, colour = family),
    shape = 21, stroke = 1.0, alpha = 0.98
  ) +
  facet_grid(rows = vars(comparison_facet), scales = "free_y", space = "free_y", drop = TRUE) +
  scale_y_discrete(labels = function(x) sub("^.*\\|\\|\\|", "", x)) +
  scale_x_continuous(
    limits = c(-xmax_go, xmax_go), labels = function(z) sprintf("%.2f", z),
    expand = expansion(mult = c(0.04, 0.04))
  ) +
  scale_size_continuous(range = c(3.5, 10.5), breaks = pretty_breaks(n = 4), name = "Leading genes\n(Count)") +
  scale_fill_gradientn(
    colours = c("#F7F7F7", "#FDE0C5", "#F4A582", "#D6604D", "#B2182B"),
    values = scales::rescale(c(0, 1.3, 3, 6, max(go_plot_data$neglog10_FDR, na.rm = TRUE))),
    oob = scales::squish, name = expression(-log[10](FDR))
  ) +
  scale_colour_manual(values = family_colors, name = "Process family", drop = TRUE) +
  labs(
    title = "Global GO landscape identifies both expected and additional cell-state programs",
    subtitle = "Representatives were selected from the complete GSEA output after significance filtering and redundancy review",
    x = "GeneRatio (leading-edge genes / pathway genes)", y = NULL,
    caption = "Negative values indicate enrichment in the second-named state and positive values the first-named state.\nGene symbols below each term are its strongest leading-edge drivers; complete GO results remain in the supplementary table."
  ) +
  theme_pub +
  theme(
    axis.text.y = element_text(size = 8.6, face = "bold", lineheight = 0.95),
    strip.text = element_text(size = 10.2),
    strip.text.y.right = element_text(angle = 0, size = 8.6, lineheight = 0.95),
    legend.position = "right",
    panel.spacing.y = unit(3.5, "mm")
  )

save_plot(p_go_global, "B_global_GO_nonredundant_bubble", 282, 238)

for (cid in comparison_order) {
  comp_label <- comparison_table |>
    filter(comparison_id == cid) |>
    pull(comparison_label)
  go_one <- go_plot_data |>
    filter(comparison_id == cid) |>
    mutate(label_key = droplevels(label_key))
  p_one <- p_go_global %+% go_one +
    facet_null() +
    labs(
      title = paste0(comp_label, ": global GO processes"),
      subtitle = NULL,
      caption = NULL
    ) +
    theme(
      legend.position = "right",
      plot.title = element_text(size = 13.2, face = "bold", margin = margin(b = 7)),
      plot.margin = margin(8, 8, 8, 8)
    )
  save_plot(p_one, paste0("B_", cid, "_global_GO_bubble"), 218, 150)
}

# -----------------------------------------------------------------------------
# C-E. KEGG pathway-leading-gene chord diagrams from the complete KEGG scan.
# Disease-labelled terms are not used as primary sectors; shared underlying
# metabolic, mitochondrial, ECM and inflammatory mechanisms are shown directly.
# -----------------------------------------------------------------------------
kegg_representatives <- tribble(
  ~comparison_id, ~enriched_in, ~ID, ~family, ~display_order,
  "TREM2_vs_Resident", "Resident Mac", "hsa04928", "Other bone-related signaling", 1,
  "TREM2_vs_Resident", "TREM2+ Mac", "hsa04512", "ECM / adhesion / migration", 2,
  "TREM2_vs_Resident", "TREM2+ Mac", "hsa04820", "ECM / adhesion / migration", 3,
  "TREM2_vs_Resident", "TREM2+ Mac", "hsa00010", "Metabolic remodeling", 4,
  "TREM2_vs_Resident", "TREM2+ Mac", "hsa04518", "ECM / adhesion / migration", 5,
  "TREM2_vs_Resident", "TREM2+ Mac", "hsa00330", "Metabolic remodeling", 6,
  "OC_vs_Resident", "Resident Mac", "hsa04380", "Osteoclast / acidification", 1,
  "OC_vs_Resident", "Resident Mac", "hsa04064", "Innate inflammation", 2,
  "OC_vs_Resident", "Resident Mac", "hsa04060", "Innate inflammation", 3,
  "OC_vs_Resident", "Resident Mac", "hsa04148", "Innate inflammation", 4,
  "OC_vs_Resident", "OC", "hsa00190", "Mitochondrial energy", 5,
  "OC_vs_Resident", "OC", "hsa04512", "ECM / adhesion / migration", 6,
  "OC_vs_Resident", "OC", "hsa04966", "Osteoclast / acidification", 7,
  "OC_vs_TREM2", "TREM2+ Mac", "hsa04060", "Innate inflammation", 1,
  "OC_vs_TREM2", "TREM2+ Mac", "hsa04064", "Innate inflammation", 2,
  "OC_vs_TREM2", "TREM2+ Mac", "hsa04668", "Innate inflammation", 3,
  "OC_vs_TREM2", "TREM2+ Mac", "hsa04148", "Innate inflammation", 4,
  "OC_vs_TREM2", "OC", "hsa00190", "Mitochondrial energy", 5,
  "OC_vs_TREM2", "OC", "hsa04512", "ECM / adhesion / migration", 6,
  "OC_vs_TREM2", "OC", "hsa04966", "Osteoclast / acidification", 7,
  "OC_vs_TREM2", "OC", "hsa01200", "Metabolic remodeling", 8
)

kegg_selected <- gsea_kegg |>
  inner_join(kegg_representatives, by = c("comparison_id", "enriched_in", "ID")) |>
  mutate(
    short_pathway = case_when(
      ID == "hsa04928" ~ "Parathyroid hormone signaling",
      ID == "hsa04512" ~ "ECM-receptor interaction",
      ID == "hsa04820" ~ "Cytoskeletal remodeling",
      ID == "hsa00010" ~ "Glycolysis / gluconeogenesis",
      ID == "hsa04518" ~ "Integrin signaling",
      ID == "hsa00330" ~ "Arginine / proline metabolism",
      ID == "hsa04380" ~ "Osteoclast differentiation",
      ID == "hsa04064" ~ "NF-kappa B signaling",
      ID == "hsa04060" ~ "Cytokine-receptor interaction",
      ID == "hsa04148" ~ "Efferocytosis",
      ID == "hsa00190" ~ "Oxidative phosphorylation",
      ID == "hsa04966" ~ "V-ATPase acidification machinery",
      ID == "hsa04668" ~ "TNF signaling",
      ID == "hsa01200" ~ "Carbon metabolism",
      TRUE ~ Description
    ),
    sector_label = paste0(short_pathway, "\n[", enriched_in, "]"),
    family = factor(family, levels = names(family_colors))
  )

make_chord_data <- function(cid, genes_per_pathway = 4L) {
  terms <- kegg_selected |>
    filter(comparison_id == cid) |>
    arrange(display_order)
  fc <- de_all |>
    filter(comparison_id == cid) |>
    select(gene, avg_log2FC, p_val_adj)

  edges <- terms |>
    mutate(gene = map2(leading_genes, comparison_id, top_leading, n = genes_per_pathway)) |>
    select(comparison_id, enriched_in, ID, Description, short_pathway, sector_label,
           family, display_order, NES, p.adjust, Count, GeneRatio_numeric, gene) |>
    unnest(gene) |>
    left_join(fc, by = "gene") |>
    mutate(value = 1)

  list(terms = terms, edges = edges)
}

draw_chord <- function(cid) {
  dat <- make_chord_data(cid)
  comp <- comparison_table |> filter(comparison_id == cid)
  edges <- dat$edges
  terms <- dat$terms
  pathway_keys <- unique(edges$sector_label)
  gene_keys <- unique(edges$gene)
  sector_order <- c(pathway_keys, gene_keys)

  path_family <- setNames(as.character(terms$family), terms$sector_label)
  gene_fc <- edges |>
    group_by(gene) |>
    summarise(avg_log2FC = first(avg_log2FC), .groups = "drop")
  gene_cols <- ifelse(
    gene_fc$avg_log2FC >= 0,
    state_colors[[comp$ident_1]],
    state_colors[[comp$ident_2]]
  )
  names(gene_cols) <- gene_fc$gene
  grid_cols <- c(family_colors[path_family[pathway_keys]], gene_cols[gene_keys])
  names(grid_cols) <- sector_order
  link_cols <- alpha(family_colors[as.character(edges$family)], 0.62)

  function() {
    # Keep an identical, centred chord field in all three standalone panels.
    # A fixed lower margin accommodates two horizontal legend rows so that the
    # panels can be aligned side by side without a right-hand legend gutter.
    par(family = "Arial", mar = c(9.2, 1.2, 3.8, 1.2), xpd = NA)
    circos.clear()
    circos.par(
      start.degree = 88,
      gap.after = c(rep(4, length(pathway_keys) - 1), 13,
                    rep(1.6, max(length(gene_keys) - 1, 0)), 13),
      track.margin = c(0.006, 0.006),
      cell.padding = c(0, 0, 0, 0),
      canvas.xlim = c(-1.42, 1.42),
      canvas.ylim = c(-1.28, 1.28)
    )
    chordDiagram(
      x = edges |> transmute(from = sector_label, to = gene, value),
      order = sector_order,
      grid.col = grid_cols,
      col = link_cols,
      transparency = 0.20,
      directional = 0,
      annotationTrack = "grid",
      preAllocateTracks = list(track.height = 0.18),
      link.sort = TRUE,
      link.largest.ontop = TRUE
    )
    circos.trackPlotRegion(
      track.index = 1, bg.border = NA,
      panel.fun = function(x, y) {
        sector <- get.cell.meta.data("sector.index")
        xlim <- get.cell.meta.data("xlim")
        ylim <- get.cell.meta.data("ylim")
        is_path <- sector %in% pathway_keys
        lab <- if (is_path) str_wrap(sector, 25) else sector
        circos.text(
          mean(xlim), ylim[1] + 0.08,
          labels = lab,
          facing = "clockwise", niceFacing = TRUE,
          adj = c(0, 0.5),
          cex = if (is_path) 0.58 else 0.66,
          font = if (is_path) 2 else 3,
          col = "black"
        )
      }
    )
    title(
      main = paste0(comp$comparison_label, ": KEGG pathways and leading genes"),
      family = "Arial", font.main = 2, cex.main = 1.16,
      line = 1.7
    )
    shown_families <- names(family_colors)[
      names(family_colors) %in% unique(as.character(terms$family))
    ]
    legend(
      x = -1.30, y = -1.37,
      legend = shown_families,
      fill = family_colors[shown_families],
      title = "Pathway family", bty = "n", cex = 0.68,
      text.font = 2, title.adj = 0, x.intersp = 0.55,
      y.intersp = 0.9, ncol = 1
    )
    legend(
      x = 0.42, y = -1.37,
      legend = c(comp$ident_1, comp$ident_2),
      fill = state_colors[c(comp$ident_1, comp$ident_2)],
      title = "Gene higher in", bty = "n", cex = 0.72,
      text.font = 2, title.adj = 0, x.intersp = 0.6,
      y.intersp = 0.9, ncol = 1
    )
    circos.clear()
  }
}

write_csv(kegg_selected, file.path(result_dir, "global_KEGG_representatives_figure_source.csv"))

chord_specs <- tribble(
  ~comparison_id, ~stem,
  "TREM2_vs_Resident", "C_KEGG_chord_TREM2_vs_Resident",
  "OC_vs_Resident", "D_KEGG_chord_OC_vs_Resident",
  "OC_vs_TREM2", "E_KEGG_chord_OC_vs_TREM2"
)

for (i in seq_len(nrow(chord_specs))) {
  cid <- chord_specs$comparison_id[i]
  write_csv(
    make_chord_data(cid)$edges,
    file.path(result_dir, paste0("KEGG_chord_edges_", cid, ".csv"))
  )
  save_base(draw_chord(cid), chord_specs$stem[i], 218, 190)
}

manifest <- tribble(
  ~panel, ~stem, ~width_mm, ~height_mm, ~role,
  "A0", "A_global_pairwise_DE_volcano", 245, 112, "Combined all-gene pairwise differential-expression landscape",
  "A1", "A_TREM2_vs_Resident_DE_volcano", 145, 118, "Standalone TREM2+ Mac versus Resident Mac volcano",
  "A2", "A_OC_vs_Resident_DE_volcano", 145, 118, "Standalone OC versus Resident Mac volcano",
  "A3", "A_OC_vs_TREM2_DE_volcano", 145, 118, "Standalone OC versus TREM2+ Mac volcano",
  "B0", "B_global_GO_nonredundant_bubble", 282, 238, "Combined whole-GO scan with non-redundant representatives and leading genes",
  "B1", "B_TREM2_vs_Resident_global_GO_bubble", 218, 150, "Standalone TREM2+ Mac versus Resident Mac GO bubble",
  "B2", "B_OC_vs_Resident_global_GO_bubble", 218, 150, "Standalone OC versus Resident Mac GO bubble",
  "B3", "B_OC_vs_TREM2_global_GO_bubble", 218, 150, "Standalone OC versus TREM2+ Mac GO bubble",
  "C", "C_KEGG_chord_TREM2_vs_Resident", 218, 190, "KEGG pathway-leading-gene chord: TREM2+ Mac vs Resident Mac",
  "D", "D_KEGG_chord_OC_vs_Resident", 218, 190, "KEGG pathway-leading-gene chord: OC vs Resident Mac",
  "E", "E_KEGG_chord_OC_vs_TREM2", 218, 190, "KEGG pathway-leading-gene chord: OC vs TREM2+ Mac"
)
write_csv(manifest, file.path(result_dir, "global_scan_figure_manifest.csv"))

writeLines(capture.output(sessionInfo()), file.path(result_dir, "sessionInfo_stage30.txt"))
message("Global pooled-state figures written to: ", normalizePath(figure_dir, winslash = "/"))
