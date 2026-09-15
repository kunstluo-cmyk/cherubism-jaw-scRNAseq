root <- Sys.getenv("CB_REVISION_ROOT")
if (!nzchar(root)) root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
config <- yaml::read_yaml(file.path(root, "config", "analysis_config.yml"))
seed <- as.integer(config$seed)
set.seed(seed)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
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
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

select <- dplyr::select

result_dir <- file.path(root, "results", "enrichment", "pooled_myeloid_pairwise")
figure_dir <- file.path(
  root, "figures", "publication_subfigures", "pooled_myeloid_pairwise_GO_KEGG_revision"
)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

input_rds <- file.path(result_dir, "pooled_myeloid_pairwise_enrichment_results.rds")
if (!file.exists(input_rds)) stop("Run stage 28 first: missing ", input_rds)
res <- readRDS(input_rds)

de_all <- bind_rows(res$de)
gsea_go <- bind_rows(map(res$gsea, "GO"))
gsea_kegg <- bind_rows(map(res$gsea, "KEGG"))
comparison_table <- res$comparison_table

state_order <- c("Resident Mac", "TREM2+ Mac", "OC")
state_colors <- c(
  "Resident Mac" = "#79B8A9",
  "TREM2+ Mac" = "#9A78B4",
  "OC" = "#527BA8"
)
process_colors <- c(
  "Osteoclast differentiation / resorption" = "#527BA8",
  "Lysosome / phagosome / lipid handling" = "#DDA85C",
  "Inflammatory / NF-kB / TNF" = "#9A78B4",
  "Mitochondrial respiration / OXPHOS" = "#79B8A9"
)
process_shapes <- c(
  "Osteoclast differentiation / resorption" = 21,
  "Lysosome / phagosome / lipid handling" = 22,
  "Inflammatory / NF-kB / TNF" = 24,
  "Mitochondrial respiration / OXPHOS" = 23
)

theme_pub <- theme_classic(base_size = 11.5, base_family = "Arial") +
  theme(
    text = element_text(family = "Arial", face = "bold", colour = "black"),
    axis.title = element_text(size = 11.5, face = "bold"),
    axis.text = element_text(size = 9.5, face = "bold", colour = "black"),
    axis.line = element_line(linewidth = 0.5, colour = "black"),
    axis.ticks = element_line(linewidth = 0.45, colour = "black"),
    strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.55),
    strip.text = element_text(size = 10.5, face = "bold", colour = "black"),
    plot.title = element_text(size = 13, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 9.4, face = "bold", colour = "grey25"),
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

save_heatmap <- function(draw_fun, stem, width_mm, height_mm) {
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

comparison_axis_labels <- c(
  "TREM2_vs_Resident" = "Resident Mac enriched  <----  |  ---->  TREM2+ Mac enriched",
  "OC_vs_Resident" = "Resident Mac enriched  <----  |  ---->  OC enriched",
  "OC_vs_TREM2" = "TREM2+ Mac enriched  <----  |  ---->  OC enriched"
)
comparison_order <- names(comparison_axis_labels)

selected_go <- tribble(
  ~ID, ~display_term, ~process_group, ~display_order,
  "GO:0030316", "Osteoclast differentiation", "Osteoclast differentiation / resorption", 1,
  "GO:0045453", "Bone resorption", "Osteoclast differentiation / resorption", 2,
  "GO:0007042", "Lysosomal lumen acidification", "Lysosome / phagosome / lipid handling", 3,
  "GO:0007041", "Lysosomal transport", "Lysosome / phagosome / lipid handling", 4,
  "GO:0071396", "Cellular response to lipid", "Lysosome / phagosome / lipid handling", 5,
  "GO:0007249", "Canonical NF-kB signaling", "Inflammatory / NF-kB / TNF", 6,
  "GO:0034612", "Response to tumor necrosis factor", "Inflammatory / NF-kB / TNF", 7,
  "GO:0006119", "Oxidative phosphorylation", "Mitochondrial respiration / OXPHOS", 8
)

selected_kegg <- tribble(
  ~ID, ~display_term, ~process_group, ~display_order,
  "hsa04380", "Osteoclast differentiation", "Osteoclast differentiation / resorption", 1,
  "hsa04142", "Lysosome biogenesis", "Lysosome / phagosome / lipid handling", 2,
  "hsa04666", "Fc gamma R-mediated\nphagosome formation", "Lysosome / phagosome / lipid handling", 3,
  "hsa00564", "Glycerophospholipid metabolism", "Lysosome / phagosome / lipid handling", 4,
  "hsa04064", "NF-kappa B signaling", "Inflammatory / NF-kB / TNF", 5,
  "hsa04668", "TNF signaling", "Inflammatory / NF-kB / TNF", 6,
  "hsa04620", "Toll-like receptor signaling", "Inflammatory / NF-kB / TNF", 7,
  "hsa00190", "Oxidative phosphorylation", "Mitochondrial respiration / OXPHOS", 8
)

parse_genes <- function(x) {
  if (length(x) == 0 || is.na(x) || !nzchar(x)) return(character())
  unique(strsplit(x, "/", fixed = TRUE)[[1]])
}

top_leading_label <- function(leading_genes, comparison_id, n = 8) {
  genes <- parse_genes(leading_genes)
  if (!length(genes)) return("")
  fc <- de_all |>
    filter(comparison_id == !!comparison_id, gene %in% genes) |>
    arrange(desc(abs(avg_log2FC)), p_val_adj, gene)
  paste(head(fc$gene, n), collapse = ", ")
}

prepare_selected <- function(x, selection, database_label) {
  out <- x |>
    inner_join(selection, by = "ID") |>
    mutate(
      database_display = database_label,
      comparison_axis = factor(comparison_axis_labels[comparison_id], levels = comparison_axis_labels[comparison_order]),
      signed_GeneRatio = if_else(enriched_in == ident_1, GeneRatio_numeric, -GeneRatio_numeric),
      neglog10_FDR = -log10(pmax(p.adjust, 1e-300)),
      FDR_lt_0_05 = p.adjust < 0.05,
      top_leading_genes = map2_chr(leading_genes, comparison_id, top_leading_label)
    ) |>
    arrange(match(comparison_id, comparison_order), display_order)
  out$display_term <- factor(out$display_term, levels = rev(selection$display_term))
  out$process_group <- factor(out$process_group, levels = names(process_colors))
  out
}

go_selected <- prepare_selected(gsea_go, selected_go, "GO Biological Process")
kegg_selected <- prepare_selected(gsea_kegg, selected_kegg, "KEGG")

write_csv(go_selected, file.path(result_dir, "selected_GO_pathways_figure_source.csv"))
write_csv(kegg_selected, file.path(result_dir, "selected_KEGG_pathways_figure_source.csv"))

plot_enrichment <- function(df, title_text, subtitle_text) {
  xmax <- max(abs(df$signed_GeneRatio), na.rm = TRUE) * 1.12
  p <- ggplot(df, aes(signed_GeneRatio, display_term)) +
    geom_vline(xintercept = 0, linewidth = 0.45, colour = "black") +
    geom_point(
      aes(size = Count, fill = neglog10_FDR, shape = process_group),
      colour = "grey35", stroke = 0.55, alpha = 0.98
    ) +
    geom_point(
      data = df |> filter(FDR_lt_0_05),
      aes(size = Count, shape = process_group),
      fill = NA, colour = "black", stroke = 1.05, show.legend = FALSE
    ) +
    facet_wrap(vars(comparison_axis), ncol = 1, scales = "free_y") +
    scale_x_continuous(
      limits = c(-xmax, xmax), labels = function(z) sprintf("%.2f", abs(z)),
      expand = expansion(mult = c(0.04, 0.04))
    ) +
    scale_size_continuous(range = c(3.3, 10), breaks = pretty_breaks(n = 4), name = "Leading genes\n(Count)") +
    scale_fill_gradientn(
      colours = c("#F7F7F7", "#FDE0C5", "#F4A582", "#D6604D", "#B2182B"),
      values = scales::rescale(c(0, 1.3, 3, 6, max(df$neglog10_FDR, na.rm = TRUE))),
      oob = scales::squish,
      name = expression(-log[10](FDR))
    ) +
    scale_shape_manual(values = process_shapes, name = "Biological process") +
    guides(
      fill = guide_colorbar(order = 1, barheight = unit(28, "mm")),
      size = guide_legend(order = 2),
      shape = guide_legend(order = 3, override.aes = list(size = 4.5, fill = "white"))
    ) +
    labs(
      title = title_text,
      subtitle = subtitle_text,
      x = "GeneRatio (leading-edge genes / pathway genes)", y = NULL,
      caption = "Point area: leading-gene count; fill: BH-adjusted P; black outer ring: FDR < 0.05.\nArrows in facet strips show enrichment direction."
    ) +
    theme_pub +
    theme(
      axis.text.y = element_text(size = 9.2, face = "bold"),
      strip.text = element_text(size = 10.2),
      legend.position = "right",
      panel.spacing.y = unit(3.5, "mm")
    )
  p
}

p_go <- plot_enrichment(
  go_selected,
  "Non-redundant GO processes and their driving genes",
  "Eight representative biological processes selected to span four shared mechanisms"
)
p_kegg <- plot_enrichment(
  kegg_selected,
  "Non-redundant KEGG pathways and their driving genes",
  "Eight representative cellular mechanisms; disease-named pathway labels are excluded from the primary display"
)

save_plot(p_go, "B_GO_nonredundant_leading_gene_enrichment", 245, 188)
save_plot(p_kegg, "C_KEGG_nonredundant_leading_gene_enrichment", 245, 188)

# Pairwise cell-state transcript effects. Adjusted P values are cell-level and exploratory.
de_plot <- de_all |>
  mutate(
    neglog10_padj = pmin(-log10(pmax(p_val_adj, 1e-300)), 80),
    status = if_else(candidate_state_DEG, higher_in, "Not selected"),
    comparison = factor(comparison, levels = comparison_table$comparison_label)
  )

manual_markers <- list(
  "TREM2_vs_Resident" = c("LCOR", "TREM2", "MSR1", "VSIG4"),
  "OC_vs_Resident" = c("ACP5", "SPP1", "ATP6V0D2", "TCIRG1", "MSR1", "VSIG4"),
  "OC_vs_TREM2" = c("ACP5", "ATP6V0D2", "TCIRG1", "SIGLEC15", "TNF", "NFKBIA", "CXCL8", "TREM2")
)
label_de <- de_plot |>
  group_by(comparison_id, higher_in) |>
  arrange(p_val_adj, desc(abs(avg_log2FC))) |>
  mutate(is_manual = map2_lgl(comparison_id, gene, ~ .y %in% manual_markers[[.x]])) |>
  filter((candidate_state_DEG & row_number() <= 4) | is_manual) |>
  ungroup() |>
  distinct(comparison_id, gene, .keep_all = TRUE)

status_colors <- c(state_colors, "Not selected" = "#D2D2D2")
p_de <- ggplot(de_plot, aes(avg_log2FC, neglog10_padj)) +
  geom_point(aes(colour = status), size = 0.85, alpha = 0.72) +
  geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", linewidth = 0.42, colour = "grey35") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.42, colour = "grey35") +
  ggrepel::geom_text_repel(
    data = label_de, aes(label = gene), family = "Arial", fontface = "bold",
    size = 3.1, colour = "black", segment.colour = "grey35", segment.size = 0.35,
    max.overlaps = Inf, box.padding = 0.35, point.padding = 0.15, min.segment.length = 0,
    seed = seed, show.legend = FALSE
  ) +
  facet_wrap(~comparison, nrow = 1, scales = "free") +
  scale_colour_manual(values = status_colors, breaks = c(state_order, "Not selected"), drop = FALSE) +
  scale_x_continuous(expand = expansion(mult = c(0.08, 0.08))) +
  labs(
    title = "Pairwise transcript differences among pooled myeloid cell states",
    subtitle = "All retained OC, TREM2+ Mac and Resident Mac cells; candidate genes require BH-adjusted P < 0.05, |avg log2FC| >= 0.25 and >=10% detection",
    x = "Average log2 fold change (first state vs second state)",
    y = expression(-log[10](cell-level~adjusted~P)), colour = "Higher in",
    caption = "Cells from CB and CTRL are pooled without condition stratification.\nCell-state and specimen effects are partly confounded; cell-level P values do not provide donor-level inference."
  ) +
  theme_pub +
  theme(legend.position = "bottom", legend.direction = "horizontal")
save_plot(p_de, "A_pooled_myeloid_pairwise_DE_volcano", 245, 112)

# Eight concept rows for driver-gene matrices: no duplicated disease labels.
driver_selection <- tribble(
  ~database, ~ID, ~row_label, ~process_group, ~row_order,
  "GO", "GO:0045453", "Bone resorption [GO]", "Osteoclast differentiation / resorption", 1,
  "KEGG", "hsa04380", "Osteoclast differentiation [KEGG]", "Osteoclast differentiation / resorption", 2,
  "GO", "GO:0007042", "Lysosomal lumen acidification [GO]", "Lysosome / phagosome / lipid handling", 3,
  "KEGG", "hsa04666", "Fc gamma R-mediated phagosome [KEGG]", "Lysosome / phagosome / lipid handling", 4,
  "GO", "GO:0071396", "Cellular response to lipid [GO]", "Lysosome / phagosome / lipid handling", 5,
  "GO", "GO:0007249", "Canonical NF-kB signaling [GO]", "Inflammatory / NF-kB / TNF", 6,
  "KEGG", "hsa04668", "TNF signaling [KEGG]", "Inflammatory / NF-kB / TNF", 7,
  "KEGG", "hsa00190", "Oxidative phosphorylation [KEGG]", "Mitochondrial respiration / OXPHOS", 8
)

driver_terms <- bind_rows(
  gsea_go |> mutate(database = "GO"),
  gsea_kegg |> mutate(database = "KEGG")
) |>
  inner_join(driver_selection, by = c("database", "ID")) |>
  mutate(row_label = factor(row_label, levels = driver_selection$row_label))

make_driver_matrix_data <- function(cid, genes_per_term = 4L, max_genes = 30L) {
  terms <- driver_terms |> filter(comparison_id == cid) |> arrange(row_order)
  fc <- de_all |>
    filter(comparison_id == cid) |>
    select(gene, avg_log2FC, p_val_adj, pct.1, pct.2)

  membership <- terms |>
    transmute(ID, row_label = as.character(row_label), process_group, row_order, leading_genes) |>
    mutate(gene = map(leading_genes, parse_genes)) |>
    unnest(gene) |>
    left_join(fc, by = "gene") |>
    filter(!is.na(avg_log2FC))

  recurrence <- membership |>
    distinct(row_label, gene) |>
    count(gene, name = "n_concepts")

  selected_genes <- membership |>
    left_join(recurrence, by = "gene") |>
    group_by(row_label) |>
    arrange(desc(n_concepts), desc(abs(avg_log2FC)), p_val_adj, gene) |>
    slice_head(n = genes_per_term) |>
    ungroup() |>
    distinct(gene) |>
    left_join(
      membership |>
        group_by(gene) |>
        summarise(
          first_row = min(row_order),
          n_concepts = n_distinct(row_label),
          abs_fc = max(abs(avg_log2FC), na.rm = TRUE),
          .groups = "drop"
        ),
      by = "gene"
    ) |>
    arrange(first_row, desc(n_concepts), desc(abs_fc), gene) |>
    slice_head(n = max_genes) |>
    pull(gene)

  plot_long <- expand_grid(
    row_label = driver_selection$row_label,
    gene = selected_genes
  ) |>
    left_join(
      membership |>
        distinct(row_label, gene, avg_log2FC, process_group, row_order),
      by = c("row_label", "gene")
    )

  list(terms = terms, membership = membership, selected_genes = selected_genes, plot_long = plot_long)
}

save_driver_matrix <- function(cid, stem) {
  dat <- make_driver_matrix_data(cid)
  comp <- comparison_table |> filter(comparison_id == cid)
  genes <- dat$selected_genes
  rows <- driver_selection$row_label

  mat <- matrix(NA_real_, nrow = length(rows), ncol = length(genes), dimnames = list(rows, genes))
  for (i in seq_len(nrow(dat$membership))) {
    rr <- as.character(dat$membership$row_label[i])
    gg <- dat$membership$gene[i]
    if (rr %in% rows && gg %in% genes) mat[rr, gg] <- dat$membership$avg_log2FC[i]
  }
  mat <- pmax(pmin(mat, 4), -4)

  higher_state <- ifelse(
    colMeans(mat, na.rm = TRUE) >= 0,
    comp$ident_1,
    comp$ident_2
  )
  higher_state[!is.finite(colMeans(mat, na.rm = TRUE))] <- comp$ident_1

  row_process <- setNames(driver_selection$process_group, driver_selection$row_label)
  row_split <- factor(row_process[rownames(mat)], levels = names(process_colors))

  col_fun <- circlize::colorRamp2(
    c(-3, 0, 3), c(state_colors[[comp$ident_2]], "#F7F7F7", state_colors[[comp$ident_1]])
  )
  top_anno <- HeatmapAnnotation(
    `Higher in` = factor(higher_state, levels = state_order),
    col = list(`Higher in` = state_colors),
    show_annotation_name = TRUE,
    annotation_name_gp = gpar(fontfamily = "Arial", fontface = "bold", fontsize = 9),
    annotation_legend_param = list(
      `Higher in` = list(
        title_gp = gpar(fontfamily = "Arial", fontface = "bold", fontsize = 10),
        labels_gp = gpar(fontfamily = "Arial", fontface = "bold", fontsize = 8.5)
      )
    )
  )
  left_anno <- rowAnnotation(
    Process = row_split,
    col = list(Process = process_colors),
    show_annotation_name = FALSE,
    annotation_legend_param = list(
      Process = list(
        title_gp = gpar(fontfamily = "Arial", fontface = "bold", fontsize = 10),
        labels_gp = gpar(fontfamily = "Arial", fontface = "bold", fontsize = 8.3)
      )
    ),
    width = unit(4.2, "mm")
  )

  draw_fun <- function() {
    ht <- Heatmap(
      mat,
      name = "avg log2FC",
      col = col_fun,
      na_col = "white",
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      top_annotation = top_anno,
      left_annotation = left_anno,
      show_row_names = TRUE,
      show_column_names = TRUE,
      row_names_side = "right",
      row_names_gp = gpar(fontfamily = "Arial", fontface = "bold", fontsize = 8.7),
      column_names_rot = 48,
      column_names_gp = gpar(fontfamily = "Arial", fontface = "bold.italic", fontsize = 8.3),
      rect_gp = gpar(col = "#E1E1E1", lwd = 0.45),
      width = unit(length(genes) * 6.2, "mm"),
      column_title = comp$comparison_label,
      column_title_gp = gpar(fontfamily = "Arial", fontface = "bold", fontsize = 11.5),
      heatmap_legend_param = list(
        title_gp = gpar(fontfamily = "Arial", fontface = "bold", fontsize = 10),
        labels_gp = gpar(fontfamily = "Arial", fontface = "bold", fontsize = 8.5),
        at = c(-3, 0, 3)
      )
    )
    draw(
      ht,
      heatmap_legend_side = "right",
      annotation_legend_side = "right",
      merge_legends = TRUE,
      padding = unit(c(7, 5, 5, 5), "mm")
    )
  }

  write_csv(dat$membership, file.path(result_dir, paste0("driver_gene_membership_", cid, ".csv")))
  save_heatmap(draw_fun, stem, 245, 173)
}

save_driver_matrix(
  "TREM2_vs_Resident",
  "D1_leading_gene_concept_matrix_TREM2_vs_Resident"
)
save_driver_matrix(
  "OC_vs_Resident",
  "D2_leading_gene_concept_matrix_OC_vs_Resident"
)
save_driver_matrix(
  "OC_vs_TREM2",
  "D3_leading_gene_concept_matrix_OC_vs_TREM2"
)

# Audit disease-named KEGG terms without using them as biological conclusions.
disease_regex <- regex(
  "Alzheimer|Parkinson|Huntington|Amyotrophic lateral sclerosis|neurodegeneration|Prion disease",
  ignore_case = TRUE
)
reference_sets <- res$kegg_term2gene |>
  filter(pathway_id %in% c("hsa00190", "hsa04142")) |>
  mutate(module = recode(pathway_id,
    "hsa00190" = "Mitochondrial / OXPHOS genes",
    "hsa04142" = "Lysosome genes"
  ))
ox_ids <- reference_sets |> filter(module == "Mitochondrial / OXPHOS genes") |> pull(ENTREZID) |> unique()
ly_ids <- reference_sets |> filter(module == "Lysosome genes") |> pull(ENTREZID) |> unique()

disease_audit <- gsea_kegg |>
  filter(str_detect(Description, disease_regex)) |>
  mutate(
    leading_ids = map(leading_entrez, parse_genes),
    n_leading = map_int(leading_ids, length),
    n_oxphos = map_int(leading_ids, ~ sum(.x %in% ox_ids)),
    n_lysosome = map_int(leading_ids, ~ sum(.x %in% setdiff(ly_ids, ox_ids))),
    n_other = pmax(n_leading - n_oxphos - n_lysosome, 0L),
    fraction_oxphos = n_oxphos / pmax(n_leading, 1),
    fraction_lysosome = n_lysosome / pmax(n_leading, 1),
    mechanistic_interpretation = case_when(
      fraction_oxphos >= 0.25 ~ "Predominantly shared mitochondrial / OXPHOS genes",
      fraction_lysosome >= 0.15 ~ "Substantial shared lysosomal genes",
      TRUE ~ "Mixed shared cellular genes; disease label not interpreted literally"
    )
  ) |>
  select(-leading_ids)
write_csv(disease_audit, file.path(result_dir, "KEGG_disease_label_shared_gene_audit.csv"))

disease_plot_data <- disease_audit |>
  filter(p.adjust < 0.05) |>
  group_by(comparison_id) |>
  arrange(p.adjust, desc(abs(NES))) |>
  slice_head(n = 5) |>
  ungroup() |>
  select(comparison_id, comparison, Description, n_oxphos, n_lysosome, n_other, n_leading) |>
  pivot_longer(c(n_oxphos, n_lysosome, n_other), names_to = "driver_class", values_to = "n_genes") |>
  mutate(
    fraction = n_genes / n_leading,
    driver_class = recode(
      driver_class,
      n_oxphos = "Mitochondrial / OXPHOS",
      n_lysosome = "Lysosome",
      n_other = "Other shared genes"
    ),
    pathway_label = str_wrap(Description, 34),
    comparison = factor(comparison, levels = comparison_table$comparison_label)
  )

if (nrow(disease_plot_data)) {
  disease_cols <- c(
    "Mitochondrial / OXPHOS" = "#79B8A9",
    "Lysosome" = "#DDA85C",
    "Other shared genes" = "#D9D9D9"
  )
  p_disease <- ggplot(disease_plot_data, aes(fraction, pathway_label, fill = driver_class)) +
    geom_col(width = 0.72, colour = "white", linewidth = 0.25) +
    facet_wrap(~comparison, ncol = 1, scales = "free_y") +
    scale_x_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1), expand = c(0, 0)) +
    scale_fill_manual(values = disease_cols) +
    labs(
      title = "Shared cellular machinery behind disease-named KEGG labels",
      subtitle = "Supplementary audit only; these labels are not interpreted as evidence of the named diseases",
      x = "Fraction of leading-edge genes", y = NULL, fill = "Actual driver class"
    ) +
    theme_pub +
    theme(legend.position = "top", axis.text.y = element_text(size = 8.8))
  save_plot(p_disease, "S1_disease_named_KEGG_shared_driver_audit", 220, 132)
}

figure_manifest <- tribble(
  ~stem, ~width_mm, ~height_mm, ~role,
  "A_pooled_myeloid_pairwise_DE_volcano", 245, 112, "Exploratory pairwise transcript differences after pooling all retained cells across both specimens",
  "B_GO_nonredundant_leading_gene_enrichment", 245, 188, "Eight representative GO processes with GeneRatio, Count, FDR and direction",
  "C_KEGG_nonredundant_leading_gene_enrichment", 245, 188, "Eight representative KEGG pathways with GeneRatio, Count, FDR and direction",
  "D1_leading_gene_concept_matrix_TREM2_vs_Resident", 245, 173, "Leading-gene overlap and effect direction",
  "D2_leading_gene_concept_matrix_OC_vs_Resident", 245, 173, "Leading-gene overlap and effect direction",
  "D3_leading_gene_concept_matrix_OC_vs_TREM2", 245, 173, "Leading-gene overlap and effect direction",
  "S1_disease_named_KEGG_shared_driver_audit", 220, 132, "Supplementary reinterpretation of disease-named KEGG terms"
)
write_csv(figure_manifest, file.path(result_dir, "figure_manifest.csv"))

writeLines(capture.output(sessionInfo()), file.path(result_dir, "sessionInfo_stage29.txt"))
message("Publication figures written to: ", figure_dir)
