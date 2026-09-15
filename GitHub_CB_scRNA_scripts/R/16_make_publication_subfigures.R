root <- Sys.getenv("CB_REVISION_ROOT")
if (!nzchar(root)) root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
Sys.setenv(CB_REVISION_ROOT = root)
project_lib <- file.path(root, "renv", "library", "R-4.4", "x86_64-w64-mingw32")
.libPaths(unique(c(project_lib, .Library, .Library.site)))
source(file.path(root, "R", "00_setup.R"), chdir = FALSE)
suppressPackageStartupMessages(library(ggplot2))

out_dir <- file.path(root, "figures", "publication_subfigures")
src_dir <- file.path(root, "results", "annotation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
rpca <- readRDS(file.path(root, "objects", "07_rpca_broad_lineages_reviewed.rds"))

pal <- readr::read_csv(file.path(root, "config", "palette.csv"), show_col_types = FALSE)
fine_colors <- stats::setNames(pal$color[pal$category == "celltype"], pal$label[pal$category == "celltype"])
if (!"pDC" %in% names(fine_colors)) fine_colors["pDC"] <- "#6A3D9A"
if (!"Plasma cells" %in% names(fine_colors)) fine_colors["Plasma cells"] <- "#F781BF"
condition_colors <- read_condition_colors()
broad_colors <- c("Stromal/osteogenic"="#8DD3C7", "Mural"="#7393B3", "Endothelial"="#43AA8B", "Myeloid"="#33A02C", "Lymphoid"="#CCEBC5", "Neural"="#FB9A99")
display_name <- c("Act OB-Stromal"="ECM-remodeling stromal state")
rename_display <- function(x) { x <- as.character(x); x[x %in% names(display_name)] <- unname(display_name[x[x %in% names(display_name)]]); x }
fine_colors_display <- fine_colors
names(fine_colors_display) <- rename_display(names(fine_colors_display))
fine_order <- c("FibroStromal","LesionStromal","Act OB-Stromal","ChondroFibroStromal","LEPR+ MSC-like","Mature OB","CNCC-OCPs","Pericytes","VSMC","Endothelial","Resident Mac","Inflam Mono","TREM2+ Mac","OC","Neutrophils","Mast cells","pDC","T/NK","B cells","Plasma cells","Schwann")
broad_order <- names(broad_colors)

theme_ref <- theme_cb(9) + theme(
  text=element_text(family="Arial",face="bold",colour="black"),
  plot.margin=margin(4,4,4,4),
  plot.title=element_text(size=10,face="bold"),
  axis.title=element_text(size=9,face="bold"),
  axis.text=element_text(size=8,face="bold",colour="black"),
  legend.title=element_text(size=8,face="bold"),
  legend.text=element_text(size=7.5,face="bold"),
  strip.text=element_text(size=9,face="bold")
)
theme_umap <- theme_ref + theme(axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank(), axis.line=element_blank(), legend.position="none")
save <- function(p, name, w=88, h=70) save_pub_r(p, file.path(out_dir,name), w, h, dpi=600)

umap <- as.data.frame(Seurat::Embeddings(rpca, reduction="umap.rpca")) |>
  tibble::rownames_to_column("cell_id")
names(umap)[2:3] <- c("UMAP_1","UMAP_2")
umap$condition <- factor(as.character(rpca$group), levels=c("CTRL","CB"))
umap$broad <- factor(as.character(rpca$broad_lineage_reviewed), levels=broad_order)
umap$celltype <- factor(rename_display(rpca$celltype_reviewed_provisional), levels=rename_display(fine_order))

label_df <- function(d, col) dplyr::summarise(dplyr::group_by(d, .data[[col]]), UMAP_1=median(UMAP_1), UMAP_2=median(UMAP_2), .groups="drop")
add_arrows <- function(p) p + annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, arrow=grid::arrow(length=unit(2.3,"mm")), linewidth=.45) + annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, arrow=grid::arrow(length=unit(2.3,"mm")), linewidth=.45) + annotate("text", x=Inf, y=-Inf, label="UMAP1", hjust=1.1, vjust=1.8, size=3.1, family="Arial", fontface="bold") + annotate("text", x=-Inf, y=Inf, label="UMAP2", angle=90, hjust=1.1, vjust=-.5, size=3.1, family="Arial", fontface="bold")

p <- ggplot(umap, aes(UMAP_1,UMAP_2,colour=celltype)) + ggrastr::geom_point_rast(size=.18, alpha=.8, raster.dpi=600) + scale_colour_manual(values=fine_colors_display, drop=FALSE) + theme_umap + guides(colour=guide_legend(ncol=2, override.aes=list(size=2.3))) + labs(colour="Cell type")
p <- p + ggrepel::geom_text_repel(data=label_df(umap,"celltype"), aes(x=UMAP_1, y=UMAP_2, label=`celltype`), size=2.75, colour="black", family="Arial", fontface="bold", seed=1, box.padding=.25, point.padding=.12, max.overlaps=Inf, inherit.aes=FALSE)
save(add_arrows(p), "UMAP_fine_cell_states", 110, 88)

p <- ggplot(umap, aes(UMAP_1,UMAP_2,colour=broad)) + ggrastr::geom_point_rast(size=.2, alpha=.82, raster.dpi=600) + scale_colour_manual(values=broad_colors, drop=FALSE) + theme_umap + guides(colour=guide_legend(ncol=2, override.aes=list(size=2.5))) + labs(colour="Major lineage")
p <- p + ggrepel::geom_text_repel(data=label_df(umap,"broad"), aes(x=UMAP_1, y=UMAP_2, label=`broad`), size=3.2, colour="black", family="Arial", fontface="bold", seed=1, box.padding=.35, inherit.aes=FALSE)
save(add_arrows(p), "UMAP_major_lineages", 110, 88)

# Absolute cell counts corresponding exactly to the six broad lineages shown in
# UMAP_major_lineages. This is descriptive because each condition has one donor.
broad_counts <- data.frame(
  condition=as.character(rpca$group),
  broad_lineage=as.character(rpca$broad_lineage_reviewed),
  stringsAsFactors=FALSE
) |>
  dplyr::filter(!is.na(condition), !is.na(broad_lineage)) |>
  dplyr::count(condition, broad_lineage, name="n_cells") |>
  tidyr::complete(
    condition=c("CTRL","CB"), broad_lineage=broad_order,
    fill=list(n_cells=0L)
  ) |>
  dplyr::group_by(condition) |>
  dplyr::mutate(total_cells=sum(n_cells), percent_cells=100*n_cells/total_cells) |>
  dplyr::ungroup()
broad_counts$condition <- factor(broad_counts$condition, levels=c("CTRL","CB"))
broad_counts$broad_lineage <- factor(broad_counts$broad_lineage, levels=rev(broad_order))
readr::write_csv(broad_counts, file.path(src_dir,"source_data_major_lineage_cell_counts.csv"))

p <- ggplot(broad_counts, aes(x=n_cells, y=broad_lineage, fill=condition)) +
  geom_col(position=position_dodge(width=.76), width=.64, colour="white", linewidth=.18) +
  geom_text(
    aes(label=scales::comma(n_cells)), position=position_dodge(width=.76),
    hjust=-.18, family="Arial", fontface="bold", size=3.4, colour="black"
  ) +
  scale_fill_manual(values=condition_colors, drop=FALSE) +
  scale_x_continuous(
    labels=scales::comma, breaks=scales::pretty_breaks(n=5),
    expand=expansion(mult=c(0,.16))
  ) +
  theme_classic(base_size=10, base_family="Arial") +
  theme(
    text=element_text(family="Arial", face="bold", colour="black"),
    axis.title.x=element_text(size=10.5, face="bold"),
    axis.title.y=element_blank(),
    axis.text.x=element_text(size=9, face="bold", colour="black"),
    axis.text.y=element_text(size=10, face="bold", colour="black", margin=margin(r=5)),
    axis.line=element_line(linewidth=.5, colour="black"),
    axis.ticks=element_line(linewidth=.45, colour="black"),
    legend.position="top",
    legend.title=element_text(size=10, face="bold"),
    legend.text=element_text(size=9.5, face="bold"),
    plot.title=element_text(size=11, face="bold", hjust=.5),
    plot.margin=margin(5,12,5,5)
  ) +
  labs(x="Number of cells", fill="Condition", title="Major lineage cell counts")
save(p, "Barplot_major_lineage_cell_counts", 145, 92)

p <- ggplot(umap, aes(UMAP_1,UMAP_2,colour=condition)) + ggrastr::geom_point_rast(size=.22, alpha=.8, raster.dpi=600) + scale_colour_manual(values=condition_colors) + facet_wrap(~condition, nrow=1) + theme_umap + theme(strip.placement="outside", strip.background=element_blank(), strip.text=element_text(size=10, face="bold")) + labs(colour="Condition")
save(add_arrows(p), "UMAP_condition_CB_CTRL_split", 160, 85)

split_labels <- umap |> dplyr::filter(!is.na(celltype)) |> dplyr::group_by(condition, celltype) |> dplyr::summarise(UMAP_1=median(UMAP_1), UMAP_2=median(UMAP_2), .groups="drop")
p <- ggplot(umap, aes(UMAP_1,UMAP_2,colour=celltype)) + ggrastr::geom_point_rast(size=.2, alpha=.84, raster.dpi=600) + scale_colour_manual(values=fine_colors_display, drop=FALSE) + facet_wrap(~condition, nrow=1) + ggrepel::geom_text_repel(data=split_labels, aes(x=UMAP_1,y=UMAP_2,label=celltype), colour="black", family="Arial", fontface="bold", size=2.35, seed=2, box.padding=.16, point.padding=.07, max.overlaps=Inf, inherit.aes=FALSE) + theme_umap + theme(strip.background=element_blank(), strip.text=element_text(size=10, face="bold"))
save(add_arrows(p), "UMAP_cell_states_CB_CTRL_split", 175, 90)

comp <- readr::read_csv(file.path(src_dir,"source_data_composition_barplots.csv"), show_col_types=FALSE)
comp$condition <- factor(comp$condition, levels=c("CTRL","CB"))
comp$celltype_display <- rename_display(comp$celltype)
comp$celltype_display <- factor(comp$celltype_display, levels=rev(rename_display(fine_order)))
comp$display_group <- factor(comp$display_group, levels=c("Mesenchymal / mural","Hematopoietic / endothelial / neural"))
comp_mes <- dplyr::filter(comp, display_group == "Mesenchymal / mural")
p <- ggplot(comp_mes, aes(percent_of_all_cells, celltype_display, fill=condition)) + geom_col(position=position_dodge(width=.72), width=.62, colour="white", linewidth=.15) + scale_fill_manual(values=condition_colors) + scale_x_continuous(expand=expansion(mult=c(0,.04))) + theme_ref + theme(legend.position="top", axis.title.y=element_blank()) + labs(x="Cells (% of all cells in each sample)", fill="Condition", title="Mesenchymal / mural")
save(p, "Barplot_mesenchymal_mural", 112, 82)
comp_hem <- dplyr::filter(comp, display_group == "Hematopoietic / endothelial / neural")
p <- ggplot(comp_hem, aes(percent_of_all_cells, celltype_display, fill=condition)) + geom_col(position=position_dodge(width=.72), width=.62, colour="white", linewidth=.15) + scale_fill_manual(values=condition_colors) + scale_x_continuous(expand=expansion(mult=c(0,.04))) + theme_ref + theme(legend.position="top", axis.title.y=element_blank()) + labs(x="Cells (% of all cells in each sample)", fill="Condition", title="Hematopoietic / endothelial / neural")
save(p, "Barplot_hematopoietic_endothelial_neural", 112, 96)

stack <- comp |> dplyr::group_by(condition, celltype_display) |> dplyr::summarise(n_cells=sum(n_cells), .groups="drop") |> dplyr::group_by(condition) |> dplyr::mutate(percent=100*n_cells/sum(n_cells))
stack$celltype_display <- factor(as.character(stack$celltype_display), levels=rename_display(fine_order))
p <- ggplot(stack, aes(condition, percent, fill=celltype_display)) + geom_col(width=.68, colour="white", linewidth=.12) + scale_fill_manual(values=fine_colors_display, drop=FALSE) + scale_y_continuous(expand=c(0,0), limits=c(0,100), breaks=seq(0,100,25)) + theme_ref + theme(legend.position="right", axis.title.x=element_blank()) + labs(y="Percentage of cells", fill="Cell type")
save(p, "Barplot_stacked_cell_state_composition", 150, 92)

make_dot <- function(df, order, colors, name, w=145, h=105, feature_map=NULL, use_scaled=TRUE) {
  df$group <- factor(rename_display(df$group), levels=rename_display(order)); df$gene <- factor(df$gene, levels=unique(df$gene))
  colour_var <- if (use_scaled) "scaled_average" else "average_expression"
  p <- ggplot(df, aes(x=gene, y=group, size=percent_expressed, colour=.data[[colour_var]])) + geom_point() + scale_size(range=c(.7,7), breaks=c(0,25,50,75,100), limits=c(0,100)) + scale_colour_gradientn(colours=c("#FCE4D6","#F46D43","#B30000"), limits=if(use_scaled)c(-1,2) else range(df$average_expression,na.rm=TRUE), oob=scales::squish) + theme_ref + theme(axis.title=element_blank(), axis.text.x=element_text(angle=55, hjust=1, vjust=1, size=6.6, face="bold"), axis.text.y=element_text(size=7.4, face="bold"), legend.position="right") + labs(size="Percent expressed", colour=if(use_scaled)"Average expression\n(scaled)" else "Average expression")
  if (!is.null(feature_map)) {
    hd <- data.frame(gene=factor(names(feature_map), levels=levels(df$gene)), feature_group=unname(feature_map), stringsAsFactors=FALSE)
    hd$feature_group <- rename_display(hd$feature_group)
    hd <- hd |> dplyr::group_by(feature_group) |> dplyr::mutate(label=ifelse(dplyr::row_number()==floor(dplyr::n()/2)+1, dplyr::first(feature_group), "")) |> dplyr::ungroup()
    short_header <- c("Stromal/osteogenic"="Stromal/\nosteogenic", "FibroStromal"="Fibro\nStromal", "LesionStromal"="Lesion\nStromal", "ECM-remodeling stromal state"="ECM-remodeling\nstromal state", "ChondroFibroStromal"="ChondroFibro\nStromal", "LEPR+ MSC-like"="LEPR+\nMSC-like", "Mature OB"="Mature\nOB", "CNCC-OCPs"="CNCC-\nOCPs", "Endothelial"="Endothelial", "Resident Mac"="Resident\nMac", "Inflam Mono"="Inflam\nMono", "TREM2+ Mac"="TREM2+\nMac", "Neutrophils"="Neutrophils", "Mast cells"="Mast\ncells", "B cells"="B\ncells", "Plasma cells"="Plasma\ncells")
    hit <- hd$label %in% names(short_header)
    hd$label[hit] <- unname(short_header[hd$label[hit]])
    header_cols <- colors; names(header_cols) <- rename_display(names(header_cols))
    ph <- ggplot(hd, aes(gene, 1, fill=feature_group)) + geom_tile(colour="white", linewidth=.35) + geom_text(aes(label=label), colour="white", family="Arial", fontface="bold", size=2.45, lineheight=.86) + scale_fill_manual(values=header_cols, drop=FALSE) + theme_void(base_family="Arial") + theme(legend.position="none", plot.margin=margin(0,4,0,4))
    p <- patchwork::wrap_plots(ph, p, ncol=1, heights=c(.13,1), guides="collect")
  }
  save(p, name, w, h)
}

# Spacious, single-panel marker dot plot used for the two main annotation figures.
# Raster exports are exactly 6300 x 2900 px; vector exports use the matching
# 600-dpi physical dimensions.  Group headers form a continuous colour ribbon
# without the previous white tile borders.
save_main_dot <- function(plot, name, width_px=6300L, height_px=2900L) {
  dpi <- 600
  width_in <- width_px / dpi
  height_in <- height_px / dpi
  path <- file.path(out_dir, name)

  svglite::svglite(paste0(path, ".svg"), width=width_in, height=height_in)
  print(plot)
  grDevices::dev.off()
  grDevices::cairo_pdf(paste0(path, ".pdf"), width=width_in, height=height_in, family="Arial")
  print(plot)
  grDevices::dev.off()
  ragg::agg_tiff(paste0(path, ".tiff"), width=width_px, height=height_px,
                 units="px", res=dpi, compression="lzw")
  print(plot)
  grDevices::dev.off()
  ragg::agg_png(paste0(path, ".png"), width=width_px, height=height_px,
                units="px", res=dpi)
  print(plot)
  grDevices::dev.off()
}

make_main_dot <- function(df, marker_groups, colors, name,
                          width_px=6300L, height_px=2900L,
                          use_scaled=TRUE) {
  group_order <- names(marker_groups)
  gene_order <- unique(unlist(marker_groups, use.names=FALSE))

  # Add deliberate between-group space in the matrix while keeping the header
  # ribbon continuous. This improves both gene-label and bubble separation.
  gap <- 1.25
  gene_layout <- list()
  group_layout <- list()
  cursor <- 1
  for (i in seq_along(marker_groups)) {
    genes_i <- marker_groups[[i]]
    xpos <- cursor + seq_along(genes_i) - 1
    gene_layout[[i]] <- data.frame(gene=genes_i, x=xpos, group=group_order[i])
    group_layout[[i]] <- data.frame(
      group=group_order[i], first=min(xpos), last=max(xpos), mid=mean(xpos)
    )
    cursor <- max(xpos) + 1 + gap
  }
  gene_layout <- dplyr::bind_rows(gene_layout)
  group_layout <- dplyr::bind_rows(group_layout)
  if (nrow(group_layout) > 1) {
    boundaries <- (group_layout$last[-nrow(group_layout)] + group_layout$first[-1]) / 2
    group_layout$xmin <- c(min(gene_layout$x)-0.5, boundaries)
    group_layout$xmax <- c(boundaries, max(gene_layout$x)+0.5)
  } else {
    group_layout$xmin <- min(gene_layout$x)-0.5
    group_layout$xmax <- max(gene_layout$x)+0.5
  }

  display_order <- rename_display(group_order)
  gene_layout$group_display <- rename_display(gene_layout$group)
  group_layout$group_display <- rename_display(group_layout$group)
  header_labels <- c(
    FibroStromal="Fibro\nStromal", LesionStromal="Lesion\nStromal",
    Pericytes="Pericytes", `LEPR+ MSC-like`="LEPR+\nMSC-like",
    Endothelial="Endo-\nthelial", VSMC="VSMC",
    ChondroFibroStromal="Chondro\nFibro\nStromal", `Mature OB`="Mature\nOB",
    `ECM-remodeling stromal state`="ECM-\nremodeling\nstromal state",
    `CNCC-OCPs`="CNCC-\nOCPs", Schwann="Schwann",
    `T/NK`="T/NK", Neutrophils="Neutrophils", `Resident Mac`="Resident\nMac",
    OC="OC", `Inflam Mono`="Inflam\nMono", `TREM2+ Mac`="TREM2+\nMac",
    `B cells`="B\ncells", `Mast cells`="Mast\ncells", pDC="pDC",
    `Plasma cells`="Plasma\ncells"
  )
  group_layout$label <- group_layout$group_display
  hit <- group_layout$group_display %in% names(header_labels)
  group_layout$label[hit] <- unname(header_labels[group_layout$group_display[hit]])
  header_cols <- colors
  names(header_cols) <- rename_display(names(header_cols))
  group_layout$header_color <- unname(header_cols[group_layout$group_display])

  plot_df <- df |>
    dplyr::mutate(group=rename_display(group)) |>
    dplyr::inner_join(dplyr::select(gene_layout, gene, x), by="gene") |>
    dplyr::mutate(y=length(display_order)-match(group, display_order)+1)
  n_groups <- length(display_order)
  if (n_groups <= 2) {
    header_ymin <- n_groups + 0.18
    header_ymax <- n_groups + 0.52
    y_upper <- n_groups + 0.60
    header_text_size <- 3.0
  } else {
    header_ymin <- n_groups + 0.56
    header_ymax <- n_groups + 1.52
    y_upper <- n_groups + 1.60
    header_text_size <- 2.95
  }
  colour_var <- if (use_scaled) "scaled_average" else "average_expression"
  colour_limits <- if (use_scaled) c(-1,2) else range(plot_df$average_expression, na.rm=TRUE)
  colour_title <- if (use_scaled) "Average expression\n(scaled)" else "Average expression"

  p <- ggplot(plot_df, aes(x=x, y=y)) +
    geom_point(aes(size=percent_expressed, colour=.data[[colour_var]]), alpha=.98) +
    geom_rect(
      data=group_layout,
      aes(xmin=xmin, xmax=xmax, ymin=header_ymin, ymax=header_ymax),
      inherit.aes=FALSE, fill=group_layout$header_color, colour=NA
    ) +
    geom_text(
      data=group_layout,
      aes(x=mid, y=(header_ymin+header_ymax)/2, label=label),
      inherit.aes=FALSE, colour="white", family="Arial", fontface="bold",
      size=header_text_size, lineheight=.82
    ) +
    scale_x_continuous(
      breaks=gene_layout$x, labels=gene_layout$gene,
      limits=c(min(group_layout$xmin), max(group_layout$xmax)),
      expand=expansion(mult=c(.005,.005))
    ) +
    scale_y_continuous(
      breaks=seq_len(n_groups), labels=rev(display_order),
      limits=c(.45, y_upper), expand=c(0,0)
    ) +
    scale_size(range=c(1.2,8.2), breaks=c(0,25,50,75,100), limits=c(0,100)) +
    scale_colour_gradientn(
      colours=c("#FCE4D6","#F46D43","#B30000"), limits=colour_limits,
      oob=scales::squish
    ) +
    theme_classic(base_size=10, base_family="Arial") +
    theme(
      text=element_text(family="Arial", face="bold", colour="black"),
      axis.title=element_blank(),
      axis.text.x=element_text(angle=55, hjust=1, vjust=1, size=9, face="bold"),
      axis.text.y=element_text(size=10.5, face="bold", margin=margin(r=5)),
      axis.ticks=element_blank(),
      axis.line.x=element_line(linewidth=.45, colour="black"),
      axis.line.y=element_line(linewidth=.45, colour="black"),
      legend.position="right",
      legend.title=element_text(size=10, face="bold"),
      legend.text=element_text(size=9.5, face="bold"),
      legend.spacing.y=unit(2.5,"mm"),
      plot.margin=margin(6,8,4,8)
    ) +
    guides(
      colour=guide_colourbar(title=colour_title, barheight=unit(27,"mm"), barwidth=unit(5,"mm"), order=1),
      size=guide_legend(title="Percent expressed", order=2, override.aes=list(colour="black"))
    )

  save_main_dot(p, name, width_px=width_px, height_px=height_px)
}
broad_groups <- list("Stromal/osteogenic"=c("COL1A1","DCN","PDGFRA","RUNX2"), Mural=c("RGS5","PDGFRB","ACTA2","MYH11"), Endothelial=c("PECAM1","EMCN","VWF","CLDN5"), Myeloid=c("LST1","TYROBP","FCER1G","LYZ"), Lymphoid=c("PTPRC","CD3D","NKG7","CD79A"), Neural=c("S100B","PLP1","MPZ","CDH19"))
broad_map <- unlist(lapply(names(broad_groups), function(n) stats::setNames(rep(n,length(broad_groups[[n]])), broad_groups[[n]])))
make_dot(readr::read_csv(file.path(src_dir,"source_data_broad_lineage_dotplot.csv"), show_col_types=FALSE), broad_order, broad_colors, "Bubbleplot_major_lineage_markers", 165, 102, broad_map, TRUE)

# Original-manuscript marker logic, expanded to the newly reviewed pDC and plasma-cell states.
# LesionStromal and the former Act OB-Stromal use full-scan-supported discriminators.
nonimmune_groups <- list(
  FibroStromal=c("ASPN","COL11A1","WNT5A"),
  LesionStromal=c("FAP","TNFSF11"),
  Pericytes=c("ABCC9","KCNJ8","RGS5"),
  `LEPR+ MSC-like`=c("LEPR","CXCL12","CFD"),
  Endothelial=c("EMCN","VWF"),
  VSMC=c("MYH11","CNN1"),
  ChondroFibroStromal=c("ACAN","COL10A1","COMP"),
  `Mature OB`=c("IFITM5","BGLAP","DMP1"),
  `ECM-remodeling stromal state`=c("POSTN","CTHRC1","IBSP"),
  `CNCC-OCPs`=c("SATB2","DLX5","BMPR1B")
)
immune_groups <- list(
  `T/NK`=c("CD3D","NKG7"),
  Neutrophils=c("FCGR3B","S100A8","CSF3R"),
  `Resident Mac`=c("CD163","C1QA","LYVE1"),
  OC=c("ACP5","SIGLEC15","TCIRG1"),
  `Inflam Mono`=c("FCN1","OLR1","S100A9"),
  `TREM2+ Mac`=c("TREM2","SLAMF8","HMOX1"),
  `B cells`=c("CD79A","MS4A1"),
  `Mast cells`=c("TPSB2","CPA3"),
  pDC=c("GZMB","LILRA4"),
  `Plasma cells`=c("JCHAIN","MZB1")
)
data_mat <- SeuratObject::LayerData(rpca, assay="RNA", layer="data")
counts_mat <- SeuratObject::LayerData(rpca, assay="RNA", layer="counts")
labels_display <- rename_display(rpca$celltype_reviewed_provisional)
summarise_dot_seurat_style <- function(marker_groups) {
  group_order <- names(marker_groups)
  gene_order <- unique(unlist(marker_groups,use.names=FALSE))
  gene_order <- intersect(gene_order,rownames(data_mat))
  out <- lapply(group_order,function(grp) {
    idx <- which(labels_display==grp)
    if (!length(idx)) return(NULL)
    avg <- log1p(Matrix::rowMeans(expm1(data_mat[gene_order,idx,drop=FALSE])))
    pct <- 100*Matrix::rowMeans(counts_mat[gene_order,idx,drop=FALSE]>0)
    data.frame(group=grp,gene=gene_order,average_expression=as.numeric(avg),percent_expressed=as.numeric(pct),stringsAsFactors=FALSE)
  }) |>
    dplyr::bind_rows() |>
    dplyr::group_by(gene) |>
    dplyr::mutate(scaled_average={s<-stats::sd(average_expression); if(is.finite(s)&&s>0) as.numeric(scale(average_expression)) else rep(0,dplyr::n())},scaled_average=pmax(-1,pmin(2,scaled_average))) |>
    dplyr::ungroup()
  out$gene <- factor(out$gene,levels=gene_order)
  out
}
make_map <- function(x) unlist(lapply(names(x),function(n) stats::setNames(rep(n,length(x[[n]])),x[[n]])))
nonimmune_dot <- summarise_dot_seurat_style(nonimmune_groups)
immune_dot <- summarise_dot_seurat_style(immune_groups)
full_groups <- c(nonimmune_groups,list(Schwann=c("PLP1","MPZ","S100B")),immune_groups)
full_dot <- summarise_dot_seurat_style(full_groups)
readr::write_csv(nonimmune_dot,file.path(src_dir,"source_data_nonimmune_marker_dotplot_original_style.csv"))
readr::write_csv(immune_dot,file.path(src_dir,"source_data_immune_marker_dotplot_original_style.csv"))
readr::write_csv(full_dot,file.path(src_dir,"source_data_fine_marker_dotplot_original_style.csv"))

# Keep legacy filenames for downstream manuscript links, while also writing explicit names.
make_main_dot(nonimmune_dot,nonimmune_groups,fine_colors_display,"Bubbleplot_mesenchymal_mural_markers",8000L,2900L)
make_main_dot(nonimmune_dot,nonimmune_groups,fine_colors_display,"Bubbleplot_mesenchymal_endothelial_other_markers",8000L,2900L)
make_main_dot(immune_dot,immune_groups,fine_colors_display,"Bubbleplot_hematopoietic_endothelial_neural_markers")
make_main_dot(immune_dot,immune_groups,fine_colors_display,"Bubbleplot_immune_cell_markers")
make_dot(full_dot,names(full_groups),fine_colors_display,"Bubbleplot_fine_cell_state_markers",300,165,NULL,TRUE)

pair_groups <- list(
  LesionStromal=c("PDZRN4","BICC1","GPC6","RUNX2"),
  `ECM-remodeling stromal state`=c("MYL9","IGFBP7","LUM","CTHRC1")
)
pair_genes <- unlist(pair_groups,use.names=FALSE)
pair_dot <- summarise_dot_seurat_style(pair_groups) |>
  dplyr::filter(group %in% names(pair_groups),as.character(gene) %in% pair_genes)
readr::write_csv(pair_dot,file.path(src_dir,"source_data_lesion_vs_ecm_remodeling_dotplot.csv"))
make_main_dot(
  pair_dot,pair_groups,fine_colors_display,"Bubbleplot_stromal_state_discriminators",
  width_px=6800L,height_px=3600L,use_scaled=FALSE
)

genes <- intersect(c("COL1A1","PDZRN4","WNT5A","CTHRC1","LUM","TREM2","ACP5","SH3BP2"), rownames(data_mat))
for (g in genes) {
  z <- data.frame(umap[,c("UMAP_1","UMAP_2")], expr=as.numeric(data_mat[g,umap$cell_id]))
  p <- ggplot(z, aes(UMAP_1,UMAP_2)) + ggrastr::geom_point_rast(colour="#78A6C8", size=.2, alpha=.45, raster.dpi=600) + ggrastr::geom_point_rast(data=z[z$expr>0,], aes(colour=expr), size=.25, alpha=.82, raster.dpi=600) + scale_colour_gradientn(colours=c("#6BAED6","#F7FBFF","#F46D43")) + theme_umap + theme(legend.position="right",plot.title=element_text(size=11,face="bold"),legend.title=element_text(size=8.5,face="bold"),legend.text=element_text(size=8,face="bold")) + labs(title=g, colour="Expression")
  save(p, paste0("FeatureUMAP_", g), 94, 78)
}

# Ten-gene localization plate matching the original manuscript's 2 x 5 layout.
# The display range is clipped to 1-2 solely for localization contrast, as in
# the reference styling; the source-data table retains the unmodified
# log-normalized expression values.
feature_grid_genes <- c(
  "COL11A1","WNT5A","CXCL12","FCGR3B","ACAN",
  "POSTN","SH3BP2","C1QA","ACP5","TREM2"
)
missing_feature_genes <- setdiff(feature_grid_genes, rownames(data_mat))
if (length(missing_feature_genes)) {
  stop("Feature-grid genes absent from RNA assay: ", paste(missing_feature_genes, collapse=", "))
}

feature_source <- lapply(feature_grid_genes, function(g) {
  data.frame(
    cell_id=umap$cell_id,
    UMAP_1=umap$UMAP_1,
    UMAP_2=umap$UMAP_2,
    gene=g,
    expression_log_normalized=as.numeric(data_mat[g,umap$cell_id]),
    stringsAsFactors=FALSE
  )
}) |>
  dplyr::bind_rows() |>
  dplyr::mutate(expression_display=pmax(1,pmin(2,expression_log_normalized)))
readr::write_csv(feature_source, file.path(src_dir,"source_data_featureplot_10gene_2x5.csv"))

feature_panels <- lapply(feature_grid_genes, function(g) {
  z <- feature_source |>
    dplyr::filter(gene == g) |>
    dplyr::arrange(expression_display)
  ggplot(z, aes(UMAP_1,UMAP_2,colour=expression_display)) +
    ggrastr::geom_point_rast(size=.16, alpha=.94, raster.dpi=600) +
    scale_colour_gradientn(
      colours=c("#78A6C8","#F7FBFF","#F47A55"),
      values=scales::rescale(c(1,1.5,2)), limits=c(1,2),
      breaks=c(1,1.25,1.5,1.75,2),
      labels=scales::label_number(accuracy=.01), oob=scales::squish
    ) +
    coord_equal(expand=TRUE) +
    theme_void(base_family="Arial") +
    theme(
      text=element_text(family="Arial",face="bold",colour="black"),
      panel.border=element_rect(colour="black",fill=NA,linewidth=.42),
      plot.title=element_text(size=8.5,face="bold",hjust=.5,margin=margin(b=2)),
      legend.position="right",
      legend.title=element_blank(),
      legend.text=element_text(size=6.8,face="bold",colour="black"),
      legend.key.height=unit(3.3,"mm"),
      legend.key.width=unit(3.2,"mm"),
      legend.margin=margin(0,0,0,1),
      legend.box.spacing=unit(.7,"mm"),
      plot.margin=margin(2.5,2.5,2.5,2.5)
    ) +
    guides(colour=guide_colourbar(barheight=unit(15,"mm"),barwidth=unit(3.2,"mm"))) +
    labs(title=g)
})
feature_grid <- patchwork::wrap_plots(feature_panels,ncol=5,nrow=2,guides="keep")
save_main_dot(feature_grid,"FeaturePlot_10gene_2x5",width_px=6107L,height_px=2036L)

manifest <- data.frame(file=sort(list.files(out_dir, pattern="\\.(svg|pdf|tiff|png)$", full.names=FALSE)), stringsAsFactors=FALSE)
readr::write_csv(manifest, file.path(src_dir,"publication_subfigure_export_manifest.csv"))
log_step("Publication subfigures exported: ", nrow(manifest), " files")
