root <- Sys.getenv("CB_REVISION_ROOT")
if (!nzchar(root)) root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
Sys.setenv(CB_REVISION_ROOT=root)
project_lib <- file.path(root,"renv","library","R-4.4","x86_64-w64-mingw32")
.libPaths(unique(c(project_lib,.Library,.Library.site)))
source(file.path(root,"R","00_setup.R"),chdir=FALSE)
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
})

out_dir <- file.path(root,"figures","publication_subfigures")
result_dir <- file.path(root,"results","annotation")
dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
dir.create(result_dir,recursive=TRUE,showWarnings=FALSE)

obj <- readRDS(file.path(root,"objects","07_rpca_broad_lineages_reviewed.rds"))
data_mat <- SeuratObject::LayerData(obj,assay="RNA",layer="data")
counts_mat <- SeuratObject::LayerData(obj,assay="RNA",layer="counts")
if (!"SH3BP2" %in% rownames(data_mat)) stop("SH3BP2 is absent from the RNA data layer")

broad_order <- c("Stromal/osteogenic","Mural","Endothelial","Myeloid","Lymphoid","Neural")
condition_order <- c("CTRL","CB")
condition_colors <- read_condition_colors()

cell_df <- data.frame(
  cell_id=colnames(obj),
  condition=as.character(obj$group),
  broad_lineage=as.character(obj$broad_lineage_reviewed),
  expression=as.numeric(data_mat["SH3BP2",colnames(obj)]),
  detected=as.numeric(counts_mat["SH3BP2",colnames(obj)] > 0),
  stringsAsFactors=FALSE
) |>
  filter(condition %in% condition_order,broad_lineage %in% broad_order) |>
  mutate(
    condition=factor(condition,levels=condition_order),
    broad_lineage=factor(broad_lineage,levels=broad_order)
  )

summary_df <- cell_df |>
  group_by(broad_lineage,condition) |>
  summarise(
    n_cells=n(),
    n_positive=sum(detected),
    mean_expression=mean(expression),
    median_expression=median(expression),
    q1_expression=quantile(expression,.25),
    q3_expression=quantile(expression,.75),
    percent_expressed=100*mean(detected),
    mean_positive_expression=ifelse(any(detected>0),mean(expression[detected>0]),NA_real_),
    .groups="drop"
  ) |>
  tidyr::complete(
    broad_lineage=factor(broad_order,levels=broad_order),
    condition=factor(condition_order,levels=condition_order)
  ) |>
  dplyr::mutate(n_cells=tidyr::replace_na(n_cells,0L))

readr::write_csv(cell_df,file.path(result_dir,"source_data_SH3BP2_cell_level_broad_lineage.csv"))
readr::write_csv(summary_df,file.path(result_dir,"source_data_SH3BP2_broad_lineage_summary.csv"))

base_theme <- theme_classic(base_size=10,base_family="Arial") +
  theme(
    text=element_text(family="Arial",face="bold",colour="black"),
    axis.title=element_text(size=10.5,face="bold"),
    axis.text=element_text(size=9.2,face="bold",colour="black"),
    axis.line=element_line(linewidth=.5,colour="black"),
    axis.ticks=element_line(linewidth=.45,colour="black"),
    plot.title=element_text(size=11,face="bold",hjust=.5),
    legend.title=element_text(size=9.5,face="bold"),
    legend.text=element_text(size=9,face="bold"),
    plot.margin=margin(5,8,5,6)
  )

# Panel A: distribution is shown descriptively. Violin widths are normalized
# within each lineage-condition group; no cell-level inferential P values.
p_dist <- ggplot(dplyr::filter(cell_df,detected>0),aes(broad_lineage,expression,fill=condition)) +
  geom_violin(
    position=position_dodge(width=.82),width=.78,scale="width",trim=TRUE,
    linewidth=.3,colour="white",alpha=.82
  ) +
  geom_boxplot(
    aes(colour=condition),position=position_dodge(width=.82),width=.16,
    outlier.shape=NA,fill="white",linewidth=.42
  ) +
  stat_summary(
    aes(colour=condition),fun=mean,geom="point",
    position=position_dodge(width=.82),shape=21,fill="white",size=2.2,stroke=.65
  ) +
  scale_fill_manual(values=condition_colors,drop=FALSE) +
  scale_colour_manual(values=condition_colors,drop=FALSE) +
  scale_y_continuous(expand=expansion(mult=c(0,.04))) +
  scale_x_discrete(limits=rev(broad_order)) +
  coord_flip() +
  base_theme +
  theme(
    legend.position="top",
    axis.text.x=element_text(size=9.2),
    axis.text.y=element_text(size=9.2),
    plot.tag=element_text(size=12,face="bold")
  ) +
  labs(
    x=NULL,y="SH3BP2 expression in positive cells\n(log-normalized)",
    fill="Condition",colour="Condition",tag="A"
  ) +
  guides(colour="none")

# Panel B: mean abundance (colour) and detection frequency (size) jointly avoid
# interpreting a dropout-heavy single-cell distribution from mean alone.
summary_dot <- summary_df |>
  mutate(
    condition=factor(condition,levels=condition_order),
    broad_lineage=factor(broad_lineage,levels=rev(broad_order))
  )
p_dot <- ggplot(summary_dot,aes(condition,broad_lineage)) +
  geom_point(aes(size=percent_expressed,colour=mean_expression),na.rm=TRUE) +
  geom_text(
    data=dplyr::filter(summary_dot,n_cells==0),label="No cells",
    family="Arial",fontface="bold",size=3.1,colour="#555555"
  ) +
  scale_size(range=c(2.2,10),limits=c(0,max(20,ceiling(max(summary_dot$percent_expressed,na.rm=TRUE)/5)*5)),breaks=scales::pretty_breaks(4)) +
  scale_colour_gradientn(colours=c("#FCE4D6","#F46D43","#B30000")) +
  base_theme +
  theme(
    axis.title=element_blank(),
    axis.text.x=element_text(size=10),
    axis.text.y=element_text(size=9.5),
    legend.position="right",
    plot.tag=element_text(size=12,face="bold")
  ) +
  labs(
    size="Percent expressed",colour="Average expression",tag="B"
  )

fig <- p_dist + p_dot +
  plot_layout(widths=c(1.72,1),guides="keep") +
  plot_annotation(
    title="SH3BP2 expression across major lineages",
    theme=theme(
      plot.title=element_text(family="Arial",face="bold",size=13,hjust=.5),
      plot.background=element_rect(fill="white",colour=NA)
    )
  )

base <- file.path(out_dir,"SH3BP2_expression_CB_CTRL_major_lineages")
save_pub_r(fig,base,210,92,dpi=600)

message("SH3BP2 broad-lineage comparison exported")
