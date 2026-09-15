root <- Sys.getenv("CB_REVISION_ROOT")
if (!nzchar(root)) root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
Sys.setenv(CB_REVISION_ROOT=root)
project_lib <- file.path(root,"renv","library","R-4.4","x86_64-w64-mingw32")
.libPaths(unique(c(project_lib,.Library,.Library.site)))
source(file.path(root,"R","00_setup.R"),chdir=FALSE)
suppressPackageStartupMessages(library(ggplot2))

fig_dir <- file.path(root,"figures","marker_reassessment")
res_dir <- file.path(root,"results","annotation")
dir.create(fig_dir,recursive=TRUE,showWarnings=FALSE)
obj <- readRDS(file.path(root,"objects","07_rpca_broad_lineages_reviewed.rds"))
mat <- SeuratObject::LayerData(obj,assay="RNA",layer="data")
counts <- SeuratObject::LayerData(obj,assay="RNA",layer="counts")
labels <- as.character(obj$celltype_reviewed_provisional)

marker_list <- list(
  FibroStromal=c("ASPN","COL11A1"),
  LesionStromal=c("WNT5A","FAP","TNFSF11"),
  `T/NK`=c("CD3D","NKG7"),
  Pericytes=c("ABCC9","KCNJ8","RGS5"),
  `LEPR+ MSC-like`=c("LEPR","CXCL12","CFD"),
  Endothelial=c("EMCN","VWF"),
  VSMC=c("MYH11","CNN1"),
  Neutrophils=c("FCGR3B","S100A8"),
  ChondroFibroStromal=c("ACAN","COL10A1","COMP"),
  `Mature OB`=c("IFITM5","BGLAP","DMP1"),
  `ECM-remodeling stromal`=c("POSTN","CTHRC1","IBSP"),
  `CNCC-OCPs`=c("SATB2","DLX5","BMPR1B"),
  `Resident Mac`=c("CD163","C1QA","LYVE1"),
  OC=c("ACP5","SIGLEC15","TCIRG1"),
  `Inflam Mono`=c("FCN1","OLR1"),
  `TREM2+ Mac`=c("TREM2","SLAMF8","HMOX1"),
  Schwann=c("PLP1","MPZ"),
  `B cells`=c("CD79A","MS4A1"),
  `Mast cells`=c("TPSB2","CPA3"),
  pDC=c("GZMB","LILRA4"),
  `Plasma cells`=c("JCHAIN","MZB1")
)
label_display <- labels
label_display[label_display=="Act OB-Stromal"] <- "ECM-remodeling stromal"
cell_order <- names(marker_list)
genes <- unique(unlist(marker_list,use.names=FALSE))
genes <- intersect(genes,rownames(mat))

# Reproduce Seurat DotPlot's mean-expression logic: average on the linear scale,
# then log1p, followed by gene-wise scaling across cell states.
dot_rows <- lapply(cell_order,function(grp){
  idx <- which(label_display==grp)
  if (!length(idx)) return(NULL)
  avg <- log1p(Matrix::rowMeans(expm1(mat[genes,idx,drop=FALSE])))
  pct <- 100*Matrix::rowMeans(counts[genes,idx,drop=FALSE]>0)
  data.frame(group=grp,gene=genes,average_expression=as.numeric(avg),percent_expressed=as.numeric(pct),stringsAsFactors=FALSE)
})
dot <- dplyr::bind_rows(dot_rows) |>
  dplyr::group_by(gene) |>
  dplyr::mutate(scaled_average={s<-stats::sd(average_expression); if(is.finite(s)&&s>0) as.numeric(scale(average_expression)) else rep(0,dplyr::n())},scaled_average=pmax(-1,pmin(2,scaled_average))) |>
  dplyr::ungroup()
feature_order <- unlist(marker_list,use.names=FALSE)
feature_order <- feature_order[feature_order %in% genes]
feature_owner <- unlist(lapply(names(marker_list),function(n) stats::setNames(rep(n,length(marker_list[[n]])),marker_list[[n]])))
dot$group <- factor(dot$group,levels=rev(cell_order))
dot$gene <- factor(dot$gene,levels=feature_order)
dot$xid <- as.numeric(dot$gene)
readr::write_csv(dot,file.path(res_dir,"source_data_original_marker_panel_recomputed.csv"))

pal <- readr::read_csv(file.path(root,"config","palette.csv"),show_col_types=FALSE)
cols <- stats::setNames(pal$color[pal$category=="celltype"],pal$label[pal$category=="celltype"])
cols["ECM-remodeling stromal"] <- cols["Act OB-Stromal"]
ann <- data.frame(gene=feature_order,owner=unname(feature_owner[feature_order]),xid=seq_along(feature_order),stringsAsFactors=FALSE) |>
  dplyr::group_by(owner) |>
  dplyr::summarise(xmin=min(xid)-.5,xmax=max(xid)+.5,xmid=mean(range(xid)),.groups="drop")
short <- c(FibroStromal="Fibro\nStromal",LesionStromal="Lesion\nStromal",`ECM-remodeling stromal`="ECM-remodeling\nstromal",ChondroFibroStromal="ChondroFibro\nStromal",`LEPR+ MSC-like`="LEPR+\nMSC-like",`Mature OB`="Mature\nOB",`CNCC-OCPs`="CNCC-\nOCPs",`Resident Mac`="Resident\nMac",`Inflam Mono`="Inflam\nMono",`TREM2+ Mac`="TREM2+\nMac",`B cells`="B\ncells",`Mast cells`="Mast\ncells",`Plasma cells`="Plasma\ncells")
ann$label <- ann$owner
ann$label[ann$label %in% names(short)] <- unname(short[ann$label[ann$label %in% names(short)]])
p_main <- ggplot(dot,aes(x=xid,y=group))+
  geom_point(aes(size=percent_expressed,colour=scaled_average))+
  scale_size(range=c(.25,5.2),limits=c(0,100),breaks=c(0,25,50,75,100))+
  scale_colour_gradientn(colours=c("#FDE0DD","#FC9272","#D7301F"),limits=c(-1,2),oob=scales::squish)+
  scale_x_continuous(breaks=seq_along(feature_order),labels=feature_order,limits=c(.5,length(feature_order)+.5),expand=c(0,0))+
  labs(x=NULL,y=NULL,size="Percent expressed",colour="Average expression\n(scaled)")+
  theme_cb(7)+theme(panel.border=element_rect(colour="black",fill=NA,linewidth=.35),axis.line=element_blank(),axis.text.x=element_text(angle=55,hjust=1,vjust=1,size=5.2),axis.text.y=element_text(size=5.7),legend.position="right",plot.margin=margin(0,3,3,3))
p_top <- ggplot(ann)+geom_rect(aes(xmin=xmin,xmax=xmax,ymin=0,ymax=1,fill=owner),colour="white",linewidth=.35)+geom_text(aes(x=xmid,y=.5,label=label),colour="white",family="Arial",fontface="bold",size=2,lineheight=.85)+scale_fill_manual(values=cols)+scale_x_continuous(limits=c(.5,length(feature_order)+.5),expand=c(0,0))+theme_void(base_family="Arial")+theme(legend.position="none",plot.margin=margin(2,3,0,3))
p_original <- patchwork::wrap_plots(p_top,p_main,ncol=1,heights=c(.13,1),guides="collect")
save_pub_r(p_original,file.path(fig_dir,"Bubbleplot_original_marker_panel_recomputed"),183,142,600)

# Pairwise reassessment: robust genes from the full scan, excluding generic
# housekeeping/ribosomal/mitochondrial genes. This is descriptive (one donor).
pair_genes <- c("PDZRN4","PRKG1","BICC1","GPC6","RORA","RUNX2","IGFBP7","MYL9","TAGLN","LUM","DCN","SFRP4","CTHRC1")
pair_genes <- intersect(pair_genes,rownames(mat))
idx_act <- which(labels=="Act OB-Stromal")
idx_les <- which(labels=="LesionStromal")
x_act <- as.matrix(mat[pair_genes,idx_act,drop=FALSE])
x_les <- as.matrix(mat[pair_genes,idx_les,drop=FALSE])
n_bal <- min(length(idx_act),length(idx_les))
set.seed(seed+18L)
B <- 500L
boot_delta <- matrix(NA_real_,nrow=length(pair_genes),ncol=B,dimnames=list(pair_genes,NULL))
for (b in seq_len(B)) {
  ia <- sample(seq_along(idx_act),n_bal,replace=TRUE)
  il <- sample(seq_along(idx_les),n_bal,replace=TRUE)
  boot_delta[,b] <- rowMeans(x_act[,ia,drop=FALSE])-rowMeans(x_les[,il,drop=FALSE])
}
auc_one <- function(a,l){
  v <- c(a,l); y <- c(rep(1,length(a)),rep(0,length(l))); r <- rank(v,ties.method="average")
  (sum(r[y==1])-length(a)*(length(a)+1)/2)/(length(a)*length(l))
}
pair <- data.frame(
  gene=pair_genes,
  mean_remodeling=rowMeans(x_act),mean_lesion=rowMeans(x_les),
  delta_remodeling_minus_lesion=rowMeans(x_act)-rowMeans(x_les),
  ci_low=apply(boot_delta,1,quantile,.025),ci_high=apply(boot_delta,1,quantile,.975),
  pct_remodeling=100*Matrix::rowMeans(counts[pair_genes,idx_act,drop=FALSE]>0),
  pct_lesion=100*Matrix::rowMeans(counts[pair_genes,idx_les,drop=FALSE]>0),
  auc_remodeling=vapply(seq_along(pair_genes),function(i) auc_one(x_act[i,],x_les[i,]),numeric(1)),
  stringsAsFactors=FALSE
)
pair$favored_state <- ifelse(pair$delta_remodeling_minus_lesion>0,"ECM-remodeling stromal","LesionStromal")
pair$directional_auc <- ifelse(pair$favored_state=="ECM-remodeling stromal",pair$auc_remodeling,1-pair$auc_remodeling)
pair$evidence_class <- ifelse(pair$ci_low>0|pair$ci_high<0,"consistent relative enrichment","not stable under balanced bootstrap")
pair <- pair[order(pair$favored_state,-pair$directional_auc),]
readr::write_csv(pair,file.path(res_dir,"lesion_vs_ecm_remodeling_marker_reassessment.csv"))

original_claims <- data.frame(gene=c("WNT5A","FAP","TNFSF11","POSTN","CTHRC1","IBSP"),claimed_state=c(rep("LesionStromal",3),rep("ECM-remodeling stromal",3)),stringsAsFactors=FALSE)
direct <- readr::read_csv(file.path(res_dir,"act_ob_vs_lesion_direct_full_gene_scan.csv.gz"),show_col_types=FALSE)
audit <- dplyr::left_join(original_claims,direct,by="gene") |>
  dplyr::mutate(abs_delta=abs(delta_act_minus_lesion),interpretation=dplyr::case_when(abs_delta<.10~"non-discriminating",abs_delta<.30~"weak relative enrichment",abs_delta<.70~"moderate relative enrichment",TRUE~"strong relative enrichment"))
readr::write_csv(audit,file.path(res_dir,"original_lesion_act_ob_marker_claim_audit.csv"))

pair_order <- c("PDZRN4","PRKG1","BICC1","GPC6","RORA","RUNX2","IGFBP7","MYL9","TAGLN","LUM","DCN","SFRP4","CTHRC1")
pair$gene <- factor(pair$gene,levels=pair_order)
abs_dot <- dplyr::bind_rows(
  data.frame(state="LesionStromal",gene=pair_genes,average_expression=rowMeans(x_les),percent_expressed=100*Matrix::rowMeans(counts[pair_genes,idx_les,drop=FALSE]>0)),
  data.frame(state="ECM-remodeling stromal",gene=pair_genes,average_expression=rowMeans(x_act),percent_expressed=100*Matrix::rowMeans(counts[pair_genes,idx_act,drop=FALSE]>0))
)
abs_dot$gene <- factor(abs_dot$gene,levels=pair_order)
abs_dot$state <- factor(abs_dot$state,levels=c("ECM-remodeling stromal","LesionStromal"))
p_abs <- ggplot(abs_dot,aes(gene,state,size=percent_expressed,colour=average_expression))+geom_point()+scale_size(range=c(.5,6),limits=c(0,100),breaks=c(0,25,50,75,100))+scale_colour_gradientn(colours=c("#FDE0DD","#FC9272","#D7301F"))+labs(x=NULL,y=NULL,size="Percent expressed",colour="Mean log-normalized\nexpression")+theme_cb(7)+theme(axis.text.x=element_text(angle=55,hjust=1,size=6),axis.text.y=element_text(size=6),legend.position="right")
save_pub_r(p_abs,file.path(fig_dir,"Bubbleplot_lesion_vs_ecm_remodeling_absolute"),150,62,600)

p_eff <- ggplot(pair,aes(delta_remodeling_minus_lesion,gene,colour=favored_state))+geom_vline(xintercept=0,linetype=2,linewidth=.35,colour="grey50")+geom_errorbar(aes(xmin=ci_low,xmax=ci_high),width=.18,linewidth=.5,orientation="y")+geom_point(size=2)+scale_colour_manual(values=c("LesionStromal"=cols[["LesionStromal"]],"ECM-remodeling stromal"=cols[["ECM-remodeling stromal"]]))+labs(x="Mean-expression difference (ECM-remodeling minus Lesion)",y=NULL,colour="Relatively enriched state")+theme_cb(7)+theme(legend.position="top",axis.text.y=element_text(face="bold"))
save_pub_r(p_eff,file.path(fig_dir,"Effectplot_lesion_vs_ecm_remodeling_balanced_bootstrap"),120,88,600)

summary <- data.frame(
  metric=c("All expressed-gene Spearman correlation","Variable-gene Spearman correlation","NMI: independent recluster vs legacy state","Act-state cells in dominant independent cluster","CTHRC1 detection: remodeling","CTHRC1 detection: lesion"),
  value=c(0.9254818582782598,0.9414827536249809,0.40720504872324575,567/616,mean(counts["CTHRC1",idx_act]>0),mean(counts["CTHRC1",idx_les]>0)),
  interpretation=c("high transcriptomic similarity","high transcriptomic similarity","moderate concordance","coherent nested subcluster","not state-exclusive","not state-exclusive")
)
readr::write_csv(summary,file.path(res_dir,"lesion_vs_ecm_remodeling_state_summary.csv"))
log_step("Marker specificity reassessment complete")
