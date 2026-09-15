args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
root <- Sys.getenv("CB_REVISION_ROOT")
if (!nzchar(root)) {
  root <- if (length(file_arg)) {
    dirname(dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/")))
  } else {
    normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  }
}

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

obj_path <- Sys.getenv(
  "CB_XIST_OBJECT",
  unset = file.path(root, "objects", "06_base_broad_lineages_reviewed.rds")
)
obj <- readRDS(obj_path)

meta <- obj[[]]
candidate_group_cols <- c("condition", "Condition", "group", "Group", "sample", "Sample", "orig.ident")
present_group_cols <- intersect(candidate_group_cols, colnames(meta))
cat("META_COLUMNS\n")
cat(paste(colnames(meta), collapse = "\t"), "\n")
cat("CANDIDATE_GROUP_VALUES\n")
for (nm in present_group_cols) {
  cat(nm, ":", paste(unique(as.character(meta[[nm]])), collapse = ","), "\n")
}

group_col <- present_group_cols[vapply(present_group_cols, function(nm) {
  vals <- toupper(unique(as.character(meta[[nm]])))
  any(grepl("CB", vals)) && any(grepl("CTRL|CONTROL", vals))
}, logical(1))][1]
if (is.na(group_col)) stop("Could not identify the CB/CTRL metadata column")

DefaultAssay(obj) <- "RNA"
genes <- c("XIST", "RPS4Y1", "DDX3Y", "KDM5D", "UTY", "EIF1AY", "ZFY")
genes <- intersect(genes, rownames(obj))
expr <- GetAssayData(obj, assay = "RNA", layer = "data")[genes, , drop = FALSE]
group <- as.character(meta[[group_col]])

summarize_gene <- function(g) {
  do.call(rbind, lapply(sort(unique(group)), function(s) {
    x <- as.numeric(expr[g, group == s])
    data.frame(
      sample = s,
      gene = g,
      cells = length(x),
      positive_cells = sum(x > 0),
      percent_positive = 100 * mean(x > 0),
      average_all_cells = mean(x),
      average_positive_cells = if (any(x > 0)) mean(x[x > 0]) else NA_real_
    )
  }))
}
summary_all <- do.call(rbind, lapply(genes, summarize_gene))
cat("GROUP_COLUMN\n", group_col, "\n", sep = "")
cat("GENE_SUMMARY\n")
print(summary_all, row.names = FALSE, digits = 5)

lineage_col <- intersect(c("broad_lineage", "Broad_lineage", "lineage", "Lineage"), colnames(meta))[1]
if (!is.na(lineage_col) && "XIST" %in% genes) {
  x <- as.numeric(expr["XIST", ])
  split_key <- interaction(group, meta[[lineage_col]], drop = TRUE, sep = " | ")
  lineage_summary <- do.call(rbind, lapply(levels(split_key), function(k) {
    idx <- split_key == k
    parts <- strsplit(k, " \\| ")[[1]]
    data.frame(sample = parts[1], lineage = parts[2], cells = sum(idx),
               positive_cells = sum(x[idx] > 0),
               percent_positive = 100 * mean(x[idx] > 0),
               average_all_cells = mean(x[idx]))
  }))
  cat("XIST_BY_LINEAGE\n")
  print(lineage_summary, row.names = FALSE, digits = 5)
}
