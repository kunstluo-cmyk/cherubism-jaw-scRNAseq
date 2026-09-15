args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg)) {
  root <- dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/"))
} else {
  root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
Sys.setenv(CB_REVISION_ROOT = root)

if (requireNamespace("renv", quietly = TRUE)) {
  try(renv::load(project = root, quiet = TRUE), silent = TRUE)
}

source(file.path(root, "R", "00_setup.R"), chdir = FALSE)
source(file.path(root, "R", "01_qc_doublets.R"), chdir = FALSE)
source(file.path(root, "R", "02_decontx.R"), chdir = FALSE)
source(file.path(root, "R", "03_build_base_objects.R"), chdir = FALSE)
source(file.path(root, "R", "04_make_qc_figures.R"), chdir = FALSE)
source(file.path(root, "R", "05_validate_outputs.R"), chdir = FALSE)

writeLines(capture.output(sessionInfo()), file.path(root, "logs", "sessionInfo.txt"))

if (requireNamespace("renv", quietly = TRUE)) {
  lock_packages <- c(
    "BiocManager", "future", "Seurat", "SeuratObject", "Matrix",
    "SingleCellExperiment", "SummarizedExperiment", "scDblFinder", "celda",
    "BiocParallel", "ggplot2", "patchwork", "dplyr", "tidyr", "readr",
    "scales", "yaml", "jsonlite", "svglite", "ragg", "renv"
  )
  try(
    renv::snapshot(
      project = root,
      lockfile = file.path(root, "renv.lock"),
      packages = lock_packages,
      prompt = FALSE
    ),
    silent = TRUE
  )
}
source(file.path(root, "R", "06_refresh_manifest.R"), chdir = FALSE)

message("Rebuild completed: ", root)
