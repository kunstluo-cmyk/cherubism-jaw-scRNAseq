root <- Sys.getenv("CB_REVISION_ROOT")
if (!nzchar(root)) stop("CB_REVISION_ROOT is not set. Run run_all.R from the repository root.")

dirs <- c(
  "objects", "results/qc", "figures/qc", "logs", "config", "docs", "R"
)
for (d in dirs) dir.create(file.path(root, d), recursive = TRUE, showWarnings = FALSE)

required <- c(
  "Seurat", "SeuratObject", "Matrix", "SingleCellExperiment", "scDblFinder",
  "celda", "BiocParallel", "ggplot2", "patchwork", "dplyr", "tidyr",
  "readr", "yaml", "jsonlite", "svglite", "ragg"
  , "RANN", "ggrastr", "ggrepel", "tibble"
)
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop(
    "Missing packages: ", paste(missing, collapse = ", "),
    ". Run 00_install_dependencies.R from the repository root first."
  )
}

config <- yaml::read_yaml(file.path(root, "config", "analysis_config.yml"))
seed <- as.integer(config$seed)
set.seed(seed)
options(future.globals.maxSize = 8 * 1024^3)
if (requireNamespace("future", quietly = TRUE)) future::plan("sequential")

log_file <- file.path(root, "logs", "pipeline.log")
log_step <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

safe_mad <- function(x) {
  out <- stats::mad(x, na.rm = TRUE)
  if (!is.finite(out) || out == 0) out <- stats::IQR(x, na.rm = TRUE) / 1.349
  if (!is.finite(out) || out == 0) out <- 1
  out
}

read_condition_colors <- function() {
  p <- readr::read_csv(file.path(root, "config", "palette.csv"), show_col_types = FALSE)
  p <- p[p$category == "condition", ]
  stats::setNames(p$color, p$label)
}

theme_cb <- function(base_size = 7) {
  ggplot2::theme_classic(base_size = base_size, base_family = "Arial") +
    ggplot2::theme(
      axis.line = ggplot2::element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = ggplot2::element_line(linewidth = 0.35, colour = "black"),
      axis.text = ggplot2::element_text(colour = "black"),
      legend.title = ggplot2::element_text(size = base_size - 0.3),
      legend.text = ggplot2::element_text(size = base_size - 0.7),
      strip.text = ggplot2::element_text(size = base_size, face = "bold"),
      plot.title = ggplot2::element_text(size = base_size + 0.5, face = "bold"),
      plot.tag = ggplot2::element_text(size = 8, face = "bold"),
      panel.grid = ggplot2::element_blank()
    )
}

save_pub_r <- function(plot, filename, width_mm, height_mm, dpi = 600) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4

  svglite::svglite(paste0(filename, ".svg"), width = width_in, height = height_in)
  print(plot)
  grDevices::dev.off()

  grDevices::cairo_pdf(paste0(filename, ".pdf"), width = width_in, height = height_in, family = "Arial")
  print(plot)
  grDevices::dev.off()

  ragg::agg_tiff(
    paste0(filename, ".tiff"), width = width_in, height = height_in,
    units = "in", res = dpi, compression = "lzw"
  )
  print(plot)
  grDevices::dev.off()

  ragg::agg_png(
    paste0(filename, ".png"), width = width_in, height = height_in,
    units = "in", res = 200
  )
  print(plot)
  grDevices::dev.off()
}

log_step("Setup complete; seed=", seed, "; future plan=sequential")
