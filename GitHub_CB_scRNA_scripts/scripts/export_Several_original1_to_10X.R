#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Rcpp)
  library(SeuratObject)
  library(Matrix)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
project_root <- Sys.getenv("CB_REVISION_ROOT")
if (!nzchar(project_root)) {
  project_root <- if (length(file_arg)) {
    dirname(dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/")))
  } else {
    normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  }
}

input_rdata <- Sys.getenv(
  "CB_INPUT_RDATA",
  unset = file.path(project_root, "data", "Several original1.Rdata")
)
output_root <- Sys.getenv(
  "CB_10X_EXPORT_DIR",
  unset = file.path(project_root, "exports", "Several_original1_10X_upload")
)
gene_info_file <- Sys.getenv(
  "CB_GENE_INFO_FILE",
  unset = file.path(project_root, "resources", "GRCh38-2020-A_geneInfo.tab")
)
gencode_gtf <- Sys.getenv(
  "CB_GENCODE_GTF",
  unset = file.path(project_root, "resources", "gencode.v32.primary_assembly.annotation.gtf.gz")
)

parse_attr <- function(x, key) {
  pattern <- paste0(".*(?:^|;[[:space:]]*)", key, " \\\"([^\\\"]+)\\\".*")
  sub(pattern, "\\1", x, perl = TRUE)
}

read_gencode_gene_map <- function(path) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  ids <- character()
  names <- character()
  repeat {
    lines <- readLines(con, n = 100000L, warn = FALSE)
    if (!length(lines)) break
    lines <- lines[!startsWith(lines, "#")]
    if (!length(lines)) next
    fields <- strsplit(lines, "\t", fixed = TRUE)
    is_gene <- vapply(fields, function(z) length(z) >= 9L && z[[3L]] == "gene", logical(1))
    if (!any(is_gene)) next
    attrs <- vapply(fields[is_gene], `[[`, character(1), 9L)
    gene_id <- sub("\\.[0-9]+$", "", parse_attr(attrs, "gene_id"))
    gene_name <- parse_attr(attrs, "gene_name")
    ids <- c(ids, gene_id)
    names <- c(names, gene_name)
  }
  keep <- !duplicated(ids)
  setNames(names[keep], ids[keep])
}

gzip_copy <- function(source, destination) {
  in_con <- file(source, open = "rb")
  out_con <- gzfile(destination, open = "wb", compression = 6)
  on.exit({
    close(in_con)
    close(out_con)
  }, add = TRUE)
  repeat {
    buffer <- readBin(in_con, what = "raw", n = 1024L * 1024L)
    if (!length(buffer)) break
    writeBin(buffer, out_con)
  }
}

write_gz_lines <- function(lines, path) {
  con <- gzfile(path, open = "wt", compression = 6)
  on.exit(close(con), add = TRUE)
  writeLines(lines, con = con, sep = "\n", useBytes = TRUE)
}

read_mtx_header <- function(path) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  format_line <- readLines(con, n = 1L, warn = FALSE)
  if (!length(format_line)) stop("Incomplete Matrix Market header: ", path)
  repeat {
    dimension_line <- readLines(con, n = 1L, warn = FALSE)
    if (!length(dimension_line)) stop("Missing Matrix Market dimensions: ", path)
    if (!startsWith(dimension_line, "%")) break
  }
  list(
    format = format_line[[1L]],
    dimensions = as.numeric(strsplit(trimws(dimension_line), "[[:space:]]+")[[1L]])
  )
}

dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

gene_info <- readLines(gene_info_file, warn = FALSE)
expected_n <- as.integer(gene_info[[1L]])
gene_ids <- gene_info[-1L]
stopifnot(length(gene_ids) == expected_n, !anyDuplicated(gene_ids))

gene_map <- read_gencode_gene_map(gencode_gtf)
if (anyNA(gene_map[gene_ids])) {
  stop("Not all GRCh38-2020-A gene IDs were found in GENCODE v32.")
}
reference_names <- unname(gene_map[gene_ids])
reference_seurat_names <- make.unique(reference_names)

objects <- new.env(parent = emptyenv())
loaded <- load(input_rdata, envir = objects)
required <- c("CB", "CTRL")
if (!all(required %in% loaded)) stop("CB and/or CTRL object missing from input RData.")

sample_map <- c(CB = "CB", CTRL = "Ctrl")
summary_rows <- list()

for (object_name in names(sample_map)) {
  sample_name <- unname(sample_map[[object_name]])
  object <- objects[[object_name]]
  counts <- LayerData(object[["RNA"]], layer = "counts")

  stopifnot(inherits(counts, "sparseMatrix"))
  stopifnot(nrow(counts) == expected_n)
  object_reference_names <- sub("\\.ENSG[0-9]+$", "", rownames(counts))
  if (!identical(object_reference_names, reference_names)) {
    mismatch <- which(object_reference_names != reference_names)
    preview <- head(mismatch, 20L)
    message("Reference-name mismatch count for ", sample_name, ": ", length(mismatch))
    message(paste(
      preview,
      rownames(counts)[preview],
      reference_names[preview],
      gene_ids[preview],
      sep = "\t",
      collapse = "\n"
    ))
    stop("Seurat feature order/names do not match the proposed 10x reference after duplicate-name normalization.")
  }
  stopifnot(!anyDuplicated(colnames(counts)))
  stopifnot(all(is.finite(counts@x)), all(counts@x >= 0), all(counts@x == round(counts@x)))

  sample_dir <- file.path(output_root, sample_name)
  dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)

  matrix_path <- file.path(sample_dir, "matrix.mtx.gz")
  barcode_path <- file.path(sample_dir, "barcodes.tsv.gz")
  feature_path <- file.path(sample_dir, "features.tsv.gz")
  matrix_tmp <- file.path(sample_dir, "matrix.mtx")

  writeMM(counts, file = matrix_tmp)
  gzip_copy(matrix_tmp, matrix_path)
  unlink(matrix_tmp)

  write_gz_lines(colnames(counts), barcode_path)
  feature_lines <- paste(gene_ids, reference_names, "Gene Expression", sep = "\t")
  write_gz_lines(feature_lines, feature_path)

  header <- read_mtx_header(matrix_path)
  stopifnot(grepl("^%%MatrixMarket matrix coordinate integer general$", header$format))
  stopifnot(all(header$dimensions == c(nrow(counts), ncol(counts), length(counts@x))))
  stopifnot(length(readLines(gzfile(barcode_path), warn = FALSE)) == ncol(counts))
  stopifnot(length(readLines(gzfile(feature_path), warn = FALSE)) == nrow(counts))

  ncount <- Matrix::colSums(counts)
  nfeature <- Matrix::colSums(counts > 0)
  metadata <- object[[]]
  count_match <- if ("nCount_RNA" %in% colnames(metadata)) all(ncount == metadata$nCount_RNA) else NA
  feature_match <- if ("nFeature_RNA" %in% colnames(metadata)) all(nfeature == metadata$nFeature_RNA) else NA
  stopifnot(isTRUE(count_match), isTRUE(feature_match))

  hashes <- tools::md5sum(c(matrix_path, barcode_path, feature_path))
  summary_rows[[sample_name]] <- data.frame(
    sample = sample_name,
    features = nrow(counts),
    cells = ncol(counts),
    nonzero_entries = length(counts@x),
    total_UMIs = sum(counts@x),
    matrix_md5 = unname(hashes[[matrix_path]]),
    barcodes_md5 = unname(hashes[[barcode_path]]),
    features_md5 = unname(hashes[[feature_path]]),
    stringsAsFactors = FALSE
  )
}

summary_table <- do.call(rbind, summary_rows)
write.table(
  summary_table,
  file = file.path(output_root, "10X_export_manifest.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

readme <- c(
  "Several original1 reconstructed 10X-compatible filtered count matrices",
  "",
  "Source: local file data/Several original1.Rdata",
  "Samples: CB and Ctrl",
  "Matrix source: integer-valued RNA counts layer stored in the original Seurat objects.",
  "Reference feature order: 10x Genomics refdata-gex-GRCh38-2020-A (GENCODE v32 / Ensembl 98).",
  "Files per sample: matrix.mtx.gz, barcodes.tsv.gz, features.tsv.gz.",
  "Cell barcode suffixes (_3 for CB and _4 for Ctrl) are preserved exactly as stored in the source objects.",
  "",
  "Important scope statement:",
  "These are reconstructed filtered feature-barcode matrices containing the cells retained in the source Seurat objects.",
  "They are not FASTQ files and are not unfiltered droplet matrices containing empty droplets.",
  "For GEO, these files should normally be described as processed 10X-compatible filtered raw-count matrices."
)
writeLines(readme, file.path(output_root, "README.txt"), useBytes = TRUE)

print(summary_table)
cat("\nExport complete: ", normalizePath(output_root, winslash = "/"), "\n", sep = "")
