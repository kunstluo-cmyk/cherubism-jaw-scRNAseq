root <- Sys.getenv("CB_REVISION_ROOT")
if (!nzchar(root)) root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

bundle_files <- list.files(
  root,
  pattern = "[.](R|csv|csv[.]gz|yml|md|rds|svg|pdf|tiff|png|txt|lock)$",
  recursive = TRUE,
  full.names = TRUE
)
bundle_files <- bundle_files[!grepl("/renv/library/", gsub("\\\\", "/", bundle_files))]
bundle_files <- bundle_files[basename(bundle_files) != "output_manifest.csv"]
info <- file.info(bundle_files)
manifest <- data.frame(
  path = substring(normalizePath(bundle_files, winslash = "/"), nchar(root) + 2),
  size_bytes = info$size,
  md5 = unname(tools::md5sum(bundle_files))
)
write.csv(
  manifest,
  file.path(root, "results", "qc", "output_manifest.csv"),
  row.names = FALSE
)
