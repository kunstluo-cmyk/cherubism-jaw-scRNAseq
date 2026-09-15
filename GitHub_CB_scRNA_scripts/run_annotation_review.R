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
source(file.path(root, "R", "11_full_gene_scan_finalize_lineages.R"), chdir = FALSE)
source(file.path(root, "R", "12_stromal_state_revalidation.R"), chdir = FALSE)
source(file.path(root, "R", "15_reviewer_targeted_expression_audit.R"), chdir = FALSE)
source(file.path(root, "R", "13_make_global_annotation_figures.R"), chdir = FALSE)
source(file.path(root, "R", "16_make_publication_subfigures.R"), chdir = FALSE)
source(file.path(root, "R", "14_validate_annotation_outputs.R"), chdir = FALSE)
source(file.path(root, "R", "17_validate_publication_subfigures.R"), chdir = FALSE)
source(file.path(root, "R", "18_marker_specificity_reassessment.R"), chdir = FALSE)
source(file.path(root, "R", "19_validate_marker_reassessment.R"), chdir = FALSE)
source(file.path(root, "R", "06_refresh_manifest.R"), chdir = FALSE)

writeLines(capture.output(sessionInfo()), file.path(root, "logs", "sessionInfo_annotation_review.txt"))
message("Annotation review completed: ", root)
