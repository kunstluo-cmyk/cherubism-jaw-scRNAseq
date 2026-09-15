args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
root <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/"))
} else {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

if (!file.exists(file.path(root, "renv.lock"))) {
  stop("renv.lock was not found. Run this script from the repository root.")
}

if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv", repos = "https://cloud.r-project.org")
}

renv::restore(project = root, lockfile = file.path(root, "renv.lock"), prompt = FALSE)
message("Dependency restoration complete for: ", root)
