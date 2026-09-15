root <- Sys.getenv("CB_REVISION_ROOT")
if (!nzchar(root)) root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
config <- yaml::read_yaml(file.path(root, "config", "analysis_config.yml"))
seed <- as.integer(config$seed)
set.seed(seed)
log_file <- file.path(root, "logs", "pipeline.log")
log_step <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(stringr)
  library(purrr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(KEGGREST)
  library(BiocParallel)
})

BiocParallel::register(BiocParallel::SerialParam())
select <- dplyr::select

log_step("Stage 28: pooled OC/TREM2+ Mac/Resident Mac pairwise DE and enrichment")

result_dir <- file.path(root, "results", "enrichment", "pooled_myeloid_pairwise")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

object_path <- file.path(root, "objects", "06_base_broad_lineages_reviewed.rds")
if (!file.exists(object_path)) stop("Missing reviewed object: ", object_path)

obj <- readRDS(object_path)
DefaultAssay(obj) <- "RNA"

celltype_col <- "celltype_reviewed_provisional"
condition_col <- "group"
states <- c("Resident Mac", "TREM2+ Mac", "OC")

if (!all(c(celltype_col, condition_col) %in% colnames(obj@meta.data))) {
  stop("Required metadata fields are absent.")
}

pooled_cells <- colnames(obj)[obj@meta.data[[celltype_col]] %in% states]
myeloid <- subset(obj, cells = pooled_cells)
myeloid[[celltype_col]] <- factor(myeloid@meta.data[[celltype_col]], levels = states)
Idents(myeloid) <- celltype_col

cell_counts <- tibble(
  cell_state = as.character(myeloid@meta.data[[celltype_col]])
) |>
  count(cell_state, name = "n_cells") |>
  mutate(
    cell_state = as.character(cell_state),
    specimen_scope = "Pooled CB and CTRL",
    biological_replicates = 2L
  ) |>
  arrange(match(cell_state, states))

write_csv(cell_counts, file.path(result_dir, "cell_counts_pooled_myeloid_states.csv"))

condition_composition <- tibble(
  cell_state = as.character(myeloid@meta.data[[celltype_col]]),
  condition = as.character(myeloid@meta.data[[condition_col]])
) |>
  count(cell_state, condition, name = "n_cells") |>
  group_by(cell_state) |>
  mutate(fraction_within_state = n_cells / sum(n_cells)) |>
  ungroup() |>
  arrange(match(cell_state, states), condition)
write_csv(condition_composition, file.path(result_dir, "condition_composition_by_cell_state.csv"))

if (any(!states %in% cell_counts$cell_state)) stop("At least one requested pooled cell state is absent.")
if (any(cell_counts$n_cells < 3)) stop("At least one requested pooled cell state has fewer than three cells.")

meta_cols <- intersect(
  c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent_mt", "decontX_contamination"),
  colnames(myeloid@meta.data)
)
if (length(meta_cols)) {
  technical_summary <- myeloid@meta.data |>
    rownames_to_column("cell") |>
    transmute(cell_state = .data[[celltype_col]], across(all_of(meta_cols))) |>
    group_by(cell_state) |>
    summarise(
      n_cells = n(),
      across(
        all_of(meta_cols),
        list(median = ~ median(.x, na.rm = TRUE), q25 = ~ quantile(.x, 0.25, na.rm = TRUE),
             q75 = ~ quantile(.x, 0.75, na.rm = TRUE)),
        .names = "{.col}_{.fn}"
      ),
      .groups = "drop"
    )
  write_csv(technical_summary, file.path(result_dir, "technical_covariate_summary.csv"))
}

comparison_table <- tribble(
  ~comparison_id, ~ident_1, ~ident_2, ~comparison_label,
  "TREM2_vs_Resident", "TREM2+ Mac", "Resident Mac", "TREM2+ Mac vs Resident Mac",
  "OC_vs_Resident", "OC", "Resident Mac", "OC vs Resident Mac",
  "OC_vs_TREM2", "OC", "TREM2+ Mac", "OC vs TREM2+ Mac"
)
write_csv(comparison_table, file.path(result_dir, "comparison_definitions.csv"))

run_de <- function(ident_1, ident_2, comparison_id, comparison_label) {
  log_step("DE: ", comparison_label)
  res <- FindMarkers(
    object = myeloid,
    ident.1 = ident_1,
    ident.2 = ident_2,
    assay = "RNA",
    slot = "data",
    test.use = "wilcox",
    logfc.threshold = 0,
    min.pct = 0.05,
    min.cells.group = 3,
    only.pos = FALSE,
    densify = FALSE,
    verbose = FALSE
  ) |>
    rownames_to_column("gene") |>
    mutate(
      comparison_id = comparison_id,
      comparison = comparison_label,
      ident_1 = ident_1,
      ident_2 = ident_2,
      higher_in = case_when(
        avg_log2FC > 0 ~ ident_1,
        avg_log2FC < 0 ~ ident_2,
        TRUE ~ "Tie"
      ),
      candidate_state_DEG =
        p_val_adj < 0.05 & abs(avg_log2FC) >= 0.25 & pmax(pct.1, pct.2) >= 0.10
    ) |>
    arrange(p_val_adj, desc(abs(avg_log2FC)), gene)
  res
}

de_list <- pmap(
  comparison_table,
  function(comparison_id, ident_1, ident_2, comparison_label) {
    run_de(ident_1, ident_2, comparison_id, comparison_label)
  }
)
names(de_list) <- comparison_table$comparison_id
de_all <- bind_rows(de_list)

walk2(
  de_list,
  names(de_list),
  ~ write_csv(.x, file.path(result_dir, paste0("DE_", .y, ".csv.gz")))
)
write_csv(de_all, file.path(result_dir, "DE_all_pairwise.csv.gz"))

de_summary <- de_all |>
  group_by(comparison_id, comparison, ident_1, ident_2) |>
  summarise(
    n_tested = n(),
    n_candidate_ident_1 = sum(candidate_state_DEG & avg_log2FC > 0, na.rm = TRUE),
    n_candidate_ident_2 = sum(candidate_state_DEG & avg_log2FC < 0, na.rm = TRUE),
    .groups = "drop"
  )
write_csv(de_summary, file.path(result_dir, "DE_summary.csv"))

# One stable SYMBOL-to-ENTREZ map is used for DE gene lists, rank vectors and pathway genes.
all_symbols <- unique(de_all$gene)
symbol_to_entrez <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = all_symbols,
  keytype = "SYMBOL",
  column = "ENTREZID",
  multiVals = "first"
)
id_map <- tibble(
  gene = names(symbol_to_entrez),
  ENTREZID = unname(symbol_to_entrez)
) |>
  filter(!is.na(ENTREZID), nzchar(ENTREZID)) |>
  distinct(gene, .keep_all = TRUE)
write_csv(id_map, file.path(result_dir, "symbol_entrez_mapping.csv"))
entrez_symbol_lookup <- id_map |>
  arrange(ENTREZID, gene) |>
  distinct(ENTREZID, .keep_all = TRUE) |>
  select(ENTREZID, gene) |>
  deframe()

# Cache the current official KEGG pathway-to-gene mapping to make the analysis auditable.
kegg_t2g_path <- file.path(result_dir, "KEGG_hsa_TERM2GENE_cached.csv.gz")
kegg_t2n_path <- file.path(result_dir, "KEGG_hsa_TERM2NAME_cached.csv")
if (file.exists(kegg_t2g_path) && file.exists(kegg_t2n_path)) {
  kegg_t2g <- read_csv(kegg_t2g_path, show_col_types = FALSE)
  kegg_t2n <- read_csv(kegg_t2n_path, show_col_types = FALSE)
} else {
  log_step("Retrieving and caching current Homo sapiens KEGG pathway definitions")
  kegg_links <- KEGGREST::keggLink("pathway", "hsa")
  kegg_names <- KEGGREST::keggList("pathway", "hsa")
  kegg_t2g <- tibble(
    ENTREZID = sub("^hsa:", "", names(kegg_links)),
    pathway_id = sub("^path:", "", unname(kegg_links))
  ) |>
    distinct(pathway_id, ENTREZID)
  kegg_t2n <- tibble(
    pathway_id = names(kegg_names),
    Description = sub(" - Homo sapiens \\(human\\)$", "", unname(kegg_names))
  ) |>
    distinct(pathway_id, .keep_all = TRUE)
  write_csv(kegg_t2g, kegg_t2g_path)
  write_csv(kegg_t2n, kegg_t2n_path)
}

split_ids <- function(x) {
  if (length(x) == 0 || is.na(x) || !nzchar(x)) return(character())
  unique(strsplit(x, "/", fixed = TRUE)[[1]])
}

entrez_to_symbols <- function(ids) {
  ids <- unique(ids[nzchar(ids)])
  if (!length(ids)) return(character())
  out <- unname(entrez_symbol_lookup[ids])
  out[is.na(out) | !nzchar(out)] <- ids[is.na(out) | !nzchar(out)]
  unique(out)
}

collapse_symbols <- function(x) paste(entrez_to_symbols(split_ids(x)), collapse = "/")

prepare_entrez_de <- function(de_df) {
  de_df |>
    inner_join(id_map, by = "gene") |>
    filter(is.finite(avg_log2FC)) |>
    arrange(ENTREZID, desc(abs(avg_log2FC)), p_val_adj, gene) |>
    distinct(ENTREZID, .keep_all = TRUE)
}

make_rank_vector <- function(de_entrez) {
  rank_df <- de_entrez |>
    filter(abs(avg_log2FC) > 1e-10) |>
    arrange(desc(avg_log2FC), ENTREZID) |>
    mutate(
      # The deterministic perturbation only resolves exact numerical ties.
      rank_metric = avg_log2FC + seq(from = 1e-12, by = 1e-12, length.out = n())
    )
  ranks <- rank_df$rank_metric
  names(ranks) <- rank_df$ENTREZID
  sort(ranks, decreasing = TRUE)
}

standardise_gsea <- function(gsea_obj, database, comp_row) {
  df <- as.data.frame(gsea_obj)
  if (!nrow(df)) return(tibble())
  core_col <- if ("core_enrichment" %in% names(df)) "core_enrichment" else "leading_edge"
  df |>
    as_tibble() |>
    mutate(
      database = database,
      comparison_id = comp_row$comparison_id,
      comparison = comp_row$comparison_label,
      ident_1 = comp_row$ident_1,
      ident_2 = comp_row$ident_2,
      enriched_in = if_else(NES >= 0, comp_row$ident_1, comp_row$ident_2),
      leading_entrez = .data[[core_col]],
      leading_genes = map_chr(.data[[core_col]], collapse_symbols),
      Count = map_int(.data[[core_col]], ~ length(split_ids(.x))),
      GeneRatio_numeric = if_else(setSize > 0, Count / setSize, NA_real_),
      GeneRatio = paste0(Count, "/", setSize),
      enrichment_method = "Ranked GSEA (fgsea backend)",
      rank_metric = "avg_log2FC with deterministic tie breaking"
    ) |>
    select(
      database, comparison_id, comparison, ident_1, ident_2, enriched_in,
      ID, Description, setSize, enrichmentScore, NES, pvalue, p.adjust, qvalue,
      Count, GeneRatio, GeneRatio_numeric, leading_genes, leading_entrez,
      enrichment_method, rank_metric, everything()
    )
}

run_ranked_enrichment <- function(de_df, comp_row) {
  log_step("Ranked GO/KEGG enrichment: ", comp_row$comparison_label)
  de_entrez <- prepare_entrez_de(de_df)
  ranks <- make_rank_vector(de_entrez)

  set.seed(seed)
  go_obj <- suppressWarnings(gseGO(
    geneList = ranks,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    minGSSize = 10,
    maxGSSize = 500,
    exponent = 1,
    eps = 1e-10,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    verbose = FALSE,
    seed = TRUE,
    by = "fgsea"
  ))

  set.seed(seed)
  kegg_obj <- suppressWarnings(GSEA(
    geneList = ranks,
    TERM2GENE = kegg_t2g[, c("pathway_id", "ENTREZID")],
    TERM2NAME = kegg_t2n[, c("pathway_id", "Description")],
    minGSSize = 10,
    maxGSSize = 500,
    exponent = 1,
    eps = 1e-10,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    verbose = FALSE,
    seed = TRUE,
    by = "fgsea"
  ))

  list(
    GO = standardise_gsea(go_obj, "GO Biological Process", comp_row),
    KEGG = standardise_gsea(kegg_obj, "KEGG", comp_row),
    ranks = tibble(ENTREZID = names(ranks), rank_metric = unname(ranks)) |>
      left_join(de_entrez |> select(ENTREZID, gene, avg_log2FC, pct.1, pct.2, p_val_adj), by = "ENTREZID")
  )
}

run_directional_ora <- function(de_df, comp_row) {
  de_entrez <- prepare_entrez_de(de_df)
  universe <- unique(de_entrez$ENTREZID)
  sig <- de_entrez |>
    filter(candidate_state_DEG)

  direction_defs <- list(
    ident_1 = sig |> filter(avg_log2FC > 0) |> pull(ENTREZID) |> unique(),
    ident_2 = sig |> filter(avg_log2FC < 0) |> pull(ENTREZID) |> unique()
  )
  direction_states <- c(ident_1 = comp_row$ident_1, ident_2 = comp_row$ident_2)

  run_one <- function(gene_ids, direction_key) {
    if (length(gene_ids) < 10) return(list(GO = tibble(), KEGG = tibble()))
    go <- suppressWarnings(enrichGO(
      gene = gene_ids,
      universe = universe,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = "BP",
      minGSSize = 10,
      maxGSSize = 500,
      pvalueCutoff = 1,
      qvalueCutoff = 1,
      pAdjustMethod = "BH",
      readable = FALSE
    ))
    kegg <- suppressWarnings(enricher(
      gene = gene_ids,
      universe = universe,
      TERM2GENE = kegg_t2g[, c("pathway_id", "ENTREZID")],
      TERM2NAME = kegg_t2n[, c("pathway_id", "Description")],
      minGSSize = 10,
      maxGSSize = 500,
      pvalueCutoff = 1,
      qvalueCutoff = 1,
      pAdjustMethod = "BH"
    ))

    standardise <- function(x, database) {
      xdf <- as.data.frame(x)
      if (!nrow(xdf)) return(tibble())
      as_tibble(xdf) |>
        mutate(
          database = database,
          comparison_id = comp_row$comparison_id,
          comparison = comp_row$comparison_label,
          ident_1 = comp_row$ident_1,
          ident_2 = comp_row$ident_2,
          enriched_in = unname(direction_states[[direction_key]]),
          contributing_entrez = geneID,
          contributing_genes = map_chr(geneID, collapse_symbols),
          GeneRatio_numeric = map_dbl(GeneRatio, ~ {
            z <- as.numeric(strsplit(.x, "/", fixed = TRUE)[[1]])
            z[1] / z[2]
          }),
          enrichment_method = "Direction-specific over-representation analysis",
          DEG_definition = "cell-level BH-adjusted P < 0.05, |avg_log2FC| >= 0.25, max detection >= 10%"
        ) |>
        select(
          database, comparison_id, comparison, ident_1, ident_2, enriched_in,
          ID, Description, GeneRatio, GeneRatio_numeric, BgRatio, RichFactor,
          FoldEnrichment, zScore, pvalue, p.adjust, qvalue, Count,
          contributing_genes, contributing_entrez, enrichment_method,
          DEG_definition, everything()
        )
    }

    list(
      GO = standardise(go, "GO Biological Process"),
      KEGG = standardise(kegg, "KEGG")
    )
  }

  imap(direction_defs, run_one)
}

gsea_list <- list()
ora_list <- list()
for (i in seq_len(nrow(comparison_table))) {
  comp_row <- comparison_table[i, ]
  cid <- comp_row$comparison_id
  gsea_list[[cid]] <- run_ranked_enrichment(de_list[[cid]], comp_row)
  ora_list[[cid]] <- run_directional_ora(de_list[[cid]], comp_row)

  write_csv(gsea_list[[cid]]$GO, file.path(result_dir, paste0("GSEA_GO_", cid, ".csv.gz")))
  write_csv(gsea_list[[cid]]$KEGG, file.path(result_dir, paste0("GSEA_KEGG_", cid, ".csv.gz")))
  write_csv(gsea_list[[cid]]$ranks, file.path(result_dir, paste0("ranked_genes_", cid, ".csv.gz")))

  ora_go <- bind_rows(map(ora_list[[cid]], "GO"))
  ora_kegg <- bind_rows(map(ora_list[[cid]], "KEGG"))
  write_csv(ora_go, file.path(result_dir, paste0("ORA_GO_", cid, ".csv.gz")))
  write_csv(ora_kegg, file.path(result_dir, paste0("ORA_KEGG_", cid, ".csv.gz")))
}

gsea_go_all <- bind_rows(map(gsea_list, "GO"))
gsea_kegg_all <- bind_rows(map(gsea_list, "KEGG"))
ora_go_all <- bind_rows(map(ora_list, ~ bind_rows(map(.x, "GO"))))
ora_kegg_all <- bind_rows(map(ora_list, ~ bind_rows(map(.x, "KEGG"))))

write_csv(gsea_go_all, file.path(result_dir, "GSEA_GO_all_pairwise.csv.gz"))
write_csv(gsea_kegg_all, file.path(result_dir, "GSEA_KEGG_all_pairwise.csv.gz"))
write_csv(ora_go_all, file.path(result_dir, "ORA_GO_all_pairwise.csv.gz"))
write_csv(ora_kegg_all, file.path(result_dir, "ORA_KEGG_all_pairwise.csv.gz"))

analysis_manifest <- tibble(
  item = c(
    "analysis_scope", "specimen", "biological_replicates", "expression_assay",
    "DE_method", "DE_interpretation", "primary_enrichment", "sensitivity_enrichment",
    "GO_domain", "KEGG_source", "KEGG_retrieval_date", "multiple_testing",
    "primary_limitation"
  ),
  value = c(
    "Pairwise comparisons among all retained Resident Mac, TREM2+ Mac and OC cells pooled across CB and CTRL specimens",
    "One CB and one CTRL specimen pooled without condition stratification", "2 specimens; not replicated within condition", "DecontX-corrected, log-normalized unintegrated RNA",
    "Seurat FindMarkers; cell-level Wilcoxon rank-sum; min.pct=0.05; logfc.threshold=0",
    "Exploratory state-associated transcript screen; cells are not biological replicates",
    "Ranked GSEA using avg_log2FC and fgsea backend",
    "Direction-specific ORA using candidate state-associated genes",
    "GO Biological Process only", "Current Homo sapiens KEGG definitions retrieved with KEGGREST",
    as.character(Sys.Date()), "Benjamini-Hochberg within each enrichment analysis",
    "Cell-state and specimen/condition effects are partly confounded because state composition differs strongly between the two specimens"
  )
)
write_csv(analysis_manifest, file.path(result_dir, "analysis_manifest.csv"))

saveRDS(
  list(
    cell_counts = cell_counts,
    comparison_table = comparison_table,
    de = de_list,
    gsea = gsea_list,
    ora = ora_list,
    kegg_term2gene = kegg_t2g,
    kegg_term2name = kegg_t2n,
    manifest = analysis_manifest
  ),
  file.path(result_dir, "pooled_myeloid_pairwise_enrichment_results.rds"),
  compress = "gzip"
)

log_step(
  "Stage 28 complete; cells=", paste(cell_counts$cell_state, cell_counts$n_cells, sep = ":", collapse = ", "),
  "; GO GSEA rows=", nrow(gsea_go_all), "; KEGG GSEA rows=", nrow(gsea_kegg_all)
)
