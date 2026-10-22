# ---- 09_pathways.R ----
# Preranked fgsea using the signed DESeq2 Wald statistic for every eligible gene.
# Gene sets are explicit, versioned GMT inputs. Missing collections are reported
# as skipped; no gene sets are downloaded or invented during a pipeline run.
source(file.path(dirname(sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value = TRUE))), '..', 'R', 'utils.R'))
cfg <- step_setup('09_pathways')
settings <- cfg$pathways
result_template <- tibble(subset = character(), contrast = character(), pathway = character(),
  pval = double(), padj = double(), log2err = double(), ES = double(), NES = double(),
  size = integer(), leading_edge = character())
coverage_template <- tibble(subset = character(), contrast = character(), pathway = character(),
  n_collection_genes = integer(), n_ranked_overlap = integer(), fraction_ranked = double(), eligible = logical())
status_template <- tibble(subset = character(), contrast = character(), status = character(), reason = character(),
                          n_genes_ranked = integer(), n_zero_statistics = integer(), n_tied_statistics = integer(), warnings = character())
results <- list(result_template); coverage_rows <- list(coverage_template); statuses <- list(status_template)
# Remove stale figures if a rerun has no configured/testable collection.
unlink(list.files(file.path(cfg$step_dir, 'figures'), pattern = '^pathway_.*\\.(png|pdf)$', full.names = TRUE))
if (is.null(settings$gmt_file)) {
  unlink(file.path(cfg$step_dir, 'gene_set_manifest.yaml'))
  statuses[[2]] <- tibble(subset = 'all', contrast = 'all', status = 'skipped',
    reason = 'No versioned GMT collection configured', n_genes_ranked = 0L, n_zero_statistics = 0L, n_tied_statistics = 0L)
} else {
  if (!requireNamespace('fgsea', quietly = TRUE)) stop('Install the pinned fgsea package in the configured R library')
  pathways <- read_gmt(settings$gmt_file)
  manifest <- list(collection = settings$collection, version = settings$version,
    gene_id_type = settings$gene_id_type, file = settings$gmt_file,
    md5 = unname(tools::md5sum(settings$gmt_file)), rank_statistic = 'signed DESeq2 Wald statistic',
    gene_universe = 'finite-statistic genes after pseudobulk DE filtering; no implicit ID conversion',
    fdr_family = 'pathways within each subset/contrast', fgsea_version = as.character(packageVersion('fgsea')))
  yaml::write_yaml(manifest, file.path(cfg$step_dir, 'gene_set_manifest.yaml'))
  input <- file.path(step_dir(cfg, '06_pseudobulk_de'), 'tables', 'pseudobulk_de_results.csv.gz')
  if (!file.exists(input)) stop('Missing DE results; run step 06 first: ', input)
  manifest$de_results_md5 <- unname(tools::md5sum(input))
  yaml::write_yaml(manifest, file.path(cfg$step_dir, 'gene_set_manifest.yaml'))
  de_results <- read_csv(input, show_col_types = FALSE)
  groups <- de_results |> distinct(subset, contrast)
  if (!nrow(groups)) statuses[[2]] <- tibble(subset = 'all', contrast = 'all', status = 'skipped',
    reason = 'DE produced no testable comparisons', n_genes_ranked = 0L, n_zero_statistics = 0L, n_tied_statistics = 0L)
  for (i in seq_len(nrow(groups))) {
    sub <- groups$subset[i]; comparison <- groups$contrast[i]
    ranks <- pathway_ranks(de_results |> filter(subset == sub, contrast == comparison))
    overlap <- vapply(pathways, function(g) sum(g %in% names(ranks)), integer(1))
    coverage <- tibble(subset = sub, contrast = comparison, pathway = names(pathways),
      n_collection_genes = lengths(pathways), n_ranked_overlap = overlap,
      fraction_ranked = overlap / lengths(pathways),
      eligible = overlap >= settings$min_size & overlap <= settings$max_size)
    coverage_rows[[length(coverage_rows) + 1L]] <- coverage
    reason <- if (length(ranks) < 2) 'Fewer than two finite ranked genes' else
      if (!any(coverage$eligible)) 'No pathways meet overlap/size limits' else
        if (!any(ranks > 0) || !any(ranks < 0)) 'Two-sided enrichment requires both positive and negative statistics' else ''
    status <- tibble(subset = sub, contrast = comparison, status = if (nzchar(reason)) 'skipped' else 'complete', reason,
      n_genes_ranked = length(ranks), n_zero_statistics = sum(ranks == 0), n_tied_statistics = sum(duplicated(ranks)))
    if (!nzchar(reason)) {
      set.seed(cfg$seed)
      captured_warnings <- character()
      fit <- withCallingHandlers(fgsea::fgseaMultilevel(pathways = pathways[coverage$eligible], stats = ranks,
        minSize = settings$min_size, maxSize = settings$max_size, eps = settings$eps, nproc = 1),
        warning = function(w) captured_warnings <<- c(captured_warnings, conditionMessage(w))) |> as_tibble()
      status$warnings <- paste(unique(captured_warnings), collapse = ' | ')
      fit <- fit |> mutate(subset = sub, contrast = comparison,
                          leading_edge = map_chr(leadingEdge, ~ paste(.x, collapse = ';'))) |>
        select(subset, contrast, pathway, pval, padj, log2err, ES, NES, size, leading_edge)
      results[[length(results) + 1L]] <- fit
      if (nrow(fit)) {
        top <- fit |> filter(is.finite(NES), !is.na(padj)) |> slice_min(padj, n = 20, with_ties = FALSE)
        if (nrow(top)) {
          p <- ggplot(top, aes(NES, reorder(pathway, NES), colour = padj < 0.05, size = size)) +
            geom_point() + scale_colour_manual(values = c('TRUE' = '#B2182B', 'FALSE' = 'grey60')) +
            theme_classic() + labs(title = paste('Pathway enrichment:', sub, comparison),
              subtitle = 'Shown by FDR rank; colours distinguish FDR < 0.05', y = NULL, colour = 'FDR < 0.05')
          stem <- paste0('pathway_', i)
          ggsave(fig_path(cfg, paste0(stem, '.png')), p, width = 12, height = 8, dpi = 100)
          ggsave(fig_path(cfg, paste0(stem, '.pdf')), p, width = 12, height = 8)
        }
      }
    }
    statuses[[length(statuses) + 1L]] <- status
  }
}
results_tbl <- bind_rows(results) |> mutate(cluster = sub('^cluster_', '', subset))
annotation_cfg <- cfg; annotation_cfg$cluster_resolution <- cfg$clustering$main_resolution
annotation_file <- file.path(step_dir(annotation_cfg, '08_annotation'), 'tables', 'cluster_annotations.csv')
if (file.exists(annotation_file)) {
  labels <- read_csv(annotation_file, col_types = cols(cluster = col_character()), show_col_types = FALSE) |>
    select(cluster, broad_type)
  results_tbl <- results_tbl |> left_join(labels, by = 'cluster')
} else results_tbl$broad_type <- NA_character_
results_tbl <- results_tbl |> mutate(broad_type = if_else(subset == 'all_bins', 'Whole tissue', replace_na(broad_type, 'Unannotated'))) |>
  select(-cluster)
write_csv(results_tbl, tab_path(cfg, 'pathway_enrichment.csv.gz'))
summary_tbl <- results_tbl |> group_by(subset, contrast, broad_type) |>
  summarise(pathways_tested = n(), significant_fdr05 = sum(padj < 0.05, na.rm = TRUE),
            positive_nes_fdr05 = sum(padj < 0.05 & NES > 0, na.rm = TRUE),
            negative_nes_fdr05 = sum(padj < 0.05 & NES < 0, na.rm = TRUE), .groups = 'drop')
write_csv(summary_tbl, tab_path(cfg, 'pathway_summary.csv'))
write_csv(bind_rows(coverage_rows), tab_path(cfg, 'pathway_gene_coverage.csv.gz'))
write_csv(bind_rows(statuses), tab_path(cfg, 'pathway_status.csv'))
step_finish(cfg)
