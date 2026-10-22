# ---- 06_pseudobulk_de.R ----
# Condition comparison with the replicate as the unit of replication: counts
# are summed per replicate (whole tissue, and per cluster), then DESeq2 tests
# condition_col against its reference level, adjusting for de$covariates
# (e.g. donor for a paired design). Bin/cell-level tests are pseudoreplicated
# and are deliberately not used.
#
# Each result also reports whether the direction agrees within every level of
# the first covariate (e.g. every donor), requiring replicated paired support.
#
# Usage: Rscript steps/06_pseudobulk_de.R config.yaml [8|16|cells]

# ---- Setup ----
source(file.path(dirname(sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value = TRUE))), '..', 'R', 'utils.R'))
cfg      <- step_setup('06_pseudobulk_de')
de       <- cfg$de
# DESeq2 used via DESeq2:: - attaching it masks dplyr verbs (S4Vectors::slice etc.)
obj      <- read_step_rds(cfg, '03_cluster')
assay    <- DefaultAssay(obj)
main_col <- paste0('cluster_res_', cfg$clustering$main_resolution)
covars   <- unlist(de$covariates)
cond     <- de$condition_col
notes <- c('Exploratory condition associations; statistical adjustment does not resolve unmeasured confounding.',
           'Whole-tissue DE can reflect composition and expression changes.',
           'Per-stratum CPM effects are descriptive and do not provide independent validation p-values.',
           unlist(de$interpretation_notes))
writeLines(notes, file.path(cfg$step_dir, 'interpretation_notes.txt'))

# Per-replicate sample table (condition, covariates)
sample_tbl <- obj[[]] |> distinct(replicate_id, across(all_of(c(cond, covars)))) |> as_tibble()
if (anyDuplicated(sample_tbl$replicate_id)) stop('Each replicate must have one condition/covariate combination')
ref_level  <- de$reference %||% sort(unique(sample_tbl[[cond]]))[1]
log_msg(cfg, 'Testing ', cond, ' vs reference "', ref_level, '"', if (length(covars)) paste0(', adjusting for ', paste(covars, collapse = ' + ')))
log_msg(cfg, 'Replicates per condition: ', paste(names(table(sample_tbl[[cond]])), table(sample_tbl[[cond]]), sep = '=', collapse = ', '))

counts_all <- LayerData(obj, assay = assay, layer = 'counts')
diagnostic_rows <- list()
status_rows <- list(); sample_rows <- list(); donor_rows <- list(); pseudobulks <- list()
# Stable schema: all-skipped runs still produce readable, header-only results.
empty_results <- tibble(subset = character(), contrast = character(), design = character(),
  dropped_covariates = character(), n_strata_compared = integer(), n_strata_expected = integer(),
  strata_support = character(), direction_agrees_available = logical(), consistent_across_strata = logical(),
  gene = character(), baseMean = double(), log2FoldChange = double(), lfcSE = double(),
  stat = double(), pvalue = double(), padj = double())


# ---- One DESeq2 test on a set of pseudobulk samples ----
# groups: per-unit pseudobulk label (replicate, or replicate x cluster); returns NULL if not testable
run_deseq <- function(cells, label) {
  rep_of <- obj$replicate_id[cells]
  n_per  <- table(rep_of)
  keep_reps <- names(n_per)[n_per >= de$min_bins]
  cells  <- cells[rep_of %in% keep_reps]
  eligibility <- sample_tbl |>
    mutate(subset = label, n_units = as.integer(n_per[replicate_id]),
           n_units = replace_na(n_units, 0L), eligible = replicate_id %in% keep_reps, .before = 1)
  sample_rows[[label]] <<- eligibility
  record_status <- function(status, reason = '', design = '', dropped = '') {
    status_rows[[label]] <<- tibble(subset = label, status, reason, design,
                                  dropped_covariates = dropped, n_replicates = length(keep_reps))
  }
  skip <- function(reason, design = '', dropped = '') {
    log_msg(cfg, label, ': ', reason, ' - skipped')
    record_status('skipped', reason, design, dropped)
    NULL
  }
  if (length(keep_reps) < 2) return(skip('fewer than two eligible replicates'))

  # Sum counts per replicate (sparse matrix product with a replicate indicator)
  ind <- Matrix::sparse.model.matrix(~ 0 + rep, data.frame(rep = factor(obj$replicate_id[cells], levels = keep_reps)))
  colnames(ind) <- keep_reps
  pb <- as.matrix(counts_all[, cells] %*% ind)
  pb <- pb[rowSums(pb) >= 10, , drop = FALSE]                     # drop near-empty genes

  if (!nrow(pb)) return(skip('no genes pass pseudobulk count filtering'))
  if (any(colSums(pb) == 0)) return(skip('zero library after gene filtering'))

  cd <- sample_tbl |> filter(replicate_id %in% keep_reps) |> slice(match(keep_reps, replicate_id)) |>
    mutate(across(all_of(c(cond, covars)), factor))
  # A subset (e.g. a cluster) may lack the reference condition or have only one condition
  if (!ref_level %in% levels(cd[[cond]]) || nlevels(cd[[cond]]) < 2) {
    log_msg(cfg, label, ': only ', paste(levels(cd[[cond]]), collapse = '/'), ' replicates with >= ', de$min_bins, ' units - skipped')
    return(skip('reference condition absent or fewer than two conditions'))
  }
  cd[[cond]] <- relevel(cd[[cond]], ref = ref_level)

  design_check <- choose_de_design(cd, cond, covars, de$covariate_policy)
  fml <- design_check$formula
  design_text <- paste(deparse(fml), collapse = '')
  dropped <- paste(design_check$dropped, collapse = ';')
  if (!design_check$ok) return(skip(design_check$reason, design_text, dropped))
  if (nzchar(dropped)) log_msg(cfg, label, ': explicit fallback removed ', dropped)
  cd <- as.data.frame(cd)
  rownames(cd) <- cd$replicate_id

  dds <- DESeq2::DESeqDataSetFromMatrix(round(pb), as.data.frame(cd), design = fml)
  dds <- DESeq2::DESeq(dds, quiet = TRUE)

  record_status('tested', design = design_text, dropped = dropped)
  pseudobulks[[label]] <<- list(counts = pb, normalized_counts = DESeq2::counts(dds, normalized = TRUE),
                               samples = cd, design = design_text)

  if (isTRUE(de$diagnostics$run)) {
    md_cols <- unique(c('replicate_id', cond, covars, de$diagnostics$metadata_cols))
    diagnostic_md <- obj[[]] |> select(any_of(md_cols)) |> distinct()
    diagnostic_status <- tryCatch({
      write_pseudobulk_diagnostics(dds, diagnostic_md, cfg, label)
      tibble(subset = label, status = 'complete', reason = '')
    }, error = function(e) {
      log_msg(cfg, label, ': diagnostics failed: ', conditionMessage(e))
      tibble(subset = label, status = 'failed', reason = conditionMessage(e))
    })
    diagnostic_rows[[label]] <<- diagnostic_status
  }

  # Within-stratum direction (e.g. per donor) from CPM, for the consistency rule
  cpm <- t(t(pb) / colSums(pb)) * 1e6
  strata <- if (length(covars)) cd[[covars[1]]] else factor(rep('all', nrow(cd)))

  out <- map_dfr(setdiff(levels(cd[[cond]]), ref_level), function(lv) {
    res <- DESeq2::results(dds, contrast = c(cond, lv, ref_level)) |> as.data.frame() |> rownames_to_column('gene') |> as_tibble()
    # Direction of the effect in each stratum that has both conditions
    expected <- if (length(covars)) sort(unique(as.character(sample_tbl[[covars[1]]]))) else character()
    effects <- map_dfr(levels(strata), function(s) {
      a <- cd[[cond]] == lv & strata == s; b <- cd[[cond]] == ref_level & strata == s
      if (!any(a) || !any(b)) return(NULL)
      tibble(subset = label, contrast = paste(lv, 'vs', ref_level), stratum = s,
             gene = rownames(pb), n_reference = sum(b), n_comparison = sum(a),
             # Descriptive CPM effect, not a donor-specific inferential test.
             log2_cpm_effect = log2(rowMeans(cpm[, a, drop = FALSE]) + 1) -
                               log2(rowMeans(cpm[, b, drop = FALSE]) + 1))
    })
    if (nrow(effects)) donor_rows[[paste(label, lv, sep = '/')]] <<- effects
    dirs <- if (nrow(effects)) map(split(effects, effects$stratum),
      \(x) sign(x$log2_cpm_effect[match(res$gene, x$gene)])) else list()
    evidence <- strata_evidence(dirs, res$log2FoldChange, length(expected), de$min_strata)
    # With no requested stratum there is no independent replication check.
    if (!length(covars)) {
      evidence$consistent[] <- NA
      evidence$support <- 'no_stratification'
    }
    res |> mutate(contrast = paste(lv, 'vs', ref_level),
                  consistent_across_strata = evidence$consistent,
                  direction_agrees_available = evidence$available,
                  n_strata_compared = if (length(covars)) length(dirs) else 0L,
                  n_strata_expected = length(expected), strata_support = evidence$support,
                  dropped_covariates = dropped, design = design_text, subset = label, .before = 1)
  })
  out
}

# ---- Whole tissue ----
results_all <- list(all_bins = run_deseq(colnames(obj), 'all_bins'))

# ---- Per cluster ----
if (isTRUE(de$by_cluster)) {
  for (cl in levels(obj[[main_col]][, 1])) {
    cells <- colnames(obj)[obj[[main_col]][, 1] == cl]
    results_all[[paste0('cluster_', cl)]] <- run_deseq(cells, paste0('cluster_', cl))
  }
}
res_tbl <- bind_rows(c(list(empty_results), results_all))
diagnostic_template <- tibble(subset = character(), status = character(), reason = character())
write_csv(bind_rows(c(list(diagnostic_template), diagnostic_rows)), tab_path(cfg, 'pseudobulk_diagnostics_status.csv'))
write_csv(bind_rows(status_rows), tab_path(cfg, 'pseudobulk_de_status.csv'))
write_csv(bind_rows(sample_rows), tab_path(cfg, 'pseudobulk_sample_eligibility.csv'))
effect_template <- tibble(subset = character(), contrast = character(), stratum = character(),
                         gene = character(), n_reference = integer(), n_comparison = integer(), log2_cpm_effect = double())
write_csv(bind_rows(c(list(effect_template), donor_rows)), tab_path(cfg, 'pseudobulk_stratum_effects.csv.gz'))
saveRDS(pseudobulks, rds_path(cfg, '06_pseudobulk_de'))
write_csv(res_tbl, tab_path(cfg, 'pseudobulk_de_results.csv.gz'))

# ---- Summary: significant genes per subset ----
summary_tbl <- res_tbl |>
  group_by(subset, contrast, design, dropped_covariates, n_strata_compared, n_strata_expected, strata_support) |>
  summarise(genes_tested = n(),
            sig_padj05     = sum(padj < 0.05, na.rm = TRUE),
            sig_and_consistent = sum(padj < 0.05 & consistent_across_strata, na.rm = TRUE),
            up = sum(padj < 0.05 & log2FoldChange > 0, na.rm = TRUE),
            down = sum(padj < 0.05 & log2FoldChange < 0, na.rm = TRUE), .groups = 'drop')
print(summary_tbl, n = Inf, width = Inf)
write_csv(summary_tbl, tab_path(cfg, 'pseudobulk_de_summary.csv'))

# ---- Volcano: whole tissue ----
if (!is.null(results_all$all_bins) && any(!is.na(results_all$all_bins$padj))) {
  vdf <- results_all$all_bins |>
    filter(!is.na(padj)) |>
    mutate(class = case_when(padj < 0.05 & consistent_across_strata ~ 'padj < 0.05, consistent',
                             padj < 0.05 ~ 'padj < 0.05', TRUE ~ 'n.s.'))
  lab <- vdf |> filter(padj < 0.05) |> group_by(sign(log2FoldChange)) |> slice_min(padj, n = 10)
  p_v <- ggplot(vdf, aes(log2FoldChange, -log10(pvalue), colour = class)) +
    ggrastr::rasterise(geom_point(size = 0.6), dpi = 300) +
    geom_text(data = lab, aes(label = gene), size = 2.5, colour = 'black', vjust = -0.6, check_overlap = TRUE) +
    scale_colour_manual(values = c('padj < 0.05, consistent' = '#B2182B', 'padj < 0.05' = '#F4A582', 'n.s.' = 'grey75')) +
    facet_wrap(~ contrast) +
    labs(title = 'Pseudobulk DE - whole tissue', subtitle = unique(vdf$design), colour = NULL) +
    theme_classic()
  ggsave(fig_path(cfg, 'volcano_all_bins.pdf'), p_v, width = 8, height = 6)
  ggsave(fig_path(cfg, 'volcano_all_bins.png'), p_v, width = 8, height = 6, dpi = 100)
} else {
  # A skipped rerun must not leave a volcano from an earlier successful run.
  unlink(c(fig_path(cfg, 'volcano_all_bins.pdf'), fig_path(cfg, 'volcano_all_bins.png')))
  log_msg(cfg, 'No testable whole-tissue results with adjusted p-values; volcano omitted')
}

step_finish(cfg)
