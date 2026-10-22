# ---- 06_pseudobulk_de.R ----
# Condition comparison with the replicate as the unit of replication: counts
# are summed per replicate (whole tissue, and per cluster), then DESeq2 tests
# condition_col against its reference level, adjusting for de$covariates
# (e.g. donor for a paired design). Bin/cell-level tests are pseudoreplicated
# and are deliberately not used.
#
# Each result also reports whether the direction agrees within every level of
# the first covariate (e.g. every donor) - the "high confidence" rule from FAQ.md.
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

# Per-replicate sample table (condition, covariates)
sample_tbl <- obj[[]] |> distinct(replicate_id, across(all_of(c(cond, covars)))) |> as_tibble()
ref_level  <- de$reference %||% sort(unique(sample_tbl[[cond]]))[1]
log_msg(cfg, 'Testing ', cond, ' vs reference "', ref_level, '"', if (length(covars)) paste0(', adjusting for ', paste(covars, collapse = ' + ')))
log_msg(cfg, 'Replicates per condition: ', paste(names(table(sample_tbl[[cond]])), table(sample_tbl[[cond]]), sep = '=', collapse = ', '))

counts_all <- LayerData(obj, assay = assay, layer = 'counts')

# ---- One DESeq2 test on a set of pseudobulk samples ----
# groups: per-unit pseudobulk label (replicate, or replicate x cluster); returns NULL if not testable
run_deseq <- function(cells, label) {
  rep_of <- obj$replicate_id[cells]
  n_per  <- table(rep_of)
  keep_reps <- names(n_per)[n_per >= de$min_bins]
  cells  <- cells[rep_of %in% keep_reps]
  if (length(keep_reps) < 2) return(NULL)

  # Sum counts per replicate (sparse matrix product with a replicate indicator)
  ind <- Matrix::sparse.model.matrix(~ 0 + rep, data.frame(rep = factor(obj$replicate_id[cells], levels = keep_reps)))
  colnames(ind) <- keep_reps
  pb <- as.matrix(counts_all[, cells] %*% ind)
  pb <- pb[rowSums(pb) >= 10, , drop = FALSE]                     # drop near-empty genes

  cd <- sample_tbl |> filter(replicate_id %in% keep_reps) |> slice(match(keep_reps, replicate_id)) |>
    mutate(across(all_of(c(cond, covars)), factor))
  # A subset (e.g. a cluster) may lack the reference condition or have only one condition
  if (!ref_level %in% levels(cd[[cond]]) || nlevels(cd[[cond]]) < 2) {
    log_msg(cfg, label, ': only ', paste(levels(cd[[cond]]), collapse = '/'), ' replicates with >= ', de$min_bins, ' units - skipped')
    return(NULL)
  }
  cd[[cond]] <- relevel(cd[[cond]], ref = ref_level)

  # Drop covariates the design can't support: one level in this subset (model.matrix
  # fails), then confounded / no residual df
  one_level <- covars[vapply(covars, \(v) nlevels(droplevels(cd[[v]])) < 2, logical(1))]
  if (length(one_level)) log_msg(cfg, label, ': dropping covariate ', paste(one_level, collapse = ', '), ' (one level in this subset)')
  use_cov <- setdiff(covars, one_level)
  cd <- cd |> mutate(across(all_of(use_cov), droplevels))
  repeat {
    fml <- as.formula(paste('~', paste(c(use_cov, cond), collapse = ' + ')))
    mm  <- model.matrix(fml, cd)
    if (qr(mm)$rank == ncol(mm) && nrow(mm) > ncol(mm)) break
    if (!length(use_cov)) { log_msg(cfg, label, ': not enough replicates for any test - skipped'); return(NULL) }
    log_msg(cfg, label, ': dropping covariate ', tail(use_cov, 1), ' (not estimable with these replicates)')
    use_cov <- head(use_cov, -1)
  }

  dds <- DESeq2::DESeqDataSetFromMatrix(round(pb), as.data.frame(cd), design = fml)
  dds <- DESeq2::DESeq(dds, quiet = TRUE)

  # Within-stratum direction (e.g. per donor) from CPM, for the consistency rule
  cpm <- t(t(pb) / colSums(pb)) * 1e6
  strata <- if (length(covars)) cd[[covars[1]]] else factor(rep('all', nrow(cd)))

  out <- map_dfr(setdiff(levels(cd[[cond]]), ref_level), function(lv) {
    res <- DESeq2::results(dds, contrast = c(cond, lv, ref_level)) |> as.data.frame() |> rownames_to_column('gene') |> as_tibble()
    # Direction of the effect in each stratum that has both conditions
    dirs <- map(levels(strata), function(s) {
      a <- cd[[cond]] == lv & strata == s; b <- cd[[cond]] == ref_level & strata == s
      if (!any(a) || !any(b)) return(NULL)
      sign(log2(rowMeans(cpm[, a, drop = FALSE]) + 1) - log2(rowMeans(cpm[, b, drop = FALSE]) + 1))
    }) |> compact()
    consistent <- if (length(dirs)) Reduce(`&`, map(dirs, \(d) d[res$gene] == sign(res$log2FoldChange))) else NA
    res |> mutate(contrast = paste(lv, 'vs', ref_level), consistent_across_strata = consistent,
                  n_strata_compared = length(dirs), design = deparse(fml), subset = label, .before = 1)
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
res_tbl <- bind_rows(results_all)
write_csv(res_tbl, tab_path(cfg, 'pseudobulk_de_results.csv.gz'))

# ---- Summary: significant genes per subset ----
summary_tbl <- res_tbl |>
  group_by(subset, contrast, design) |>
  summarise(genes_tested = n(),
            sig_padj05     = sum(padj < 0.05, na.rm = TRUE),
            sig_and_consistent = sum(padj < 0.05 & consistent_across_strata, na.rm = TRUE),
            up = sum(padj < 0.05 & log2FoldChange > 0, na.rm = TRUE),
            down = sum(padj < 0.05 & log2FoldChange < 0, na.rm = TRUE), .groups = 'drop')
print(summary_tbl, n = Inf, width = Inf)
write_csv(summary_tbl, tab_path(cfg, 'pseudobulk_de_summary.csv'))

# ---- Volcano: whole tissue ----
if (!is.null(results_all$all_bins)) {
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
}

step_finish(cfg)
