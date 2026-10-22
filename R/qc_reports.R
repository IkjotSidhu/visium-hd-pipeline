# ---- Readable QC summaries ----
# All means are calculated over individual bins/cells, never means of means.
# Retention starts at assigned, test-selected units entering the pipeline's QC.
qc_report_settings <- function(cfg, sheet) {
  sample <- cfg$qc$report$sample_col
  groups <- cfg$qc$report$group_cols
  if (is.null(groups)) groups <- intersect(unique(c(cfg$de$condition_col, 'donor', 'tissue')), names(sheet))
  list(sample_col = sample, group_cols = as.character(unlist(groups)))
}

qc_level_summary <- function(data, entities, keys) {
  stat <- function(x, fun) if (!length(x) || all(is.na(x))) NA_real_ else fun(x, na.rm = TRUE)
  summarize_stage <- function(x, stage) {
    values <- x |> group_by(across(all_of(keys))) |>
      summarise(n_units = n(), n_samples = n_distinct(.sample), n_libraries = n_distinct(library_id),
        n_replicates = n_distinct(replicate_id),
        n_donors = if ('.donor' %in% names(x)) n_distinct(.donor, na.rm = TRUE) else NA_integer_,
        total_umi = sum(.umi), mean_umi = stat(.umi, mean), median_umi = stat(.umi, median),
        mean_genes = stat(.genes, mean), median_genes = stat(.genes, median),
        mean_mt_pct = stat(.mt, mean), median_mt_pct = stat(.mt, median), n_mt_missing = sum(is.na(.mt)), .groups = 'drop')
    entities |> left_join(values, by = keys) |>
      mutate(across(any_of(c('n_units','n_samples','n_libraries','n_replicates','total_umi','n_mt_missing')), ~ replace_na(.x, 0)),
             stage = stage, .before = 1)
  }
  before <- summarize_stage(data, 'before_qc')
  after <- summarize_stage(data |> filter(pass), 'after_qc')
  filters <- data |> group_by(across(all_of(keys))) |>
    summarise(fail_umi = sum(.fail_umi), fail_genes = sum(.fail_genes), fail_mt = sum(.fail_mt),
              fail_any = sum(!pass), .groups = 'drop')
  retention <- before |> select(all_of(keys), n_before = n_units) |>
    left_join(after |> select(all_of(keys), n_after = n_units), by = keys) |>
    mutate(n_removed = n_before - n_after, pct_kept = if_else(n_before > 0, 100*n_after/n_before, NA_real_)) |>
    left_join(filters, by = keys) |>
    mutate(across(starts_with('fail_'), ~ replace_na(.x, 0)))
  summary <- retention |> left_join(after |> select(-stage, -n_units), by = keys)
  list(summary = summary, metrics = bind_rows(before,after), retention = retention, keys = keys)
}

build_qc_reports <- function(qc_df, sheet, cfg, count_col, feature_col, flow = tibble()) {
  settings <- qc_report_settings(cfg, sheet)
  for (v in c(count_col, feature_col))
    if (any(!is.finite(qc_df[[v]]) | qc_df[[v]] < 0)) stop('QC counts/features must be finite and nonnegative: ', v)
  data <- qc_df |> mutate(.sample = as.character(.data[[settings$sample_col]]),
    .umi = .data[[count_col]], .genes = .data[[feature_col]], .mt = percent.mt,
    .fail_umi = .umi < qc_value(cfg,'min_counts'), .fail_genes = .genes < qc_value(cfg,'min_features'),
    .fail_mt = replace_na(.mt,0) > qc_value(cfg,'max_mt_pct'))
  if ('donor' %in% names(data)) data$.donor <- data$donor
  group_parts <- lapply(settings$group_cols, function(column) {
    x <- data |> mutate(group_variable = column, group = as.character(.data[[column]]))
    entities <- sheet |> transmute(group_variable = column, group = as.character(.data[[column]])) |> distinct()
    qc_level_summary(x, entities, c('group_variable','group'))
  })
  if (length(group_parts)) {
    group <- list(keys = c('group_variable','group'))
    for (field in c('summary','metrics','retention')) group[[field]] <- bind_rows(lapply(group_parts, `[[`, field))
  } else group <- qc_level_summary(data[0,] |> mutate(group_variable = character(), group = character()),
                                   tibble(group_variable = character(), group = character()), c('group_variable','group'))
  sample <- qc_level_summary(data |> mutate(sample_id = .sample),
    sheet |> transmute(sample_id = as.character(.data[[settings$sample_col]])) |> distinct(), 'sample_id')
  replicate <- qc_level_summary(data, sheet |> distinct(replicate_id), 'replicate_id')
  context <- sheet |> select(any_of(unique(c('replicate_id',settings$sample_col,'library_id','donor',cfg$de$condition_col,'tissue')))) |>
    mutate(sample_id = as.character(.data[[settings$sample_col]]))
  for (field in c('summary','metrics','retention')) replicate[[field]] <- replicate[[field]] |> left_join(context, by = 'replicate_id')
  if (nrow(flow)) flow <- flow |> left_join(data |> group_by(library_id) |> summarise(n_after_qc = sum(pass), .groups = 'drop'), by = 'library_id')
  list(group = group, sample = sample, replicate = replicate, flow = flow, settings = settings,
       unit_name = if (cfg$bin_size == 'cells') 'cells' else paste0(as.integer(cfg$bin_size), ' µm bins/spots'),
       res_label = cfg$res_label, test_fraction = cfg$test_fraction,
       thresholds = list(min_counts = qc_value(cfg,'min_counts'), min_features = qc_value(cfg,'min_features'), max_mt_pct = qc_value(cfg,'max_mt_pct')))
}

qc_pretty_table <- function(data, columns, max_rows = 100) {
  columns <- intersect(columns, names(data)); data <- data[,columns,drop = FALSE]
  labels <- c(group_variable = 'Group by', group = 'Group', sample_id = 'Sample', replicate_id = 'Replicate',
    n_before = 'Before QC', n_after = 'Kept', n_removed = 'Removed', pct_kept = 'Kept (%)', n_units = 'Bins/cells',
    n_samples = 'Samples', n_replicates = 'Replicates', mean_umi = 'Mean UMI/unit', median_umi = 'Median UMI/unit',
    mean_genes = 'Mean genes/unit', median_genes = 'Median genes/unit', total_umi = 'Total UMIs',
    mean_mt_pct = 'Mean MT (%)', median_mt_pct = 'Median MT (%)', fail_umi = 'Fail UMI', fail_genes = 'Fail genes',
    fail_mt = 'Fail MT', fail_any = 'Any QC failure', library_id = 'Capture area',
    n_loaded = 'Loaded', n_assigned = 'Assigned', n_unassigned = 'Unassigned', n_selected = 'Entering QC', n_after_qc = 'Kept')
  for (name in names(data)) if (is.numeric(data[[name]])) {
    digits <- if (name %in% c('pct_kept','mean_umi','median_umi','mean_genes','median_genes','mean_mt_pct','median_mt_pct')) 1 else 0
    data[[name]] <- ifelse(is.na(data[[name]]), '—', formatC(data[[name]], format = 'f', digits = digits, big.mark = ','))
  } else data[[name]][is.na(data[[name]])] <- '(Missing)'
  names(data) <- ifelse(names(data) %in% names(labels), labels[names(data)], names(data))
  html_table(data, max_rows)
}

render_qc_sections <- function(reports, table_prefix = 'tables/', max_rows = 100, id_prefix = '') {
  id <- function(level) paste0(id_prefix,'qc-',level)
  unit <- html_escape(reports$unit_name)
  out <- c(paste0('<h2>QC summary — ',unit,'</h2>'),
    paste0('<p>Sample = <code>',html_escape(reports$settings$sample_col),'</code>. Group tables summarize each configured metadata column separately.</p>'),
    paste0('<p>QC thresholds: UMIs ≥ ',reports$thresholds$min_counts,'; genes ≥ ',reports$thresholds$min_features,
           '; mitochondrial UMIs ≤ ',reports$thresholds$max_mt_pct,'%. Means/medians refer to individual ',unit,'.</p>'),
    '<p>Before QC means assigned units selected for this run. Test sampling and unassigned units are separate from QC failures. Gene metrics count detected genes per unit, not distinct genes across a group.</p>')
  if (!is.null(reports$test_fraction)) out <- c(out, paste0('<p class="notice">Test run: ',100*reports$test_fraction,'% of each replicate was selected before QC. Counts are not full-dataset totals.</p>'))
  for (level in c('group','sample','replicate')) {
    x <- reports[[level]]; keys <- x$keys
    out <- c(out, paste0('<section id="',id(level),'"><h3>',switch(level,group='Groups',sample='Samples',replicate='Replicates'),'</h3>'),
      '<h4>Retention</h4>', qc_pretty_table(x$retention,c(keys,'n_before','n_after','n_removed','pct_kept'),max_rows),
      '<h4>Quality after QC</h4>', qc_pretty_table(x$metrics |> filter(stage == 'after_qc'),
        c(keys,'n_units','mean_umi','median_umi','mean_genes','median_genes','mean_mt_pct','n_samples','n_replicates'),max_rows),
      '<details><summary>Before-QC quality and individual filters</summary>',
      qc_pretty_table(x$metrics |> filter(stage == 'before_qc'),
        c(keys,'n_units','mean_umi','median_umi','mean_genes','median_genes','mean_mt_pct'),max_rows),
      '<p>Filter counts overlap: one unit may fail more than one criterion. “Any QC failure” counts each removed unit once.</p>',
      qc_pretty_table(x$retention,c(keys,'fail_umi','fail_genes','fail_mt','fail_any'),max_rows), '</details>',
      paste0('<p>Download CSV: <a href="',html_escape(paste0(table_prefix,'qc_by_',level,'.csv')),'">Summary</a> · ',
        '<a href="',html_escape(paste0(table_prefix,'qc_metrics_by_',level,'.csv')),'">Before/after metrics</a> · ',
        '<a href="',html_escape(paste0(table_prefix,'qc_retention_by_',level,'.csv')),'">Retention and filters</a></p></section>'))
  }
  if (nrow(reports$flow)) out <- c(out,'<h3>Loading and assignment flow</h3>',
    '<p>Loaded counts come from Space Ranger filtered matrices; they are not all 2 µm squares or all image-detected nuclei. No capture areas are pooled across different bin sizes.</p>',
    qc_pretty_table(reports$flow,c('library_id','n_loaded','n_assigned','n_unassigned','n_selected','n_after_qc'),max_rows))
  out
}

write_qc_reports <- function(reports, cfg) {
  for (level in c('group','sample','replicate')) {
    write_csv(reports[[level]]$summary,tab_path(cfg,paste0('qc_by_',level,'.csv')))
    write_csv(reports[[level]]$metrics,tab_path(cfg,paste0('qc_metrics_by_',level,'.csv')))
    write_csv(reports[[level]]$retention,tab_path(cfg,paste0('qc_retention_by_',level,'.csv')))
  }
  write_csv(reports$flow,tab_path(cfg,'qc_loading_flow.csv'))
  saveRDS(reports,file.path(cfg$step_dir,'qc_reports.rds'))
  document <- c('<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">',
    '<title>QC report</title><style>body{font:16px/1.5 system-ui,sans-serif;max-width:1250px;margin:35px auto;padding:0 20px;color:#243247;background:#f8fafc}section{padding:20px;background:white;border:1px solid #dce2e8;border-radius:8px;margin:20px 0}table{border-collapse:collapse;width:100%;font-size:14px}th,td{padding:9px;border-bottom:1px solid #ddd;text-align:left;white-space:nowrap}th{background:#edf3f8}.table-scroll{overflow:auto}a{color:#165c9c}summary{cursor:pointer}.notice{background:#fff4db;padding:12px}</style></head><body>',
    paste0('<h1>',html_escape(cfg$project_name),' — quality control</h1>'),
    '<nav><a href="#qc-group">Groups</a> · <a href="#qc-sample">Samples</a> · <a href="#qc-replicate">Replicates</a></nav>',
    render_qc_sections(reports,max_rows = cfg$qc$report$max_rows),
    '<p>Raw CSVs retain numeric precision; displayed values are rounded. Post-QC statistics for empty groups are unavailable (—), not zero.</p></body></html>')
  writeLines(document,file.path(cfg$step_dir,'qc_report.html'),useBytes = TRUE)
}
