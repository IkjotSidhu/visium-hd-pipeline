# ---- Source Space Ranger overview (full libraries, before pipeline sampling/QC) ----
read_spaceranger_overview <- function(sheet) {
  bins <- list(); sequencing <- list(); cells <- list(); status <- list()
  for (lib in unique(sheet$library_id)) {
    path <- file.path(sheet$spaceranger_dir[match(lib,sheet$library_id)], 'metrics_summary.csv')
    if (!file.exists(path)) {
      status[[lib]] <- tibble(library_id=lib, source=path, status='Missing metrics_summary.csv')
      next
    }
    source <- tryCatch(read_csv(path, col_types=cols(.default=col_character()), show_col_types=FALSE), error=function(e) NULL)
    if (is.null(source) || nrow(source) != 1) {
      status[[lib]] <- tibble(library_id=lib, source=path, status='Unreadable or not a single-library metrics file')
      next
    }
    value <- function(name) {
      if (!name %in% names(source)) return(NA_real_)
      x <- source[[name]][1]; v <- suppressWarnings(as.numeric(gsub('[,%]', '', x)))
      if (grepl('%',x,fixed=TRUE)) v <- v/100
      v
    }
    status[[lib]] <- tibble(library_id=lib, source=path, status='Loaded')
    for (size in c(2,8,16)) {
      unit <- if (size == 2) 'Square' else 'Bin'
      plural <- if (size == 2) 'Squares' else 'Bins'
      suffix <- paste0(' ',size,' µm')
      row <- tibble(library_id=lib, bin_size_um=size,
        bins_under_tissue=value(paste0('Number of ',plural,' Under Tissue',suffix)),
        mean_reads=value(paste0('Mean Reads Under Tissue per ',unit,suffix)),
        mean_umi=value(paste0('Mean UMIs Under Tissue per ',unit,suffix)),
        mean_genes=value(paste0('Mean Genes Under Tissue per ',unit,suffix)),
        tissue_bins_pct=100*value(paste0('Fraction of ',plural,' Under Tissue',suffix)),
        genes_detected=value(paste0('Total Genes Detected Under Tissue',suffix)))
      if (any(!is.na(unlist(row[,-c(1,2)])))) bins[[paste(lib,size)]] <- row
    }
    sequencing[[lib]] <- tibble(library_id=lib, reads=value('Number of Reads'),
      saturation_pct=100*value('Sequencing Saturation'), valid_barcodes_pct=100*value('Valid Barcodes'),
      mapped_probe_pct=100*value('Reads Mapped Confidently to Probe Set'),
      reads_tissue_pct=100*value('Fraction Reads in Squares Under Tissue'),
      umi_per_mm2=value('UMIs per sq mm of Tissue'), reads_per_mm2=value('Reads per sq mm of Tissue'))
    cells[[lib]] <- tibble(library_id=lib, n_cells=value('Number of Cells'),
      reads_cells_pct=100*value('Reads in Cells'), umi_cells_pct=100*value('UMIs in Cells'),
      mean_reads=value('Mean Reads per Cell'), median_umi=value('Median UMIs per Cell'),
      median_genes=value('Median Genes per Cell'), median_cell_area_um2=value('Median Cell Area (μm²)'),
      median_nucleus_area_um2=value('Median Nucleus Area (μm²)'))
  }
  list(bins=bind_rows(bins), sequencing=bind_rows(sequencing), cells=bind_rows(cells), status=bind_rows(status))
}

render_spaceranger_overview <- function(x, table_prefix, max_rows, id_prefix) {
  if (is.null(x)) return(character())
  labels <- c(bin_size_um='Bin size (µm)', bins_under_tissue='Bins under tissue', mean_reads='Mean reads/unit',
    mean_umi='Mean UMIs/unit',mean_genes='Mean genes/unit',tissue_bins_pct='Bins under tissue (%)',
    genes_detected='Distinct genes detected',reads='Sequenced reads',saturation_pct='Sequencing saturation (%)',
    valid_barcodes_pct='Valid barcodes (%)',mapped_probe_pct='Confident probe mapping (%)',
    reads_tissue_pct='Reads under tissue (%)',umi_per_mm2='UMIs/mm² tissue',reads_per_mm2='Reads/mm² tissue',
    n_cells='Cells',reads_cells_pct='Reads in cells (%)',umi_cells_pct='UMIs in cells (%)',
    median_umi='Median UMIs/cell',median_genes='Median genes/cell',median_cell_area_um2='Median cell area (µm²)',
    median_nucleus_area_um2='Median nucleus area (µm²)',library_id='Capture area')
  table <- function(data) {
    for (n in names(data)) if (is.numeric(data[[n]])) {
      digits <- if (grepl('pct$|^mean_|^median_',n)) 1 else 0
      data[[n]] <- ifelse(is.na(data[[n]]),'—',formatC(data[[n]],format='f',digits=digits,big.mark=','))
    }
    names(data) <- ifelse(names(data) %in% names(labels),labels[names(data)],names(data))
    html_table(data,max_rows)
  }
  out <- c(paste0('<section id="',id_prefix,'qc-bin-overview"><h3>Bin Metrics Overview — Space Ranger</h3>'),
    '<p>Full capture-area metrics imported from Space Ranger, before pipeline replicate assignment, test sampling and QC. These values remain full-library totals even in a test run. 2 µm squares are shown for context; the pipeline analyzes the configured 8 µm, 16 µm and cell outputs.</p>',
    '<p>Mean genes refers to genes per bin; distinct genes is the union detected across tissue. The percentage of bins under tissue describes capture-area coverage, not the fraction passing pipeline QC.</p>')
  for (field in c('bins','sequencing','cells')) {
    data <- x[[field]]
    out <- c(out,paste0('<h4>',switch(field,bins='Bins under tissue',sequencing='Sequencing and capture',cells='Cell segmentation metrics'),'</h4>'))
    out <- c(out,if (nrow(data)) table(data) else '<p>Unavailable in supplied Space Ranger metrics.</p>',
      paste0('<p><a href="',html_escape(paste0(table_prefix,'spaceranger_',field,'_overview.csv')),'">Download CSV</a></p>'))
  }
  out <- c(out,'<details><summary>Source files and availability</summary>',table(x$status),'</details></section>')
  out
}

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
       spaceranger = if ('spaceranger_dir' %in% names(sheet)) read_spaceranger_overview(sheet) else NULL,
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
  out <- c(out, render_spaceranger_overview(reports$spaceranger,table_prefix,max_rows,id_prefix))
  out <- c(out,render_qc_relationships(reports$relationships,table_prefix))
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
  if (!is.null(reports$spaceranger)) for (field in names(reports$spaceranger))
    write_csv(reports$spaceranger[[field]],tab_path(cfg,paste0('spaceranger_',field,'_overview.csv')))
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
    '<nav><a href="#qc-bin-overview">Bin metrics</a> · <a href="#qc-group">Groups</a> · <a href="#qc-sample">Samples</a> · <a href="#qc-replicate">Replicates</a></nav>',
    render_qc_sections(reports,max_rows = cfg$qc$report$max_rows),
    '<p>Raw CSVs retain numeric precision; displayed values are rounded. Post-QC statistics for empty groups are unavailable (—), not zero.</p></body></html>')
  writeLines(document,file.path(cfg$step_dir,'qc_report.html'),useBytes = TRUE)
}

# Fit every unit; cap only the number of displayed points for readable exports.
write_qc_relationship_plots <- function(qc_df, cfg, count_col, feature_col) {
  sample_col <- cfg$qc$report$sample_col
  dir.create(fig_path(cfg,'umi_genes'),recursive=TRUE,showWarnings=FALSE)
  results <- list()
  for (level in c('sample','replicate')) {
    column <- if (level == 'sample') sample_col else 'replicate_id'
    groups <- unique(as.character(qc_df[[column]]))
    for (i in seq_along(groups)) for (stage in c('before_qc','after_qc')) {
      id <- groups[i]
      selected <- if (is.na(id)) is.na(qc_df[[column]]) else !is.na(qc_df[[column]]) & qc_df[[column]] == id
      if (stage == 'after_qc') selected <- selected & qc_df$pass
      data <- tibble(umi=qc_df[[count_col]][selected],genes=qc_df[[feature_col]][selected])
      data <- data |> filter(is.finite(umi),is.finite(genes))
      n_units <- nrow(data)
      status <- if (n_units < 3) 'Fewer than three units' else if (n_distinct(data$umi)<2 || n_distinct(data$genes)<2) 'Constant UMI or gene counts' else 'ok'
      fit <- if (status == 'ok') lm(genes ~ umi,data=data) else NULL
      r2 <- if (is.null(fit)) NA_real_ else summary(fit)$r.squared
      slope <- if (is.null(fit)) NA_real_ else unname(coef(fit)[2])
      intercept <- if (is.null(fit)) NA_real_ else unname(coef(fit)[1])
      # Evenly spaced row indices are deterministic and do not change the RNG.
      display <- data[unique(round(seq(1,max(1,n_units),length.out=min(n_units,20000)))),,drop=FALSE]
      label <- if (is.na(id)) '(Missing)' else id
      subtitle <- paste0(if (stage=='before_qc') 'Before QC' else 'After QC',
        ' | n = ',format(n_units,big.mark=','),' | R² = ',if (is.na(r2)) 'unavailable' else sprintf('%.3f',r2))
      caption <- paste0('OLS: detected genes ~ UMI, with intercept; raw counts.\nFit uses all ',n_units,
        ' units; ',nrow(display),' points displayed.',if (!is.null(cfg$test_fraction)) paste0('\nTest run: ',100*cfg$test_fraction,'% selected before QC.') else '',
        if (status!='ok') paste0('\n',status,'.') else '')
      p <- ggplot(display,aes(umi,genes)) + geom_point(size=0.65,alpha=0.22,colour='#225f9c') +
        labs(title=paste(label,'—',if (cfg$bin_size=='cells') 'cells' else paste0(as.integer(cfg$bin_size),' µm bins')),
             subtitle=subtitle,x='UMIs per bin/cell',y='Detected genes per bin/cell',caption=caption) +
        theme_classic(base_size=12) + theme(plot.caption=element_text(hjust=0))
      if (!is.null(fit)) p <- p + geom_abline(intercept=intercept,slope=slope,colour='#c54627',linewidth=0.8)
      stem <- paste0('umi_genes/',level,'_',sprintf('%03d',i),'_',stage)
      for (extension in c('png','pdf')) ggsave(fig_path(cfg,paste0(stem,'.',extension)),p,
        width=7,height=5.5,dpi=140,bg='white')
      results[[length(results)+1]] <- tibble(level=level,id=id,stage=stage,n_units=n_units,
        n_displayed=nrow(display),r_squared=r2,slope=slope,intercept=intercept,status=status,
        png=paste0(stem,'.png'),pdf=paste0(stem,'.pdf'))
    }
  }
  result <- bind_rows(results)
  write_csv(result,tab_path(cfg,'umi_genes_fit_statistics.csv'))
  result
}

render_qc_relationships <- function(x, table_prefix) {
  if (is.null(x) || !nrow(x)) return(character())
  figures <- sub('tables/$','figures/',table_prefix)
  out <- c('<h3>UMIs versus detected genes</h3>',
    '<p>Separate sample and replicate plots before and after QC. R² describes an ordinary least-squares fit of detected genes against UMIs on raw counts, with an intercept; it is not a stand-alone quality score. Fits use all units, even when displayed points are capped at 20,000. Gene-count saturation can produce curvature that a straight line does not capture.</p>',
    paste0('<p><a href="',html_escape(paste0(table_prefix,'umi_genes_fit_statistics.csv')),'">Download fit statistics</a></p>'))
  for (level in c('sample','replicate')) {
    out <- c(out,paste0('<details><summary>',if (level=='sample') 'Sample plots' else 'Replicate plots','</summary>'))
    for (i in which(x$level==level)) out <- c(out,
      paste0('<h4>',html_escape(paste(x$id[i],x$stage[i])),'</h4><a href="',html_escape(paste0(figures,x$pdf[i])),
        '">PDF</a><br><img loading="lazy" style="max-width:100%;width:760px" alt="UMIs versus detected genes" src="',
        html_escape(paste0(figures,x$png[i])),'">'))
    out <- c(out,'</details>')
  }
  out
}
