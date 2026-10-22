# ---- 10_report.R ----
# Lightweight static HTML report. No pandoc, browser, server or network needed.
# The launcher submits this once after every selected resolution chain finishes.
source(file.path(dirname(sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value = TRUE))), '..', 'R', 'utils.R'))
cfg <- step_setup('10_report')
filter <- Sys.getenv('VHD_REPORT_RESOLUTIONS', '')
resolutions <- as.character(unlist(cfg$bin_sizes))
if (nzchar(filter)) resolutions <- intersect(resolutions, strsplit(filter, ':', fixed = TRUE)[[1]])
relative <- function(path) substring(path, nchar(cfg$output_dir) + 2L)
href <- function(path) html_escape(URLencode(relative(path), reserved = FALSE))
link <- function(path, label = basename(path)) paste0('<a href="', href(path), '">', html_escape(label), '</a>')
parts <- c('<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">',
  paste0('<title>', html_escape(cfg$project_name), ' — analysis report</title>'),
  '<style>body{font:16px/1.5 system-ui,sans-serif;max-width:1150px;margin:40px auto;padding:0 20px;color:#243247;background:#fafbfc}h1,h2,h3{line-height:1.2}a{color:#165c9c}section{background:white;padding:22px;margin:20px 0;border:1px solid #dce2e8;border-radius:8px}.table-scroll{overflow:auto}table{border-collapse:collapse;width:100%;font-size:13px}th,td{padding:7px;border-bottom:1px solid #ddd;text-align:left;vertical-align:top}th{background:#eff4f8}img{max-width:100%;height:auto}figure{margin:20px 0}small{color:#556}code{overflow-wrap:anywhere}.notice{border-left:4px solid #c28a2c;padding:12px;background:#fff6e5}</style></head><body>',
  paste0('<h1>', html_escape(cfg$project_name), '</h1><p>Visium HD analysis overview · Generated ', html_escape(format(Sys.time())), '</p>'),
  '<p class="notice">Exploratory analysis. Missing, disabled and skipped outputs are shown explicitly. Existing artifacts may come from different partial runs; inspect step timestamps, settings and provenance before combining results.</p>',
  paste0('<p>Resolutions: ', html_escape(paste(resolutions, collapse = ', ')), '. ',
         'Pipeline commit: <code>', html_escape(cfg$pipeline_commit %||% 'unavailable'), '</code>. ',
         'Working tree: ', html_escape(cfg$pipeline_worktree %||% 'unavailable'), '.</p>'),
  '<nav><a href="#study">Study</a> · <a href="#provenance">Provenance</a></nav>')
ss <- read_sample_sheet(cfg)
parts <- c(parts, '<section id="study"><h2>Study design and interpretation</h2>',
           html_table(ss |> select(-spaceranger_dir), cfg$report$max_table_rows))
notes <- c('Tissue strips and bins/cells do not increase the number of independent donors.',
           'Whole-tissue differences can reflect cell composition and expression.', unlist(cfg$de$interpretation_notes))
parts <- c(parts, paste0('<ul>', paste0('<li>', html_escape(notes), '</li>', collapse = ''), '</ul></section>'))
assignment <- step_dir(cfg, '01_assign_replicates')
if (dir.exists(assignment)) {
  images <- list.files(file.path(assignment, 'figures'), pattern = '\\.png$', full.names = TRUE)
  parts <- c(parts, '<section><h2>Replicate assignment</h2>',
             paste0('<p>', vapply(images, link, character(1)), '</p>'), '</section>')
}
summary_files <- c('02_load_qc/tables/qc_summary.csv',
  '06_pseudobulk_de/tables/pseudobulk_de_status.csv',
  '06_pseudobulk_de/tables/pseudobulk_de_summary.csv',
  '06_pseudobulk_de/tables/pseudobulk_diagnostics_status.csv',
  '08_annotation/tables/annotation_status.csv', '09_pathways/tables/pathway_status.csv',
  '09_pathways/tables/pathway_summary.csv')
artifact_rows <- list()
for (res in resolutions) {
  label <- if (res == 'cells') 'cells' else sprintf('%03dum', as.integer(res))
  root <- file.path(cfg$output_dir, label)
  parts <- c(parts, paste0('<section><h2>', html_escape(label), '</h2>'))
  for (step in c('02_load_qc', '03_cluster', '04_markers', '05_banksy', '06_pseudobulk_de', '07_rctd', '08_annotation', '09_pathways')) {
    params <- file.path(root, step, 'params_used.yaml')
    complete <- file.exists(params) && file.exists(file.path(root, step, 'sessionInfo.txt'))
    enabled <- switch(step, '05_banksy' = cfg$banksy$run, '06_pseudobulk_de' = cfg$de$run,
      '07_rctd' = cfg$rctd$run, '08_annotation' = cfg$annotation$run,
      '09_pathways' = cfg$pathways$run && cfg$de$run, TRUE)
    state <- if (!isTRUE(enabled)) 'disabled in current config' else if (complete) 'artifacts present' else 'missing/incomplete'
    parts <- c(parts, paste0('<p><strong>', html_escape(step), '</strong>: ', html_escape(state),
      if (complete) paste0(' · last completed ', html_escape(as.character(file.info(params)$mtime)), ' · ', link(params, 'settings')) else '', '</p>'))
  }
  for (file in summary_files) {
    path <- file.path(root, file)
    if (file.exists(path)) {
      tbl <- tryCatch(read_csv(path, show_col_types = FALSE), error = function(e) NULL)
      if (!is.null(tbl)) parts <- c(parts, paste0('<h3>', html_escape(basename(file)), '</h3><p>', link(path, 'Download table'), '</p>'),
                                    html_table(tbl, cfg$report$max_table_rows))
    }
  }
  all_files <- list.files(root, recursive = TRUE, full.names = TRUE)
  all_files <- all_files[!file.info(all_files)$isdir]
  if (length(all_files)) {
    artifact_rows[[label]] <- tibble(resolution = label, path = vapply(all_files, relative, character(1)),
                                    bytes = file.info(all_files)$size, modified = as.character(file.info(all_files)$mtime))
    images <- all_files[grepl('\\.png$', all_files) &
      grepl('^(all_bins_pca_|all_bins_sample_distances|annotation_composition|umap_cluster_res_|pathway_)', basename(all_files))]
    for (path in head(images, 12)) parts <- c(parts, paste0('<figure><a href="', href(path), '"><img loading="lazy" src="',
      href(path), '" alt="', html_escape(basename(path)), '"></a><figcaption>', html_escape(relative(path)), '</figcaption></figure>'))
    parts <- c(parts, '<details><summary>All artifacts, tables and figures</summary><ul>',
               paste0('<li>', vapply(all_files, link, character(1)), '</li>'), '</ul></details>')
  }
  parts <- c(parts, '</section>')
}
write_csv(bind_rows(c(list(tibble(resolution = character(), path = character(), bytes = double(), modified = character())), artifact_rows)),
          tab_path(cfg, 'artifact_index.csv'))
step_finish(cfg)
parts <- c(parts, '<section id="provenance"><h2>Report provenance</h2>',
  paste0('<p>', link(file.path(cfg$step_dir, 'params_used.yaml'), 'Resolved configuration'), ' · ',
         link(file.path(cfg$step_dir, 'sessionInfo.txt'), 'R session'), ' · ',
         link(tab_path(cfg, 'artifact_index.csv'), 'Artifact index and timestamps'), '</p>'),
  '<p>Relative links keep this report usable when the entire output directory is copied. Large RDS files remain external.</p></section></body></html>')
writeLines(parts, file.path(cfg$output_dir, 'index.html'), useBytes = TRUE)
log_msg(cfg, 'Report: ', file.path(cfg$output_dir, 'index.html'))
