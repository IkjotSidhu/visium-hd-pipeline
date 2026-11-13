# ---- utils.R ----
# Shared helpers for every pipeline step: config loading + validation, step
# setup/teardown (library path, seed, logging, provenance), output paths,
# sample-sheet handling, colour palettes and hires spatial plotting.
# Each step sources this file and calls step_setup() first.

# ---- Default configuration ----
# Every key a step reads has a default here; the YAML config only needs to set
# what differs. See config_template.yaml for documentation of each key.
default_config <- list(
  project_name = 'VisiumHD_project',
  output_dir   = 'output',
  r_library    = '/sc/arion/projects/naiklab/ikjot/R',
  seed         = 1234,
  species      = 'human',
  bin_sizes    = list(8, 16, 'cells'),    # resolutions to analyse: 8, 16 and/or "cells" (Space Ranger 4 segmentation)
  sample_sheet = NULL,                    # CSV: one row per replicate (required)
  replicates = list(
    layout          = 'none',             # none | left_to_right | top_to_bottom | manual
    manual_dir      = NULL,               # <library_id>.csv files (barcode, replicate_id) for layout = manual
    block_um        = 100,                # strip detection: block size
    min_fill        = 0.25,               # strip detection: fraction of a block that must be tissue
    link_blocks     = 2,                  # strip detection: dense blocks this close join one strip
    attach_blocks   = 1                   # strip detection: edge blocks this close attach to a strip
  ),
  qc = list(                              # each value: one number for all resolutions, or a map
    min_counts   = 50,                    #   per resolution, e.g. {"8": 50, "16": 100, cells: 50}
    min_features = 0,                     # min genes per bin/cell
    max_mt_pct   = 100,                   # max mitochondrial %; 100 = no filter (default)
    report = list(sample_col = 'library_id', group_cols = NULL, max_rows = 100)
  ),
  clustering = list(
    n_variable_features = 2000,
    sketch_cells        = 25000,          # cells sketched PER LAYER (= per sample/integration group)
    pca_dims            = 30,
    resolutions         = (1:20) / 10,
    main_resolution     = 1,
    integration         = 'auto',         # auto (Harmony when >1 donor/sample) | harmony | none
    integration_by      = 'donor'         # sample-sheet column Harmony corrects for (auto falls back to replicate_id)
  ),
  annotation = list(run = FALSE, mapping_files = list(), require_complete = FALSE),   # step 08: turn on to apply reviewed labels
  pathways = list(run = TRUE, gmt_file = NULL, collection = NULL, version = NULL,
                  gene_id_type = 'symbol', min_size = 15, max_size = 500, eps = 1e-10),
  report = list(run = TRUE, max_table_rows = 20),
  markers = list(top_n = 10, min_pct = 0.1, logfc_threshold = 0.25),
  spotsweeper = list(run = TRUE, metrics = c('umi', 'genes', 'mito'),    # step 02: label (never remove) local outliers, as a check
                     n_neighbors = 36, cutoff = 3),          # SpotSweeper package defaults
  banksy = list(run = TRUE, lambda = 0.2, k_geom = 15, resolutions = (1:20) / 10, pca_dims = 30),
  de = list(
    run            = TRUE,
    condition_col  = 'condition',
    reference      = NULL,                # reference level of condition_col
    covariates     = list(),              # e.g. donor (paired design)
    by_cluster     = TRUE,                # also test within each cluster at main_resolution
    min_bins       = 50,                  # pseudobulk sample needs >= this many bins to be used
    covariate_policy = 'strict',          # strict skips unsupported designs; drop explicitly permits fallback
    min_strata     = 2,                   # minimum paired strata for the consistency flag
    diagnostics = list(run = TRUE, top_variable_genes = 500,
                       metadata_cols = c('donor', 'condition', 'tissue', 'region', 'library_id')),
    interpretation_notes = list()         # study-specific limits copied alongside results
  ),
  rctd = list(
    run            = FALSE,
    reference_rds  = NULL,                # Seurat object of a single-cell reference
    celltype_col   = NULL,                # metadata column with cell-type labels
    reference_assay = 'RNA',              # raw, uncorrected counts; never inferred from DefaultAssay
    force_recompute = FALSE,             # bypass even a matching cache
    min_cells_type = 25,                  # drop reference types with fewer cells
    umi_min        = 50,
    mode           = 'doublet',
    max_cores      = 4
  ),
  plots = list(image_alpha = 1,   width = 8, height = 8, dpi = 300,   # PNG resolution
               point_size = NULL, cell_diameter_um = 10,     # NULL = dots sized to one bin / cell
               strip_um_per_px = 0.8, strip_height = 12,     # per-strip maps: H&E resolution, page height (in)
               gdal_setup = 'module load gdal/3.9.1'),       # makes gdal_translate available (reads the original scan)
  lsf = list(
    project = 'acc_naiklab', queue = 'premium',
    resources = list(                     # cores / memory per core (MB) / walltime per step
      '01_assign_replicates' = list(cores = 1, mem_mb = 32000,  walltime = '02:00'),
      '02_load_qc'           = list(cores = 1, mem_mb = 64000,  walltime = '04:00'),
      '03_cluster'           = list(cores = 4, mem_mb = 32000,  walltime = '12:00'),
      '04_markers'           = list(cores = 4, mem_mb = 16000,  walltime = '06:00'),
      '05_banksy'            = list(cores = 2, mem_mb = 64000,  walltime = '12:00'),
      '06_pseudobulk_de'     = list(cores = 1, mem_mb = 64000,  walltime = '04:00'),
      '07_rctd'              = list(cores = 4, mem_mb = 75000,  walltime = '48:00'),
      '08_annotation'        = list(cores = 1, mem_mb = 16000, walltime = '02:00'),
      '09_pathways'          = list(cores = 1, mem_mb = 16000, walltime = '02:00'),
      '10_report'            = list(cores = 1, mem_mb = 8000, walltime = '01:00')
    )
  )
)

# ---- Config ----
# Read the YAML config, merge onto the defaults, resolve relative paths against
# the config file's directory, and fail early on anything missing.
load_config <- function(config_path, bin_size = NULL) {
  if (!file.exists(config_path)) stop('Config file not found: ', config_path)
  user_cfg <- yaml::read_yaml(config_path)
  cfg <- utils::modifyList(default_config, user_cfg, keep.null = TRUE)

  cfg$clustering$resolutions <- unlist(cfg$clustering$resolutions)
  cfg$banksy$resolutions <- unlist(cfg$banksy$resolutions)
  cfg$spotsweeper$metrics <- unlist(cfg$spotsweeper$metrics)

  # The resolution this step runs at (passed by run_pipeline.sh); defaults to the first listed
  cfg$bin_size <- as.character(bin_size %||% cfg$bin_sizes[[1]])

  # Relative paths are relative to the directory holding the config
  cfg$config_path <- normalizePath(config_path)
  cfg$base_dir    <- dirname(cfg$config_path)
  resolve <- function(p) if (is.null(p) || grepl('^/', p)) p else file.path(cfg$base_dir, p)
  cfg$output_dir            <- resolve(cfg$output_dir)

  # Test mode (run_pipeline.sh -t): random fraction of each replicate, separate output folder
  tf <- Sys.getenv('VHD_TEST_FRACTION', '')
  cfg$test_fraction <- if (nzchar(tf)) as.numeric(tf) else NULL
  if (!is.null(cfg$test_fraction)) cfg$output_dir <- paste0(cfg$output_dir, '_test')
  cfg$sample_sheet          <- resolve(cfg$sample_sheet)
  cfg$replicates$manual_dir <- resolve(cfg$replicates$manual_dir)
  cfg$rctd$reference_rds    <- resolve(cfg$rctd$reference_rds)
  resolve_mapping <- function(x) if (is.list(x)) lapply(x, resolve_mapping) else resolve(x)
  cfg$annotation$mapping_files <- resolve_mapping(cfg$annotation$mapping_files)
  cfg$pathways$gmt_file <- resolve(cfg$pathways$gmt_file)

  # Resolution label used in file names and output folders
  cfg$res_label <- if (cfg$bin_size == 'cells') 'cells' else sprintf('%03dum', as.integer(cfg$bin_size))

  validate_config(cfg)
  cfg
}

# Checks that would otherwise only fail hours into a job
validate_config <- function(cfg) {
  errs <- character()
  if (is.null(cfg$sample_sheet) || !file.exists(cfg$sample_sheet))
    errs <- c(errs, paste('sample_sheet not found:', cfg$sample_sheet))
  if (!all(as.character(unlist(cfg$bin_sizes)) %in% c('8', '16', 'cells')))
    errs <- c(errs, 'bin_sizes may only contain 8, 16 and "cells"')
  if (!cfg$bin_size %in% c('8', '16', 'cells'))
    errs <- c(errs, 'bin size must be 8, 16 or "cells"')
  if (!cfg$replicates$layout %in% c('none', 'left_to_right', 'top_to_bottom', 'manual'))
    errs <- c(errs, 'replicates$layout must be none, left_to_right, top_to_bottom or manual')
  if (cfg$replicates$layout == 'manual' && (is.null(cfg$replicates$manual_dir) || !dir.exists(cfg$replicates$manual_dir)))
    errs <- c(errs, 'replicates$manual_dir must exist when layout = manual')
  for (res in list(cfg$clustering$resolutions, cfg$banksy$resolutions)) {
    values <- unlist(res)
    if (!is.numeric(values) || !length(values) || any(!is.finite(values) | values <= 0) || anyDuplicated(values))
      errs <- c(errs, 'Clustering/BANKSY resolutions must be unique positive finite numbers')
  }
  if (!cfg$clustering$main_resolution %in% cfg$clustering$resolutions)
    errs <- c(errs, 'clustering$main_resolution must be one of clustering$resolutions')
  if (!cfg$clustering$integration %in% c('auto', 'harmony', 'none'))
    errs <- c(errs, 'clustering$integration must be auto, harmony or none')
  if (isTRUE(cfg$rctd$run) && (is.null(cfg$rctd$reference_rds) || !file.exists(cfg$rctd$reference_rds)))
    errs <- c(errs, 'rctd$reference_rds must exist when rctd$run is true')

  if (!cfg$de$covariate_policy %in% c('strict', 'drop'))
    errs <- c(errs, 'de$covariate_policy must be strict or drop')
  if (length(cfg$de$min_strata) != 1 || !is.numeric(cfg$de$min_strata) ||
      !is.finite(cfg$de$min_strata) || cfg$de$min_strata < 2 || cfg$de$min_strata %% 1 != 0)
    errs <- c(errs, 'de$min_strata must be an integer >= 2')
  if (isTRUE(cfg$rctd$run)) {
    for (key in c('celltype_col', 'reference_assay'))
      if (!is.character(cfg$rctd[[key]]) || length(cfg$rctd[[key]]) != 1 || !nzchar(cfg$rctd[[key]]))
        errs <- c(errs, paste0('rctd$', key, ' must be a non-empty string'))
    if (!cfg$rctd$mode %in% c('doublet', 'full', 'multi'))
      errs <- c(errs, 'rctd$mode must be doublet, full or multi')
    if (!is.logical(cfg$rctd$force_recompute) || length(cfg$rctd$force_recompute) != 1 || is.na(cfg$rctd$force_recompute))
      errs <- c(errs, 'rctd$force_recompute must be true or false')
  }

  # QC report settings
  sample_col <- cfg$qc$report$sample_col
  if (!is.character(sample_col) || length(sample_col) != 1 || is.na(sample_col) || !nzchar(sample_col))
    errs <- c(errs, 'qc$report$sample_col must name one sample-sheet column')
  if (!is.numeric(cfg$qc$report$max_rows) || length(cfg$qc$report$max_rows) != 1 ||
      !is.finite(cfg$qc$report$max_rows) || cfg$qc$report$max_rows < 1 || cfg$qc$report$max_rows %% 1 != 0)
    errs <- c(errs, 'qc$report$max_rows must be a positive integer')
  # Sample sheet columns and Space Ranger folders
  if (!is.null(cfg$sample_sheet) && file.exists(cfg$sample_sheet)) {
    ss <- utils::read.csv(cfg$sample_sheet, check.names = FALSE)
    miss <- setdiff(c('replicate_id', 'library_id', 'spaceranger_dir'), colnames(ss))
    if (length(miss)) errs <- c(errs, paste('sample_sheet is missing columns:', paste(miss, collapse = ', ')))
    if (anyDuplicated(ss$replicate_id)) errs <- c(errs, 'sample_sheet replicate_id values must be unique')
    if ('spaceranger_dir' %in% colnames(ss)) {
      dirs <- unique(ifelse(grepl('^/', ss$spaceranger_dir), ss$spaceranger_dir, file.path(cfg$base_dir, ss$spaceranger_dir)))
      bad <- dirs[!dir.exists(file.path(dirs, 'binned_outputs'))]
      if (length(bad)) errs <- c(errs, paste('no binned_outputs/ in:', paste(bad, collapse = ', ')))
    }
    if (layout_splits(cfg) == FALSE && any(duplicated(ss$library_id)))
      errs <- c(errs, 'several replicates share a library_id but replicates$layout is none')
    if (isTRUE(cfg$de$run) && !cfg$de$condition_col %in% colnames(ss))
      errs <- c(errs, paste('de$condition_col not in sample_sheet:', cfg$de$condition_col))
    if (isTRUE(cfg$de$run)) {
      de_cols <- unique(c(cfg$de$condition_col, unlist(cfg$de$covariates)))
      missing_de <- setdiff(de_cols, colnames(ss))
      if (length(missing_de)) errs <- c(errs, paste('DE columns missing:', paste(missing_de, collapse = ', ')))
      for (v in intersect(de_cols, colnames(ss)))
        if (any(is.na(ss[[v]]) | trimws(as.character(ss[[v]])) == ''))
          errs <- c(errs, paste('Missing DE metadata in column:', v))
      if (!is.null(cfg$de$reference) && cfg$de$condition_col %in% colnames(ss) &&
          !cfg$de$reference %in% ss[[cfg$de$condition_col]])
        errs <- c(errs, 'de$reference is not present in the sample sheet')
      if (cfg$de$condition_col %in% unlist(cfg$de$covariates))
        errs <- c(errs, 'de$condition_col cannot also be a covariate')
    }
    # QC report columns must exist, and every row needs a sample identifier
    qc_cols <- unique(c(sample_col, unlist(cfg$qc$report$group_cols)))
    missing_qc <- setdiff(qc_cols, colnames(ss))
    if (length(missing_qc)) errs <- c(errs, paste('QC reporting columns missing:', paste(missing_qc, collapse = ', ')))
    if (length(sample_col) == 1 && !is.na(sample_col) && sample_col %in% colnames(ss) &&
        any(is.na(ss[[sample_col]]) | trimws(as.character(ss[[sample_col]])) == ''))
      errs <- c(errs, 'QC sample identifiers must not be missing')
    if (cfg$clustering$integration == 'harmony' && !cfg$clustering$integration_by %in% colnames(ss))
      errs <- c(errs, paste('clustering$integration_by not in sample_sheet:', cfg$clustering$integration_by))
  }
  for (p in unlist(cfg$annotation$mapping_files))
    if (!file.exists(p)) errs <- c(errs, paste('annotation mapping not found:', p))
  if (!is.null(cfg$pathways$gmt_file)) {
    if (!file.exists(cfg$pathways$gmt_file)) errs <- c(errs, 'pathways$gmt_file not found')
    for (key in c('collection', 'version', 'gene_id_type'))
      if (!is.character(cfg$pathways[[key]]) || length(cfg$pathways[[key]]) != 1 || !nzchar(cfg$pathways[[key]]))
        errs <- c(errs, paste0('pathways$', key, ' is required with a GMT file'))
  }
  positive_integer <- function(x) is.numeric(x) && length(x) == 1 && is.finite(x) && x >= 1 && x %% 1 == 0
  for (x in list(cfg$pathways$min_size, cfg$pathways$max_size, cfg$de$diagnostics$top_variable_genes, cfg$report$max_table_rows))
    if (!positive_integer(x)) errs <- c(errs, 'Pathway sizes, PCA gene count and report row limit must be positive integers')
  if (cfg$pathways$min_size > cfg$pathways$max_size) errs <- c(errs, 'pathways min_size exceeds max_size')
  if (!is.numeric(cfg$pathways$eps) || length(cfg$pathways$eps) != 1 || !is.finite(cfg$pathways$eps) || cfg$pathways$eps < 0)
    errs <- c(errs, 'pathways$eps must be nonnegative')
  if (length(errs)) stop('Config problems:\n  - ', paste(errs, collapse = '\n  - '), call. = FALSE)
  invisible(TRUE)
}

# TRUE when capture areas are split into several replicates
layout_splits <- function(cfg) cfg$replicates$layout != 'none'

# Harmony decision: returns the sample-sheet column to correct for, or NULL for none.
# auto/harmony correct for integration_by (or replicate_id when that column is absent)
# whenever it has more than one value; none never corrects.
resolve_integration <- function(cfg, sheet) {
  cl <- cfg$clustering
  if (cl$integration == 'none') return(NULL)
  by <- if (cl$integration_by %in% colnames(sheet)) cl$integration_by else 'replicate_id'
  if (by != cl$integration_by) log_msg(cfg, 'No "', cl$integration_by, '" column in the sample sheet - integrating by replicate_id')
  if (n_distinct(sheet[[by]]) < 2) {
    log_msg(cfg, 'Only one ', by, ' - no integration needed')
    return(NULL)
  }
  by
}

# A QC threshold for the current resolution: either one value or a per-resolution map
qc_value <- function(cfg, key) {
  v <- cfg$qc[[key]]
  if (is.list(v)) {
    if (is.null(v[[cfg$bin_size]])) stop('qc$', key, ' has no value for resolution ', cfg$bin_size, call. = FALSE)
    v[[cfg$bin_size]]
  } else v
}

# ---- Step setup / teardown ----
# Every step starts with: cfg <- step_setup('03_cluster')
step_setup <- function(step) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 1) stop('Usage: Rscript steps/', step, '.R <config.yaml> [8|16|cells]', call. = FALSE)
  cfg <- load_config(args[1], bin_size = if (length(args) > 1 && step != '10_report') args[2] else NULL)

  # Project R library (compiled against the Rocky9 R/4.4.1 module)
  if (!is.null(cfg$r_library)) .libPaths(c(cfg$r_library, .libPaths()))
  suppressPackageStartupMessages({
    library(Seurat)
    library(tidyverse)
    library(patchwork)
  })

  # Reproducibility: same seed for every step
  set.seed(cfg$seed)

  # Parallel workers = cores LSF gave this job (1 outside LSF); BLAS stays at 1 thread
  n_cores <- as.integer(Sys.getenv('LSB_DJOB_NUMPROC', '1'))
  if (n_cores > 1) future::plan('multicore', workers = n_cores) else future::plan('sequential')
  options(future.globals.maxSize = 64 * 1024^3, future.seed = TRUE)

  cfg$pipeline_dir <- .pipeline_dir
  git_value <- function(arguments) tryCatch(paste(system2('git', c('-C', shQuote(.pipeline_dir), arguments),
      stdout = TRUE, stderr = FALSE), collapse = '\n'), error = function(e) 'unavailable')
  cfg$pipeline_commit <- git_value(c('rev-parse', 'HEAD'))
  cfg$pipeline_worktree <- if (nzchar(git_value(c('status', '--porcelain')))) 'modified' else 'clean'
  cfg$step       <- step
  cfg$n_cores    <- n_cores
  cfg$step_dir   <- step_dir(cfg, step)
  cfg$start_time <- Sys.time()
  for (d in file.path(cfg$step_dir, c('figures', 'tables'))) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  dir.create(rds_dir(cfg), recursive = TRUE, showWarnings = FALSE)

  log_msg(cfg, '==== ', step, ' | ', cfg$project_name, ' | ', cfg$res_label, ' | cores: ', n_cores, ' ====')
  cfg
}

# Record what ran: parameters + session info next to the step's outputs
step_finish <- function(cfg) {
  yaml::write_yaml(cfg[setdiff(names(cfg), 'start_time')], file.path(cfg$step_dir, 'params_used.yaml'))
  writeLines(capture.output(sessionInfo()), file.path(cfg$step_dir, 'sessionInfo.txt'))
  log_msg(cfg, 'Done in ', format(round(difftime(Sys.time(), cfg$start_time, units = 'mins'), 1)))
}

log_msg <- function(cfg, ...) message(format(Sys.time(), '%H:%M:%S'), ' [', cfg$step %||% 'setup', '] ', ...)

# ---- Job plan for run_pipeline.sh ----
# One line per job: step res cores mem_mb walltime queue project logdir.
# Steps 01 and 10 run once; 02-09 run per resolution; disabled steps are skipped.
plan_jobs <- function(cfg, first = '01', last = '10', res_filter = '') {
  steps <- c('01_assign_replicates', '02_load_qc', '03_cluster', '04_markers', '05_banksy', '06_pseudobulk_de', '07_rctd', '08_annotation', '09_pathways', '10_report')
  keep  <- c(layout_splits(cfg), TRUE, TRUE, TRUE, isTRUE(cfg$banksy$run), isTRUE(cfg$de$run), isTRUE(cfg$rctd$run),
             isTRUE(cfg$annotation$run), isTRUE(cfg$pathways$run) && isTRUE(cfg$de$run), isTRUE(cfg$report$run))
  num   <- as.integer(substr(steps, 1, 2))
  steps <- steps[keep & num >= as.integer(first) & num <= as.integer(last)]

  res <- as.character(unlist(cfg$bin_sizes))
  if (nzchar(res_filter)) res <- intersect(res, strsplit(res_filter, ',')[[1]])
  if (!length(res)) stop('No configured resolutions match the requested resolution filter')

  logdir <- file.path(cfg$output_dir, 'logs')
  line <- function(step, r) {
    rs <- cfg$lsf$resources[[step]]
    paste(step, r, rs$cores, rs$mem_mb, rs$walltime, cfg$lsf$queue, cfg$lsf$project, logdir)
  }
  c(if ('01_assign_replicates' %in% steps) line('01_assign_replicates', res[1]),
    unlist(lapply(res, \(r) vapply(setdiff(steps, c('01_assign_replicates', '10_report')), line, character(1), r = r))),
    if ('10_report' %in% steps && length(res)) line('10_report', 'all'))
}

# ---- Paths ----
# Replicate assignment is shared by every resolution; everything else is per resolution
step_dir <- function(cfg, step) {
  root <- if (step %in% c('01_assign_replicates', '10_report')) file.path(cfg$output_dir, step) else file.path(cfg$output_dir, cfg$res_label, step)
  if (!is.null(cfg$cluster_resolution) && step %in% c('03_cluster', '04_markers', '05_banksy', '08_annotation'))
    root <- file.path(root, paste0('res_', resolution_label(cfg$cluster_resolution)))
  root
}
rds_dir  <- function(cfg) file.path(cfg$output_dir, cfg$res_label, 'rds')
rds_path <- function(cfg, step) {
  suffix <- if (step == '08_annotation' && !is.null(cfg$cluster_resolution)) paste0('_res_', resolution_label(cfg$cluster_resolution)) else ''
  file.path(rds_dir(cfg), paste0(step, suffix, '.rds'))
}
# SpotSweeper flag labels, in plotting order (same colours in step 02 maps and step 03 bars)
spotsweeper_labels <- c(umi = 'low UMIs', genes = 'low genes', mito = 'high mito %')

# Secondary QC (SpotSweeper) output, kept apart from the main step folders: <res>/secondary_qc/
secondary_qc_cfg <- function(cfg) {
  cfg$step_dir <- file.path(cfg$output_dir, cfg$res_label, 'secondary_qc')
  for (d in c('figures', 'tables')) dir.create(file.path(cfg$step_dir, d), recursive = TRUE, showWarnings = FALSE)
  cfg
}
# ---- Per-cluster QC summary (secondary QC) ----
# Table + box plots of UMIs, genes and mito % for each cluster of one resolution,
# written to <res>/secondary_qc/. Linear scale; boxes
# are drawn from precomputed quartiles (Tukey whiskers, outlier points not drawn).
write_cluster_qc <- function(md, column, r, cfg) {
  qc2 <- secondary_qc_cfg(cfg)
  suffix <- paste0('_res_', resolution_label(r))
  assay <- paste0('Spatial.', cfg$res_label)                         # e.g. Spatial.008um
  metrics <- c(UMIs = paste0('nCount_', assay), Genes = paste0('nFeature_', assay), `Mito %` = 'percent.mt')
  long <- md |>
    select(cluster = all_of(column), all_of(unname(metrics))) |>
    pivot_longer(-cluster, names_to = 'column', values_to = 'value') |>
    mutate(metric = factor(names(metrics)[match(column, metrics)], levels = names(metrics)),
           value = replace_na(value, 0))                             # mito % is NA for bins with no counts
  # Table: n, median and quartiles per cluster and metric
  tbl <- long |>
    group_by(cluster, metric) |>
    summarise(median = median(value), q25 = quantile(value, 0.25), q75 = quantile(value, 0.75), .groups = 'drop') |>
    mutate(metric = recode(metric, UMIs = 'umi', Genes = 'genes', `Mito %` = 'mito_pct')) |>
    pivot_wider(names_from = metric, values_from = c(median, q25, q75), names_glue = '{metric}_{.value}') |>
    left_join(md |> count(cluster = .data[[column]], name = 'n_units'), by = 'cluster') |>
    relocate(n_units, .after = cluster) |>
    mutate(across(where(is.double), \(x) round(x, 2)))
  write_csv(tbl, tab_path(qc2, paste0('cluster_qc', suffix, '.csv')))
  # Box statistics per cluster (whiskers = most extreme value within 1.5 IQR of the box)
  box <- long |>
    group_by(metric, cluster) |>
    summarise(q25 = quantile(value, 0.25), q50 = median(value), q75 = quantile(value, 0.75),
              lo = min(value[value >= q25 - 1.5 * (q75 - q25)]), hi = max(value[value <= q75 + 1.5 * (q75 - q25)]),
              .groups = 'drop')
  p <- ggplot(box, aes(cluster, fill = cluster)) +
    geom_boxplot(aes(ymin = lo, lower = q25, middle = q50, upper = q75, ymax = hi), stat = 'identity') +
    facet_wrap(~metric, ncol = 1, scales = 'free_y') +
    scale_fill_manual(values = make_palette(levels(factor(box$cluster))), guide = 'none') +
    labs(title = paste('QC metrics per cluster - resolution', resolution_label(r)),
         x = 'Cluster', y = NULL) +
    theme_classic()
  ggsave(fig_path(qc2, paste0('cluster_qc', suffix, '.pdf')), p, width = 12, height = 9)
  ggsave(fig_path(qc2, paste0('cluster_qc', suffix, '.png')), p, width = 12, height = 9, dpi = 100)
  invisible(tbl)
}

# ---- SpotSweeper check (secondary QC; labels only, never filters) ----
# Output folder: <res>/secondary_qc/spotsweeper/{figures,tables}
spotsweeper_cfg <- function(cfg) {
  cfg$step_dir <- file.path(cfg$output_dir, cfg$res_label, 'secondary_qc', 'spotsweeper')
  for (d in c('figures', 'tables')) dir.create(file.path(cfg$step_dir, d), recursive = TRUE, showWarnings = FALSE)
  cfg
}

# Each unit is compared with its n_neighbors nearest units in the same strip
# (robust z-score of log1p values); |z| beyond the cutoff in the bad direction
# flags it. md needs: cell, replicate_id, library_id, x_fullres, y_fullres,
# nCount_/nFeature_<assay> and percent.mt. Writes flags, a per-strip summary,
# % flagged per strip, z-score histograms and H&E maps; returns the flag table.
run_spotsweeper <- function(md, cfg, sheet) {
  ss <- cfg$spotsweeper
  out <- spotsweeper_cfg(cfg)
  assay <- paste0('Spatial.', cfg$res_label)
  tests <- tibble(key = c('umi', 'genes', 'mito'), column = c(paste0('nCount_', assay), paste0('nFeature_', assay), 'percent.mt'),
                  direction = c('lower', 'lower', 'higher'), label = unname(spotsweeper_labels[c('umi', 'genes', 'mito')])) |>
    filter(key %in% ss$metrics)
  md <- md |> mutate(percent.mt = replace_na(percent.mt, 0))            # mito % is NA for bins with no counts
  spe <- SpatialExperiment::SpatialExperiment(
    colData = S4Vectors::DataFrame(md |> select(cell, replicate_id, x_fullres, y_fullres, all_of(tests$column)), row.names = md$cell),
    spatialCoordsNames = c('x_fullres', 'y_fullres'))
  for (i in seq_len(nrow(tests)))
    spe <- SpotSweeper::localOutliers(spe, metric = tests$column[i], direction = tests$direction[i], n_neighbors = ss$n_neighbors,
                                      samples = 'replicate_id', log = TRUE, cutoff = ss$cutoff)
  cd <- as.data.frame(SummarizedExperiment::colData(spe))[md$cell, ]     # localOutliers reorders by strip

  # One z-score and one flag column per metric, a combined flag and a single label for maps
  flags <- md |> select(cell, replicate_id, library_id)
  for (i in seq_len(nrow(tests))) {
    flags[[paste0('spotsweeper_', tests$key[i], '_z')]]    <- cd[[paste0(tests$column[i], '_z')]]
    flags[[paste0('spotsweeper_', tests$key[i], '_flag')]] <- cd[[paste0(tests$column[i], '_outliers')]]
  }
  flag_cols <- paste0('spotsweeper_', tests$key, '_flag')
  n_flags <- rowSums(as.matrix(flags[flag_cols]))
  flags$spotsweeper_any <- n_flags > 0
  flags$spotsweeper <- case_when(n_flags == 0 ~ 'not flagged', n_flags > 1 ~ 'multiple',
                                 TRUE ~ tests$label[max.col(as.matrix(flags[flag_cols]), ties.method = 'first')])
  write_csv(flags, tab_path(out, 'spotsweeper_flags.csv.gz'))

  # Per-strip summary: units and % flagged by each metric and by any
  summary_ss <- flags |> group_by(replicate_id) |>
    summarise(n_units = n(), across(all_of(c(flag_cols, 'spotsweeper_any')), list(n = sum, pct = \(x) round(100 * mean(x), 2))))
  write_csv(summary_ss, tab_path(out, 'spotsweeper_summary.csv'))

  # % flagged per strip, per metric and any metric
  metric_pal <- make_palette(c(tests$label, 'any metric'))
  strip_long <- flags |> pivot_longer(all_of(c(flag_cols, 'spotsweeper_any')), names_to = 'metric', values_to = 'flagged') |>
    mutate(metric = factor(c(setNames(tests$label, flag_cols), spotsweeper_any = 'any metric')[metric], levels = names(metric_pal))) |>
    group_by(replicate_id, metric) |> summarise(pct = 100 * mean(flagged), .groups = 'drop')
  p_strip <- ggplot(strip_long, aes(replicate_id, pct, fill = metric)) +
    geom_col(position = position_dodge(preserve = 'single')) +
    scale_fill_manual(values = metric_pal) +
    labs(title = paste0('SpotSweeper local outliers per strip (cutoff ', ss$cutoff, ', ', ss$n_neighbors, ' neighbours)'),
         x = NULL, y = '% of units flagged', fill = 'Flagged for') +
    theme_classic()
  ggsave(fig_path(out, 'spotsweeper_flagged_by_strip.pdf'), p_strip, width = 10, height = 4)
  ggsave(fig_path(out, 'spotsweeper_flagged_by_strip.png'), p_strip, width = 10, height = 4, dpi = 150)

  # z-score histograms per strip and metric, with the cutoff
  zl <- flags |> pivot_longer(ends_with('_z'), names_to = 'key', values_to = 'z') |>
    mutate(key = sub('spotsweeper_(.*)_z', '\\1', key)) |>
    left_join(tests |> select(key, label), by = 'key') |> mutate(label = factor(label, levels = tests$label))
  cut_lines <- tests |> mutate(z = if_else(direction == 'lower', -ss$cutoff, ss$cutoff), label = factor(label, levels = tests$label))
  p_z <- ggplot(zl, aes(z)) +
    geom_histogram(bins = 100, fill = 'grey40') +
    geom_vline(data = cut_lines, aes(xintercept = z), colour = highlight[['sig']]) +
    scale_y_continuous(trans = 'log1p', breaks = c(0, 10, 100, 1000, 1e4, 1e5)) +
    facet_grid(replicate_id ~ label, scales = 'free_x') +
    labs(title = paste0('SpotSweeper robust z-scores vs ', ss$n_neighbors, ' nearest neighbours (red line = cutoff ', ss$cutoff, ')'),
         x = 'Robust z-score of log1p(metric)', y = 'Units (log scale)') +
    theme_classic()
  ggsave(fig_path(out, 'spotsweeper_zscores.pdf'), p_z, width = 11, height = 9)
  ggsave(fig_path(out, 'spotsweeper_zscores.png'), p_z, width = 11, height = 9, dpi = 150)

  # Maps of the flags over the H&E (per capture area and per strip)
  levels_ss <- c('not flagged', tests$label, 'multiple')
  pal_ss <- c('not flagged' = highlight[['ns']], make_palette(c(tests$label, 'multiple')))
  map_df <- md |> select(cell, replicate_id, library_id, x_fullres, y_fullres) |>
    mutate(spotsweeper = factor(flags$spotsweeper, levels = levels_ss))
  save_spatial(map_df, 'spotsweeper', out, sheet, 'spotsweeper_flags', palette = pal_ss[levels_ss], label = 'SpotSweeper')
  flags
}

# Share of SpotSweeper-flagged units per cluster at one resolution (table + bars)
write_spotsweeper_by_cluster <- function(md, column, r, cfg) {
  out <- spotsweeper_cfg(cfg)
  suffix <- paste0('_res_', resolution_label(r))
  flag_cols <- grep('^spotsweeper_.*_flag$|^spotsweeper_any$', colnames(md), value = TRUE)
  ss_tbl <- md |> group_by(cluster = .data[[column]]) |>
    summarise(n_units = n(), across(all_of(flag_cols), list(n = sum, pct = \(x) round(100 * mean(x), 2))))
  write_csv(ss_tbl, tab_path(out, paste0('spotsweeper_by_cluster', suffix, '.csv')))
  ss_levels <- c(unname(spotsweeper_labels[cfg$spotsweeper$metrics]), 'multiple')
  ss_comp <- md |> count(cluster = .data[[column]], spotsweeper) |>
    group_by(cluster) |> mutate(pct = 100 * n / sum(n)) |> ungroup() |> filter(spotsweeper != 'not flagged') |>
    mutate(spotsweeper = factor(spotsweeper, levels = ss_levels))
  p_ss <- ggplot(ss_comp, aes(cluster, pct, fill = spotsweeper)) +
    geom_col() +
    scale_x_discrete(drop = FALSE) +                                   # clusters with no flags shown at 0
    scale_fill_manual(values = make_palette(ss_levels), drop = FALSE) +
    labs(title = paste('SpotSweeper-flagged units per cluster - resolution', resolution_label(r)),
         x = 'Cluster', y = '% of cluster units flagged', fill = 'SpotSweeper') +
    theme_classic()
  ggsave(fig_path(out, paste0('spotsweeper_by_cluster', suffix, '.pdf')), p_ss, width = 12, height = 5)
  ggsave(fig_path(out, paste0('spotsweeper_by_cluster', suffix, '.png')), p_ss, width = 12, height = 5, dpi = 100)
  invisible(ss_tbl)
}

fig_path <- function(cfg, file) file.path(cfg$step_dir, 'figures', file)
tab_path <- function(cfg, file) file.path(cfg$step_dir, 'tables', file)

# Read the previous step's object, with a clear error if it hasn't run
read_step_rds <- function(cfg, step) {
  p <- rds_path(cfg, step)
  if (!file.exists(p)) stop('Missing ', p, ' - run step ', step, ' first', call. = FALSE)
  log_msg(cfg, 'Reading ', basename(p))
  readRDS(p)
}

# ---- Sample sheet ----
# One row per replicate; spaceranger_dir resolved to an absolute path
read_sample_sheet <- function(cfg) {
  utils::read.csv(cfg$sample_sheet, check.names = FALSE) |>
    as_tibble() |>
    mutate(spaceranger_dir = ifelse(grepl('^/', spaceranger_dir), spaceranger_dir, file.path(cfg$base_dir, spaceranger_dir)))
}

# Spatial folder holding the hires image + scale factors for a library
spatial_dir <- function(sr_dir) file.path(sr_dir, 'binned_outputs/square_008um/spatial')

# ---- Palettes ----
# Lab standard swatches (Development/ColorPalettes.pdf). make_palette() uses the
# smallest swatch that fits; past 18 groups it continues with categorical_palette.
lab_swatches <- list(
  s3 = c('#237D26', '#FCAF17', '#DC217D', '#6786C4'),                                   # 4, colourblind friendly
  s4 = c('#682E88', '#BFA58A', '#F28500', '#27B14B', '#025A48', '#8DCBEF', '#C355A0', '#791316'),
  s1 = c('#B6DEF6', '#6FAB43', '#FDD700', '#8151A1', '#FAA51A', '#AFD243', '#595A5A', '#D3B081',
         '#4B5DAA', '#ED2224', '#A8703D', '#A290AA'),
  s6 = c('#38A12B', '#FFC81C', '#C355A0', '#5570B6', '#8DCBEF', '#AD96C7', '#CC843D', '#682E88',
         '#DF6FA1', '#F27549', '#7F8080', '#CBDF7A', '#DCE3E2'),
  s2 = c('#8DCBEF', '#5570B6', '#CC843D', '#AD96C7', '#EB2627', '#7F8080', '#DEBEB4', '#682E88', '#791316',
         '#C355A0', '#FFC81C', '#27B14B', '#8B5B2C', '#838B8A', '#025A48', '#BFCCCA', '#F27549', '#CBDF7A'))
two_colours  <- c('#5570B6', '#EB2627')   # swatch 5 "bright": any two groups
ctrl_treated <- c('#999999', '#C2000B')   # swatch 5 "ctrl vs treated": reference first
highlight    <- c(ns = '#999999', sig = '#C2000B', weak = '#FEA621')   # swatch 7 (significance highlights)

# Fallback past the lab swatches: Okabe-Ito, then tab10/tab20 hues
categorical_palette <- c('#0072B2', '#D55E00', '#009E73', '#CC79A7', '#E69F00', '#56B4E9', '#F0E442', '#332288',
                         '#8C564B', '#7F7F7F', '#BCBD22', '#17BECF', '#AEC7E8', '#FFBB78', '#98DF8A', '#FF9896',
                         '#C5B0D5', '#C49C94', '#F7B6D2', '#DBDB8D', '#9EDAE5', '#393B79', '#637939', '#8C6D31',
                         '#843C39', '#7B4173', '#3182BD', '#E6550D', '#31A354', '#756BB1', '#636363', '#FD8D3C',
                         '#74C476', '#9E9AC8', '#FDAE6B', '#A1D99B', '#BCBDDC', '#D6616B', '#E7969C', '#CE6DBD')

# Named colour vector in level order. With a reference level (e.g. the DE
# control) two groups get grey for the reference and red for the other.
make_palette <- function(levels, reference = NULL) {
  levels <- as.character(levels)
  n <- length(levels)
  if (n == 2 && !is.null(reference) && reference %in% levels)
    return(setNames(ctrl_treated, c(reference, setdiff(levels, reference))))
  pool <- if (n == 2) two_colours else
    Find(\(s) length(s) >= n, lab_swatches[c('s3', 's4', 's1', 's6', 's2')]) %||%
    c(lab_swatches$s2, setdiff(categorical_palette, c(lab_swatches$s2, '#7F7F7F', '#636363')))   # past 18: lab swatch, then fallback (minus greys that duplicate the swatch grey)
  # Beyond every listed colour, interpolate - neighbours get less distinct
  cols <- if (n <= length(pool)) pool[seq_len(n)] else grDevices::colorRampPalette(pool)(n)
  setNames(cols, levels)
}

# Cluster IDs sorted numerically so palettes and legends read 0, 1, 2, ...
sort_levels <- function(x) {
  u <- unique(as.character(x))
  if (all(grepl('^[0-9]+$', u))) u[order(as.integer(u))] else sort(u)
}

# ---- Spatial plotting ----
# Every bin/cell is drawn as a vector dot over the H&E, so spots stay sharp at
# any zoom. Positions come from x_fullres / y_fullres / library_id, so bins and
# segmented cells work the same. Two layouts:
#   plot_spatial()        - one page per capture area, over Space Ranger's hires image
#   plot_spatial_strips() - one page per replicate strip, over a crop of the
#                           original scan written by step 01 (write_strip_images)
.image_cache <- new.env()
load_hires <- function(sr_dir) {
  key <- normalizePath(sr_dir)
  if (is.null(.image_cache[[key]])) {
    sp <- spatial_dir(sr_dir)
    .image_cache[[key]] <- list(img = png::readPNG(file.path(sp, 'tissue_hires_image.png')),   # full resolution, no downsampling
                                sf  = jsonlite::read_json(file.path(sp, 'scalefactors_json.json')))
  }
  .image_cache[[key]]
}

# Colour limits for continuous maps: 0 to the 99th percentile on a linear scale,
# so a few extreme bins don't flatten the rest; values above are shown at the cap
spatial_limits <- function(x) c(0, unname(quantile(x, 0.99, na.rm = TRUE)))

# One panel: img spans x0..x1, y0..y1 in image pixels (y negated); df has px/py in
# the same pixels. unit_px = one bin (or cell) in image pixels sets the dot size.
spatial_panel <- function(df, img, x0, x1, y0, y1, unit_px, color_by, cfg, width, height, title,
                          palette = NULL, continuous = FALSE, limits = NULL, label = color_by) {
  img <- img * cfg$plots$image_alpha + (1 - cfg$plots$image_alpha)   # 1 = H&E as scanned; lower fades toward white

  # Dot diameter = one bin, converted to mm on the page. The panel is roughly the
  # page minus legend/title; ggplot point size ~ diameter / 0.75 mm.
  mm_per_px <- min((width * 25.4 - 40) / (x1 - x0), (height * 25.4 - 15) / (y1 - y0))
  size <- cfg$plots$point_size %||% (0.9 * unit_px * mm_per_px / 0.75)

  p <- ggplot(df, aes(px, py, colour = .data[[color_by]])) +
    annotation_raster(img, xmin = x0, xmax = x1, ymin = -y1, ymax = -y0) +
    geom_point(shape = 16, size = size, stroke = 0) +
    coord_fixed(xlim = c(x0, x1), ylim = c(-y1, -y0), expand = FALSE) +
    labs(title = title, x = NULL, y = NULL) +
    theme_void() +
    theme(plot.background = element_rect(fill = 'white', colour = NA))

  if (continuous) {
    # Sequential, light-low (Blues), linear; capped values noted in the legend
    p + scale_colour_gradientn(colours = RColorBrewer::brewer.pal(9, 'Blues'), limits = limits, oob = scales::squish,
                               name = paste0(label, '\n(capped at ', signif(limits[2], 3), ')'))
  } else {
    p + scale_colour_manual(values = palette, breaks = names(palette), name = label,
                            guide = guide_legend(override.aes = list(size = 3)))
  }
}

# Unit (bin / cell) size in micrometres for dot sizing
unit_um <- function(cfg) if (cfg$res_label == 'cells') cfg$plots$cell_diameter_um else as.numeric(cfg$bin_size)

plot_spatial <- function(md, color_by, cfg, sheet, palette = NULL, continuous = FALSE,
                         limits = NULL, title_prefix = NULL, label = color_by,
                         width = cfg$plots$width, height = cfg$plots$height) {
  if (continuous && is.null(limits)) limits <- spatial_limits(md[[color_by]])
  lib_dirs <- sheet |> distinct(library_id, spaceranger_dir)
  map(set_names(lib_dirs$library_id), function(lib) {
    hi <- load_hires(lib_dirs$spaceranger_dir[lib_dirs$library_id == lib])
    sc <- hi$sf$tissue_hires_scalef

    # Bin/cell positions in hires pixels (image y runs downward, so negate it)
    df <- md |>
      filter(library_id == lib) |>
      mutate(px = x_fullres * sc, py = -y_fullres * sc)

    # Crop the H&E to the plotted tissue (+ margin) so only that part is embedded
    m  <- 50
    x0 <- max(1, floor(min(df$px) - m));  x1 <- min(dim(hi$img)[2], ceiling(max(df$px) + m))
    y0 <- max(1, floor(-max(df$py) - m)); y1 <- min(dim(hi$img)[1], ceiling(-min(df$py) + m))
    spatial_panel(df, hi$img[y0:y1, x0:x1, 1:3], x0 - 1, x1, y0 - 1, y1,
                  unit_um(cfg) / hi$sf$microns_per_pixel * sc, color_by, cfg, width, height,
                  paste(c(title_prefix, lib), collapse = ' - '), palette, continuous, limits, label)
  })
}

# ---- Per-strip H&E images (written once by step 01) ----
strip_image_dir <- function(cfg) file.path(cfg$output_dir, '01_assign_replicates', 'strip_images')

# Original microscope scan given to spaceranger --image, if a copy is inside the
# Space Ranger folder (the original path is often not readable)
find_scan <- function(sr_dir) {
  cmd <- file.path(sr_dir, '_cmdline')
  if (!file.exists(cmd)) return(NULL)
  image <- regmatches(readLines(cmd, warn = FALSE), regexpr('--image=[^ ]+', readLines(cmd, warn = FALSE)))
  if (!length(image)) return(NULL)
  hits <- list.files(sr_dir, pattern = paste0('^', basename(sub('--image=', '', image[1])), '$'),
                     recursive = TRUE, full.names = TRUE)
  if (length(hits)) hits[1] else NULL
}

# Cut each strip (bounding box of its bins + margin) out of the original scan at
# plots$strip_um_per_px with gdal_translate; without a scan (or GDAL), crop the
# hires image instead. Records where each crop sits in fullres pixels.
write_strip_images <- function(cfg, extents, sheet) {
  out <- strip_image_dir(cfg)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  rows <- pmap(extents, function(replicate_id, library_id, xmin, xmax, ymin, ymax) {
    sr_dir <- sheet$spaceranger_dir[sheet$library_id == library_id][1]
    hi     <- load_hires(sr_dir)
    m      <- 150                                                        # margin, fullres pixels
    x0 <- max(0, floor(xmin - m)); y0 <- max(0, floor(ymin - m))
    w  <- ceiling(xmax + m) - x0;  h  <- ceiling(ymax + m) - y0
    path <- file.path(out, paste0(replicate_id, '.png'))
    scan <- find_scan(sr_dir)

    # Original scan via GDAL (reads only this window), downsampled to the target resolution
    f <- cfg$plots$strip_um_per_px / hi$sf$microns_per_pixel             # fullres pixels per output pixel
    ok <- !is.null(scan) && system2('bash', c('-lc', shQuote(paste(cfg$plots$gdal_setup, '&& GDAL_PAM_ENABLED=NO gdal_translate -q -of PNG -b 1 -b 2 -b 3 -r average',
      '-srcwin', x0, y0, w, h, '-outsize', round(w / f), round(h / f), shQuote(scan), shQuote(path))))) == 0
    if (ok) {
      source <- scan
    } else {
      # Fallback: crop of Space Ranger's hires image
      log_msg(cfg, replicate_id, ': original scan not available, using the hires image')
      sc  <- hi$sf$tissue_hires_scalef
      cx0 <- max(1, floor(x0 * sc) + 1); cy0 <- max(1, floor(y0 * sc) + 1)
      cx1 <- min(dim(hi$img)[2], ceiling((x0 + w) * sc)); cy1 <- min(dim(hi$img)[1], ceiling((y0 + h) * sc))
      png::writePNG(hi$img[cy0:cy1, cx0:cx1, 1:3], path)
      x0 <- (cx0 - 1) / sc; y0 <- (cy0 - 1) / sc; w <- (cx1 - cx0 + 1) / sc
      source <- 'tissue_hires_image.png'
    }
    img_w <- dim(png::readPNG(path, info = FALSE))[2]
    tibble(replicate_id, library_id, x0, y0, fullres_per_px = w / img_w, image = basename(path), source)
  })
  write_csv(bind_rows(rows), file.path(out, 'strip_images.csv'))
}

plot_spatial_strips <- function(md, color_by, cfg, sheet, palette = NULL, continuous = FALSE,
                                limits = NULL, title_prefix = NULL, label = color_by, height = cfg$plots$strip_height) {
  index <- file.path(strip_image_dir(cfg), 'strip_images.csv')
  if (!file.exists(index)) return(list())
  strips <- read_csv(index, show_col_types = FALSE) |>
    filter(replicate_id %in% md$replicate_id, !replicate_id %in% sheet$library_id)   # no name clash with capture-area pages
  if (continuous && is.null(limits)) limits <- spatial_limits(md[[color_by]])
  plots <- pmap(strips, function(replicate_id, library_id, x0, y0, fullres_per_px, image, source) {
    img <- png::readPNG(file.path(strip_image_dir(cfg), image))
    mpp <- load_hires(sheet$spaceranger_dir[sheet$library_id == library_id][1])$sf$microns_per_pixel
    df  <- md[md$replicate_id == replicate_id, ] |>
      mutate(px = (x_fullres - x0) / fullres_per_px, py = -(y_fullres - y0) / fullres_per_px)
    # Page shaped like the strip, plus room for the legend
    width <- min(14, max(5, height * dim(img)[2] / dim(img)[1] + 2.5))
    p <- spatial_panel(df, img, 0, dim(img)[2], 0, dim(img)[1], unit_um(cfg) / mpp / fullres_per_px,
                       color_by, cfg, width, height, paste(c(title_prefix, replicate_id), collapse = ' - '),
                       palette, continuous, limits, label)
    attr(p, 'page') <- c(width, height)                                  # save_plot_list uses this page size
    p
  })
  set_names(plots, strips$replicate_id)
}

# Capture-area pages plus one page per replicate strip, as <stem>_<library or replicate>.pdf/.png
save_spatial <- function(md, color_by, cfg, sheet, stem, ..., width = cfg$plots$width, height = cfg$plots$height) {
  save_plot_list(plot_spatial(md, color_by, cfg, sheet, ..., width = width, height = height), cfg, stem, width, height)
  save_plot_list(plot_spatial_strips(md, color_by, cfg, sheet, ...), cfg, stem)
}

# Save each plot in a list (e.g. one per sample/library) as its own PDF and PNG: <stem>_<name>.pdf/.png.
# A plot carrying a 'page' attribute (width, height) is saved at that size.
save_plot_list <- function(plots, cfg, stem, width = cfg$plots$width, height = cfg$plots$height) {
  iwalk(plots, function(p, nm) {
    page <- attr(p, 'page') %||% c(width, height)
    ggsave(fig_path(cfg, paste0(stem, '_', nm, '.pdf')), p, width = page[1], height = page[2])
    ggsave(fig_path(cfg, paste0(stem, '_', nm, '.png')), p, width = page[1], height = page[2], dpi = cfg$plots$dpi)
  })
}

# ---- Analysis safeguards ----
# A requested adjustment is preserved unless the user explicitly opts into a
# different design. Returning a reason lets step 06 audit every skipped subset.
choose_de_design <- function(cd, cond, covars, policy = 'strict') {
  use_cov <- covars
  repeat {
    single <- use_cov[vapply(use_cov, \(v) nlevels(factor(cd[[v]])) < 2, logical(1))]
    fml <- reformulate(c(use_cov, cond))
    if (length(single)) {
      reason <- paste('single-level covariate:', paste(single, collapse = ', '))
    } else {
      mm <- model.matrix(fml, cd)
      reason <- if (qr(mm)$rank < ncol(mm)) 'confounded design' else
        if (nrow(mm) <= ncol(mm)) 'no residual degrees of freedom' else NULL
    }
    if (is.null(reason)) return(list(ok = TRUE, formula = fml, dropped = setdiff(covars, use_cov), reason = ''))
    if (policy == 'strict' || !length(use_cov))
      return(list(ok = FALSE, formula = fml, dropped = setdiff(covars, use_cov), reason = reason))
    use_cov <- if (length(single)) setdiff(use_cov, single) else head(use_cov, -1)
  }
}

# One agreeing donor is reported as available evidence, never as replication.
# All expected strata must contribute, with at least min_strata paired strata.
strata_evidence <- function(directions, effect, n_expected, min_strata = 2L) {
  n <- length(directions)
  available <- if (n) Reduce(`&`, lapply(directions, \(d) d == sign(effect))) else rep(NA, length(effect))
  supported <- n >= min_strata && n == n_expected
  list(available = available, consistent = if (supported) available else rep(NA, length(effect)),
       support = if (!n) 'none' else if (n == 1) 'single_stratum' else
         if (n < n_expected) 'partial_strata' else 'all_strata')
}

# RCTD requires raw nonnegative integer counts. Rounding corrected counts does
# not reconstruct raw RNA and can disguise an unsuitable reference assay.
validate_raw_counts <- function(counts, label = 'Reference') {
  values <- if (inherits(counts, 'sparseMatrix')) counts@x else as.vector(counts)
  if (!nrow(counts) || !ncol(counts) || any(!is.finite(values)) || any(values < 0) ||
      any(abs(values - round(values)) > 1e-8))
    stop(label, ': expected finite, nonnegative integer raw counts', call. = FALSE)
  if (is.null(rownames(counts)) || is.null(colnames(counts)) ||
      anyDuplicated(rownames(counts)) || anyDuplicated(colnames(counts)))
    stop(label, ': gene and cell names must be present and unique', call. = FALSE)
  if (any(Matrix::colSums(counts) <= 0))
    stop(label, ': zero-count cells must be removed', call. = FALSE)
  invisible(TRUE)
}

rctd_cache_key <- function(reference_key, counts, coords, settings) {
  digest::digest(list(schema = 1L, reference = reference_key, counts = counts,
                      coords = coords, settings = settings), algo = 'sha256')
}

read_rctd_cache <- function(path, key, force = FALSE) {
  if (force || !file.exists(path)) return(NULL)
  cached <- tryCatch(readRDS(path), error = function(e) NULL)
  if (!is.list(cached) || !identical(cached$key, key)) return(NULL)
  cached$result
}

# Rename on the same filesystem prevents an interrupted write from looking valid.
write_rctd_cache <- function(result, path, key) {
  tmp <- tempfile(pattern = '.rctd-', tmpdir = dirname(path))
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(list(key = key, result = result), tmp)
  if (!file.rename(tmp, path)) stop('Could not save RCTD cache: ', path)
}

# Extract explicitly selected raw counts and align annotations by cell name.
# The object's active assay is deliberately irrelevant to this operation.
rctd_reference_inputs <- function(ref_obj, assay, celltype_col) {
  if (!assay %in% SeuratObject::Assays(ref_obj)) stop('Reference assay not found: ', assay)
  if (inherits(ref_obj[[assay]], 'SCTAssay') || tolower(assay) %in% c('sct', 'integrated'))
    stop('RCTD reference must use original raw RNA counts, not SCT/integrated counts')
  if (!celltype_col %in% colnames(ref_obj[[]])) stop('Reference celltype_col not found: ', celltype_col)
  if (inherits(ref_obj[[assay]], 'Assay5'))
    ref_obj[[assay]] <- SeuratObject::JoinLayers(ref_obj[[assay]])
  counts <- SeuratObject::LayerData(ref_obj, assay = assay, layer = 'counts')
  validate_raw_counts(counts)
  if (!all(colnames(counts) %in% rownames(ref_obj[[]]))) stop('Reference counts and metadata do not align')
  list(counts = as(counts, 'dgCMatrix'),
       types = as.character(ref_obj[[]][colnames(counts), celltype_col]))
}

# ---- Output modules ----
.pipeline_dir <- normalizePath(file.path(dirname(sys.frame(1)$ofile), '..'))
source(file.path(.pipeline_dir, 'R', 'outputs.R'))

source(file.path(.pipeline_dir, 'R', 'qc_reports.R'))
