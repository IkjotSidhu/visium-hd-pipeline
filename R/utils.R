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
    max_mt_pct   = 50                     # max mitochondrial %; 100 disables
  ),
  clustering = list(
    n_variable_features = 2000,
    sketch_cells        = 25000,          # cells sketched PER LAYER (= per sample/integration group)
    pca_dims            = 30,
    resolutions         = c(0.2, 0.5, 1, 2),
    main_resolution     = 1,
    integration         = 'auto',         # auto (Harmony when >1 donor/sample) | harmony | none
    integration_by      = 'donor'         # sample-sheet column Harmony corrects for (auto falls back to replicate_id)
  ),
  markers = list(top_n = 10, min_pct = 0.1, logfc_threshold = 0.25),
  banksy = list(run = TRUE, lambda = 0.2, k_geom = 15, resolutions = c(0.5, 1), pca_dims = 30),
  de = list(
    run            = TRUE,
    condition_col  = 'condition',
    reference      = NULL,                # reference level of condition_col
    covariates     = list(),              # e.g. donor (paired design)
    by_cluster     = TRUE,                # also test within each cluster at main_resolution
    min_bins       = 50                   # pseudobulk sample needs >= this many bins to be used
  ),
  rctd = list(
    run            = FALSE,
    reference_rds  = NULL,                # Seurat object of a single-cell reference
    celltype_col   = NULL,                # metadata column with cell-type labels
    min_cells_type = 25,                  # drop reference types with fewer cells
    umi_min        = 50,
    mode           = 'doublet',
    max_cores      = 4
  ),
  plots = list(image_alpha = 0.4, width = 8, height = 8),
  lsf = list(
    project = 'acc_naiklab', queue = 'premium',
    resources = list(                     # cores / memory per core (MB) / walltime per step
      '01_assign_replicates' = list(cores = 1, mem_mb = 32000,  walltime = '02:00'),
      '02_load_qc'           = list(cores = 1, mem_mb = 64000,  walltime = '04:00'),
      '03_cluster'           = list(cores = 4, mem_mb = 32000,  walltime = '12:00'),
      '04_markers'           = list(cores = 4, mem_mb = 16000,  walltime = '06:00'),
      '05_banksy'            = list(cores = 2, mem_mb = 64000,  walltime = '12:00'),
      '06_pseudobulk_de'     = list(cores = 1, mem_mb = 64000,  walltime = '04:00'),
      '07_rctd'              = list(cores = 4, mem_mb = 75000,  walltime = '48:00')
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
  if (!cfg$clustering$main_resolution %in% cfg$clustering$resolutions)
    errs <- c(errs, 'clustering$main_resolution must be one of clustering$resolutions')
  if (!cfg$clustering$integration %in% c('auto', 'harmony', 'none'))
    errs <- c(errs, 'clustering$integration must be auto, harmony or none')
  if (isTRUE(cfg$rctd$run) && (is.null(cfg$rctd$reference_rds) || !file.exists(cfg$rctd$reference_rds)))
    errs <- c(errs, 'rctd$reference_rds must exist when rctd$run is true')

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
    if (cfg$clustering$integration == 'harmony' && !cfg$clustering$integration_by %in% colnames(ss))
      errs <- c(errs, paste('clustering$integration_by not in sample_sheet:', cfg$clustering$integration_by))
  }
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
  cfg <- load_config(args[1], bin_size = if (length(args) > 1) args[2] else NULL)

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
# Step 01 runs once; 02-07 once per resolution; switched-off steps are skipped.
plan_jobs <- function(cfg, first = '01', last = '07', res_filter = '') {
  steps <- c('01_assign_replicates', '02_load_qc', '03_cluster', '04_markers', '05_banksy', '06_pseudobulk_de', '07_rctd')
  keep  <- c(layout_splits(cfg), TRUE, TRUE, TRUE, isTRUE(cfg$banksy$run), isTRUE(cfg$de$run), isTRUE(cfg$rctd$run))
  num   <- as.integer(substr(steps, 1, 2))
  steps <- steps[keep & num >= as.integer(first) & num <= as.integer(last)]

  res <- as.character(unlist(cfg$bin_sizes))
  if (nzchar(res_filter)) res <- intersect(res, strsplit(res_filter, ',')[[1]])

  logdir <- file.path(cfg$output_dir, 'logs')
  line <- function(step, r) {
    rs <- cfg$lsf$resources[[step]]
    paste(step, r, rs$cores, rs$mem_mb, rs$walltime, cfg$lsf$queue, cfg$lsf$project, logdir)
  }
  c(if ('01_assign_replicates' %in% steps) line('01_assign_replicates', res[1]),
    unlist(lapply(res, \(r) vapply(setdiff(steps, '01_assign_replicates'), line, character(1), r = r))))
}

# ---- Paths ----
# Replicate assignment is shared by every resolution; everything else is per resolution
step_dir <- function(cfg, step) {
  if (step == '01_assign_replicates') file.path(cfg$output_dir, step) else file.path(cfg$output_dir, cfg$res_label, step)
}
rds_dir  <- function(cfg) file.path(cfg$output_dir, cfg$res_label, 'rds')
rds_path <- function(cfg, step) file.path(rds_dir(cfg), paste0(step, '.rds'))
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
# Okabe-Ito first (colourblind-safe), then extra distinct hues for larger sets
categorical_palette <- c('#0072B2', '#D55E00', '#009E73', '#CC79A7', '#E69F00', '#56B4E9', '#F0E442', '#332288',
                         '#8C564B', '#7F7F7F', '#BCBD22', '#17BECF', '#AEC7E8', '#FFBB78', '#98DF8A', '#FF9896',
                         '#C5B0D5', '#C49C94', '#F7B6D2', '#DBDB8D', '#9EDAE5', '#393B79', '#637939', '#8C6D31',
                         '#843C39', '#7B4173', '#3182BD', '#E6550D', '#31A354', '#756BB1', '#636363', '#FD8D3C',
                         '#74C476', '#9E9AC8', '#FDAE6B', '#A1D99B', '#BCBDDC', '#D6616B', '#E7969C', '#CE6DBD')

# Named colour vector in a fixed order (levels keep their position across plots)
make_palette <- function(levels) {
  levels <- as.character(levels)
  # Beyond the hand-picked set (e.g. resolution 2 clusterings), interpolate - neighbours get less distinct
  cols <- if (length(levels) <= length(categorical_palette)) categorical_palette[seq_along(levels)] else
    grDevices::colorRampPalette(categorical_palette)(length(levels))
  setNames(cols, levels)
}

# Cluster IDs sorted numerically so palettes and legends read 0, 1, 2, ...
sort_levels <- function(x) {
  u <- unique(as.character(x))
  if (all(grepl('^[0-9]+$', u))) u[order(as.integer(u))] else sort(u)
}

# ---- Spatial plotting ----
# Hires H&E (dimmed) with each bin/cell drawn on top, rasterised so PDFs stay
# small. Uses x_fullres / y_fullres / library_id from metadata, so it works the
# same for bins and segmented cells. Returns one ggplot per library.
.image_cache <- new.env()
load_hires <- function(sr_dir) {
  key <- normalizePath(sr_dir)
  if (is.null(.image_cache[[key]])) {
    sp  <- spatial_dir(sr_dir)
    img <- png::readPNG(file.path(sp, 'tissue_hires_image.png'))
    sf  <- jsonlite::read_json(file.path(sp, 'scalefactors_json.json'))
    # Downsample to <= 2500 px on the long side: embedding the full ~6000 px hires
    # raster made every page ~25 MB; scale factor adjusted to match
    f <- max(1, ceiling(max(dim(img)[1:2]) / 2500))
    img <- img[seq(1, dim(img)[1], by = f), seq(1, dim(img)[2], by = f), , drop = FALSE]
    sf$tissue_hires_scalef <- sf$tissue_hires_scalef / f
    .image_cache[[key]] <- list(img = img, sf = sf)
  }
  .image_cache[[key]]
}

plot_spatial <- function(md, color_by, cfg, sheet, palette = NULL, continuous = FALSE,
                         limits = NULL, title_prefix = NULL) {
  lib_dirs <- sheet |> distinct(library_id, spaceranger_dir)
  map(set_names(lib_dirs$library_id), function(lib) {
    hi  <- load_hires(lib_dirs$spaceranger_dir[lib_dirs$library_id == lib])
    img <- hi$img[, , 1:3]
    img <- img * cfg$plots$image_alpha + (1 - cfg$plots$image_alpha)   # blend toward white = dim
    sc  <- hi$sf$tissue_hires_scalef

    # Bin/cell positions in hires pixels (image y runs downward, so negate it)
    df <- md |>
      filter(library_id == lib) |>
      mutate(px = x_fullres * sc, py = -y_fullres * sc)

    # Bins drawn as tiles of their true width; cells as points
    if (cfg$res_label == 'cells') {
      layer <- ggrastr::rasterise(geom_point(size = 0.5, stroke = 0), dpi = 300)
    } else {
      w <- as.numeric(cfg$bin_size) / hi$sf$microns_per_pixel * sc
      layer <- ggrastr::rasterise(geom_tile(width = w, height = w), dpi = 300)
    }
    aes_fill <- if (cfg$res_label == 'cells') aes(colour = .data[[color_by]]) else aes(fill = .data[[color_by]])

    p <- ggplot(df, aes(px, py)) +
      annotation_raster(img, xmin = 0, xmax = dim(img)[2], ymin = -dim(img)[1], ymax = 0) +
      layer + aes_fill +
      coord_fixed(xlim = range(df$px) + c(-50, 50), ylim = range(df$py) + c(-50, 50)) +
      labs(title = paste(c(title_prefix, lib), collapse = ' - '), x = NULL, y = NULL) +
      theme_void() +
      theme(plot.background = element_rect(fill = 'white', colour = NA))

    is_cells <- cfg$res_label == 'cells'
    if (continuous) {
      # Sequential, light-low (Blues): low values recede into the tissue
      blues <- RColorBrewer::brewer.pal(9, 'Blues')
      p + if (is_cells) scale_colour_gradientn(colours = blues, limits = limits, oob = scales::squish, name = color_by)
          else          scale_fill_gradientn(colours = blues, limits = limits, oob = scales::squish, name = color_by)
    } else {
      p + if (is_cells) scale_colour_manual(values = palette, breaks = names(palette), name = color_by, guide = guide_legend(override.aes = list(size = 3)))
          else          scale_fill_manual(values = palette, breaks = names(palette), name = color_by)
    }
  })
}

# Save a list of plots one per page of a PDF plus a PNG of each for quick viewing
save_plot_list <- function(plots, cfg, stem, width = cfg$plots$width, height = cfg$plots$height) {
  pdf(fig_path(cfg, paste0(stem, '.pdf')), width = width, height = height)
  walk(plots, print)
  invisible(dev.off())
  iwalk(plots, \(p, nm) ggsave(fig_path(cfg, paste0(stem, '_', nm, '.png')), p, width = width, height = height, dpi = 100))
}
