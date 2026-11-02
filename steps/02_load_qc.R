# ---- 02_load_qc.R ----
# Load one resolution (8 um bins, 16 um bins or Space Ranger 4 segmented cells)
# for every capture area, keep bins/cells assigned to a replicate, attach the
# sample sheet, compute QC metrics, plot them, filter, and save one merged object.
#
# Usage: Rscript steps/02_load_qc.R config.yaml [8|16|cells]

# ---- Setup ----
source(file.path(dirname(sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value = TRUE))), '..', 'R', 'utils.R'))
cfg   <- step_setup('02_load_qc')
sheet <- read_sample_sheet(cfg)

# QC thresholds for this resolution
min_counts   <- qc_value(cfg, 'min_counts')
min_features <- qc_value(cfg, 'min_features')
max_mt_pct   <- qc_value(cfg, 'max_mt_pct')
mt_pattern   <- if (cfg$species == 'mouse') '^mt-' else '^MT-'
unit_label   <- if (cfg$res_label == 'cells') 'cells' else paste(cfg$res_label, 'bins')
log_msg(cfg, 'QC: UMI >= ', min_counts, ', genes >= ', min_features, ', MT% <= ', max_mt_pct)

# ---- Loaders ----
# Binned data: Seurat's loader with the hires image (default is lowres);
# slice = library so capture areas get distinct image names
load_bins <- function(sr_dir, lib) {
  obj <- Load10X_Spatial(sr_dir, bin.size = as.integer(cfg$bin_size), slice = lib,
                         image.name = 'tissue_hires_image.png')
  xy  <- GetTissueCoordinates(obj)
  obj$x_fullres <- xy$x[match(colnames(obj), xy$cell)]
  obj$y_fullres <- xy$y[match(colnames(obj), xy$cell)]
  obj
}

# Segmented cells: counts from segmented_outputs/, centroid = mean position of the
# 2 um squares inside each cell (Seurat's loader doesn't read SR4 cells)
load_cells <- function(sr_dir, lib) {
  counts <- Read10X_h5(file.path(sr_dir, 'segmented_outputs/filtered_feature_cell_matrix.h5'))
  pos_2um <- nanoparquet::read_parquet(file.path(sr_dir, 'binned_outputs/square_002um/spatial/tissue_positions.parquet'),
                                       col_select = c('barcode', 'pxl_row_in_fullres', 'pxl_col_in_fullres'))
  centroids <- nanoparquet::read_parquet(file.path(sr_dir, 'barcode_mappings.parquet'),
                                         col_select = c('square_002um', 'cell_id', 'in_cell')) |>
    filter(in_cell, cell_id %in% colnames(counts)) |>
    inner_join(pos_2um, by = c(square_002um = 'barcode')) |>
    group_by(cell_id) |>
    summarise(x_fullres = mean(pxl_col_in_fullres), y_fullres = mean(pxl_row_in_fullres), area_2um_squares = n())
  obj <- CreateSeuratObject(counts, assay = 'Spatial.cells')
  idx <- match(colnames(obj), centroids$cell_id)
  obj$x_fullres <- centroids$x_fullres[idx]
  obj$y_fullres <- centroids$y_fullres[idx]
  obj$cell_area_um2 <- centroids$area_2um_squares[idx] * 4    # each 2 um square = 4 um^2
  obj
}

# ---- Load each capture area ----
loading_rows <- list()
obj_list <- map(set_names(unique(sheet$library_id)), function(lib) {
  lib_rows <- sheet |> filter(library_id == lib)
  sr_dir   <- lib_rows$spaceranger_dir[1]
  obj <- if (cfg$res_label == 'cells') load_cells(sr_dir, lib) else load_bins(sr_dir, lib)

  # Replicate per bin/cell: from step 01, or the library's single sample-sheet row
  if (layout_splits(cfg)) {
    assign_file <- file.path(step_dir(cfg, '01_assign_replicates'), 'tables', sprintf('%s_%s.csv.gz', lib, cfg$res_label))
    if (!file.exists(assign_file)) stop('Missing ', assign_file, ' - run step 01_assign_replicates first')
    assign <- read_csv(assign_file, show_col_types = FALSE) |> filter(!is.na(replicate_id))
    rep_id <- assign$replicate_id[match(colnames(obj), assign$barcode)]
  } else {
    rep_id <- rep(lib_rows$replicate_id, ncol(obj))
  }
  n_before <- ncol(obj)
  keep <- !is.na(rep_id)
  n_assigned <- sum(keep)

  # Test mode: keep a random fraction of each replicate
  if (!is.null(cfg$test_fraction)) {
    idx  <- which(keep)
    sel  <- unlist(lapply(split(idx, rep_id[idx]), \(i) i[sample.int(length(i), ceiling(length(i) * cfg$test_fraction))]))
    keep <- seq_along(rep_id) %in% sel
    log_msg(cfg, lib, ': test mode - keeping ', round(100 * cfg$test_fraction), '% of each replicate')
  }
  # Loading flow for the QC report: loaded -> assigned to a replicate -> selected (test mode)
  loading_rows[[lib]] <<- tibble(library_id = lib, n_loaded = n_before, n_assigned = n_assigned,
                                 n_unassigned = n_before - n_assigned, n_selected = sum(keep))
  obj <- subset(obj, cells = colnames(obj)[keep])
  rep_id <- rep_id[keep]
  log_msg(cfg, lib, ': ', ncol(obj), ' of ', n_before, ' ', unit_label, ' assigned to a replicate')

  # Sample-sheet columns (donor, condition, ...) onto every bin/cell
  md <- sheet |> select(-spaceranger_dir)
  md <- md[match(rep_id, md$replicate_id), ]
  for (col in colnames(md)) obj[[col]] <- md[[col]]

  # Mitochondrial %
  obj[['percent.mt']] <- PercentageFeatureSet(obj, pattern = mt_pattern)
  obj
})

# Count columns are named after the assay (e.g. nCount_Spatial.008um)
assay     <- DefaultAssay(obj_list[[1]])
count_col <- paste0('nCount_', assay)
feat_col  <- paste0('nFeature_', assay)

# ---- QC table + removal report ----
qc_df <- map_dfr(obj_list, \(o) o[[]] |> rownames_to_column('barcode')) |>
  mutate(pass = .data[[count_col]] >= min_counts & .data[[feat_col]] >= min_features &
                replace_na(percent.mt, 0) <= max_mt_pct)

# Per replicate: what each filter removes
qc_summary <- qc_df |>
  group_by(across(any_of(c('replicate_id', 'library_id', 'donor', cfg$de$condition_col)))) |>
  summarise(units          = n(),
            median_umi     = median(.data[[count_col]]),
            median_genes   = median(.data[[feat_col]]),
            median_mt_pct  = round(median(percent.mt, na.rm = TRUE), 2),
            fail_umi       = sum(.data[[count_col]] < min_counts),
            fail_mt        = sum(replace_na(percent.mt, 0) > max_mt_pct),
            pct_kept       = round(100 * mean(pass), 1),
            .groups = 'drop')
print(qc_summary, width = Inf)
write_csv(qc_summary, tab_path(cfg, 'qc_summary.csv'))

# Per-unit metrics and pass/fail, then the group/sample/replicate QC report (R/qc_reports.R)
write_csv(qc_df |> mutate(unit_id = paste(library_id, barcode, sep = '_'), .before = 1),
          tab_path(cfg, 'qc_unit_metrics.csv.gz'))
reports <- build_qc_reports(qc_df, sheet, cfg, count_col, feat_col, bind_rows(loading_rows))
reports$relationships <- write_qc_relationship_plots(qc_df, cfg, count_col, feat_col)   # UMI vs genes plots
write_qc_reports(reports, cfg)

# Flag uneven loss between conditions (biases composition comparisons)
if (cfg$de$condition_col %in% colnames(qc_df)) {
  by_cond <- qc_df |> group_by(cond = .data[[cfg$de$condition_col]]) |> summarise(pct_kept = 100 * mean(pass))
  if (diff(range(by_cond$pct_kept)) > 15)
    warning('QC keeps unevenly across conditions: ',
            paste0(by_cond$cond, ' ', round(by_cond$pct_kept, 1), '%', collapse = ', '),
            ' - keep this in mind for composition comparisons', call. = FALSE)
}

# ---- QC distributions ----
rep_pal <- make_palette(sheet$replicate_id)
qc_violin <- function(y, title, hline = NULL, log = FALSE) {
  p <- ggplot(qc_df, aes(replicate_id, .data[[y]], fill = replicate_id)) +
    geom_violin(scale = 'width', linewidth = 0.2) +
    scale_fill_manual(values = rep_pal, guide = 'none') +
    labs(title = title, x = NULL, y = NULL) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  if (log) p <- p + scale_y_continuous(transform = scales::pseudo_log_trans(base = 10), breaks = c(0, 10, 100, 1000, 10000))  # handles 0 counts
  if (!is.null(hline)) p <- p + geom_hline(yintercept = hline, linetype = 2, colour = two_colours[2])
  p
}
p_dist <- qc_violin(count_col, 'UMIs', hline = min_counts, log = TRUE) |
  qc_violin(feat_col, 'Genes', hline = if (min_features > 0) min_features, log = TRUE) |
  qc_violin('percent.mt', 'Mitochondrial %', hline = if (max_mt_pct < 100) max_mt_pct)
ggsave(fig_path(cfg, 'qc_violins.pdf'), p_dist, width = 13, height = 5)
ggsave(fig_path(cfg, 'qc_violins.png'), p_dist, width = 13, height = 5, dpi = 100)

# ---- Spatial QC (before filtering) ----
# Raw values on a linear colour scale from 0 to the 99th percentile (higher values shown at the cap)
qc_df$qc_status <- ifelse(qc_df$pass, 'pass', 'fail')
save_spatial(qc_df, count_col, cfg, sheet, 'qc_spatial_umi', continuous = TRUE, label = 'UMIs')
save_spatial(qc_df, 'percent.mt', cfg, sheet, 'qc_spatial_percent_mt', continuous = TRUE, label = 'Mitochondrial %')
save_spatial(qc_df, 'qc_status', cfg, sheet, 'qc_spatial_status', palette = setNames(two_colours, c('pass', 'fail')))

# ---- Filter + merge ----
obj_list <- map(obj_list, function(o) {
  keep <- o[[count_col]][, 1] >= min_counts & o[[feat_col]][, 1] >= min_features &
          replace_na(o$percent.mt, 0) <= max_mt_pct
  subset(o, cells = colnames(o)[keep])
})

# Bin barcodes / cell IDs repeat across capture areas, so prefix with the library
merged <- if (length(obj_list) > 1) merge(obj_list[[1]], obj_list[-1], add.cell.ids = names(obj_list)) else
  RenameCells(obj_list[[1]], add.cell.id = names(obj_list))
log_msg(cfg, 'Merged: ', ncol(merged), ' ', unit_label, ' x ', nrow(merged), ' genes')

# ---- SpotSweeper: label local outliers (nothing is removed) ----
# Each retained unit is compared with its n_neighbors nearest units in the same
# strip (robust z-score of log1p values); |z| beyond the cutoff in the bad
# direction flags it. Flags and z-scores go into the metadata for later steps.
if (isTRUE(cfg$spotsweeper$run)) {
  ss <- cfg$spotsweeper
  qc2 <- secondary_qc_cfg(cfg)                                          # tables/figures go to <res>/secondary_qc/
  tests <- tibble(key = c('umi', 'genes', 'mito'), column = c(count_col, feat_col, 'percent.mt'),
                  direction = c('lower', 'lower', 'higher'), label = unname(spotsweeper_labels[c('umi', 'genes', 'mito')])) |>
    filter(key %in% ss$metrics)
  md <- merged[[]] |> rownames_to_column('cell') |> mutate(percent.mt = replace_na(percent.mt, 0))
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
  for (col in setdiff(names(flags), c('cell', 'replicate_id', 'library_id'))) merged[[col]] <- flags[[col]]
  write_csv(flags, tab_path(qc2, 'spotsweeper_flags.csv.gz'))

  # Per-strip summary: units and % flagged by each metric and by any
  summary_ss <- flags |> group_by(replicate_id) |>
    summarise(n_units = n(), across(all_of(c(flag_cols, 'spotsweeper_any')), list(n = sum, pct = \(x) round(100 * mean(x), 2))))
  write_csv(summary_ss, tab_path(qc2, 'spotsweeper_summary.csv'))
  log_msg(cfg, 'SpotSweeper: ', sum(flags$spotsweeper_any), ' of ', nrow(flags), ' ', unit_label, ' flagged (cutoff ',
          ss$cutoff, ', ', ss$n_neighbors, ' neighbours)')

  # Maps of the flags over the H&E (per capture area and per strip)
  levels_ss <- c('not flagged', tests$label, 'multiple')
  pal_ss <- c('not flagged' = highlight[['ns']], make_palette(c(tests$label, 'multiple')))
  map_df <- md |> select(cell, replicate_id, library_id, x_fullres, y_fullres) |>
    mutate(spotsweeper = factor(flags$spotsweeper, levels = levels_ss))
  save_spatial(map_df, 'spotsweeper', qc2, sheet, 'spotsweeper_flags', palette = pal_ss[levels_ss], label = 'SpotSweeper')
}

saveRDS(merged, rds_path(cfg, '02_load_qc'))
step_finish(cfg)
