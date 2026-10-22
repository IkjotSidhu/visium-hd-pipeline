# ---- Regression checks for statistical and cache safeguards ----
# Run from the repository root with the lab R/4.4.1 module.
source('R/utils.R')
.libPaths(c(default_config$r_library, .libPaths()))
suppressPackageStartupMessages({ library(Seurat); library(tidyverse) })
expect_error <- function(expr) stopifnot(inherits(tryCatch(force(expr), error = identity), 'error'))

# ---- Requested paired designs must remain paired ----
paired <- data.frame(donor = factor(rep(c('D1', 'D2'), each = 2)),
                     condition = factor(rep(c('Control', 'Treated'), 2)))
stopifnot(choose_de_design(paired, 'condition', 'donor')$ok)
confounded <- data.frame(donor = factor(c('D1', 'D1', 'D2', 'D2')),
                         condition = factor(c('Control', 'Control', 'Treated', 'Treated')))
stopifnot(!choose_de_design(confounded, 'condition', 'donor')$ok)
fallback <- choose_de_design(confounded, 'condition', 'donor', 'drop')
stopifnot(fallback$ok, identical(fallback$dropped, 'donor'))
stopifnot(!choose_de_design(paired[1:3, ], 'condition', 'donor')$ok)
single <- paired; single$donor <- factor('D1')
stopifnot(!choose_de_design(single, 'condition', 'donor')$ok)
stopifnot(choose_de_design(single, 'condition', 'donor', 'drop')$ok)

# ---- Independent support, including disagreement and missing strata ----
one <- strata_evidence(list(c(1, -1)), c(2, -2), 2)
stopifnot(all(one$available), all(is.na(one$consistent)), one$support == 'single_stratum')
two <- strata_evidence(list(c(1, -1), c(1, 1)), c(2, -2), 2)
stopifnot(identical(two$consistent, c(TRUE, FALSE)))
partial <- strata_evidence(list(c(1, -1), c(1, -1)), c(2, -2), 3)
stopifnot(all(is.na(partial$consistent)), partial$support == 'partial_strata')
stopifnot(all(is.na(strata_evidence(list(), c(2, -2), 2)$consistent)))

# ---- Raw counts must not be reconstructed by rounding ----
counts <- Matrix::Matrix(matrix(c(1, 2, 3, 4), 2,
  dimnames = list(c('G1', 'G2'), c('C1', 'C2'))), sparse = TRUE)
validate_raw_counts(counts)
bad <- counts; bad[1, 1] <- 0.5; expect_error(validate_raw_counts(bad))
bad[1, 1] <- -1; expect_error(validate_raw_counts(bad))
bad[1, 1] <- NA_real_; expect_error(validate_raw_counts(bad))
bad <- counts; bad[, 1] <- 0; expect_error(validate_raw_counts(bad))

# Explicit RNA extraction works even when another assay is active.
ref <- CreateSeuratObject(counts)
ref$cell_type <- c('Type1', 'Type2')
ref[['other']] <- CreateAssay5Object(counts = counts * 2)
DefaultAssay(ref) <- 'other'
ref_input <- rctd_reference_inputs(ref, 'RNA', 'cell_type')
stopifnot(identical(as.matrix(ref_input$counts), as.matrix(counts)),
          identical(ref_input$types, c('Type1', 'Type2')))
expect_error(rctd_reference_inputs(ref, 'missing', 'cell_type'))
expect_error(rctd_reference_inputs(ref, 'RNA', 'missing'))
ref[['integrated']] <- CreateAssay5Object(counts = counts)
expect_error(rctd_reference_inputs(ref, 'integrated', 'cell_type'))

# ---- Cache hits require identical effective inputs ----
coords <- data.frame(x = c(1, 2), y = c(3, 4), row.names = colnames(counts))
settings <- list(mode = 'doublet', umi_min = 50, version = 'test')
key <- rctd_cache_key('reference-A', counts, coords, settings)
path <- tempfile(fileext = '.rds')
write_rctd_cache(list(fitted = TRUE), path, key)
stopifnot(isTRUE(read_rctd_cache(path, key)$fitted), is.null(read_rctd_cache(path, key, TRUE)))
changed <- counts; changed[1, 1] <- 2
stopifnot(is.null(read_rctd_cache(path, rctd_cache_key('reference-A', changed, coords, settings))))
stopifnot(is.null(read_rctd_cache(path, rctd_cache_key('reference-B', counts, coords, settings))))
changed_xy <- coords; changed_xy$x[1] <- 10
stopifnot(is.null(read_rctd_cache(path, rctd_cache_key('reference-A', counts, changed_xy, settings))))
settings$mode <- 'full'
stopifnot(is.null(read_rctd_cache(path, rctd_cache_key('reference-A', counts, coords, settings))))
saveRDS('legacy unvalidated object', path)
stopifnot(is.null(read_rctd_cache(path, key)))
writeLines('corrupt cache', path)
stopifnot(is.null(suppressWarnings(read_rctd_cache(path, key))))
unlink(path)

# ---- All-skipped results retain usable schemas ----
# Evaluate the actual step using an isolated synthetic object; no DESeq2 fit is
# needed when one condition is present. This also exercises the output writer.
step <- normalizePath('steps/06_pseudobulk_de.R')
expressions <- parse(step)
fixture <- CreateSeuratObject(counts * 10)
fixture$replicate_id <- c('R1', 'R2')
fixture$condition <- 'Control'
fixture$donor <- c('D1', 'D2')
fixture$cluster_res_1 <- factor(c('0', '0'))
test_dir <- tempfile('vhd-empty-'); dir.create(test_dir)
cfg <- default_config; cfg$de$reference <- 'Control'; cfg$de$min_bins <- 1
cfg$de$covariates <- list('donor'); cfg$step_dir <- test_dir
cfg$output_dir <- test_dir; cfg$res_label <- '008um'; cfg$step <- '06_pseudobulk_de'
cfg$start_time <- Sys.time(); cfg$bin_size <- '8'
dir.create(file.path(test_dir, 'tables')); dir.create(file.path(test_dir, 'figures'))
dir.create(rds_dir(cfg), recursive = TRUE)
env <- new.env(parent = globalenv())
env$step_setup <- function(...) cfg
env$read_step_rds <- function(...) fixture
# First expression sources utils; omit it so fixture I/O overrides remain in scope.
for (expr in as.list(expressions)[-1]) eval(expr, env)
result <- readr::read_csv(tab_path(cfg, 'pseudobulk_de_results.csv.gz'), show_col_types = FALSE)
status <- readr::read_csv(tab_path(cfg, 'pseudobulk_de_status.csv'), show_col_types = FALSE)
stopifnot(nrow(result) == 0, all(c('gene', 'padj', 'n_strata_compared') %in% names(result)),
          all(status$status == 'skipped'), nrow(status) == 2,
          file.exists(file.path(test_dir, 'sessionInfo.txt')))
# A single replicate also finishes normally and removes stale volcano figures.
fixture$replicate_id <- 'R1'; fixture$donor <- 'D1'
file.create(fig_path(cfg, 'volcano_all_bins.pdf'), fig_path(cfg, 'volcano_all_bins.png'))
for (expr in as.list(expressions)[-1]) eval(expr, env)
status <- readr::read_csv(tab_path(cfg, 'pseudobulk_de_status.csv'), show_col_types = FALSE)
stopifnot(all(status$reason == 'fewer than two eligible replicates'),
          !file.exists(fig_path(cfg, 'volcano_all_bins.pdf')),
          !file.exists(fig_path(cfg, 'volcano_all_bins.png')))
unlink(test_dir, recursive = TRUE)
# ---- RCTD wiring with a mocked fitter (not biological validation) ----
# Exercise the actual step, S4 result caching and full-mode output without an
# expensive deconvolution or inventing a biological single-cell reference.
methods::setClass('VhdTestRCTD', slots = c(results = 'list'))
test_dir <- tempfile('vhd-rctd-'); dir.create(test_dir)
cfg <- default_config; cfg$step_dir <- test_dir; cfg$output_dir <- test_dir
cfg$res_label <- '008um'; cfg$step <- '07_rctd'; cfg$bin_size <- '8'; cfg$start_time <- Sys.time()
cfg$rctd$celltype_col <- 'cell_type'; cfg$rctd$min_cells_type <- 1
cfg$rctd$reference_rds <- file.path(test_dir, 'reference.rds')
ref <- CreateSeuratObject(counts); ref$cell_type <- c('Type1', 'Type2')
saveRDS(ref, cfg$rctd$reference_rds)
dir.create(file.path(test_dir, 'tables')); dir.create(file.path(test_dir, 'figures'))
dir.create(rds_dir(cfg), recursive = TRUE)
spatial_fixture <- CreateSeuratObject(counts * 10)
spatial_fixture$replicate_id <- 'R1'; spatial_fixture$library_id <- 'L1'
spatial_fixture$x_fullres <- c(1, 2); spatial_fixture$y_fullres <- c(3, 4)
env <- new.env(parent = globalenv()); env$fit_calls <- 0L
env$step_setup <- function(...) cfg
env$read_step_rds <- function(...) spatial_fixture
env$read_sample_sheet <- function(...) data.frame()
env$Reference <- function(...) list()
env$SpatialRNA <- function(coords, counts, ...) list(cells = colnames(counts))
env$create.RCTD <- function(puck, ...) {
  env$fit_calls <- env$fit_calls + 1L
  w <- matrix(rep(c(0.7, 0.3), each = length(puck$cells)), ncol = 2,
              dimnames = list(puck$cells, c('Type1', 'Type2')))
  methods::new('VhdTestRCTD', results = list(weights = w))
}
env$run.RCTD <- function(rctd, ...) rctd
env$normalize_weights <- identity
env$plot_spatial <- function(...) list()
env$save_plot_list <- function(...) invisible(NULL)
expressions <- as.list(parse('steps/07_rctd.R'))[-1]
run_fixture <- function() for (expr in expressions) eval(expr, env)
run_fixture(); stopifnot(env$fit_calls == 1L)
run_fixture(); stopifnot(env$fit_calls == 1L)
manifest <- read_csv(tab_path(cfg, 'rctd_cache_manifest.csv'), show_col_types = FALSE)
stopifnot(all(manifest$cache_reused))
# Changing reference content at the same path must invalidate a fitted object.
ref[['RNA']] <- CreateAssay5Object(counts = counts * 2)
saveRDS(ref, cfg$rctd$reference_rds)
run_fixture(); stopifnot(env$fit_calls == 2L)
cfg$rctd$force_recompute <- TRUE
run_fixture(); stopifnot(env$fit_calls == 3L)
cfg$rctd$force_recompute <- FALSE; cfg$rctd$mode <- 'full'
run_fixture(); stopifnot(env$fit_calls == 4L)
weights <- read_csv(tab_path(cfg, 'rctd_weights.csv.gz'), show_col_types = FALSE)
stopifnot(nrow(weights) == ncol(spatial_fixture), all(is.na(weights$rctd_first_type)))
unlink(test_dir, recursive = TRUE)
cat('All regression checks passed\n')
