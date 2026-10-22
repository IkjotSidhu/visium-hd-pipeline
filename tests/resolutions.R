# ---- Actual marker sweep, including a one-cluster resolution ----
source('R/utils.R')
.libPaths(c(default_config$r_library, .libPaths()))
suppressPackageStartupMessages({ library(Seurat); library(tidyverse) })
set.seed(17)
counts <- matrix(rpois(80 * 120, 1), 80, dimnames = list(paste0('G', 1:80), paste0('C', 1:120)))
counts[1:10, 1:60] <- counts[1:10, 1:60] + 10
counts[11:20, 61:120] <- counts[11:20, 61:120] + 10
fixture <- CreateSeuratObject(counts)
fixture[['sketch']] <- CreateAssay5Object(counts = counts)
fixture <- NormalizeData(fixture, assay = 'sketch', verbose = FALSE)
fixture$cluster_res_0.1 <- factor(rep('0', 120))
fixture$cluster_res_1 <- factor(rep(c('0','1'), each = 60))
root <- tempfile('vhd-marker-sweep-'); dir.create(root)
cfg <- default_config; cfg$output_dir <- root; cfg$res_label <- '008um'; cfg$bin_size <- '8'
cfg$clustering$resolutions <- c(0.1, 1); cfg$step <- '04_markers'
cfg$step_dir <- step_dir(cfg, cfg$step); cfg$start_time <- Sys.time()
dir.create(file.path(cfg$step_dir, 'tables'), recursive = TRUE); dir.create(file.path(cfg$step_dir, 'figures'))
env <- new.env(parent = globalenv()); env$step_setup <- function(...) cfg; env$read_step_rds <- function(...) fixture
for (expr in as.list(parse('steps/04_markers.R'))[-1]) eval(expr, env)
for (r in c(0.1, 1)) {
  current <- resolution_cfg(cfg, r, FALSE)
  stopifnot(length(list.files(file.path(current$step_dir, 'figures'), pattern = '\\.png$')) == 2,
            length(list.files(file.path(current$step_dir, 'figures'), pattern = '\\.pdf$')) == 2,
            file.exists(file.path(current$step_dir, 'params_used.yaml')))
  for (path in list.files(file.path(current$step_dir, 'figures'), pattern = '\\.png$', full.names = TRUE)) {
    img <- png::readPNG(path)
    stopifnot(all(img[1, 1, 1:3] == 1), dim(img)[3] == 3 || all(img[, , 4] == 1))
  }
}
low <- read_csv(tab_path(resolution_cfg(cfg, 0.1, FALSE), 'marker_status.csv'), show_col_types = FALSE)
high <- read_csv(tab_path(resolution_cfg(cfg, 1, FALSE), 'marker_status.csv'), show_col_types = FALSE)
stopifnot(low$status == 'no_markers', high$status == 'complete')
unlink(root, recursive = TRUE)
cat('Multi-resolution marker plots passed\n')
