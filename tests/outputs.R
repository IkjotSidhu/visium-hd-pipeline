# ---- Output and planning regression checks ----
source('R/utils.R')
.libPaths(c(default_config$r_library, .libPaths()))
suppressPackageStartupMessages({ library(Seurat); library(tidyverse) })
expect_error <- function(expr) stopifnot(inherits(tryCatch(force(expr), error = identity), 'error'))

# ---- Annotation mapping identity and explicit evidence ----
path <- tempfile(fileext = '.csv')
template <- read_annotation_map(NULL, c('0', '1'), 'key-A')
write.csv(template, path, row.names = FALSE)
stopifnot(nrow(read_annotation_map(path, c('0', '1'), 'key-A')) == 2)
expect_error(read_annotation_map(path, c('0', '1'), 'key-B'))
template$broad_type[1] <- 'Example type'
write.csv(template, path, row.names = FALSE)
expect_error(read_annotation_map(path, c('0', '1'), 'key-A'))
template$confidence[1] <- 'low'; template$evidence[1] <- 'Example marker evidence'
write.csv(template, path, row.names = FALSE)
stopifnot(read_annotation_map(path, c('0', '1'), 'key-A')$broad_type[1] == 'Example type')
write.csv(rbind(template, template[1, ]), path, row.names = FALSE)
expect_error(read_annotation_map(path, c('0', '1'), 'key-A'))
unlink(path)

# ---- GMT/rank semantics: retain negative and nonsignificant genes ----
path <- tempfile(fileext = '.gmt')
writeLines(c('Path1\tdescription\tG1\tG2\tG2', 'Path2\tdescription\tG2\tG3'), path)
sets <- read_gmt(path)
stopifnot(identical(sets$Path1, c('G1', 'G2')))
ranks <- pathway_ranks(tibble(gene = c('G1', 'G2', 'G3', 'G4'), stat = c(-3, 2, NA, 0), padj = c(0.5, 0.01, NA, 0.9)))
stopifnot(identical(names(ranks), c('G2', 'G4', 'G1')), ranks[['G1']] == -3)
expect_error(pathway_ranks(tibble(gene = c('G1', 'G1'), stat = c(1, -1))))
writeLines('malformed', path); expect_error(read_gmt(path)); unlink(path)
stopifnot(!grepl('<script>', html_escape('<script>"&'), fixed = TRUE))
stopifnot(grepl('&lt;script&gt;', html_table(data.frame(label = '<script>')), fixed = TRUE))

# ---- Report scheduled once, after all per-resolution work ----
cfg <- default_config; cfg$output_dir <- tempdir()
plan <- plan_jobs(cfg)
stopifnot(sum(grepl('^10_report ', plan)) == 1, startsWith(tail(plan, 1), '10_report all '),
          sum(grepl('^08_annotation ', plan)) == 3, sum(grepl('^09_pathways ', plan)) == 3)
cfg$de$run <- FALSE
stopifnot(!any(grepl('^09_pathways ', plan_jobs(cfg))))
cfg$report$run <- FALSE
stopifnot(!any(grepl('^10_report ', plan_jobs(cfg))))
cat('Output helper checks passed\n')

# ---- Actual annotation step with a partial reviewed mapping ----
counts <- matrix(1:8, nrow = 2, dimnames = list(c('G1', 'G2'), paste0('C', 1:4)))
fixture <- CreateSeuratObject(counts)
fixture$cluster_res_1 <- factor(c('0', '0', '1', '1'))
fixture$replicate_id <- c('R1', 'R1', 'R2', 'R2'); fixture$library_id <- 'L1'
fixture$x_fullres <- 1:4; fixture$y_fullres <- 1:4
root <- tempfile('vhd-annotation-'); dir.create(root)
cfg <- default_config; cfg$output_dir <- root; cfg$res_label <- '008um'; cfg$bin_size <- '8'
cfg$clustering$resolutions <- 1
cfg$step <- '08_annotation'; cfg$step_dir <- step_dir(cfg, cfg$step); cfg$start_time <- Sys.time()
dir.create(file.path(cfg$step_dir, 'tables'), recursive = TRUE)
dir.create(file.path(cfg$step_dir, 'figures')); dir.create(rds_dir(cfg), recursive = TRUE)
key <- clustering_key(fixture, 'cluster_res_1')
mapping <- read_annotation_map(NULL, c('0', '1'), key)
mapping$broad_type[1] <- 'Example type'; mapping$confidence[1] <- 'low'; mapping$evidence[1] <- 'Test evidence'
path <- file.path(root, 'mapping.csv'); write.csv(mapping, path, row.names = FALSE)
cfg$annotation$mapping_files <- list('8' = path)
env <- new.env(parent = globalenv()); env$step_setup <- function(...) cfg
env$read_step_rds <- function(...) fixture
env$read_sample_sheet <- function(...) tibble(replicate_id = c('R1', 'R2'), library_id = 'L1', spaceranger_dir = root)
env$plot_spatial <- function(...) list(); env$save_plot_list <- function(...) invisible(NULL)
for (expr in as.list(parse('steps/08_annotation.R'))[-1]) eval(expr, env)
current <- resolution_cfg(cfg, 1, create = FALSE)
units <- read_csv(tab_path(current, 'unit_annotations.csv.gz'), show_col_types = FALSE)
composition <- read_csv(tab_path(current, 'annotation_composition.csv'), show_col_types = FALSE)
stopifnot(sum(units$annotation_status == 'reviewed') == 2, sum(units$broad_type == 'Unannotated') == 2,
          all(tapply(composition$fraction, composition$replicate_id, sum) == 1))
cfg$annotation$require_complete <- TRUE
expect_error(for (expr in as.list(parse('steps/08_annotation.R'))[-1]) eval(expr, env))
unlink(root, recursive = TRUE)

# ---- Report with missing/disabled steps and hostile metadata ----
root <- tempfile('vhd-report-'); dir.create(root)
cfg <- default_config; cfg$output_dir <- root; cfg$res_label <- '008um'; cfg$bin_size <- '8'; cfg$bin_sizes <- list(8)
cfg$step <- '10_report'; cfg$step_dir <- step_dir(cfg, cfg$step); cfg$start_time <- Sys.time()
dir.create(file.path(cfg$step_dir, 'tables'), recursive = TRUE); dir.create(file.path(cfg$step_dir, 'figures'))
env <- new.env(parent = globalenv()); env$step_setup <- function(...) cfg
env$read_sample_sheet <- function(...) tibble(replicate_id = '<script>alert(1)</script>', spaceranger_dir = root)
for (expr in as.list(parse('steps/10_report.R'))[-1]) eval(expr, env)
html <- paste(readLines(file.path(root, 'index.html')), collapse = '\n')
stopifnot(grepl('missing/incomplete', html, fixed = TRUE), grepl('disabled in current config', html, fixed = TRUE),
          !grepl('<script>', html, fixed = TRUE), grepl('&lt;script&gt;', html, fixed = TRUE))
unlink(root, recursive = TRUE)
cat('Annotation/report integration fixtures passed\n')

# ---- Stable folder names, distinct bin sizes and resolution-specific labels ----
stopifnot(identical(default_config$clustering$resolutions, (1:20)/10),
          identical(vapply(c(0.1, 1, 2), resolution_label, character(1)), c('0.1','1.0','2.0')))
cfg <- default_config; cfg$output_dir <- tempdir(); cfg$res_label <- '008um'; cfg$bin_size <- '8'; cfg$step <- '08_annotation'
cfg$annotation$mapping_files <- list('8' = list('0.1' = 'low.csv', '1.0' = 'primary.csv'))
stopifnot(annotation_mapping(resolution_cfg(cfg, 0.1, FALSE)) == 'low.csv',
          annotation_mapping(resolution_cfg(cfg, 1, FALSE)) == 'primary.csv')
cfg$annotation$mapping_files <- list('8' = 'legacy.csv')
stopifnot(is.null(annotation_mapping(resolution_cfg(cfg, 0.1, FALSE))),
          annotation_mapping(resolution_cfg(cfg, 1, FALSE)) == 'legacy.csv')
stopifnot(grepl('/008um/08_annotation/res_0.1$', resolution_cfg(cfg, 0.1, FALSE)$step_dir),
          identical(rds_path(resolution_cfg(cfg, 0.1, FALSE), '03_cluster'), rds_path(cfg, '03_cluster')),
          rds_path(resolution_cfg(cfg, 0.1, FALSE), '08_annotation') != rds_path(resolution_cfg(cfg, 1, FALSE), '08_annotation'))
expect_error(assert_cluster_resolutions(fixture, c(0.1, 1)))
cat('Resolution path and mapping checks passed\n')
