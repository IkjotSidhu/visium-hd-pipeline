# ---- 07_rctd.R ----
# RCTD deconvolution (spacexr) of every replicate against a single-cell
# reference. Per-replicate results are cached in rds/rctd/ and reused, so a
# rerun only computes what's missing. Writes per-bin weights, the dominant
# type, composition per replicate and hires dominant-type maps.
#
# Reference: a Seurat object (rctd$reference_rds) with raw counts and a
# cell-type column (rctd$celltype_col). Use a tissue-matched reference: types
# absent from the tissue can be assigned spuriously (see FAQ.md).
#
# Usage: Rscript steps/07_rctd.R config.yaml [8|16|cells]

# ---- Setup ----
source(file.path(dirname(sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value = TRUE))), '..', 'R', 'utils.R'))
cfg   <- step_setup('07_rctd')
rc    <- cfg$rctd
sheet <- read_sample_sheet(cfg)
suppressPackageStartupMessages(library(spacexr))
cache_dir <- file.path(rds_dir(cfg), 'rctd')
dir.create(cache_dir, showWarnings = FALSE)

# ---- Reference ----
# Labels sanitised: RCTD rejects "/" and stray whitespace/Unicode splits types
sanitize_label <- function(x) {
  x |> iconv(to = 'ASCII//TRANSLIT') |> str_replace_all('/', '-') |> str_squish() |> str_replace_all(' ', '_')
}
ref_obj <- readRDS(rc$reference_rds)
ref_assay <- rc$reference_assay
ref_input <- rctd_reference_inputs(ref_obj, ref_assay, rc$celltype_col)
ref_counts <- ref_input$counts
ref_types <- sanitize_label(ref_input$types)
rm(ref_input)
names(ref_types) <- colnames(ref_counts)

# Drop rare and unlabelled types
type_n <- table(ref_types)
drop   <- names(type_n)[type_n < rc$min_cells_type | names(type_n) %in% c('', 'NA', 'Unknown')]
if (length(drop)) log_msg(cfg, 'Dropping reference types: ', paste(drop, collapse = ', '))
keep <- !ref_types %in% drop & !is.na(ref_types)
if (nlevels(factor(ref_types[keep])) < 2) stop('RCTD needs at least two retained reference types')
reference_key <- digest::digest(list(counts = ref_counts[, keep], types = ref_types[keep],
                                    assay = ref_assay, min_cells_type = rc$min_cells_type), algo = 'sha256')
reference <- Reference(ref_counts[, keep], factor(ref_types[keep]), colSums(ref_counts[, keep]))
write_csv(tibble(cell = colnames(ref_counts), type = ref_types, retained = keep),
          tab_path(cfg, 'rctd_reference_cells.csv.gz'))
rm(ref_obj, ref_counts); gc()
log_msg(cfg, 'Reference: ', sum(keep), ' cells, ', nlevels(factor(ref_types[keep])), ' types')

# ---- Spatial data ----
obj   <- read_step_rds(cfg, '03_cluster')
assay <- DefaultAssay(obj)
counts_all <- LayerData(obj, assay = assay, layer = 'counts')

# ---- RCTD per replicate (cached) ----
weights_list <- list(); first_type <- list(); cache_rows <- list()
cache_settings <- list(mode = rc$mode, umi_min = rc$umi_min, seed = cfg$seed,
                       max_cores = rc$max_cores,
                       spacexr = as.character(packageVersion('spacexr')),
                       Seurat = as.character(packageVersion('Seurat')))
for (rep_id in unique(obj$replicate_id)) {
  # Hash the actual effective inputs, not file names or modification times.
  cells <- colnames(obj)[obj$replicate_id == rep_id]
  coords <- data.frame(x = obj$x_fullres[cells], y = obj$y_fullres[cells], row.names = cells)
  if (any(!is.finite(as.matrix(coords)))) stop(rep_id, ': non-finite spatial coordinates')
  spatial_counts <- counts_all[, cells, drop = FALSE]
  validate_raw_counts(spatial_counts, rep_id)
  cache <- file.path(cache_dir, paste0(rep_id, '.rds'))
  key <- rctd_cache_key(reference_key, spatial_counts, coords, cache_settings)
  rctd <- read_rctd_cache(cache, key, rc$force_recompute)
  reused <- !is.null(rctd)
  if (reused) {
    log_msg(cfg, rep_id, ': using validated cached RCTD result')
  } else {
    log_msg(cfg, rep_id, ': running RCTD (', rc$mode, ' mode; missing/stale/forced cache)')
    # Per-replicate seed makes cache misses independent of earlier cache hits.
    set.seed(cfg$seed)
    puck <- SpatialRNA(coords, spatial_counts, colSums(spatial_counts))
    rctd <- create.RCTD(puck, reference, max_cores = rc$max_cores, UMI_min = rc$umi_min)
    rctd <- run.RCTD(rctd, doublet_mode = rc$mode)
    write_rctd_cache(rctd, cache, key)
  }
  cache_rows[[rep_id]] <- tibble(replicate_id = rep_id, cache_reused = reused, cache_key = key,
                                reference_key = reference_key, n_input_units = length(cells))
  # Per-bin type weights (rows sum to 1) and RCTD's first-choice type
  w <- as.matrix(normalize_weights(rctd@results$weights))
  weights_list[[rep_id]] <- w
  if (!is.null(rctd@results$results_df)) {
    first_type[[rep_id]] <- tibble(cell = rownames(rctd@results$results_df),
                                   rctd_first_type = as.character(rctd@results$results_df$first_type),
                                   rctd_spot_class = as.character(rctd@results$results_df$spot_class))
  }
}

# ---- Combine ----
types <- sort(unique(unlist(map(weights_list, colnames))))
weights <- map_dfr(weights_list, \(w) as_tibble(w, rownames = 'cell')) |>
  mutate(across(all_of(types), \(x) replace_na(x, 0)))
dominant <- weights |>
  pivot_longer(all_of(types), names_to = 'type', values_to = 'w') |>
  slice_max(w, n = 1, with_ties = FALSE, by = cell) |>
  transmute(cell, rctd_dominant_type = type, rctd_dominant_weight = w)
first_template <- tibble(cell = character(), rctd_first_type = character(), rctd_spot_class = character())
rctd_tbl <- weights |> left_join(dominant, by = 'cell') |>
  left_join(bind_rows(c(list(first_template), first_type)), by = 'cell')
write_csv(bind_rows(cache_rows), tab_path(cfg, 'rctd_cache_manifest.csv'))
write_csv(rctd_tbl, tab_path(cfg, 'rctd_weights.csv.gz'))
saveRDS(rctd_tbl, rds_path(cfg, '07_rctd'))
log_msg(cfg, nrow(rctd_tbl), ' of ', ncol(obj), ' units deconvolved (others below UMI_min = ', rc$umi_min, ')')

# ---- Composition per replicate (mean weight) ----
md <- obj[[]] |> rownames_to_column('cell') |> select(cell, replicate_id, library_id, x_fullres, y_fullres)
type_pal <- make_palette(types)
comp <- weights |> left_join(md, by = 'cell') |>
  pivot_longer(all_of(types), names_to = 'type', values_to = 'w') |>
  group_by(replicate_id, type) |> summarise(mean_weight = mean(w), .groups = 'drop') |>
  mutate(type = factor(type, levels = types))
write_csv(comp, tab_path(cfg, 'rctd_composition_per_replicate.csv'))
p_comp <- ggplot(comp, aes(replicate_id, mean_weight, fill = type)) +
  geom_col() +
  scale_fill_manual(values = type_pal) +
  labs(x = NULL, y = 'Mean RCTD weight', fill = 'Cell type', title = 'RCTD composition per replicate') +
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(fig_path(cfg, 'rctd_composition.pdf'), p_comp, width = 9, height = 6)
ggsave(fig_path(cfg, 'rctd_composition.png'), p_comp, width = 9, height = 6, dpi = 100)

# ---- Dominant type maps (hires) ----
md_dom <- md |> inner_join(dominant, by = 'cell')
save_spatial(md_dom, 'rctd_dominant_type', cfg, sheet, 'spatial_rctd_dominant_type',
             palette = type_pal[sort(unique(md_dom$rctd_dominant_type))], width = 10, height = 9)

step_finish(cfg)
