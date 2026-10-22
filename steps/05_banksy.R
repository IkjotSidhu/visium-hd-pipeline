# ---- 05_banksy.R ----
# BANKSY spatial clustering: each bin/cell's expression is augmented with its
# spatial neighbourhood, so clusters follow tissue domains rather than only
# cell identity. lambda (0-1) sets the neighbourhood weight (~0.2 for cell-type
# like clusters, ~0.8 for domains); k_geom the neighbourhood size.
# Saves cluster labels (not the large BANKSY assay) + hires spatial maps.
#
# Usage: Rscript steps/05_banksy.R config.yaml [8|16|cells]

# ---- Setup ----
source(file.path(dirname(sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value = TRUE))), '..', 'R', 'utils.R'))
cfg   <- step_setup('05_banksy')
bk    <- cfg$banksy
sheet <- read_sample_sheet(cfg)
suppressPackageStartupMessages({ library(Banksy); library(SeuratWrappers) })

obj   <- read_step_rds(cfg, '03_cluster')
assay <- DefaultAssay(obj)
dims  <- seq_len(bk$pca_dims)

# Drop the sketch assay and projected reductions (not needed, slow RunBanksy down);
# keep the variable features found on the full data in step 03
vf  <- VariableFeatures(obj[[assay]])
obj <- DietSeurat(obj, assays = assay, layers = c('counts', 'data'))
VariableFeatures(obj) <- vf

# ---- Run each lambda x k_geom combination ----
labels <- tibble(cell = colnames(obj))
for (lambda in unlist(bk$lambda)) for (k in unlist(bk$k_geom)) {
  tag <- paste0('banksy_l', lambda, '_k', k)
  log_msg(cfg, 'BANKSY lambda = ', lambda, ', k_geom = ', k)

  # Neighbourhoods are built within each replicate (group), on full-res pixel coordinates
  b <- RunBanksy(obj, lambda = lambda, assay = assay, slot = 'data', features = 'variable',
                 group = 'replicate_id', dimx = 'x_fullres', dimy = 'y_fullres', k_geom = k,
                 split.scale = TRUE, assay_name = 'BANKSY', verbose = FALSE)
  b <- RunPCA(b, assay = 'BANKSY', reduction.name = 'pca.banksy', features = rownames(b[['BANKSY']]),
              npcs = max(dims), verbose = FALSE)
  b <- FindNeighbors(b, reduction = 'pca.banksy', dims = dims, verbose = FALSE)
  for (r in unlist(bk$resolutions)) {
    b <- FindClusters(b, resolution = r, cluster.name = paste0(tag, '_res_', r), random.seed = cfg$seed, verbose = FALSE)
  }
  labels <- labels |> left_join(b[[]] |> select(starts_with(tag)) |> rownames_to_column('cell'), by = 'cell')
  rm(b); gc()
}

# ---- Save labels ----
write_csv(labels, tab_path(cfg, 'banksy_clusters.csv.gz'))
saveRDS(labels, rds_path(cfg, '05_banksy'))

# ---- Spatial maps ----
md <- obj[[]] |> rownames_to_column('cell') |> select(cell, library_id, replicate_id, x_fullres, y_fullres) |>
  left_join(labels, by = 'cell')
for (col in setdiff(colnames(labels), 'cell')) {
  md[[col]] <- as.character(md[[col]])
  save_plot_list(plot_spatial(md, col, cfg, sheet, palette = make_palette(sort_levels(md[[col]]))),
                 cfg, paste0('spatial_', col), width = 10, height = 9)
}

step_finish(cfg)
