# ---- 03_cluster.R ----
# Normalise once, sketch a balanced subset per sample, cluster the sketch
# (optionally after Harmony integration), and project PCA/UMAP/clusters back to
# every bin/cell. Writes UMAPs, a clustree, cluster composition per replicate and
# hires spatial cluster maps.
#
# Usage: Rscript steps/03_cluster.R config.yaml [8|16|cells]

# ---- Setup ----
source(file.path(dirname(sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value = TRUE))), '..', 'R', 'utils.R'))
cfg   <- step_setup('03_cluster')
cl    <- cfg$clustering
sheet <- read_sample_sheet(cfg)
obj   <- read_step_rds(cfg, '02_load_qc')
assay <- DefaultAssay(obj)
dims  <- seq_len(cl$pca_dims)
integrate_by <- resolve_integration(cfg, sheet)                 # NULL = no integration
use_harmony  <- !is.null(integrate_by)
log_msg(cfg, 'Integration: ', if (use_harmony) paste('Harmony by', integrate_by) else 'none')

# ---- Layers ----
# Sketching samples ncells from EACH layer, and Harmony integrates across layers:
# split by the integration variable, otherwise by replicate (balanced sketch)
split_col <- if (use_harmony) integrate_by else 'replicate_id'
obj[[assay]] <- JoinLayers(obj[[assay]])
obj[[assay]] <- split(obj[[assay]], f = obj[[split_col]][, 1])
log_msg(cfg, 'Layers split by ', split_col, ': ', paste(unique(obj[[split_col]][, 1]), collapse = ', '))

# ---- Normalise + variable features (once, on the full object) ----
obj <- NormalizeData(obj, verbose = FALSE)
obj <- FindVariableFeatures(obj, nfeatures = cl$n_variable_features, verbose = FALSE)

# ---- Sketch ----
# Leverage-score sketch keeps rare populations better than uniform sampling;
# layers smaller than sketch_cells are kept whole
obj <- SketchData(obj, ncells = cl$sketch_cells, method = 'LeverageScore', sketched.assay = 'sketch',
                  seed = cfg$seed, verbose = FALSE)
DefaultAssay(obj) <- 'sketch'
log_msg(cfg, 'Sketched ', ncol(obj[['sketch']]), ' of ', ncol(obj), ' units')

# ---- PCA (+ Harmony) on the sketch ----
obj <- FindVariableFeatures(obj, nfeatures = cl$n_variable_features, verbose = FALSE)
obj <- ScaleData(obj, verbose = FALSE)
obj <- RunPCA(obj, npcs = max(dims), reduction.name = 'pca', verbose = FALSE)
if (use_harmony) {
  obj <- IntegrateLayers(obj, method = HarmonyIntegration, orig.reduction = 'pca',
                         new.reduction = 'harmony', verbose = FALSE)
}
red <- if (use_harmony) 'harmony' else 'pca'

# ---- Cluster + UMAP on the sketch ----
obj <- FindNeighbors(obj, reduction = red, dims = dims, verbose = FALSE)
obj <- FindClusters(obj, resolution = cl$resolutions, random.seed = cfg$seed, verbose = FALSE)
obj <- RunUMAP(obj, reduction = red, dims = dims, return.model = TRUE, reduction.name = 'umap.sketch',
               seed.use = cfg$seed, verbose = FALSE)

# ---- Project to every bin/cell ----
res_cols <- paste0('sketch_snn_res.', cl$resolutions)
if (use_harmony) {
  obj <- ProjectIntegration(obj, sketched.assay = 'sketch', assay = assay, reduction = 'harmony', verbose = FALSE)
  full_red <- 'harmony.full'
  obj <- ProjectData(obj, assay = assay, sketched.assay = 'sketch', sketched.reduction = full_red,
                     full.reduction = full_red, dims = dims, umap.model = 'umap.sketch',
                     refdata = setNames(as.list(res_cols), paste0('cluster_res_', cl$resolutions)), verbose = FALSE)
} else {
  obj <- ProjectData(obj, assay = assay, sketched.assay = 'sketch', sketched.reduction = 'pca',
                     full.reduction = 'pca.full', dims = dims, umap.model = 'umap.sketch',
                     refdata = setNames(as.list(res_cols), paste0('cluster_res_', cl$resolutions)), verbose = FALSE)
}
DefaultAssay(obj) <- assay
obj[[assay]] <- JoinLayers(obj[[assay]])      # downstream steps read counts from one layer

# Projected UMAP name differs by Seurat version - find it
umap_full <- setdiff(grep('umap', Reductions(obj), value = TRUE), 'umap.sketch')[1]
log_msg(cfg, 'Projected reductions: ', paste(Reductions(obj), collapse = ', '))

# Clusters as factors in numeric order; main resolution as the identity
for (r in cl$resolutions) {
  col <- paste0('cluster_res_', r)
  obj[[col]] <- factor(obj[[col]][, 1], levels = sort_levels(obj[[col]][, 1]))
}
main_col <- paste0('cluster_res_', cl$main_resolution)
Idents(obj) <- main_col

# Save before plotting, so a plotting error never costs the clustering
saveRDS(obj, rds_path(cfg, '03_cluster'))

# ---- Clustree (how sketch clusters split across resolutions) ----
sketch_md <- obj[[]][colnames(obj[['sketch']]), ]
# Edge legends hidden: ggraph's edge_colourbar guide is unknown to ggplot2 >= 4.0
p_tree <- clustree::clustree(sketch_md, prefix = 'sketch_snn_res.') +
  ggraph::scale_edge_colour_gradientn(colours = c('grey80', 'black'), guide = 'none') +
  ggraph::scale_edge_alpha(guide = 'none')
ggsave(fig_path(cfg, 'clustree.pdf'), p_tree, width = 10, height = max(10, 0.65 * length(cl$resolutions)))
ggsave(fig_path(cfg, 'clustree.png'), p_tree, width = 10, height = max(10, 0.65 * length(cl$resolutions)), dpi = 100)

# ---- UMAPs (all bins/cells, rasterised) ----
md <- obj[[]] |> mutate(umap_1 = Embeddings(obj, umap_full)[, 1], umap_2 = Embeddings(obj, umap_full)[, 2])
umap_plot <- function(col, label = FALSE) {
  pal <- make_palette(levels(factor(md[[col]], levels = sort_levels(md[[col]]))))
  p <- ggplot(md, aes(umap_1, umap_2, colour = .data[[col]])) +
    ggrastr::rasterise(geom_point(size = 0.1, stroke = 0), dpi = 300) +
    scale_colour_manual(values = pal) +
    guides(colour = guide_legend(override.aes = list(size = 3), ncol = if (length(pal) > 20) 2 else 1)) +
    labs(title = col, x = 'UMAP 1', y = 'UMAP 2') +
    theme_classic() + coord_fixed()
  if (label) {
    centers <- md |> group_by(.data[[col]]) |> summarise(umap_1 = median(umap_1), umap_2 = median(umap_2))
    p <- p + geom_text(data = centers, aes(label = .data[[col]]), colour = 'black', size = 3)
  }
  p
}
# Metadata views and clustree are shared; every clustering gets its own gallery.
meta_cols <- intersect(c('replicate_id', 'donor', cfg$de$condition_col, 'library_id'), colnames(md))
save_plot_list(map(set_names(meta_cols), umap_plot), cfg, 'umap_metadata', width = 9, height = 8)
summary_rows <- list()
for (r in cl$resolutions) {
  current <- resolution_cfg(cfg, r)
  column <- cluster_column(r)
  log_msg(cfg, 'Writing plots for clustering resolution ', resolution_label(r))
  save_plot_list(list(clusters = umap_plot(column, label = TRUE)), current, 'umap', width = 9, height = 8)
  write_csv(md |> rownames_to_column('cell') |> select(cell, replicate_id, library_id, all_of(column)),
            tab_path(current, 'cluster_labels.csv.gz'))

  comp <- md |> count(cluster = .data[[column]], replicate_id) |>
    group_by(cluster) |> mutate(frac = n / sum(n)) |> ungroup()
  write_csv(md |> count(replicate_id, cluster = .data[[column]]) |>
    group_by(replicate_id) |> mutate(pct_of_replicate = round(100 * n / sum(n), 2)) |> ungroup(),
    tab_path(current, 'cluster_composition.csv'))
  p_comp <- ggplot(comp, aes(cluster, frac, fill = replicate_id)) + geom_col() +
    scale_fill_manual(values = make_palette(sheet$replicate_id)) +
    labs(title = paste('Replicate make-up of clusters - resolution', resolution_label(r)),
         x = 'Cluster', y = 'Fraction of units') + theme_classic()
  ggsave(fig_path(current, 'cluster_by_replicate.pdf'), p_comp, width = 12, height = 5)
  ggsave(fig_path(current, 'cluster_by_replicate.png'), p_comp, width = 12, height = 5, dpi = 100)
  n_single <- NA_integer_
  if ('donor' %in% colnames(md)) {
    donor_mix <- md |> count(cluster = .data[[column]], donor) |>
      group_by(cluster) |> summarise(max_donor_frac = round(max(n) / sum(n), 3), n = sum(n)) |>
      arrange(desc(max_donor_frac))
    write_csv(donor_mix, tab_path(current, 'donor_mixing_per_cluster.csv'))
    n_single <- sum(donor_mix$max_donor_frac > 0.95)
  }
  md_sp <- md |> mutate(across(all_of(column), as.character))
  save_plot_list(plot_spatial(md_sp, column, current, sheet,
      palette = make_palette(levels(obj[[column]][, 1])), title_prefix = paste('Resolution', resolution_label(r)), width = 10, height = 9),
    current, 'spatial_clusters', width = 10, height = 9)
  summary_rows[[resolution_label(r)]] <- tibble(clustering_resolution = r,
    n_clusters = nlevels(obj[[column]][, 1]), n_units = nrow(md), clusters_gt95pct_one_donor = n_single,
    folder = paste0('res_', resolution_label(r)))
  step_finish(current)
}
write_csv(bind_rows(summary_rows), tab_path(cfg, 'clustering_resolution_summary.csv'))
step_finish(cfg)
