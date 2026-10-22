# ---- 04_markers.R ----
# Cluster markers at clustering$main_resolution, computed on the sketch (a
# balanced subset - fast, and avoids letting the largest samples dominate).
# Writes the full marker table, top-N per cluster, a dot plot and a z-scored
# average-expression heatmap.
#
# Usage: Rscript steps/04_markers.R config.yaml [8|16|cells]

# ---- Setup ----
source(file.path(dirname(sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value = TRUE))), '..', 'R', 'utils.R'))
cfg      <- step_setup('04_markers')
mk       <- cfg$markers
obj      <- read_step_rds(cfg, '03_cluster')
main_col <- paste0('cluster_res_', cfg$clustering$main_resolution)

# ---- Sketch subset with joined layers ----
# Seurat v5: markers need a single counts/data layer
sk <- subset(obj, cells = colnames(obj[['sketch']]))
DefaultAssay(sk) <- 'sketch'
sk[['sketch']] <- JoinLayers(sk[['sketch']])
Idents(sk) <- main_col
rm(obj); gc()
log_msg(cfg, 'Markers on ', ncol(sk), ' sketched units, ', nlevels(Idents(sk)), ' clusters (', main_col, ')')

# ---- Find markers (Wilcoxon; Seurat uses presto when installed) ----
markers <- FindAllMarkers(sk, only.pos = TRUE, min.pct = mk$min_pct, logfc.threshold = mk$logfc_threshold, verbose = FALSE) |>
  as_tibble() |>
  arrange(cluster, p_val_adj, desc(avg_log2FC))
write_csv(markers, tab_path(cfg, paste0('markers_', main_col, '.csv')))

top <- markers |> filter(p_val_adj < 0.05) |> group_by(cluster) |> slice_max(avg_log2FC, n = mk$top_n, with_ties = FALSE) |> ungroup()
write_csv(top, tab_path(cfg, paste0('top', mk$top_n, '_markers_', main_col, '.csv')))
log_msg(cfg, nrow(markers), ' markers; ', n_distinct(top$cluster), ' clusters with significant markers')

# ---- Dot plot: top 3 per cluster ----
dot_genes <- top |> group_by(cluster) |> slice_head(n = 3) |> pull(gene) |> unique()
p_dot <- DotPlot(sk, features = dot_genes, cluster.idents = FALSE) +
  scale_colour_gradientn(colours = RColorBrewer::brewer.pal(9, 'Blues')) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)) +
  labs(x = NULL, y = 'Cluster', title = paste('Top 3 markers per cluster -', main_col))
w <- max(8, length(dot_genes) * 0.18)
ggsave(fig_path(cfg, paste0('dotplot_top3_', main_col, '.pdf')), p_dot, width = w, height = 7, limitsize = FALSE)
ggsave(fig_path(cfg, paste0('dotplot_top3_', main_col, '.png')), p_dot, width = w, height = 7, dpi = 100, limitsize = FALSE)

# ---- Heatmap: z-scored mean expression of top-N markers ----
sk_assay <- DefaultAssay(sk)                                    # resolve at point of use
avg <- as.matrix(AverageExpression(sk, assays = sk_assay, features = unique(top$gene), group.by = main_col)[[sk_assay]])
avg <- avg[apply(avg, 1, var) > 0, , drop = FALSE]              # zero-variance rows would give NaN z-scores
z   <- t(scale(t(avg)))
z[!is.finite(z)] <- 0

# Genes in cluster order (each gene placed at the cluster it marks first)
gene_order <- top |> filter(gene %in% rownames(z)) |> distinct(gene, .keep_all = TRUE) |> pull(gene)
heat_df <- as_tibble(z, rownames = 'gene') |>
  pivot_longer(-gene, names_to = 'cluster', values_to = 'z') |>
  mutate(gene = factor(gene, levels = rev(gene_order)),
         cluster = factor(sub('^g', '', cluster), levels = sort_levels(sub('^g', '', cluster))))
p_heat <- ggplot(heat_df, aes(cluster, gene, fill = z)) +
  geom_tile() +
  # Diverging z-scores: blue-grey-red with a neutral midpoint at 0
  scale_fill_gradient2(low = '#2166AC', mid = 'grey92', high = '#B2182B', midpoint = 0, limits = c(-2.5, 2.5),
                       oob = scales::squish, name = 'z-score') +
  labs(x = 'Cluster', y = NULL, title = paste0('Top ', mk$top_n, ' markers - ', main_col)) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 5), panel.grid = element_blank())
h <- max(6, length(gene_order) * 0.07)
ggsave(fig_path(cfg, paste0('heatmap_top', mk$top_n, '_', main_col, '.pdf')), p_heat, width = 9, height = h, limitsize = FALSE)
ggsave(fig_path(cfg, paste0('heatmap_top', mk$top_n, '_', main_col, '.png')), p_heat, width = 9, height = h, dpi = 100, limitsize = FALSE)

step_finish(cfg)
