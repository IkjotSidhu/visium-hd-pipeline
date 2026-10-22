# ---- 08_annotation.R ----
# Reviewable cluster annotations. No cell-type labels are inferred automatically.
# A mapping template includes a fingerprint so labels cannot silently follow a
# different clustering. Original cluster IDs are retained in every cell table.
source(file.path(dirname(sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value = TRUE))), '..', 'R', 'utils.R'))
cfg <- step_setup('08_annotation')
obj <- read_step_rds(cfg, '03_cluster')
sheet <- read_sample_sheet(cfg)
assert_cluster_resolutions(obj, cfg$clustering$resolutions)
run_annotation <- function(cfg) {
column <- cluster_column(cfg$cluster_resolution)
key <- clustering_key(obj, column)
cfg$clustering_key <- key
clusters <- sort_levels(obj[[column]][, 1])
template <- read_annotation_map(NULL, clusters, key)
write_csv(template, tab_path(cfg, 'annotation_template.csv'))
path <- annotation_mapping(cfg)
cfg$annotation_mapping_md5 <- if (!is.null(path)) unname(tools::md5sum(path)) else NULL
annotation <- read_annotation_map(path, clusters, key) |>
  mutate(annotation_status = if_else(!is.na(broad_type) & nzchar(trimws(broad_type)), 'reviewed', 'unannotated'),
         broad_type = if_else(annotation_status == 'reviewed', trimws(broad_type), 'Unannotated'),
         cell_state = if_else(is.na(cell_state) | !nzchar(trimws(cell_state)), 'Unspecified', trimws(cell_state)))
if (isTRUE(cfg$annotation$require_complete) && any(annotation$annotation_status != 'reviewed'))
  stop('Annotation is incomplete; fill annotation_template.csv and configure the reviewed mapping')
write_csv(annotation, tab_path(cfg, 'cluster_annotations.csv'))

# ---- Marker evidence for manual review ----
marker_path <- file.path(step_dir(cfg, '04_markers'), 'tables', paste0('top', cfg$markers$top_n, '_markers_', column, '.csv'))
if (file.exists(marker_path)) {
  markers <- read_csv(marker_path, col_types = cols(cluster = col_character()), show_col_types = FALSE)
  marker_params <- file.path(step_dir(cfg, '04_markers'), 'params_used.yaml')
  marker_key <- if (file.exists(marker_params)) yaml::read_yaml(marker_params)$clustering_key else NULL
  if (!is.null(marker_key) && !identical(marker_key, key)) stop('Marker evidence belongs to a different clustering; rerun step 04')
  evidence <- markers |> left_join(annotation, by = 'cluster') |>
    mutate(marker_provenance = if (is.null(marker_key)) 'legacy_unverified' else 'matching_clustering')
  write_csv(evidence, tab_path(cfg, 'annotation_marker_evidence.csv'))
} else {
  log_msg(cfg, 'Marker evidence not found; run step 04 for the matching clustering')
}

# ---- Per-unit annotations and replicate composition ----
md <- obj[[]] |> rownames_to_column('cell') |>
  mutate(cluster = as.character(.data[[column]])) |>
  select(cell, cluster, replicate_id, library_id, x_fullres, y_fullres) |>
  left_join(annotation, by = 'cluster')
write_csv(md, tab_path(cfg, 'unit_annotations.csv.gz'))
saveRDS(md, rds_path(cfg, '08_annotation'))
composition <- md |> count(replicate_id, broad_type, name = 'n_units') |>
  complete(replicate_id = unique(md$replicate_id), broad_type = unique(md$broad_type), fill = list(n_units = 0L)) |>
  group_by(replicate_id) |> mutate(total_units = sum(n_units), fraction = n_units / total_units) |> ungroup() |>
  left_join(sheet |> select(-spaceranger_dir), by = 'replicate_id')
write_csv(composition, tab_path(cfg, 'annotation_composition.csv'))
status <- md |> count(cluster, broad_type, annotation_status, name = 'n_units')
write_csv(status, tab_path(cfg, 'annotation_status.csv'))
palette <- make_palette(sort(unique(md$broad_type)))
if ('Unannotated' %in% names(palette)) palette['Unannotated'] <- 'grey70'
p <- ggplot(composition, aes(replicate_id, fraction, fill = broad_type)) + geom_col() +
  scale_fill_manual(values = palette) + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = paste('Annotated unit composition:', cfg$res_label, 'resolution', resolution_label(cfg$cluster_resolution)),
       subtitle = 'Denominator: QC-retained bins/cells; unannotated units retained',
       x = NULL, y = 'Fraction of retained units', fill = 'Broad type')
ggsave(fig_path(cfg, 'annotation_composition.png'), p, width = 10, height = 7, dpi = 100)
ggsave(fig_path(cfg, 'annotation_composition.pdf'), p, width = 10, height = 7)
save_plot_list(plot_spatial(md, 'broad_type', cfg, sheet, palette = palette, title_prefix = paste('Resolution', resolution_label(cfg$cluster_resolution))), cfg, 'spatial_broad_type')
step_finish(cfg)
}
for (r in cfg$clustering$resolutions) run_annotation(resolution_cfg(cfg, r))
step_finish(cfg)
