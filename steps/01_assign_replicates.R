# ---- 01_assign_replicates.R ----
# Split capture areas that hold several replicate tissue pieces into one
# replicate per piece. Runs once (resolution independent): pieces are found on
# 8 um bins and the labels are propagated to 16 um bins and segmented cells via
# Space Ranger's barcode_mappings.parquet.
#
# replicates$layout:
#   left_to_right / top_to_bottom - automatic: tissue pieces are ordered along the
#     image axis and matched to the sample-sheet rows of that library in order
#   manual - <manual_dir>/<library_id>.csv with columns barcode (8 um) + replicate_id
#     (e.g. a Loupe Browser export)
#   none - nothing to do (one replicate per capture area)
#
# Usage: Rscript steps/01_assign_replicates.R config.yaml

# ---- Setup ----
source(file.path(dirname(sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value = TRUE))), '..', 'R', 'utils.R'))
cfg <- step_setup('01_assign_replicates')
rp  <- cfg$replicates
if (!layout_splits(cfg)) {
  log_msg(cfg, 'replicates$layout is none - nothing to assign')
  step_finish(cfg)
  quit(save = 'no')
}
sheet <- read_sample_sheet(cfg)

# ---- Automatic piece detection on 8 um bins ----
# Coarsen tissue onto blocks, keep densely filled blocks (Space Ranger also calls
# scattered specks of "tissue" on bare glass, which would chain pieces together),
# link nearby dense blocks into pieces, keep the N largest, order them along the
# layout axis, then attach sparse edge blocks next to each piece.
detect_pieces <- function(pos, n_reps, lib) {
  bins_per_block <- (rp$block_um / 8)^2
  blocks <- pos |> count(bx, by, name = 'n_tissue') |> mutate(dense = n_tissue / bins_per_block >= rp$min_fill)
  dense  <- blocks |> filter(dense) |> mutate(block_id = row_number())

  # Chebyshev-neighbourhood offsets of radius r (in blocks)
  offsets <- function(r) expand_grid(dx = -r:r, dy = -r:r) |> filter(!(dx == 0 & dy == 0))

  # Connected components of dense blocks
  edges <- dense |>
    cross_join(offsets(rp$link_blocks)) |>
    mutate(nbx = bx + dx, nby = by + dy) |>
    inner_join(dense |> select(bx, by, block_id), by = c(nbx = 'bx', nby = 'by'), suffix = c('', '_nb')) |>
    select(from = block_id, to = block_id_nb)
  g <- igraph::graph_from_data_frame(edges, directed = FALSE, vertices = dense |> select(block_id))
  dense$piece <- igraph::components(g)$membership[as.character(dense$block_id)]

  # Piece sizes and position along the layout axis
  axis_col <- if (rp$layout == 'left_to_right') 'bx' else 'by'
  pieces <- dense |>
    group_by(piece) |>
    summarise(n_bins = sum(n_tissue), axis_pos = weighted.mean(.data[[axis_col]], n_tissue)) |>
    arrange(desc(n_bins))
  log_msg(cfg, lib, ': largest tissue pieces (8 um bins): ', paste(head(pieces$n_bins, n_reps + 2), collapse = ', '))

  # Warn when the split is ambiguous: next piece almost as big as the smallest kept one
  if (nrow(pieces) > n_reps && pieces$n_bins[n_reps + 1] > 0.25 * pieces$n_bins[n_reps])
    warning(lib, ': piece ', n_reps + 1, ' (', pieces$n_bins[n_reps + 1], ' bins) is close in size to piece ', n_reps,
            ' (', pieces$n_bins[n_reps], ') - check the QC overlay; consider layout: manual', call. = FALSE)
  if (nrow(pieces) < n_reps) stop(lib, ': found ', nrow(pieces), ' tissue pieces but the sample sheet lists ', n_reps)

  # N largest pieces, ordered along the axis
  strips <- pieces |> slice_head(n = n_reps) |> arrange(axis_pos)

  # Each occupied block takes the replicate of the nearest strip-core block within attach_blocks
  list(strips = strips, dense = dense, blocks = blocks,
       attach = function(rep_ids) {
         strips$replicate_id <- rep_ids
         core <- dense |> inner_join(strips |> select(piece, replicate_id), by = 'piece') |> select(cbx = bx, cby = by, replicate_id)
         blocks |>
           cross_join(bind_rows(tibble(dx = 0, dy = 0), offsets(rp$attach_blocks))) |>
           mutate(cbx = bx + dx, cby = by + dy, d = pmax(abs(dx), abs(dy))) |>
           inner_join(core, by = c('cbx', 'cby')) |>
           slice_min(d, n = 1, with_ties = FALSE, by = c(bx, by)) |>
           select(bx, by, replicate_id)
       })
}

# ---- Per library ----
summary_rows  <- list()
overlay_plots <- list()
for (lib in unique(sheet$library_id)) {
  lib_rows <- sheet |> filter(library_id == lib)          # sample-sheet order = layout order
  sr_dir   <- lib_rows$spaceranger_dir[1]
  sf       <- jsonlite::read_json(file.path(spatial_dir(sr_dir), 'scalefactors_json.json'))

  # 8 um tissue bins in image coordinates (um); array_row/col are rotated vs the H&E, so use pixels
  pos <- nanoparquet::read_parquet(file.path(sr_dir, 'binned_outputs/square_008um/spatial/tissue_positions.parquet')) |>
    as_tibble() |>
    filter(in_tissue == 1) |>
    mutate(x_um = pxl_col_in_fullres * sf$microns_per_pixel,
           y_um = pxl_row_in_fullres * sf$microns_per_pixel,
           bx = floor(x_um / rp$block_um), by = floor(y_um / rp$block_um))

  if (rp$layout == 'manual') {
    # Manual labels (8 um barcodes); bins not listed stay unassigned
    manual <- read_csv(file.path(rp$manual_dir, paste0(lib, '.csv')), show_col_types = FALSE)
    unknown <- setdiff(unique(manual$replicate_id), lib_rows$replicate_id)
    if (length(unknown)) stop(lib, ': manual replicate_id not in sample sheet: ', paste(unknown, collapse = ', '))
    bins_008 <- pos |> left_join(manual |> select(barcode, replicate_id), by = 'barcode')
  } else {
    det <- detect_pieces(pos, nrow(lib_rows), lib)
    bins_008 <- pos |> left_join(det$attach(lib_rows$replicate_id), by = c('bx', 'by'))
  }
  bins_008 <- bins_008 |> mutate(library_id = lib)
  log_msg(cfg, lib, ': ', sum(is.na(bins_008$replicate_id)), ' of ', nrow(bins_008), ' tissue bins left unassigned')

  # ---- Propagate to 16 um bins and cells (majority of their 2 um squares) ----
  map2 <- nanoparquet::read_parquet(file.path(sr_dir, 'barcode_mappings.parquet')) |>
    as_tibble() |>
    inner_join(bins_008 |> filter(!is.na(replicate_id)) |> select(square_008um = barcode, replicate_id), by = 'square_008um')
  majority <- function(df, id_col) {
    df |> count(barcode = .data[[id_col]], replicate_id) |>
      slice_max(n, n = 1, with_ties = FALSE, by = barcode) |>
      transmute(library_id = lib, barcode, replicate_id)
  }
  bins_016 <- majority(map2, 'square_016um')
  cells    <- majority(map2 |> filter(in_cell, !is.na(cell_id)), 'cell_id')

  # ---- Write ----
  write_csv(bins_008 |> select(library_id, barcode, replicate_id), tab_path(cfg, paste0(lib, '_008um.csv.gz')))
  write_csv(bins_016, tab_path(cfg, paste0(lib, '_016um.csv.gz')))
  write_csv(cells,    tab_path(cfg, paste0(lib, '_cells.csv.gz')))

  summary_rows[[lib]] <- full_join(count(bins_008, replicate_id, name = 'bins_008um'),
                                   count(bins_016, replicate_id, name = 'bins_016um'), by = 'replicate_id') |>
    full_join(count(cells, replicate_id, name = 'cells'), by = 'replicate_id') |>
    mutate(library_id = lib, .before = 1)

  # ---- QC overlay: assignment over the H&E ----
  md <- bins_008 |>
    transmute(library_id, x_fullres = pxl_col_in_fullres, y_fullres = pxl_row_in_fullres,
              replicate = replace_na(replicate_id, 'unassigned'))
  pal <- c(make_palette(lib_rows$replicate_id), unassigned = 'grey60')
  cfg_8 <- modifyList(cfg, list(bin_size = '8', res_label = '008um'))     # tiles drawn at 8 um
  overlay_plots[[lib]] <- plot_spatial(md, 'replicate', cfg_8, lib_rows, palette = pal,
                                       title_prefix = 'Replicate assignment')[[lib]]
}

# One PDF page per library (+ replicate_assignment_<library>.png)
save_plot_list(overlay_plots, cfg, 'replicate_assignment')

# ---- Summary ----
summary_tbl <- bind_rows(summary_rows)
print(summary_tbl)
write_csv(summary_tbl, tab_path(cfg, 'replicate_summary.csv'))
step_finish(cfg)
