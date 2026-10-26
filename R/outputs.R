# ---- Output helpers ----
# Shared by annotation, pseudobulk diagnostics, enrichment and the static report.
# No packages are attached here: configuration planning must remain lightweight.

clustering_key <- function(obj, column) {
  digest::digest(list(column = column, cells = colnames(obj), labels = as.character(obj[[column]][, 1])),
                 algo = 'sha256')
}

read_annotation_map <- function(path, clusters, key) {
  template <- data.frame(cluster = as.character(clusters), broad_type = '', cell_state = '',
                         confidence = 'unreviewed', evidence = '', clustering_key = key)
  if (is.null(path)) return(template)
  mapping <- utils::read.csv(path, colClasses = 'character', check.names = FALSE, na.strings = c('', 'NA'))
  required <- names(template)
  if (!all(required %in% names(mapping))) stop('Annotation mapping needs columns: ', paste(required, collapse = ', '))
  if (anyNA(mapping$cluster) || anyDuplicated(mapping$cluster)) stop('Annotation cluster IDs must be present and unique')
  if (any(!mapping$cluster %in% clusters)) stop('Annotation contains cluster IDs absent from this object')
  if (any(is.na(mapping$clustering_key) | mapping$clustering_key != key))
    stop('Annotation mapping belongs to a different clustering; regenerate and review the template')
  assigned <- !is.na(mapping$broad_type) & nzchar(trimws(mapping$broad_type))
  if (any(assigned & (is.na(mapping$confidence) | !mapping$confidence %in% c('low', 'medium', 'high'))) ||
      any(assigned & (is.na(mapping$evidence) | !nzchar(trimws(mapping$evidence)))))
    stop('Assigned labels require low/medium/high confidence and nonempty marker/histology evidence')
  idx <- match(template$cluster, mapping$cluster)
  present <- !is.na(idx)
  template[present, ] <- mapping[idx[present], required]
  template
}

# Use a full variance-stabilizing transformation rather than the fast vst()
# approximation, which requires enough genes above its mean-count threshold.
write_pseudobulk_diagnostics <- function(dds, metadata, cfg, label) {
  if (anyDuplicated(metadata$replicate_id)) stop('Diagnostic metadata must have one row per replicate')
  vst <- SummarizedExperiment::assay(DESeq2::varianceStabilizingTransformation(dds, blind = FALSE))
  metadata <- metadata[match(colnames(dds), metadata$replicate_id), , drop = FALSE]
  stopifnot(identical(as.character(metadata$replicate_id), colnames(dds)))
  metadata$library_size <- colSums(DESeq2::counts(dds))
  metadata$size_factor <- DESeq2::sizeFactors(dds)
  write_csv(metadata, tab_path(cfg, paste0(label, '_sample_diagnostics.csv')))
  design <- model.matrix(DESeq2::design(dds), as.data.frame(SummarizedExperiment::colData(dds)))
  write_csv(as_tibble(design, rownames = 'replicate_id'), tab_path(cfg, paste0(label, '_model_matrix.csv')))
  saveRDS(vst, file.path(cfg$step_dir, paste0(label, '_vst.rds')))
  # PCA on the high-variance genes, saved as <label>_pca[_suffix]_<group> for each metadata column
  plot_pca <- function(mat, suffix, subtitle) {
    variance <- apply(mat, 1, var)
    selected <- names(sort(variance[is.finite(variance) & variance > 0], decreasing = TRUE))
    selected <- head(selected, cfg$de$diagnostics$top_variable_genes)
    if (length(selected) < 2 || ncol(mat) < 3) stop('PCA requires >=2 variable genes and >=3 samples')
    pc <- prcomp(t(mat[selected, , drop = FALSE]), center = TRUE, scale. = FALSE)
    percent <- 100 * pc$sdev^2 / sum(pc$sdev^2)
    scores <- as_tibble(pc$x[, 1:2, drop = FALSE], rownames = 'replicate_id') |>
      left_join(metadata, by = 'replicate_id')
    prefix <- paste0(label, '_pca', suffix)
    write_csv(scores, tab_path(cfg, paste0(prefix, '_scores.csv')))
    write_csv(tibble(gene = selected, variance = variance[selected]), tab_path(cfg, paste0(prefix, '_genes.csv')))
    write_csv(tibble(component = seq_along(percent), percent_variance = percent),
              tab_path(cfg, paste0(prefix, '_variance.csv')))
    groups <- intersect(unique(c(cfg$de$condition_col, cfg$de$diagnostics$metadata_cols)), colnames(scores))
    for (group in groups) {
      p <- ggplot(scores, aes(PC1, PC2, colour = factor(.data[[group]]))) +
        geom_point(size = 3) +
        ggrepel::geom_text_repel(aes(label = replicate_id), seed = cfg$seed, max.overlaps = Inf,
                                 box.padding = 0.5, point.padding = 0.4, size = 3, show.legend = FALSE) +
        scale_x_continuous(expand = expansion(mult = 0.2)) +
        scale_y_continuous(expand = expansion(mult = 0.2)) +
        labs(title = paste('Pseudobulk PCA:', label), colour = group, subtitle = subtitle,
             x = sprintf('PC1 (%.1f%%)', percent[1]), y = sprintf('PC2 (%.1f%%)', percent[2])) + theme_classic()
      stem <- paste0(prefix, '_', gsub('[^[:alnum:]_-]', '_', group))
      ggsave(fig_path(cfg, paste0(stem, '.png')), p, width = 9, height = 7, dpi = 100)
      ggsave(fig_path(cfg, paste0(stem, '.pdf')), p, width = 9, height = 7)
    }
  }
  plot_pca(vst, '', 'VST; no batch correction')

  # Same PCA with the DE covariates (e.g. donor) regressed out while keeping the
  # condition effect, as in the DESeq2 vignette. For display only; DE itself
  # adjusts for covariates through its design formula.
  cd <- as.data.frame(SummarizedExperiment::colData(dds))
  covars <- setdiff(all.vars(DESeq2::design(dds)), cfg$de$condition_col)
  if (length(covars)) {
    corrected <- limma::removeBatchEffect(vst, covariates = model.matrix(reformulate(covars), cd)[, -1, drop = FALSE],
                                          design = model.matrix(reformulate(cfg$de$condition_col), cd))
    plot_pca(corrected, '_corrected', paste('VST;', paste(covars, collapse = ' + '), 'removed (display only)'))
  }

  # Distances use all retained genes; PCA uses the selected high-variance genes.
  distances <- as.matrix(dist(t(vst)))
  write_csv(as_tibble(distances, rownames = 'replicate_id'), tab_path(cfg, paste0(label, '_sample_distances.csv')))
  order <- rownames(distances)[hclust(as.dist(distances))$order]
  long <- as_tibble(distances, rownames = 'sample_a') |>
    pivot_longer(-sample_a, names_to = 'sample_b', values_to = 'distance') |>
    mutate(sample_a = factor(sample_a, levels = order), sample_b = factor(sample_b, levels = rev(order)))
  p <- ggplot(long, aes(sample_a, sample_b, fill = distance)) + geom_tile() +
    scale_fill_gradient(low = 'white', high = '#2166AC') + coord_equal() +
    theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = paste('VST sample distances:', label), x = NULL, y = NULL)
  ggsave(fig_path(cfg, paste0(label, '_sample_distances.png')), p, width = 8, height = 7, dpi = 100)
  ggsave(fig_path(cfg, paste0(label, '_sample_distances.pdf')), p, width = 8, height = 7)
  invisible(vst)
}

read_gmt <- function(path) {
  fields <- strsplit(readLines(path, warn = FALSE), '\t', fixed = TRUE)
  if (!length(fields) || any(lengths(fields) < 3)) stop('GMT requires pathway, description and at least one gene per line')
  names <- vapply(fields, `[`, character(1), 1)
  if (any(!nzchar(names)) || anyDuplicated(names)) stop('GMT pathway names must be nonempty and unique')
  setNames(lapply(fields, function(x) unique(x[-c(1, 2)][nzchar(x[-c(1, 2)])])), names)
}

# Keep signed Wald statistics, including non-significant genes. Duplicate gene
# identifiers are an input problem, not a reason to pick the most favourable row.
pathway_ranks <- function(tbl) {
  if (anyDuplicated(tbl$gene)) stop('Duplicate gene identifiers in pathway ranking')
  keep <- !is.na(tbl$gene) & nzchar(tbl$gene) & is.finite(tbl$stat)
  stats <- setNames(tbl$stat[keep], tbl$gene[keep])
  stats[order(-stats, names(stats))]
}

html_escape <- function(x) {
  x <- gsub('&', '&amp;', as.character(x), fixed = TRUE)
  x <- gsub('<', '&lt;', x, fixed = TRUE)
  x <- gsub('>', '&gt;', x, fixed = TRUE)
  x <- gsub('"', '&quot;', x, fixed = TRUE)
  gsub("'", '&#39;', x, fixed = TRUE)
}

html_table <- function(tbl, max_rows = 20L) {
  if (!nrow(tbl)) return('<p>No rows.</p>')
  shown <- head(tbl, max_rows)
  header <- paste0('<th>', html_escape(names(shown)), '</th>', collapse = '')
  rows <- apply(as.data.frame(shown), 1, function(r)
    paste0('<tr>', paste0('<td>', html_escape(r), '</td>', collapse = ''), '</tr>'))
  paste0('<div class="table-scroll"><table><thead><tr>', header, '</tr></thead><tbody>',
         paste(rows, collapse = ''), '</tbody></table></div><p>Showing ', nrow(shown), ' of ', nrow(tbl), ' rows.</p>')
}

# ---- Clustering-resolution output contexts ----
# Bin sizes (8/16/cells) and clustering resolutions are separate dimensions.
resolution_label <- function(resolution) format(as.numeric(resolution), scientific = FALSE, trim = TRUE, nsmall = 1, digits = 12)
cluster_column <- function(resolution) paste0('cluster_res_', format(as.numeric(resolution), scientific = FALSE, trim = TRUE, digits = 12))
resolution_cfg <- function(cfg, resolution, create = TRUE) {
  cfg$cluster_resolution <- as.numeric(resolution)
  cfg$step_dir <- step_dir(cfg, cfg$step)
  if (create) for (name in c('figures', 'tables')) dir.create(file.path(cfg$step_dir, name), recursive = TRUE, showWarnings = FALSE)
  cfg
}
assert_cluster_resolutions <- function(obj, resolutions) {
  columns <- vapply(as.numeric(unlist(resolutions)), cluster_column, character(1))
  missing <- setdiff(columns, colnames(obj[[]]))
  if (length(missing)) stop('Clustering columns missing: ', paste(missing, collapse = ', '), '; rerun step 03 with the requested resolution sweep')
  invisible(TRUE)
}
annotation_mapping <- function(cfg) {
  entry <- cfg$annotation$mapping_files[[cfg$bin_size]]
  resolution <- cfg$cluster_resolution %||% cfg$clustering$main_resolution
  if (is.list(entry)) return(entry[[resolution_label(resolution)]] %||% entry[[as.character(resolution)]])
  # A legacy single-file mapping belongs only to the primary resolution.
  if (isTRUE(all.equal(as.numeric(resolution), as.numeric(cfg$clustering$main_resolution)))) entry else NULL
}
