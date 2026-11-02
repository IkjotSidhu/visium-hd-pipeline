# Visium HD Analysis Pipeline

Takes 10x Visium HD Space Ranger 4 output from raw counts to QC, clusters, markers,
spatial domains, condition comparisons and pathways — for **8 µm bins, 16 µm bins and
segmented cells** in one run. Capture areas that hold several tissue pieces (e.g.
replicates or conditions side by side) are split automatically.

Runs on Minerva as a chain of LSF jobs. One config file, one command.

---

## Quick start

**1. Make a sample sheet** — a CSV with one row per tissue piece:

| replicate_id | library_id | spaceranger_dir | donor | condition |
|---|---|---|---|---|
| D1_Control | Run1 | data/Run1 | D1 | Control |
| D1_Treated | Run1 | data/Run1 | D1 | Treated |
| D2_Control | Run2 | data/Run2 | D2 | Control |
| D2_Treated | Run2 | data/Run2 | D2 | Treated |

- `replicate_id` — a unique name for the piece
- `library_id` — the Space Ranger run (capture area) the piece is on
- `spaceranger_dir` — that run's output folder (the one containing `binned_outputs/`)
- Add any other columns you want carried onto every bin/cell (`donor`, `condition`, `tissue`, ...)

**2. Copy and edit the config:**

```bash
cp /path/to/visium-hd-pipeline/config_template.yaml config.yaml
```

At minimum set `sample_sheet`, and — if several pieces share a capture area —
`replicates: layout: left_to_right` (pieces listed in the sheet in the order they sit on
the slide). Every option is explained in the template.

**3. Run it:**

```bash
bash /path/to/visium-hd-pipeline/run_pipeline.sh config.yaml
```

Follow the jobs with `bjobs`; logs are in `output/logs/`. When the jobs finish,
open **`output/index.html`** — one page with QC, clustering, DE and pathway summaries
and every figure, for all bin sizes.

---

## Check these first

| Look at | To check |
|---|---|
| `output/01_assign_replicates/figures/replicate_assignment_*.png` | Each tissue piece got the right replicate name |
| `output/008um/02_load_qc/qc_report.html` | How many bins each group, sample and replicate kept, and why bins were removed |
| `output/008um/03_cluster/figures/umap_metadata_donor.png` | Donors mix after batch correction |
| `output/008um/03_cluster/res_1.0/figures/spatial_clusters_*.png` | Clusters on the H&E, per capture area and per strip |
| `output/008um/04_markers/res_1.0/tables/top10_markers_*.csv` | Marker genes for naming clusters |
| `output/008um/06_pseudobulk_de/tables/pseudobulk_de_status.csv` | Which comparisons were tested, and why any were skipped |

The same folders exist for `016um/` and `cells/`.

---

## Common tasks

**Try a new dataset quickly** — run on 10% of each replicate, written to `<output_dir>_test/`:
```bash
bash run_pipeline.sh config.yaml -t 0.1
```

**Preview the jobs without submitting:**
```bash
bash run_pipeline.sh config.yaml -n
```

**Rerun from a given step** (earlier results are reused):
```bash
bash run_pipeline.sh config.yaml -s 03-10        # steps 03 to 10
bash run_pipeline.sh config.yaml -s 04-04        # just step 04
bash run_pipeline.sh config.yaml -s 10-10        # just rebuild index.html
```

**Run only some bin sizes:**
```bash
bash run_pipeline.sh config.yaml -r 8            # 8 µm only
bash run_pipeline.sh config.yaml -r 16,cells
```

**Run one step interactively:**
```bash
module purge && module load R/4.4.1
Rscript /path/to/visium-hd-pipeline/steps/03_cluster.R config.yaml 8
```

**Change QC thresholds** — one value for all bin sizes, or one per bin size
(`max_mt_pct: 100` turns the mitochondrial filter off):
```yaml
qc:
  min_counts: {"8": 50, "16": 100, cells: 50}
  max_mt_pct: 50
```

**Choose how the QC report groups samples:**
```yaml
qc:
  report:
    sample_col: visium_sample_id           # default: library_id
    group_cols: [condition, donor, tissue] # each summarised separately
```

**Flag low-quality spots (SpotSweeper)** — after QC filtering, each bin/cell is compared with
its nearest neighbours in the same strip and labelled if its UMIs or genes are far lower, or
its mitochondrial % far higher. Nothing is removed. Flags are stored in the metadata
(`spotsweeper`, `spotsweeper_any`, one flag and z-score per metric). Maps, tables and the
share of flagged units per cluster at every resolution are in `<bin size>/secondary_qc/`.
```yaml
spotsweeper:
  metrics: [umi, genes, mito]
  n_neighbors: 36
  cutoff: 2.5
```

**Choose clustering resolutions** — the default sweeps 0.1 to 2.0 in steps of 0.1; every
resolution gets its own UMAP, spatial maps and markers. List fewer to save time and space:
```yaml
clustering:
  resolutions: [0.2, 0.5, 1, 2]
  main_resolution: 1            # used for DE and pathways
```

**Spatial plots** — every bin or cell is drawn as a dot (vector in the PDFs) over the H&E,
once per capture area (`<plot>_<library_id>`) and once per tissue strip (`<plot>_<replicate_id>`).
Strip pages use the original microscope scan when Space Ranger's `--image` file is inside the
Space Ranger folder, otherwise the hires image. Continuous maps (UMIs, mitochondrial %) use a
linear colour scale from 0 to the 99th percentile. To change the look:
```yaml
plots:
  image_alpha: 1                # H&E strength (1 = full, lower fades toward white)
  dpi: 300                      # PNG resolution
  point_size: null              # null = one bin wide; or a fixed size
  strip_um_per_px: 0.8          # H&E resolution of strip pages (smaller = sharper, larger files)
```

**Batch correction** — Harmony runs automatically when the sample sheet has more than one
donor (or, without a `donor` column, more than one sample). To correct for a different
column, or to switch it off, set and rerun from step 03:
```yaml
clustering:
  integration_by: batch          # any sample-sheet column
  integration: none              # disable
```

**Compare conditions** — pseudobulk per replicate, optionally paired by donor:
```yaml
de:
  run: true
  condition_col: condition
  reference: Control
  covariates: [donor]
```
A subset (e.g. a cluster) that cannot support the requested design is skipped, with the
reason in `pseudobulk_de_status.csv`. A gene counts as `consistent_across_strata` only
when its direction agrees in every donor.

**Add pathway enrichment** — give a GMT file whose gene IDs match the DE results
(gene symbols for the 10x probe sets), then rerun `-s 09-10`:
```yaml
pathways:
  gmt_file: metadata/gene_sets/reactome_human_v97.gmt
  collection: Reactome human pathways
  version: "97"
```

**Apply cell-type labels** — turn on step 08 and run it once to get an
`annotation_template.csv` per resolution (`08_annotation/res_<r>/tables/`). Copy it
out of `output/`, fill `broad_type`, `confidence` (low/medium/high) and `evidence`,
keep the `cluster` and `clustering_key` columns, then point the config at it and rerun
`-s 08-10`:
```yaml
annotation:
  run: true
  mapping_files:
    "8":
      "1.0": metadata/annotations_8_res1.0.csv
```
A label file made for a different clustering is rejected, so rerun step 08 for a fresh
template after reclustering.

**Deconvolve with a single-cell reference (RCTD):**
```yaml
rctd:
  run: true
  reference_rds: /path/to/reference_seurat.rds
  celltype_col: cell_type
  reference_assay: RNA           # the assay holding raw counts
```
Results are cached per replicate and recomputed automatically when the inputs or
settings change (`force_recompute: true` forces a fresh run).

**Split tissue pieces by hand** (irregular layouts) — export one CSV per capture area
from Loupe Browser with columns `barcode` (8 µm) and `replicate_id`, save them as
`<library_id>.csv` in one folder, and set:
```yaml
replicates:
  layout: manual
  manual_dir: replicate_labels/
```

**Give a step more memory or time:**
```yaml
lsf:
  resources:
    03_cluster: {cores: 8, mem_mb: 48000, walltime: "24:00"}
```

---

## What each step does

| Step | Does | Key outputs |
|---|---|---|
| **01** assign replicates | Finds each tissue piece and names it from the sample sheet | Assignment tables, overlay on the H&E, one H&E image per strip |
| **02** load & QC | Loads bins/cells with the high-resolution H&E, adds sample-sheet columns, filters, labels SpotSweeper outliers | `qc_report.html`, QC tables by group/sample/replicate, Space Ranger metrics, UMI-vs-genes plots, QC maps |
| **03** cluster | Normalises, sketches, clusters at each resolution, projects to every bin | UMAPs (also split by condition), clustree, cluster make-up per replicate and condition, cluster maps |
| **04** markers | Marker genes for each cluster, at each resolution | Marker tables, dot plot, heatmap |
| **05** BANKSY | Spatial domains (expression + neighbourhood) | Domain labels and maps |
| **06** pseudobulk DE | Condition comparison with replicates as the unit, whole tissue and per cluster | Results, summary, skip reasons, volcano, sample PCA (raw and donor-corrected) and distance plots |
| **07** RCTD | Cell-type proportions from a single-cell reference | Weights, dominant type maps, composition |
| **08** annotation | Applies your reviewed cluster labels | Label template, labelled maps, composition per replicate |
| **09** pathways | Preranked enrichment (fgsea) from the DE statistics | Enrichment table, leading-edge genes, summary plots |
| **10** report | Collects everything into one page | `output/index.html` |

Steps 05–10 can be switched on or off with `run:` in the config. Step 07 (needs a
reference) and step 08 (needs label files) are off by default; step 09 is skipped
when no GMT file is given.

## Output layout

```text
output/
├── 01_assign_replicates/        strip_images/ holds the per-strip H&E
├── 008um/                       also 016um/ and cells/
│   ├── 02_load_qc/              qc_report.html, figures/, tables/
│   ├── 03_cluster/              shared UMAPs + clustree; res_<r>/ per resolution
│   ├── 04_markers/res_<r>/
│   ├── 05_banksy/res_<r>/
│   ├── 06_pseudobulk_de/
│   ├── 07_rctd/  08_annotation/res_<r>/  09_pathways/
│   ├── secondary_qc/            SpotSweeper flags: maps, tables, flagged share per cluster
│   └── rds/                     Seurat objects passed between steps
├── index.html                   overview of the whole run
└── logs/
```

Every step folder holds `figures/`, `tables/`, `params_used.yaml` (the exact settings
used) and `sessionInfo.txt`.

---

## Requirements

Minerva with the `R/4.4.1` and `gdal/3.9.1` modules (GDAL reads the original scan for
strip pages). The R packages are already in the lab R library (`r_library` in the config). To rebuild the environment elsewhere:
`renv::restore()` with the included `renv.lock`.

## More

- [FAQ](FAQ.md) — QC choices, replicate splitting, statistics with few replicates, RCTD and HPC troubleshooting.
- `config_template.yaml` — every option, with defaults.
- `tests/` — quick checks after changing the code: `Rscript tests/<name>.R` (with `module load R/4.4.1`).
