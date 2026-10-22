# Visium HD Analysis Pipeline

Takes 10x Visium HD Space Ranger 4 output from raw counts to QC, clusters, markers,
spatial domains and condition comparisons — for **8 µm bins, 16 µm bins and segmented
cells** in one run. Capture areas that hold several tissue pieces (e.g. replicates or
conditions side by side) are split automatically.

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

The default clustering sweep is **0.1, 0.2, …, 2.0** for each bin/cell size.
Every resolution gets its own UMAP, spatial cluster maps, composition plots and
marker dot plots/heatmaps. BANKSY spatial domains use the same default sweep.

At minimum set `sample_sheet`, and — if several pieces share a capture area —
`replicates: layout: left_to_right` (pieces listed in the sheet in the order they sit on
the slide). Every option is explained in the template.

**3. Run it:**

```bash
bash /path/to/visium-hd-pipeline/run_pipeline.sh config.yaml
```

Results appear in `output/`. Open `output/index.html` for the overview after completion.
Follow the jobs with `bjobs`; logs are in `output/logs/`.

---

## Check these first

| Look at | To check |
|---|---|
| `output/01_assign_replicates/figures/replicate_assignment_*.png` | Each tissue piece got the right replicate name |
| `output/008um/02_load_qc/tables/qc_summary.csv` | How many bins each replicate kept, and why bins were removed |
| `output/008um/03_cluster/figures/umap_metadata_donor.png` | Donors mix after batch correction |
| `output/008um/03_cluster/res_1.0/figures/spatial_clusters_*.png` | Clusters on the H&E |
| `output/008um/04_markers/res_1.0/tables/top10_markers_*.csv` | Marker genes for naming clusters |

The same folders exist for `016um/` and `cells/`.

---

## Common tasks

**Try a new dataset quickly** — run on 10% of each replicate, written to `output_test/`:
```bash
bash run_pipeline.sh config.yaml -t 0.1
```

**Preview the jobs without submitting:**
```bash
bash run_pipeline.sh config.yaml -n
```

**Rerun from a given step** (earlier results are reused):
```bash
bash run_pipeline.sh config.yaml -s 03-07        # steps 03 to 07
bash run_pipeline.sh config.yaml -s 04-04        # just step 04
```

**Run only some resolutions:**
```bash
bash run_pipeline.sh config.yaml -r 8            # 8 µm only
bash run_pipeline.sh config.yaml -r 16,cells
```

**Run one step interactively:**
```bash
module purge && module load R/4.4.1
Rscript /path/to/visium-hd-pipeline/steps/03_cluster.R config.yaml 8
```

**Change QC thresholds** — in the config; one value for all resolutions or one per resolution:
```yaml
qc:
  min_counts: {"8": 50, "16": 100, cells: 50}
  max_mt_pct: 50
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

DE preserves requested covariates by default and writes a status table for skipped
comparisons. The summary includes donor/stratum support; agreement in one donor
does not count as consistency across donors.

**Deconvolve with a single-cell reference (RCTD):**
```yaml
rctd:
  run: true
  reference_rds: /path/to/reference_seurat.rds
  celltype_col: cell_type
  reference_assay: RNA
```

RCTD checks for raw integer counts in the selected assay. Changed inputs or settings
invalidate its cache; `force_recompute: true` forces a fresh fit.

**Review and apply cell-type labels:**

Step 08 produces `08_annotation/res_0.1/tables/annotation_template.csv` through
`res_2.0/` for each bin/cell size. Each clustering has its own fingerprint.
Copy each template outside the output folder, fill `broad_type`, optional `cell_state`,
`confidence` (`low`, `medium`, `high`) and `evidence`, retaining `cluster` and
`clustering_key`. Configure one reviewed file per resolution:

```yaml
annotation:
  mapping_files:
    "8":
      "0.1": metadata/annotations_8_res0.1.csv
      "1.0": metadata/annotations_8_res1.0.csv
    "16":
      "1.0": metadata/annotations_16_res1.0.csv
  require_complete: false
```

A legacy single-file mapping per bin size is used only for `main_resolution`.
Rerun `-s 08-10`. Missing labels remain `Unannotated`; the original cluster IDs
are retained. Mappings from a different clustering are rejected. Marker tests are
annotation evidence, not condition-level inference. Legacy marker tables without
a matching fingerprint are labelled `legacy_unverified`; rerun step 04 to verify them.

**Add a pathway collection:**

```yaml
pathways:
  gmt_file: metadata/gene_sets/reactome_human_v97.gmt
  collection: Reactome human pathways
  version: "97"
  gene_id_type: symbol
```

Rerun `-s 09-10`. GMT gene identifiers must match the DE result identifiers; no ID
conversion is performed. Enrichment uses all finite signed Wald statistics,
including nonsignificant genes. Coverage tables explain which sets meet size
limits; FDR is within each subset/contrast. Ranking ties and numerical warnings are
recorded in the pathway status table; inspect them before interpreting results.
The collection version and checksum
are saved. With no GMT, step 09 reports a skip and writes empty results.

**Rebuild only the overview:**

```bash
bash run_pipeline.sh config.yaml -s 10-10
```

The report needs no pandoc or web server. Copy the whole output folder to preserve
its relative links. It shows available/disabled/missing artifacts and timestamps;
partial reruns can contain outputs from different runs, so inspect provenance.
A `-r 8` filter limits the overview to the selected resolution.

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
| **01** assign replicates | Finds each tissue piece and names it from the sample sheet | Assignment tables, overlay on the H&E |
| **02** load & QC | Loads bins/cells with the high-resolution H&E, adds sample-sheet columns, filters | QC summary, violin plots, UMI / MT% / pass-fail maps |
| **03** cluster | Normalises, sketches, clusters at several resolutions, projects to every bin | UMAPs, clustree, cluster make-up per replicate, cluster maps |
| **04** markers | Marker genes for each cluster | Marker tables, dot plot, heatmap |
| **05** BANKSY | Spatial domains (expression + neighbourhood) | Domain labels and maps |
| **06** pseudobulk DE | Condition comparison with replicates as the unit, whole tissue and per cluster | Results table, summary, volcano plot |
| **07** RCTD | Cell-type proportions from a single-cell reference | Weights, dominant type maps, composition |
| **08** annotation | Applies reviewed cluster labels and gathers marker evidence | Mapping template, annotations, spatial maps, replicate composition |
| **09** pathways | Preranked enrichment from pseudobulk Wald statistics | Enrichment/FDR, leading-edge genes, overlap coverage, collection provenance |
| **10** report | Collects summaries and figures across resolutions | Static `output/index.html`, artifact index, settings and timestamps |

Steps 05–10 can be switched off in the config (`run: false`). Pathways also require
DE to be enabled. The report is submitted once, after all selected resolution chains.

## Output layout

```text
output/
├── 01_assign_replicates/
├── 008um/                       also 016um/ and cells/
│   ├── 03_cluster/
│   │   ├── figures/             shared metadata UMAPs and clustree
│   │   ├── tables/              comparison across all 20 resolutions
│   │   ├── res_0.1/             figures/, tables/, settings and session
│   │   ├── res_0.2/
│   │   └── ... res_2.0/
│   ├── 04_markers/res_0.1/ ... res_2.0/
│   ├── 05_banksy/res_0.1/ ... res_2.0/
│   ├── 08_annotation/res_0.1/ ... res_2.0/
│   ├── 02_load_qc/  06_pseudobulk_de/  07_rctd/  09_pathways/
│   └── rds/                     one shared clustering object, per-resolution annotations
├── 10_report/
├── index.html                   expandable gallery for every resolution
└── logs/
```

`main_resolution` selects the default identity and the clustering for DE/pathway
comparisons. Plotting, marker discovery and annotation run at **every** configured
clustering resolution. The shared Seurat object stores all cluster labels, without
copying the full expression matrix 20 times. Missing requested cluster columns
cause an explicit error directing you to rerun step 03.

Each clustering folder contains UMAP, per-library spatial maps, replicate
composition and donor-mixing tables. Marker folders contain full/top marker tables
and PNG/PDF dot plots and heatmaps. A no-marker resolution has explanatory figures
and a status table. The report includes every resolution's PNGs without truncation.

`params_used.yaml` records the exact settings each step ran with.

Step 06 also writes VST sample PCA plots coloured by available metadata, sample-distance
heatmaps, library sizes, size factors, PCA coordinates/genes and model matrices.
`pseudobulk_diagnostics_status.csv` records success or failure per tested subset.

Step 06 also writes:

- `pseudobulk_de_status.csv`: tested/skipped subsets, reasons and fitted designs.
- `pseudobulk_sample_eligibility.csv`: units retained per replicate and subset.
- `pseudobulk_stratum_effects.csv.gz`: descriptive CPM effects within each donor/stratum.
- `rds/06_pseudobulk_de.rds`: raw and DESeq2-normalized counts, sample metadata and designs for tested subsets.
- `interpretation_notes.txt`: interpretation guidance and configured study notes.

Step 07 writes `rctd_reference_cells.csv.gz` and `rctd_cache_manifest.csv` for
reference inclusion and cache provenance.

---

## Requirements

Minerva with the `R/4.4.1` module. The packages are already in the lab R library
(`r_library` in the config). To rebuild the environment elsewhere:
`renv::restore()` with the included `renv.lock`.

## More

- [FAQ](FAQ.md) — QC choices, replicate splitting, statistics with few replicates, RCTD and HPC troubleshooting.
- `config_template.yaml` — every option, with defaults.

Focused checks (from the repository root):

```bash
module purge && module load R/4.4.1
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 Rscript tests/regression.R
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 Rscript tests/outputs.R
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 Rscript tests/resolutions.R
```
