# FAQ — Lessons & Limitations from Real Analyses

[← Back to README](README.md)

Non-obvious issues and design decisions distilled from production Visium HD runs
(the intestinal Duodenum/Ileum project). Each answer includes a concrete example.

## Environment & HPC (Minerva / LSF)

**Q: My R packages won't load on the cluster — I get ABI / symbol errors. What's wrong?**
Always `module purge && module load R/4.4.1` and set
`.libPaths("/sc/arion/projects/naiklab/ikjot/R")`. That library is compiled
against the Rocky9 R/4.4.1 module; **conda R or any other R build causes ABI
mismatches and package load failures.** *Example:* switching to a conda R env to
"save time" broke `Seurat`/`spacexr` loading with cryptic symbol errors — reverting
to the module fixed it instantly. Never use conda R for project scripts.

**Q: My single-core job (e.g. RCTD) segfaults with no R error message. Why?**
OpenBLAS spawns one thread **per physical core** by default. When a step requests a
single core (`create.RCTD(max_cores = 1)`), the extra threads exhaust `RLIMIT_NPROC`
and the process segfaults silently. **Fix:** pin all BLAS/OMP backends to one thread
(`export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLAS_NUM_THREADS=1`).
*Example:* RCTD on a 372k-spot 8µm sample died with no traceback until these caps
were added — `submit_step.sh` sets them for every step. This is the #1 "job died mysteriously" cause.

**Q: How much memory / walltime do I actually need?**
Profiles that worked in practice: **8µm RCTD** (~372k spots) needed **300 GB/core × 2
(600 GB)** after OOM-ing at 200 GB; **16µm RCTD** ran at ~200 GB/core; both on the
`premium` queue with **144 h** walltime and `span[hosts=1]`. `rusage[mem=...]` is **per
core**, not total. Start high on the largest sample and trim down.

**Q: `.Rmd` or `.R` on the cluster?**
Run heavy compute as a plain **`.R` via `bsub`** — which is why every pipeline step is
an `.R` script. Compute nodes don't put `pandoc` on the R module's PATH, so rendering an
`.Rmd` to HTML fails there (`pandoc version 1.12.3 or higher is required`); use
`module load pandoc`, or `rmarkdown::render(..., run_pandoc = FALSE)` when you only need
the chunks to run. *Caution:* keep any `.R`/`.Rmd` pair in sync — a stale `.R` diverging
from its notebook silently ran old logic in a past project.

**Q: A step failed halfway — do I rerun everything?**
No. Each step reads the previous step's RDS from `output/<res>/rds/`, so rerun just that
step: `bash run_pipeline.sh config.yaml -s 04-04 -r 8`. Steps save their object *before*
plotting, so a plotting error never costs the compute.

## RCTD Deconvolution

**Q: Can I use one shared reference for all samples/tissues?**
Prefer **tissue/group-specific references**, and understand the trade-off — this is a
real **limitation**, not just a setting. If a cell type is structurally absent from a
tissue but left in that tissue's reference, RCTD can assign it spuriously; but if you
*exclude* it, any cross-tissue comparison of that type becomes **"present in reference"
vs "structurally forced to zero"** — an artifact, not biology. *Example:*
`Mid_villus_proximal_enterocytes` is excluded from the Ileum reference and
`Mid_villus_distal_enterocytes` from the Duodenum reference; we therefore **drop those
types entirely** from cross-tissue composition plots and tests rather than compare them.

**Q: RCTD errors on my labels or counts.**
Three recurring fixes: (a) RCTD **prohibits `/` in cell-type labels** → sanitize to `-`;
(b) merged **SCT objects store non-integer RNA counts** → `round()` and coerce to
`dgCMatrix` before building the reference; (c) drop rare types (**min ~25 cells**) and
special labels (`Unknown`, `Mix cluster remove`). *Example:* a `TA_cells / Early_progenitors`
label crashed `Reference()` until the slash was replaced.

**Q: Do I have to recompute RCTD every run?**
No — **cache per-sample RCTD objects to RDS and reload**, recomputing only when the file
is missing. A full re-run of a 4-sample 8µm deconvolution is many hours; caching makes
iterating on the downstream plots cheap.

## Reference & Annotation Hygiene

**Q: Cell types duplicate, split, or lose their color in plots. Why?**
Label whitespace/Unicode inconsistencies in the annotation sheet. Sanitize consistently:
Greek α (U+03B1) → `aa`, em/en-dashes → `-`, collapse internal whitespace; RCTD column
names replace non-alphanumerics with `_`. **Don't trust the sheet — fix it in code.**
*Examples from one sheet:* a broad label `"Cd8aa _IEL"` (stray space) silently split into
a second cell type; and `Mid_villus_proximal_enterocytes` carried **two different hex
codes** across its clusters, so first-occurrence dedup picked the wrong color. Define the
palette explicitly in code rather than relying on sheet row order.

**Q: How are clusters mapped to names?**
Via an external CSV: cluster `#` → `Population` (fine, ~41 labels) → `Broad` (~24 labels).
Keep a two-level scheme so you can deconvolve at either granularity without re-annotating.

## Statistical Limitations (read this before claiming significance)

**Q: Why can't my n=2 vs n=2 comparison ever reach significance at the sample level?**
With **2 vs 2** replicates, a two-sided Wilcoxon's **minimum possible p-value is 1/3** —
it *cannot* be significant no matter how large the effect. Meanwhile **spot-level** tests
pool hundreds of thousands of spots and are **pseudoreplicated** (massively inflated
significance). **Neither alone is trustworthy.** *Approach used:* a two-tier rule — call a
result "high-confidence" only if it is **spot-level significant AND consistent in
direction across both replicates** in each group, and corroborate with **pseudobulk
DESeq2** (the correct replication unit). *Example:* Ileum-vs-Duodenum RCTD-proportion and
DEG comparisons.

**Q: What does `06_pseudobulk_de` do with so few replicates?**
It sums counts per replicate (whole tissue and per cluster) and fits DESeq2 with
`de$covariates` (e.g. donor, for a paired design). Covariates the replicates can't
support are dropped with a log message, and subsets with too few replicates are skipped.
Every gene also gets `consistent_across_strata`: whether the direction agrees within
each level of the first covariate (e.g. every donor) — report genes that are significant
**and** consistent. With 1–2 replicates per condition, treat results as hypotheses.

## Batch Correction

**Q: When does the pipeline run Harmony?**
By default (`integration: auto`) whenever the sample sheet has more than one donor — or,
without a `donor` column, more than one sample. Multi-donor Visium HD data clusters by
donor otherwise. *Example (IBD colon, 2 donors):* without correction, 8 of 19 clusters at
8µm were >95% from a single donor. Use `integration_by` to correct for another column
(e.g. batch), or `integration: none` to switch it off. `03_cluster` writes
`donor_mixing_per_cluster.csv` so you can check the mixing after correction.

## Seurat v5 Gotchas

**Q: `AverageExpression` / assay operations behave oddly in Seurat v5.**
v5 returns a **sparse matrix** from `AverageExpression` → wrap in `as.matrix()`. Set
`DefaultAssay()` on the *object you pass* (not a copy). Call `JoinLayers()` before pulling
`counts` from a merged object, or you'll silently get one layer.

## Visualization Conventions

**Q: My composition bars order cell types randomly and colors drift between plots.**
ggplot orders a character fill **alphabetically** by default. Define an **explicit,
ordered palette** (a named color vector) and apply it as **factor levels** to *both* the
fill variable and the dominant-cell-type column so every plot matches. *Example:* the 24
broad types are ordered by compartment (Epithelia → Immune → Stroma, crypt→villus within),
and that same `broad_levels` vector drives legends, stacked bars, and dominant-type maps.

## Bin Size & Segmented Cells

**Q: Do 8µm, 16µm and segmented cells need different handling?**
Yes — the pipeline runs all three (`bin_sizes: [8, 16, cells]`), each in its own
`output/<res>/` folder and job chain. A 16µm bin pools four 8µm bins, so it carries
roughly 4× the UMIs: set QC per resolution (`qc: min_counts: {"8": 50, "16": 100, cells: 50}`)
rather than reusing the 8µm value. 2µm is not supported — it is far heavier and too
sparse for this workflow. *Example (IBD colon):* mean UMIs per tissue bin were ~165 at
8µm vs ~624 at 16µm in the same capture area.

**Q: Can I trust segmented-cell composition across conditions?**
Check first. Space Ranger 4 segments nuclei on the H&E and expands them; tissue where
nuclei are hard to call gets few cells even when the bins there have plenty of UMIs.
*Example (IBD colon):* non-inflamed strips had 2.5–3.4k cells vs 8–19k on inflamed
strips of similar area, with crypt epithelium largely unsegmented. Compare the
`02_load_qc` spatial plots for cells against the 8µm bins before reading anything into
cell-level composition differences.

## QC

**Q: What QC thresholds should I use?**
Defaults: **UMI ≥ 50** per 8µm bin (lab standard from past analyses) and **MT ≤ 50%**.
`02_load_qc` reports what each filter removes per replicate and warns when conditions
lose very different fractions (that biases composition comparisons).

**Q: Is a mitochondrial cutoff appropriate for Visium HD? For gut tissue?**
Not the usual 10–20% from single-cell work. In scRNA-seq, high MT% flags cells damaged
during dissociation; tissue sections aren't dissociated and bins aren't cells, and
colonocytes/enterocytes are mitochondria-rich. High MT% in bins is mostly a low-UMI
artefact (2 of 8 reads being MT = 25%) plus damaged tissue edges. **50%** only removes
near-pure-MT bins, which carry little else. *Example (IBD colon):* above UMI 50, a 20%
cutoff would have removed 5% of one non-inflamed strip (its torn edges) and almost
nothing elsewhere; 50% removed 54 bins in total.

**Q: Should I run ambient RNA removal (SoupX, CellBender, DecontX)?**
No. Those model free mRNA in a dissociated cell suspension, estimated from empty
droplets — there is no suspension in a tissue section. The spatial analogue is lateral
**bleed**: abundant transcripts drift a few µm into neighbouring bins and onto bare
glass. It is small per bin and local, so a global "soup" correction doesn't fit, and
spot-swapping tools (e.g. SpotClean) were built for 55µm Visium, not HD. *Example
(IBD colon):* off-tissue bins averaged 0.25–0.84 UMIs (vs ~175 on tissue) and the top
off-tissue genes were the tissue's own top genes (IGKC, PIGR, FABP1, MT-). Practical
consequence: a little IGKC/IGHG1 in epithelial clusters is expected bleed — don't call
plasma-cell regions from immunoglobulins alone.

**Q: Space Ranger marked bare glass as tissue. Does that matter?**
It can. Its tissue detection calls scattered specks (and debris, folds, fibres) as
"in tissue". They have ~0 UMIs, so the UMI filter removes them — but they also chain
separate tissue pieces together, which is why `01_assign_replicates` finds pieces from
*densely* filled blocks only. Always check the overlay in `output/01_assign_replicates/figures/`.

## Multiple Replicates per Capture Area

**Q: My capture area holds several tissue pieces from different conditions. How do I split them?**
List one sample-sheet row per piece (same `library_id`), in the order they sit on the
slide, and set `replicates: layout: left_to_right` (or `top_to_bottom`). Step 01 detects
the pieces, orders them along that axis and labels 8µm bins, 16µm bins and segmented
cells. Torn pieces are bridged if the gap is under ~200 µm. If pieces touch or the
order is irregular, export labels from Loupe Browser and use `layout: manual`.
*Gotcha:* Visium HD `array_row`/`array_col` are rotated relative to the H&E, so
"left" must come from image pixel coordinates — step 01 does this.

## Plots

**Q: Why use the hires image, and why are my spatial PDFs so big?**
Seurat's loader defaults to the lowres image; the pipeline always loads
`tissue_hires_image.png`. Large PDFs come from tens of thousands of vector points plus
the embedded full-size image: the pipeline rasterises the bins/cells (ggrastr) and
downsamples the image to ≤2500 px. *Example (IBD colon):* spatial QC PDFs drawn with
vector points over the full hires image were 26–103 MB; rasterised with the downsampled
image, a two-library PDF is ~4 MB.

**Q: `clustree` fails with "Unknown guide: edge_colourbar".**
ggraph's edge guides aren't recognised by ggplot2 ≥ 4.0. `03_cluster` hides those
legends (`ggraph::scale_edge_colour_gradientn(guide = "none")`).
