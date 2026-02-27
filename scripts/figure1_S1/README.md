# Figure 1C–D — CRISPR Screen

## Panels

| Panel | Description |
|-------|-------------|
| Fig 1C | β-score rank plot (GDT201 1:1 vs NC 1:1) |
| Fig 1D | Top positively / negatively selected gene tables |

## Input data

Place input in **`demo/figure1_S1/`**.

**MAGeCK MLE gene summary** — `demo/figure1_S1/gene_summary_demo.txt` (included), or the pipeline output at `crispr_analysis/output/combined_res_all_replicates/gadeta_mle_all_reps.gene_summary.txt`. TSV with columns including:

```
Gene      sgRNA  ctrl_05_1|beta  gdt_05_1|beta  ctrl_1_1|beta  gdt_1_1|beta  ...
SEMA4D    4      0.43176         0.28362        0.27451        0.56219       ...
```

`diff = gdt_1_1|beta − ctrl_1_1|beta`. The script uses `ReadBeta()` (or MAGeCKFlute) to keep `Gene` and all `|beta` columns.

## Demo data

**`demo/figure1_S1/gene_summary_demo.txt`** — included demo; same column layout, fewer genes. Used when the full pipeline output is not available.

## How to run

```bash
cd scripts/figure1_S1/
Rscript figure1_S1_crispr_screen.R
```

## Dependencies

R packages: `tidyverse`, `ggrepel`, `MAGeCKFlute` (optional).

### Environment / package versions

Tested with **R 4.5.2**. Example versions:

| Package     | Version |
|------------|---------|
| tidyverse  | 2.0.0   |
| ggrepel    | 0.9.7   |
| MAGeCKFlute| (optional) |

To capture your R environment: run `sessionInfo()` and `packageVersion("pkgname")` in R, or use `renv::snapshot()`.
