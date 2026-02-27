# Nature submission — figure reproducibility

Scripts to reproduce the manuscript figures. All input data lives in **`demo/`**: use the included demo files, or put your own data there (same folder structure and formats below).

## Layout

```
.
├── README.md           (this file)
├── demo/               Input data — use included files or add your own here
│   ├── figure1_S1/
│   ├── figure3_S3/
│   └── figure4_S4/
└── scripts/            Figure scripts and panel mapping
    ├── README.md
    ├── figure1_S1/     CRISPR screen (Fig 1C–D)
    ├── figure3_S3/     scRNA-seq + bulk RNA (Fig 3D–G, Sup 3A–D)
    └── figure4_S4/     ML classifier + DESeq2 (Fig 4A–D, Sup 4A–D)
```

## Input data (demo/)

Place your input files in the matching `demo/figureN_SN/` folder. Scripts use the included demo files when your own files are not present. Expected files and formats:

### demo/figure1_S1/

| File | Description |
|------|-------------|
| Any `*gene_summary*.txt` | MAGeCK MLE gene summary (TSV: Gene, sgRNA, ctrl_1_1\|beta, gdt_1_1\|beta, …). If none, script uses `gene_summary_demo.txt`. |

### demo/figure3_S3/

| File | Description |
|------|-------------|
| `coad_merged_harmony.rds` | **OR** build it by placing the five Synapse RDS files below in this folder and running `figure3_S3_scrna_coad.R`. |
| `syn26844071_normal_epithelial_Seuratv3_SC00483.rds` | From Synapse syn26844071 (Joanito et al. *Nat Genet* 2022). |
| `syn26844071_normal_stroma_Seuratv3_SC00486.rds` | |
| `syn26844071_tumour_epithelial_Seuratv3_SC00481.rds` | |
| `syn26844071_tumour_epithelial_Seuratv3_SC00482.rds` | |
| `syn26844071_tumor_stroma_Seuratv3_SC00485.rds` | |
| `RNAseq_TCGA_GTEx.csv` or `RNAseq_TCGA_GTEx_demo.csv` | Bulk RNA TPM (TCGA + GTEx). Columns: gene_name, tpm, tpm_log2, omicsoft_land, primary_indication, sample_id, gross_anatomical_region. |

### demo/figure4_S4/

Figure 4 scripts use the included `*_demo.csv` files by default. To use your own data, add files with the same names as in the table (with or without `_demo` suffix). Scripts look for the non-demo name first, then the demo name.

| File | Description |
|------|-------------|
| `ccle_gadeta_expression.csv` | Cell lines × genes + resp_status |
| `ccle_indications.csv` | cell_line, indication |
| `normal_gadeta_expression.csv` | Normal cell line expression |
| `normal_indications.csv` | cell_line, sub_anatomical_region1 |
| `ccle_dge_counts.csv` | Genes (ENSEMBL) × cell lines, integer counts |
| `ccle_metadata.csv` | cell_line, condition (non_responder / responder) |
| `crispr_hit_genes.csv` | gene_name (optional) |
| `pos_genes_gdt_vs_nc_1_1_mle_combined.csv` | Gene, diff (CRISPR β-score) |

## Quick start

Run from the repo root (this directory):

```bash
# Figure 1
Rscript scripts/figure1_S1/figure1_S1_crispr_screen.R

# Figure 3 bulk RNA
Rscript scripts/figure3_S3/figure3D_S3D_bulk_rna.R

# Figure 3 scRNA (needs demo/figure3_S3/coad_merged_harmony.rds or the 5 Synapse RDS files in demo/figure3_S3/)
Rscript scripts/figure3_S3/figure3_S3_scrna_coad.R

# Figure 4
conda activate gdt201
python scripts/figure4_S4/figure4_ml_classifier.py
Rscript scripts/figure4_S4/supfig4_AB_deseq2.R
python scripts/figure4_S4/supfig4_CD_ml_features.py
```

For panel-to-script mapping, per-figure details, and **R/Python environment and package versions**, see **`scripts/README.md`** and the README in each `scripts/figureN_SN/` folder.

## Figure 3 scRNA RDS

There is no pre-built demo Seurat object. To build `demo/figure3_S3/coad_merged_harmony.rds`:

1. Download the five pre-processed Seurat v3 objects from **Synapse (syn26844071)** — Joanito et al. *Nat Genet* 54, 963–975 (2022).
2. Place the five RDS files in **`demo/figure3_S3/`** (exact filenames in the table above and in `scripts/figure3_S3/README.md`).
3. Run `Rscript scripts/figure3_S3/figure3_S3_scrna_coad.R`. It will merge compartments, run Harmony, and save the RDS in the same folder.
