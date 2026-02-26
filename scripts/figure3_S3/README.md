# Figure 3D–G & Supplementary Figure 3A–D — scRNA-seq + Bulk RNA

## Scripts

| Script | Panels |
|--------|--------|
| `figure3_S3_scrna_coad.R` | Fig 3D–G, Sup Fig 3C–D |
| `figure3D_S3D_bulk_rna.R` | Sup Fig 3A–B |

## Panels

| Panel | Description | Script |
|-------|-------------|--------|
| Fig 3D | UMAP by tissue type | `figure3_S3_scrna_coad.R` |
| Fig 3E | Expression bins (FUT4, MGAT5, PROCR) by tissue type | `figure3_S3_scrna_coad.R` |
| Fig 3F | UMAP by cell type | `figure3_S3_scrna_coad.R` |
| Fig 3G | Expression bins (FUT4, MGAT5, PROCR) by cell type (tumor / normal) | `figure3_S3_scrna_coad.R` |
| Sup 3A | mRNA expression heatmap [Log2(TPM+0.1)] — TCGA + GTEx, all CRISPR hit genes | `figure3D_S3D_bulk_rna.R` |
| Sup 3B | Gene signature ssGSEA score distribution (PROCR, MGAT5, FUT4) | `figure3D_S3D_bulk_rna.R` |
| Sup 3C | Expression bins — all CRISPR hit + FUT genes by tissue type | `figure3_S3_scrna_coad.R` |
| Sup 3D | Expression bins (FUT4, MGAT5, PROCR) by CMS, MSI vs MSS, and iCMS | `figure3_S3_scrna_coad.R` |

## Input data

### scRNA-seq

**Seurat RDS:** `demo/figure3_S3/coad_merged_harmony.rds`

If this file does not exist, the script builds it from Synapse data:

1. **Download** the five pre-processed Seurat v3 objects from Synapse (**syn26844071**, Joanito et al. *Nat Genet* 54, 963–975, 2022).
2. Place the five RDS files in **`demo/figure3_S3/`** with these names:
   - `syn26844071_normal_epithelial_Seuratv3_SC00483.rds`
   - `syn26844071_normal_stroma_Seuratv3_SC00486.rds`
   - `syn26844071_tumour_epithelial_Seuratv3_SC00481.rds`
   - `syn26844071_tumour_epithelial_Seuratv3_SC00482.rds`
   - `syn26844071_tumor_stroma_Seuratv3_SC00485.rds`
3. Run `figure3_S3_scrna_coad.R`. It will merge the compartments, run Harmony batch correction, and save `coad_merged_harmony.rds` in `demo/figure3_S3/`.

Requires ~60 GB RAM for the full object. No demo RDS is provided; use the Synapse pipeline above.

### Bulk RNA (Sup 3A–B)

**`demo/figure3_S3/RNAseq_TCGA_GTEx.csv`** or **`RNAseq_TCGA_GTEx_demo.csv`** — per-sample, per-gene TPM (TCGA + GTEx).

## How to run

```bash
cd scripts/figure3_S3/

# scRNA (requires ~60 GB RAM; build RDS from Synapse files in demo/figure3_S3/ first if needed)
Rscript figure3_S3_scrna_coad.R

# Bulk RNA + ssGSEA
Rscript figure3D_S3D_bulk_rna.R
```

## Dependencies

- **scRNA:** `Seurat` (≥5.0), `harmony`, `ggplot2`, `patchwork`, `dplyr`, `tidyr`, `data.table`, `ggpubr`, `future`
- **Bulk RNA:** `ggplot2`, `dplyr`, `tidyr`, `pheatmap`, `tibble`, `scales`, `GSVA`
