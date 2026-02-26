# Figure 4A–D & Supplementary Figure 4A–D — ML Classifier + DESeq2

## Scripts

| Script | Language | Panels |
|--------|----------|--------|
| `figure4_ml_classifier.py` | Python | Fig 4A–D |
| `supfig4_AB_deseq2.R` | R | Sup Fig 4A–B |
| `supfig4_CD_ml_features.py` | Python | Sup Fig 4C–D |

## Panels

| Panel | Description | Script |
|-------|-------------|--------|
| Fig 4A | CCLE responder histogram (# models predicting responder) | `figure4_ml_classifier.py` |
| Fig 4B | Among CCLE with ≥6 models: count by indication | `figure4_ml_classifier.py` |
| Fig 4C | Normal cell line responder histogram | `figure4_ml_classifier.py` |
| Fig 4D | Among Normal with ≥3 models: count by sub_anatomical_region1 | `figure4_ml_classifier.py` |
| Sup Fig 4A | Compiled volcano plot — DGE with top 20, fuc, EPCR, CRISPR hit annotations | `supfig4_AB_deseq2.R` |
| Sup Fig 4B | Top 20 DEG heatmap (responders vs non-responders) | `supfig4_AB_deseq2.R` |
| Sup Fig 4C | ML feature importance (normalised 0–1): Adaboost, RF, GB Trees, XGBoost, LogReg | `supfig4_CD_ml_features.py` |
| Sup Fig 4D | Sum of normalised ML feature importances vs CRISPR β-score | `supfig4_CD_ml_features.py` |

## Input data

All inputs: **`demo/figure4_S4/`**. Scripts use the included *_demo.csv files, or your own files with the same names.

| File | Used by | Description |
|------|---------|-------------|
| `ccle_gadeta_expression.csv` | Fig 4, Sup 4C–D | Cell lines × genes + `resp_status` |
| `ccle_indications.csv` | Fig 4 | Column `indication` |
| `normal_gadeta_expression.csv` | Fig 4 | Normal cell line expression |
| `normal_indications.csv` | Fig 4 | Column `sub_anatomical_region1` |
| `ccle_dge_counts.csv` | Sup 4A–B | Genes × cell lines, integer counts |
| `ccle_metadata.csv` | Sup 4A–B | `cell_line`, `condition` |
| `crispr_hit_genes.csv` | Sup 4A–B | `gene_name` (optional) |
| `pos_genes_gdt_vs_nc_1_1_mle_combined.csv` | Sup 4D | `Gene`, `diff` (CRISPR β-score) |

## How to run

```bash
cd scripts/figure4_S4/

# Figure 4 (Python)
conda activate gdt201
python figure4_ml_classifier.py

# Supplementary 4A–B (R)
Rscript supfig4_AB_deseq2.R

# Supplementary 4C–D (Python)
python supfig4_CD_ml_features.py
```

## Dependencies

- **Python:** `scikit-learn`, `xgboost`, `pandas`, `numpy`, `matplotlib`, `seaborn`, `scipy`
- **R:** `DESeq2`, `tidyverse`, `pheatmap`, `annotables`, `EnhancedVolcano`, `ggrepel`, `patchwork`, `ashr`
