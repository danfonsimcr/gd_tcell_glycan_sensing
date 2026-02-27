# Publication Figure Scripts

Annotated, self-contained scripts that reproduce all computationally-generated
figures in the manuscript. Each subfolder maps to a main figure and its
supplementary panels.

## Structure

```
scripts/
├── figure1_S1/                        CRISPR screen (R)
│   ├── figure1_S1_crispr_screen.R     Fig 1C–D
│   └── README.md
├── figure3_S3/                        scRNA-seq + bulk RNA-seq (R)
│   ├── figure3_S3_scrna_coad.R        Fig 3D–G, Sup Fig 3C–D
│   ├── figure3D_S3D_bulk_rna.R        Sup Fig 3A–B
│   └── README.md
└── figure4_S4/                        ML classifier (Python) + DESeq2 (R)
    ├── figure4_ml_classifier.py       Fig 4A–D
    ├── supfig4_AB_deseq2.R            Sup Fig 4A–B
    ├── supfig4_CD_ml_features.py      Sup Fig 4C–D
    └── README.md
```

## Panel mapping

| Figure | Panel | Description | Script |
|--------|-------|-------------|--------|
| 1 | C | CRISPR β-score rank plot | `figure1_S1/figure1_S1_crispr_screen.R` |
| 1 | D | Top CRISPR hit genes table | `figure1_S1/figure1_S1_crispr_screen.R` |
| 3 | D | UMAP by tissue type | `figure3_S3/figure3_S3_scrna_coad.R` |
| 3 | E | Expression bins (FUT4, MGAT5, PROCR) by tissue type | `figure3_S3/figure3_S3_scrna_coad.R` |
| 3 | F | UMAP by cell type | `figure3_S3/figure3_S3_scrna_coad.R` |
| 3 | G | Expression bins by cell type (tumor / normal) | `figure3_S3/figure3_S3_scrna_coad.R` |
| S3 | A | mRNA expression heatmap (TCGA + GTEx, CRISPR hit genes) | `figure3_S3/figure3D_S3D_bulk_rna.R` |
| S3 | B | Composite Z-score violin (PROCR, MGAT5, FUT4 signature) | `figure3_S3/figure3D_S3D_bulk_rna.R` |
| S3 | C | Expression bins — all CRISPR hit + FUT genes by tissue type | `figure3_S3/figure3_S3_scrna_coad.R` |
| S3 | D | Expression bins (FUT4, MGAT5, PROCR) by CMS, MSI vs MSS, iCMS | `figure3_S3/figure3_S3_scrna_coad.R` |
| 4 | A | CCLE responder histogram | `figure4_S4/figure4_ml_classifier.py` |
| 4 | B | CCLE ≥6 models: count by indication | `figure4_S4/figure4_ml_classifier.py` |
| 4 | C | Normal cell line histogram | `figure4_S4/figure4_ml_classifier.py` |
| 4 | D | Normal ≥3 models: count by sub_anatomical_region1 | `figure4_S4/figure4_ml_classifier.py` |
| S4 | A | Compiled volcano plot (DGE) | `figure4_S4/supfig4_AB_deseq2.R` |
| S4 | B | Top 20 DEG heatmap | `figure4_S4/supfig4_AB_deseq2.R` |
| S4 | C | ML feature importance (normalised 0–1) | `figure4_S4/supfig4_CD_ml_features.py` |
| S4 | D | Sum norm features vs CRISPR β-score | `figure4_S4/supfig4_CD_ml_features.py` |

## How to run

```bash
# R scripts
Rscript figure1_S1/figure1_S1_crispr_screen.R

# Python
conda activate gdt201
python figure4_S4/figure4_ml_classifier.py
python figure4_S4/supfig4_CD_ml_features.py
```

Each script header contains the input data description and a tidy head representation of the expected data format. See the per-folder README for detailed input file tables and **environment / package versions** (R `sessionInfo()` and Python `conda list` / `pip list`).

## Data (demo/)

All input data lives in **`demo/`**. Use the included demo files or put your own data in the same folder structure. See the main **README.md** at the repo root for the full list of expected files per figure.

| Figure | Data location |
|--------|----------------|
| 1 | `demo/figure1_S1/` (*gene_summary*.txt or gene_summary_demo.txt) |
| 3 | `demo/figure3_S3/` (bulk CSV; scRNA: build RDS from Synapse files, see figure3_S3/README.md) |
| 4 | `demo/figure4_S4/` (included *_demo.csv files or your own) |
