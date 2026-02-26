#!/usr/bin/env Rscript
# =============================================================================
# Figure 3D–G & Supplementary Figure 3C–D — scRNA-seq COAD Analysis
#
#   Fig 3D  — UMAP by tissue type
#   Fig 3E  — Expression bins (FUT4, MGAT5, PROCR) by tissue type
#   Fig 3F  — UMAP by cell type
#   Fig 3G  — Expression bins (FUT4, MGAT5, PROCR) by cell type
#             (tumor vs normal cells)
#   Sup 3C  — Expression bins for all CRISPR hit + FUT genes by tissue type
#   Sup 3D  — Expression bins (FUT4, MGAT5, PROCR) by CMS, MSI vs MSS,
#             and iCMS subtype
#
# Requires ~60 GB RAM for the Seurat object.
#
# --- Input data ---------------------------------------------------------------
#
# Seurat RDS: demo/figure3_S3/coad_merged_harmony.rds
#   If missing, script builds it from Synapse RDS files in demo/figure3_S3/
#   (Synapse syn26844071, Joanito et al. Nat Genet 2022).
#
#   The Synapse objects were pre-processed by Joanito et al.:
#     - CellRanger v3.1, QC, DoubletFinder v2.0.3, RCA2 v1.2 cell-type
#       annotation, CMS / iCMS / MSI classification
#   Our processing:
#     1. Load 5 Synapse Seurat v3 objects
#     2. Merge compartments
#     3. NormalizeData → FindVariableFeatures → ScaleData → PCA
#     4. Harmony (v0.1.1) batch correction by sample_type
#     5. UMAP + FindClusters
#
#   Key metadata columns:
#     sample_type  — normal_epithelial, normal_stroma, tumor_epithelial_1/2,
#                    tumor_stroma
#     cell.type    — RCA2-assigned cell type labels
#     CMS          — Consensus Molecular Subtype
#     iCMS         — intrinsic CMS (iCMS2 / iCMS3)
#     msi          — microsatellite instability status (MSI / MSS)
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(ggpubr)
  library(future)
})

Sys.setenv("R_MAX_VSIZE" = 60 * 1024^3)
plan("multicore", workers = 16)
options(future.globals.maxSize = 60 * 1024^3)
options(future.rng.onMisuse = "ignore")
setDTthreads(threads = 16)

# =============================================================================
# Build coad_merged_harmony.rds  (skip if the file already exists)
#
# Place the five Synapse RDS files in demo/figure3_S3/, then run this script.
# It merges the 5 compartments and applies Harmony batch correction.
# =============================================================================

data_dir   <- "../../demo/figure3_S3"
harmony_rds <- file.path(data_dir, "coad_merged_harmony.rds")

if (!file.exists(harmony_rds)) {

  synapse_files <- list(
    normal_epithelial  = file.path(data_dir, "syn26844071_normal_epithelial_Seuratv3_SC00483.rds"),
    normal_stroma      = file.path(data_dir, "syn26844071_normal_stroma_Seuratv3_SC00486.rds"),
    tumor_epithelial_1 = file.path(data_dir, "syn26844071_tumour_epithelial_Seuratv3_SC00481.rds"),
    tumor_epithelial_2 = file.path(data_dir, "syn26844071_tumour_epithelial_Seuratv3_SC00482.rds"),
    tumor_stroma       = file.path(data_dir, "syn26844071_tumor_stroma_Seuratv3_SC00485.rds")
  )

  if (!all(file.exists(unlist(synapse_files)))) {
    stop(
      "coad_merged_harmony.rds not found and not all Synapse RDS files present in ", data_dir, ".\n",
      "Download the five Seurat v3 objects from Synapse syn26844071 (Joanito et al. Nat Genet 2022),\n",
      "place them in demo/figure3_S3/, then re-run this script. See scripts/figure3_S3/README.md."
    )
  }

  # Step 1 — Load each compartment and tag with sample_type
  seurat_list <- lapply(names(synapse_files), function(nm) {
    obj <- readRDS(synapse_files[[nm]])
    DefaultAssay(obj) <- "RNA"
    obj$sample_type <- nm
    obj$orig.ident  <- nm
    obj
  })
  names(seurat_list) <- names(synapse_files)

  # Step 2 — Merge
  merged <- merge(seurat_list[[1]],
                  y = seurat_list[2:length(seurat_list)],
                  add.cell.ids = names(seurat_list))
  rm(seurat_list); gc()

  # Step 3 — Normalize → Variable features → Scale → PCA
  merged <- NormalizeData(merged, verbose = FALSE)
  merged <- FindVariableFeatures(merged, selection.method = "vst",
                                  nfeatures = 2000, verbose = FALSE)
  merged <- ScaleData(merged, verbose = FALSE)
  merged <- RunPCA(merged, npcs = 50, verbose = FALSE)

  # Step 4 — Harmony batch correction by sample_type
  merged <- RunHarmony(merged, group.by.vars = "sample_type",
                       dims.use = 1:50, max.iter.harmony = 20,
                       verbose = TRUE)

  # Step 5 — UMAP + clustering
  merged <- RunUMAP(merged, reduction = "harmony", dims = 1:50, verbose = FALSE)
  merged <- FindNeighbors(merged, reduction = "harmony", dims = 1:50, verbose = FALSE)
  merged <- FindClusters(merged, resolution = 0.5, verbose = FALSE)

  # Save
  dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(merged, harmony_rds)
  rm(merged); gc()
}

# -- Theme + colours ----------------------------------------------------------

theme_publication <- function(base_size = 7) {
  theme_bw(base_size = base_size) +
    theme(
      panel.background  = element_rect(fill = "white", colour = NA),
      plot.background   = element_rect(fill = "white", colour = NA),
      panel.grid.major  = element_line(colour = "grey92", linewidth = 0.25),
      panel.grid.minor  = element_blank(),
      panel.border      = element_rect(colour = "black", fill = NA, linewidth = 0.5),
      axis.line         = element_line(colour = "black", linewidth = 0.5),
      axis.ticks        = element_line(colour = "black", linewidth = 0.5),
      axis.text         = element_text(size = base_size, colour = "black"),
      axis.title        = element_text(size = base_size + 1, colour = "black", face = "bold"),
      legend.background = element_rect(fill = "white", colour = NA),
      legend.key        = element_rect(fill = "white", colour = NA),
      legend.key.size   = unit(3, "mm"),
      legend.text       = element_text(size = base_size - 0.5),
      legend.title      = element_blank(),
      plot.title        = element_text(size = base_size + 1, face = "bold", hjust = 0),
      strip.background  = element_rect(fill = "grey90", colour = "black", linewidth = 0.5),
      strip.text        = element_text(size = base_size, face = "bold", colour = "black")
    )
}

colors_tissue <- c(
  "Normal Epithelial" = "#31a354",
  "Normal Stroma"     = "#a1d99b",
  "Tumor Epithelial"  = "#de2d26",
  "Tumor Stroma"      = "#fc9272"
)

colors_red <- c(
  "Zero"              = "grey90",
  "Low (0-0.5)"       = "#fee5d9",
  "Medium (0.5-1.5)"  = "#fc9272",
  "High (>1.5)"       = "#de2d26"
)

# -- Expression bins: stacked bar of expression levels per group --------------

make_expression_bins <- function(obj, genes, group_col) {
  if (!group_col %in% colnames(obj@meta.data)) return(NULL)
  genes <- intersect(genes, rownames(obj))
  if (length(genes) == 0) return(NULL)

  data_list <- lapply(genes, function(gene) {
    expr <- FetchData(obj, vars = c(gene, group_col))
    colnames(expr) <- c("expression", "group")
    expr <- expr[!is.na(expr$group), , drop = FALSE]
    expr$gene <- gene
    expr$group <- factor(expr$group)
    expr$expr_bin <- cut(expr$expression, breaks = c(-Inf, 0, 0.5, 1.5, Inf),
                         labels = c("Zero", "Low (0-0.5)", "Medium (0.5-1.5)", "High (>1.5)"))
    expr
  })

  plot_data <- bind_rows(data_list) %>%
    group_by(gene, group, expr_bin, .drop = FALSE) %>%
    summarise(count = n(), .groups = "drop_last") %>%
    mutate(pct = 100 * count / sum(count)) %>%
    ungroup()

  ggplot(plot_data, aes(x = group, y = pct, fill = expr_bin)) +
    geom_bar(stat = "identity", position = "stack", width = 0.7) +
    facet_wrap(~ gene, ncol = 4) +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          legend.position = "bottom", legend.key.size = unit(2.5, "mm"),
          strip.text = element_text(size = 9, face = "bold")) +
    labs(title = "", x = "", y = "% of cells", fill = "")
}

# -- Load Seurat object -------------------------------------------------------
merged <- readRDS(harmony_rds)
DefaultAssay(merged) <- "RNA"

merged$tissue_group <- dplyr::case_when(
  merged$sample_type == "normal_epithelial" ~ "Normal Epithelial",
  merged$sample_type == "normal_stroma"     ~ "Normal Stroma",
  merged$sample_type %in% c("tumor_epithelial_1", "tumor_epithelial_2") ~ "Tumor Epithelial",
  merged$sample_type == "tumor_stroma"      ~ "Tumor Stroma",
  TRUE ~ NA_character_
)
merged$tissue_group <- factor(merged$tissue_group, levels = names(colors_tissue))
merged$tumor_normal <- ifelse(grepl("normal", merged$sample_type, ignore.case = TRUE),
                              "Normal", "Tumor")

merged_clean <- subset(merged, cells = colnames(merged)[!is.na(merged$tissue_group)])
set.seed(42)
merged_clean <- merged_clean[, sample(colnames(merged_clean))]
gc()

# Gene sets (matching coad_analysis.rmd)
genes_sparse <- c("FUT3", "FUT4", "FUT5", "FUT6", "FUT9", "GMDS", "TSTA3",
                   "MGAT5", "PROCR", "SLC35C1", "TM9SF1", "TM9SF3", "MAN2A1",
                   "MGAT1", "MGAT2", "MAN1A2")
genes_subset <- c("PROCR", "FUT4", "MGAT5")

genes_sparse_present <- intersect(genes_sparse, rownames(merged_clean))
genes_subset_present <- intersect(genes_subset, rownames(merged_clean))

# =============================================================================
# Figure 3D — UMAP coloured by tissue type
# =============================================================================
DimPlot(merged_clean, reduction = "umap", group.by = "tissue_group",
        cols = colors_tissue, pt.size = 0.1, shuffle = TRUE) +
  theme_publication() + labs(title = "", x = "UMAP 1", y = "UMAP 2")

# =============================================================================
# Figure 3E — Expression bins (FUT4, MGAT5, PROCR) by tissue type
# =============================================================================
make_expression_bins(merged_clean, genes_subset_present, "tissue_group") +
  scale_fill_manual(values = colors_red) +
  coord_cartesian(ylim = c(0, 62))

# =============================================================================
# Figure 3F — UMAP coloured by cell type
# =============================================================================
DimPlot(merged, reduction = "umap", group.by = "cell.type",
        pt.size = 0.1, label = TRUE, repel = TRUE) +
  theme_publication() + labs(title = "", x = "UMAP 1", y = "UMAP 2")

# =============================================================================
# Figure 3G — Expression bins (FUT4, MGAT5, PROCR) by cell type
#             Split into tumor and normal cells
# =============================================================================
merged_tumor  <- subset(merged_clean, tumor_normal == "Tumor")
merged_normal <- subset(merged_clean, tumor_normal == "Normal")

p_tumor <- make_expression_bins(merged_tumor, genes_subset_present, "cell.type") +
  scale_fill_manual(values = colors_red) +
  labs(title = "Tumor")

p_normal <- make_expression_bins(merged_normal, genes_subset_present, "cell.type") +
  scale_fill_manual(values = colors_red) +
  labs(title = "Normal")

p_tumor / p_normal + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

# =============================================================================
# Supplementary Figure 3C — Expression bins for all CRISPR hit + FUT genes
#                            by tissue type
# =============================================================================
make_expression_bins(merged_clean, genes_sparse_present, "tissue_group") +
  scale_fill_manual(values = colors_red)

# =============================================================================
# Supplementary Figure 3D — Expression bins (FUT4, MGAT5, PROCR) by
#                            CMS, MSI vs MSS, and iCMS subtype
# =============================================================================

# CMS subtype
if ("CMS" %in% colnames(merged_clean@meta.data)) {
  make_expression_bins(merged_clean, genes_subset_present, "CMS") +
    scale_fill_manual(values = colors_red) +
    labs(title = "CMS subtype")
}

# MSI vs MSS
if ("msi" %in% colnames(merged_clean@meta.data)) {
  make_expression_bins(merged_clean, genes_subset_present, "msi") +
    scale_fill_manual(values = colors_red) +
    labs(title = "MSI status")
}

# iCMS subtype
if ("iCMS" %in% colnames(merged_clean@meta.data)) {
  make_expression_bins(merged_clean, genes_subset_present, "iCMS") +
    scale_fill_manual(values = colors_red) +
    labs(title = "iCMS subtype")
}
