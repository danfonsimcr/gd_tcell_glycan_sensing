#!/usr/bin/env Rscript
# =============================================================================
# Supplementary Figure 3A–B — Bulk RNA-seq (TCGA + GTEx)
#
#   Sup 3A  — mRNA expression [Log2(TPM+0.1)] heatmap for TCGA + GTEx
#             across all CRISPR hit genes
#   Sup 3B  — ssGSEA score distribution (PROCR, MGAT5, FUT4) via GSVA
#
# --- Input data ---------------------------------------------------------------
#
# Script looks in demo/figure3_S3/ (RNAseq_TCGA_GTEx.csv or RNAseq_TCGA_GTEx_demo.csv).
# Per-sample, per-gene TPM (TCGA + GTEx). Columns: gene_name, tpm, tpm_log2,
#   omicsoft_land, primary_indication, sample_id, gross_anatomical_region.
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(pheatmap)
  library(tibble)
  library(scales)
  library(GSVA)
})

# -- Load full expression data -------------------------------------------------

data_dir   <- "../../demo/figure3_S3"
bulk_file  <- file.path(data_dir, "RNAseq_TCGA_GTEx.csv")
bulk_demo  <- file.path(data_dir, "RNAseq_TCGA_GTEx_demo.csv")
if (file.exists(bulk_file)) {
  bulk_rna_all <- read.csv(bulk_file)
} else if (file.exists(bulk_demo)) {
  message("Using demo data: ", bulk_demo)
  bulk_rna_all <- read.csv(bulk_demo)
} else {
  stop("No bulk RNA file found in ", data_dir, " (expect RNAseq_TCGA_GTEx.csv or RNAseq_TCGA_GTEx_demo.csv)")
}

crispr_hit_genes <- c("BTNL3", "BTNL8", "FUT3", "FUT4", "FUT5", "FUT9",
                       "GMDS", "MAN1A1", "MAN1A2", "MGAT1", "MGAT2",
                       "MGAT5", "PROCR", "SLC35C1", "TM9SF1", "TM9SF3")

# Subset for heatmap only — ssGSEA needs the full transcriptome
bulk_rna_data <- bulk_rna_all %>% filter(gene_name %in% crispr_hit_genes)

anatomical_order <- c(
  "Brain cancer", "CNS", "Head and Neck cancer", "Head and Neck", "PNS",
  "Thyroid cancer", "Lung cancer", "Trachea", "Heart", "Vascular system",
  "Breast cancer", "Esophageal cancer", "Gastroesophageal", "Gastric cancer",
  "Stomach cancer", "Small intestine", "Colorectal cancer", "Large intestine",
  "Liver cancer", "Pancreatic cancer", "Biliary cancer", "Renal cancer",
  "Urothelial cancer", "Urinary tract", "Prostate cancer", "Testicular cancer",
  "Ovarian cancer", "Uterine cancer", "Cervical cancer", "Reproductive female",
  "Reproductive male", "Melanoma cancer", "Skin", "Soft/Connective Tissues",
  "Sarcoma cancer", "Hemo-lymphocytic", "Leukemia", "Lymphoma", "Blood cancer",
  "Endocrine cancer", "Adrenal gland cancer", "Mesothelioma", "Thymic cancer",
  "Thymus cancer"
)

# =============================================================================
# Supplementary Figure 3A — mRNA expression heatmap [Log2(TPM+0.1)]
# TCGA (left) + GTEx (right), anatomically ordered, all CRISPR hit genes
# =============================================================================

tcga_data <- bulk_rna_data %>%
  filter(omicsoft_land == "TCGA", !is.na(primary_indication)) %>%
  group_by(gene_name, primary_indication) %>%
  summarise(mean_tpm_log2 = mean(tpm_log2, na.rm = TRUE), .groups = "drop") %>%
  mutate(region = primary_indication, dataset = "TCGA")

gtex_data <- bulk_rna_data %>%
  filter(omicsoft_land == "GTEX", !is.na(gross_anatomical_region)) %>%
  group_by(gene_name, gross_anatomical_region) %>%
  summarise(mean_tpm_log2 = mean(tpm_log2, na.rm = TRUE), .groups = "drop") %>%
  mutate(region = gross_anatomical_region, dataset = "GTEX")

combined_data <- bind_rows(
  tcga_data %>% select(gene_name, region, mean_tpm_log2, dataset),
  gtex_data %>% select(gene_name, region, mean_tpm_log2, dataset)
)

tcga_regions <- unique(tcga_data$region)
gtex_regions <- unique(gtex_data$region)
tcga_final_order <- c(anatomical_order[anatomical_order %in% tcga_regions],
                      sort(setdiff(tcga_regions, anatomical_order)))
gtex_final_order <- c(anatomical_order[anatomical_order %in% gtex_regions],
                      sort(setdiff(gtex_regions, anatomical_order)))
final_order <- c(tcga_final_order, gtex_final_order)

combined_matrix <- combined_data %>%
  select(gene_name, region, mean_tpm_log2) %>%
  pivot_wider(names_from = region, values_from = mean_tpm_log2, values_fill = 0) %>%
  column_to_rownames("gene_name") %>%
  as.matrix()
combined_matrix <- combined_matrix[, final_order[final_order %in% colnames(combined_matrix)]]

gap_col <- matrix(NA, nrow = nrow(combined_matrix), ncol = 1)
colnames(gap_col) <- " "
rownames(gap_col) <- rownames(combined_matrix)
tcga_count <- length(tcga_final_order)
combined_matrix_with_gap <- cbind(
  combined_matrix[, 1:tcga_count], gap_col,
  combined_matrix[, (tcga_count + 1):ncol(combined_matrix)]
)

annotation_col <- data.frame(
  Dataset = sapply(colnames(combined_matrix_with_gap), function(x) {
    if (x == " ") NA
    else unique(combined_data$dataset[combined_data$region == x])
  }),
  row.names = colnames(combined_matrix_with_gap)
)
ann_colors <- list(Dataset = c(TCGA = "#003B7F", GTEX = "#FDBD13"))

pheatmap(combined_matrix_with_gap,
         color = hcl.colors(30, "Emrld"),
         cluster_rows = FALSE, cluster_cols = FALSE, scale = "none",
         main = "mRNA Expression [Log2(TPM+0.1)]: TCGA and GTEx",
         fontsize = 7, fontsize_row = 9, fontsize_col = 7, angle_col = 45,
         border_color = NA, cellwidth = 12, cellheight = 15,
         annotation_col = annotation_col, annotation_colors = ann_colors,
         annotation_legend = TRUE, na_col = "white")

# =============================================================================
# Supplementary Figure 3B — ssGSEA score distribution (PROCR, MGAT5, FUT4)
#
# ssGSEA ranks ALL genes in each sample's expression profile and calculates
# an enrichment score for the gene set against that full ranking.  We must
# pass the complete genes × samples TPM matrix, not just the 3 signature genes.
# =============================================================================

sig_genes <- c("PROCR", "MGAT5", "FUT4")
gene_set  <- list(PROCR_MGAT5_FUT4 = sig_genes)

# Build a full genes × samples TPM matrix (all genes, TCGA only)
tcga_tpm <- bulk_rna_all %>%
  filter(omicsoft_land == "TCGA", !is.na(sample_id), !is.na(tpm)) %>%
  select(gene_name, sample_id, tpm) %>%
  pivot_wider(names_from = sample_id, values_from = tpm,
              values_fn = mean, values_fill = 0) %>%
  column_to_rownames("gene_name") %>%
  as.matrix()

ssgsea_scores <- gsva(
  tcga_tpm,
  gene_set,
  method = "ssgsea",
  kcdf   = "Gaussian",
  verbose = FALSE
)

# Map sample_id → primary_indication
sample_meta <- bulk_rna_all %>%
  filter(omicsoft_land == "TCGA", !is.na(primary_indication)) %>%
  select(sample_id, primary_indication) %>%
  distinct()

ssgsea_df <- data.frame(
  sample_id    = colnames(ssgsea_scores),
  ssgsea_score = as.numeric(ssgsea_scores["PROCR_MGAT5_FUT4", ])
) %>%
  inner_join(sample_meta, by = "sample_id")

ssgsea_stats <- ssgsea_df %>%
  group_by(primary_indication) %>%
  summarise(
    median_score = median(ssgsea_score, na.rm = TRUE),
    sd_score     = sd(ssgsea_score, na.rm = TRUE),
    n = n(), .groups = "drop"
  ) %>%
  arrange(desc(median_score)) %>%
  mutate(primary_indication = factor(primary_indication, levels = primary_indication))

ggplot(ssgsea_stats, aes(x = primary_indication, y = median_score)) +
  geom_col(fill = "#de2d26", alpha = 0.8) +
  geom_errorbar(aes(ymin = median_score - sd_score / sqrt(n),
                     ymax = median_score + sd_score / sqrt(n)),
                width = 0.4, color = "black", linewidth = 0.5) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 11),
        axis.title  = element_text(size = 12, face = "bold"),
        plot.title  = element_text(size = 14, face = "bold", hjust = 0.5)) +
  labs(title = "ssGSEA Enrichment Score (PROCR, MGAT5, FUT4)",
       subtitle = "TCGA by indication — calculated via GSVA",
       x = "Primary Indication", y = "Median ssGSEA Score (normalised)")
