#!/usr/bin/env Rscript
# =============================================================================
# Supplementary Figure 3A–B — Bulk RNA-seq (TCGA + GTEx)
#
#   Sup 3A — mRNA expression [Log2(TPM+0.1)] heatmap (all CRISPR hit genes).
#   Sup 3B — Composite Z-score violin (PROCR, MGAT5, FUT4 signature).
#
# Z-score method:
#   Per sample, each gene is z-scored against the full transcriptome:
#       z_g = (log2(TPM+0.1)_g - mean_all_genes) / sd_all_genes
#   Composite = mean(z_PROCR, z_MGAT5, z_FUT4).
#   Full transcriptome background (~19–20k genes) ensures a stable, unbiased
#   reference (signature = 3/20k ≈ 0.015% of background).
#
# Input:
#   RNAseq_TCGA_GTEx.csv (or _demo.csv) in demo/figure3_S3/.
#   Columns: gene_name, tpm, tpm_log2, omicsoft_land, primary_indication,
#            sample_id, gross_anatomical_region. Optional: sample_category.
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(pheatmap)
  library(tibble)
  library(scales)
})

# -- Paths ---------------------------------------------------------------------
data_dir  <- "../../demo/figure3_S3"
bulk_file <- file.path(data_dir, "RNAseq_TCGA_GTEx.csv")
bulk_demo <- file.path(data_dir, "RNAseq_TCGA_GTEx_demo.csv")
if (file.exists(bulk_file)) {
  bulk_rna_all <- read.csv(bulk_file)
} else if (file.exists(bulk_demo)) {
  message("Using demo data: ", bulk_demo)
  bulk_rna_all <- read.csv(bulk_demo)
} else {
  stop("No bulk RNA file found in ", data_dir,
       " (expect RNAseq_TCGA_GTEx.csv or RNAseq_TCGA_GTEx_demo.csv)")
}

# Optional: restrict TCGA to Primary Tumor when column exists
if ("sample_category" %in% names(bulk_rna_all)) {
  bulk_rna_all <- bulk_rna_all %>%
    filter(
      (omicsoft_land == "TCGA" & sample_category == "Primary Tumor") |
        omicsoft_land == "GTEX"
    )
  message("Filtered TCGA to Primary Tumor; GTEx unchanged.")
}

# -- Gene sets -----------------------------------------------------------------
# 16 genes from the CRISPR screen that are reliably quantified in bulk
# RNA-seq (TCGA/GTEx). Includes butyrophilins (BTNL3/8) and MAN1A1 that
# are present in bulk but sparse in scRNA-seq. Compare with the scRNA
# gene set in figure3_S3_scrna_coad.R, which swaps these for MAN2A1,
# TSTA3, and FUT6 (detectable in the single-cell data).
crispr_hit_genes <- c("BTNL3", "BTNL8", "FUT3", "FUT4", "FUT5", "FUT9",
                      "GMDS", "MAN1A1", "MAN1A2", "MGAT1", "MGAT2",
                      "MGAT5", "PROCR", "SLC35C1", "TM9SF1", "TM9SF3")
sig_genes <- c("PROCR", "MGAT5", "FUT4")

bulk_rna_data <- bulk_rna_all %>% filter(gene_name %in% crispr_hit_genes)

# Anatomical ordering — same as land_rna_sig_v2
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
order_anatomical <- function(regions) {
  known <- anatomical_order[anatomical_order %in% regions]
  c(known, sort(setdiff(regions, known)))
}

# -- Colors --------------------------------------------------------------------
TCGA_COL <- "#003B7F"
GTEX_COL <- "#FDBD13"
DATASET_COLORS <- c(TCGA = TCGA_COL, GTEX = GTEX_COL)

# =============================================================================
# Supplementary Figure 3A — mRNA expression heatmap [Log2(TPM+0.1)]
# TCGA (left) + GTEx (right), anatomically ordered, all 16 CRISPR hit genes
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

tcga_regions <- order_anatomical(unique(tcga_data$region))
gtex_regions <- order_anatomical(unique(gtex_data$region))
final_order <- c(tcga_regions, gtex_regions)

combined_matrix <- combined_data %>%
  select(gene_name, region, mean_tpm_log2) %>%
  pivot_wider(names_from = region, values_from = mean_tpm_log2, values_fill = NA) %>%
  column_to_rownames("gene_name") %>%
  as.matrix()
combined_matrix <- combined_matrix[, final_order[final_order %in% colnames(combined_matrix)]]

# Gap column between TCGA and GTEx
tcga_count <- length(tcga_regions[tcga_regions %in% colnames(combined_matrix)])
n_col <- ncol(combined_matrix)
gap_col <- matrix(NA, nrow = nrow(combined_matrix), ncol = 1)
colnames(gap_col) <- " "
rownames(gap_col) <- rownames(combined_matrix)
if (tcga_count == 0) {
  combined_matrix_with_gap <- cbind(gap_col, combined_matrix)
} else if (tcga_count >= n_col) {
  combined_matrix_with_gap <- cbind(combined_matrix, gap_col)
} else {
  combined_matrix_with_gap <- cbind(
    combined_matrix[, seq_len(tcga_count), drop = FALSE], gap_col,
    combined_matrix[, (tcga_count + 1):n_col, drop = FALSE]
  )
}

annotation_col <- data.frame(
  Dataset = sapply(colnames(combined_matrix_with_gap), function(x) {
    if (x == " ") NA else combined_data$dataset[combined_data$region == x][1]
  }),
  row.names = colnames(combined_matrix_with_gap)
)

pheatmap(combined_matrix_with_gap,
         color = hcl.colors(30, "Emrld"),
         cluster_rows = FALSE, cluster_cols = FALSE, scale = "none",
         main = "mRNA Expression [Log2(TPM+0.1)]: TCGA and GTEx",
         fontsize = 7, fontsize_row = 9, fontsize_col = 7, angle_col = 45,
         border_color = NA, cellwidth = 12, cellheight = 15,
         annotation_col = annotation_col,
         annotation_colors = list(Dataset = DATASET_COLORS),
         annotation_legend = TRUE, na_col = "white")

# =============================================================================
# Supplementary Figure 3B — Composite Z-score (PROCR, MGAT5, FUT4)
# =============================================================================

expr_long <- bulk_rna_all %>%
  filter(!is.na(tpm_log2), !is.na(sample_id)) %>%
  select(gene_name, sample_id, tpm_log2)

expr_mat <- expr_long %>%
  pivot_wider(names_from = sample_id, values_from = tpm_log2,
              values_fn = mean, values_fill = NA) %>%
  column_to_rownames("gene_name") %>%
  as.matrix()

sig_present <- sig_genes[sig_genes %in% rownames(expr_mat)]

sample_mean <- apply(expr_mat, 2L, mean, na.rm = TRUE)
sample_sd   <- apply(expr_mat, 2L, sd,   na.rm = TRUE)
sample_sd[sample_sd == 0 | !is.finite(sample_sd)] <- NA

z_mat <- sweep(expr_mat, 2L, sample_mean, "-")
z_mat <- sweep(z_mat, 2L, sample_sd, "/")
z_composite <- colMeans(z_mat[sig_present, , drop = FALSE], na.rm = TRUE)

sample_meta <- bulk_rna_all %>%
  select(sample_id, omicsoft_land, primary_indication, gross_anatomical_region) %>%
  distinct(sample_id, .keep_all = TRUE)

z_df <- data.frame(
  sample_id = names(z_composite),
  zscore_composite = as.numeric(z_composite)
) %>%
  inner_join(sample_meta, by = "sample_id") %>%
  mutate(
    region  = ifelse(omicsoft_land == "TCGA", primary_indication, gross_anatomical_region),
    dataset = ifelse(omicsoft_land == "TCGA", "TCGA", "GTEX")
  ) %>%
  filter(!is.na(region))

z_stats <- z_df %>%
  group_by(dataset, region) %>%
  summarise(
    median_z = median(zscore_composite, na.rm = TRUE),
    mean_z   = mean(zscore_composite, na.rm = TRUE),
    sd_z     = sd(zscore_composite, na.rm = TRUE),
    n        = n(), .groups = "drop"
  ) %>%
  arrange(desc(mean_z))

top_n <- 25
top_tcga <- z_stats %>% filter(dataset == "TCGA") %>% head(top_n) %>% pull(region)
top_gtex <- z_stats %>% filter(dataset == "GTEX") %>% head(top_n) %>% pull(region)
top_regions <- union(top_tcga, top_gtex)

region_order_z <- z_df %>%
  filter(region %in% top_regions) %>%
  group_by(region) %>%
  summarise(avg = mean(zscore_composite, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(avg)) %>%
  pull(region)

z_plot_data <- z_df %>%
  filter(region %in% top_regions) %>%
  mutate(region = factor(region, levels = region_order_z))

p_violin <- ggplot(z_plot_data,
                   aes(x = region, y = zscore_composite, fill = dataset)) +
  geom_violin(trim = FALSE, alpha = 0.7, scale = "width") +
  geom_boxplot(width = 0.2, position = position_dodge(0.9),
               outlier.shape = NA, alpha = 0.5) +
  stat_summary(fun = median, geom = "point", shape = 23, size = 2,
               position = position_dodge(0.9), fill = "white", color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = DATASET_COLORS, name = "Dataset") +
  scale_x_discrete(drop = FALSE) +
  theme_minimal() +
  theme(
    axis.text.x   = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y   = element_text(size = 10),
    axis.title     = element_text(size = 12, face = "bold"),
    plot.title     = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle  = element_text(size = 11, hjust = 0.5),
    legend.position = "top",
    legend.title   = element_text(face = "bold")
  ) +
  labs(
    title    = "Composite Z-Score — PROCR, MGAT5, FUT4",
    subtitle = paste0("Full transcriptome background (", nrow(expr_mat), " genes) | TCGA + GTEx"),
    x = "Tissue / Cancer Type", y = "Composite Signature Z-Score"
  )
print(p_violin)
