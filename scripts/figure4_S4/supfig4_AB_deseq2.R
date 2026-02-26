#!/usr/bin/env Rscript
# =============================================================================
# Supplementary Figure 4A–B — Differential Gene Expression (DESeq2)
#
#   Sup 4A  — Compiled volcano plot (top 20 genes, fucosylation, EPCR,
#             CRISPR hit genes annotated)
#   Sup 4B  — Top 20 DEG heatmap (responders vs non-responders)
#
# --- Input data ---------------------------------------------------------------
#
# Script looks in demo/figure4_S4/.
#   ccle_dge_counts.csv  — genes (ENSEMBL) × cell lines, integer counts
#   ccle_metadata.csv    — cell_line, condition (non_responder / responder)
#   crispr_hit_genes.csv — optional, gene_name for volcano annotation
# DESeq2: design ~ condition. LFC shrinkage via ashr.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
  library(pheatmap)
  library(annotables)
  library(ggrepel)
  library(patchwork)
})

# -- Load counts and metadata -------------------------------------------------

data_dir   <- "../../demo/figure4_S4"
counts_file <- if (file.exists(file.path(data_dir, "ccle_dge_counts.csv"))) file.path(data_dir, "ccle_dge_counts.csv") else file.path(data_dir, "ccle_dge_counts_demo.csv")
meta_file  <- if (file.exists(file.path(data_dir, "ccle_metadata.csv"))) file.path(data_dir, "ccle_metadata.csv") else file.path(data_dir, "ccle_metadata_demo.csv")
if (!file.exists(counts_file)) stop("Counts not found in ", data_dir)
if (!file.exists(meta_file))  stop("Metadata not found in ", data_dir)
if (grepl("_demo\\.csv", counts_file)) message("Using demo data: ", counts_file)

counts_data <- read.csv(counts_file, row.names = 1, check.names = FALSE)
counts_data <- as.matrix(round(counts_data))

col_data <- read.csv(meta_file, row.names = 1)
col_data$condition <- factor(col_data$condition,
                              levels = c("non_responder", "responder"))

# -- DESeq2 -------------------------------------------------------------------

dds <- DESeqDataSetFromMatrix(countData = counts_data,
                               colData = col_data,
                               design = ~ condition)
dds <- DESeq(dds)

results_ashr <- lfcShrink(dds,
                           coef = "condition_responder_vs_non_responder",
                           type = "ashr")

# -- Annotate ENSEMBL → gene symbols -----------------------------------------

df_res <- as.data.frame(results_ashr) %>%
  mutate(ENSEMBL = rownames(results_ashr), .before = 1) %>%
  inner_join(grch38, by = c("ENSEMBL" = "ensgene")) %>%
  arrange(padj)

df_res <- df_res %>%
  mutate(gene_category = case_when(
    symbol %in% paste0("FUT", 1:11)                        ~ "Fucosyltransferases",
    symbol %in% c("SLC35C1", "GMDS", "TSTA3")              ~ "Fucosylation",
    symbol %in% c("MAN2A1", "MGAT2", "MGAT5")              ~ "Glycosylation",
    symbol %in% c("TM9SF1", "TM9SF2", "TM9SF3", "TM9SF4") ~ "Golgi Regulation",
    symbol == "PROCR"                                        ~ "EPCR",
    TRUE                                                     ~ "Other"
  ))

# -- Heatmap annotation -------------------------------------------------------

ann_df <- as.data.frame(colData(dds)[, "condition", drop = FALSE])
annot_colours <- list(condition = c(non_responder = "#003B7F",
                                     responder = "#FDBD13"))

vst_dds <- vst(dds, blind = FALSE)

# =============================================================================
# Supplementary Figure 4A — Compiled volcano plot
# =============================================================================

volcano_plot <- function(dds_results, gene_labels, graph_subtitle = "") {
  dds_results <- dds_results %>%
    mutate(gene_category = replace_na(gene_category, "Other"),
           padj = pmax(padj, 1e-300, na.rm = TRUE))

  label_df <- dds_results %>%
    semi_join(gene_labels %>% distinct(symbol), by = "symbol")

  cat_palette <- c(
    "Fucosyltransferases" = "#8B2323", "Fucosylation" = "#00008B",
    "Glycosylation" = "#800080", "O-linked fucosylation" = "#FF8C00",
    "Golgi Regulation" = "#228B22", "EPCR" = "#DC143C",
    "Top 20 Significantly expressed" = "#104E8B", "Other" = "#BEBEBE"
  )

  ggplot(dds_results, aes(x = log2FoldChange, y = -log10(padj), colour = gene_category)) +
    geom_point(size = 1.2, alpha = 0.5) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.5, alpha = 0.6) +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", linewidth = 0.5, alpha = 0.6) +
    geom_text_repel(data = label_df, aes(label = symbol), size = 3,
                    box.padding = 0.4, max.overlaps = 25, segment.size = 0.4) +
    scale_colour_manual(values = cat_palette, drop = FALSE, name = "Gene Category") +
    labs(subtitle = graph_subtitle,
         x = expression(log[2] ~ "Fold Change"),
         y = expression(-log[10] ~ "(Adjusted " * italic(P) * "-value)")) +
    theme_classic(base_size = 12) +
    theme(plot.subtitle = element_text(face = "bold", size = 11, hjust = 0.5),
          axis.title = element_text(size = 10, face = "bold"),
          axis.text = element_text(size = 9, color = "black"))
}

df_tbl <- df_res %>%
  filter(!is.na(padj), !is.na(log2FoldChange), !is.na(symbol)) %>%
  mutate(padj = pmax(padj, 1e-300))

# CRISPR hit genes for volcano annotation (optional)
crispr_hit_path <- file.path(data_dir, "crispr_hit_genes.csv")
if (!file.exists(crispr_hit_path)) crispr_hit_path <- file.path(data_dir, "crispr_hit_genes_demo.csv")
if (file.exists(crispr_hit_path)) {
  crispr_hits <- read_csv(crispr_hit_path, show_col_types = FALSE)
  crispr_hit_df <- df_tbl %>% filter(symbol %in% crispr_hits$gene_name)
} else {
  crispr_hit_df <- df_tbl %>% head(0)
}

futx_genes     <- df_tbl %>% filter(symbol %in% paste0("FUT", 1:11))
top20_for_volc <- df_tbl %>% arrange(padj) %>% slice_head(n = 20)
df_top20_cat   <- df_tbl %>%
  mutate(gene_category = if_else(symbol %in% top20_for_volc$symbol,
                                  "Top 20 Significantly expressed", gene_category))
o_fucos <- df_tbl %>% filter(symbol %in% c("POFUT1", "POFUT2"))

p1 <- volcano_plot(df_tbl, crispr_hit_df, "CRISPR Hit Genes")
p2 <- volcano_plot(df_tbl, futx_genes, "Fucosyltransferases (FUT1-11)")
p3 <- volcano_plot(df_top20_cat, top20_for_volc, "Top 20 Genes by Adjusted P-value")
p4 <- volcano_plot(df_tbl, o_fucos, "O-Fucosylation Genes")

(p1 | p2) / (p3 | p4) +
  plot_annotation(
    title = "Differential Gene Expression: Responder vs Non-Responder",
    theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
  ) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

# =============================================================================
# Supplementary Figure 4B — Top 20 DEG heatmap
# =============================================================================

top20 <- df_res %>%
  filter(!is.na(padj), !is.na(symbol)) %>%
  arrange(padj) %>%
  slice_head(n = 20)

top20_assay <- assay(vst_dds)[top20$ENSEMBL, ]
rownames(top20_assay) <- top20$symbol

pheatmap(top20_assay,
         cluster_rows = TRUE, cluster_cols = TRUE,
         annotation_col = ann_df,
         annotation_colors = annot_colours,
         annotation_names_col = FALSE,
         color = hcl.colors(30, "Emrld"),
         fontsize = 10, border_color = "black",
         main = "Top 20 Most Significantly Expressed Genes")
