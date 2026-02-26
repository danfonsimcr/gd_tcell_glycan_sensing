#!/usr/bin/env Rscript
# =============================================================================
# Figure 1C — CRISPR beta-score rank plot (GDT201 1:1 vs NC 1:1, all reps)
# Figure 1D — Curated table of top positively selected genes with functions
#
# --- Pipeline that produced the input data -----------------------------------
#
# See crispr_analysis/ for the full pipeline:
#   0. TrimGalore on R1 FASTQs
#   1. MAGeCK count (brunello_library.csv, trimmed reads)
#   2. CRISPRcleanR (16 samples, ncontrols=3) -> corrected counts
#   3. MAGeCK MLE (--norm-method none) on CRISPRcleanR-corrected counts
#      with combined_design_matrix_all_replicates.txt
#   4. Gene summary -> beta-scores per condition
#
# --- Input -------------------------------------------------------------------
#
# MAGeCK MLE gene summary (TSV) from the all-replicates analysis:
#   crispr_analysis/output/combined_res_all_replicates/
#       gadeta_mle_all_reps.gene_summary.txt
#
# Falls back to demo/figure1_S1/gene_summary_demo.txt when pipeline output
# is not available.
# =============================================================================

library(tidyverse)
library(ggrepel)

if (!requireNamespace("MAGeCKFlute", quietly = TRUE)) {
  ReadBeta <- function(x) {
    df <- tibble::as_tibble(x)
    beta_cols <- grep("\\|beta$", names(df), value = TRUE)
    out <- dplyr::select(df, Gene, dplyr::all_of(beta_cols))
    names(out) <- sub("\\|beta$", "", names(out))
    out
  }
  NormalizeBeta <- function(beta, ...) beta
} else {
  library(MAGeCKFlute)
}

BETA_CUTOFF <- 0.5

# -- Locate gene summary ------------------------------------------------------

script_dir <- tryCatch(
  dirname(sys.frame(1)$ofile),
  error = function(e) "."
)
project_root <- normalizePath(file.path(script_dir, "../../../crispr_analysis"),
                              mustWork = FALSE)

primary_file <- file.path(project_root, "output", "combined_res_all_replicates",
                          "gadeta_mle_all_reps.gene_summary.txt")
demo_file <- normalizePath(
  file.path(script_dir, "../../demo/figure1_S1/gene_summary_demo.txt"),
  mustWork = FALSE)

if (file.exists(primary_file)) {
  gene_sum_file <- primary_file
  message("Using pipeline output: ", gene_sum_file)
} else if (file.exists(demo_file)) {
  gene_sum_file <- demo_file
  message("Pipeline output not found; using demo: ", demo_file)
} else {
  stop("No gene summary found.\n",
       "  Expected: ", primary_file, "\n",
       "  Or demo:  ", demo_file)
}

gene_sum <- read_tsv(gene_sum_file, show_col_types = FALSE)
gdata    <- ReadBeta(gene_sum)
gdata    <- NormalizeBeta(beta = gdata, id = 1, method = "loess")

# -- Compute beta-score diff (GDT201 1:1 minus NC 1:1) ------------------------

df_1_1 <- gdata %>%
  dplyr::select(Gene, ctrl_1_1, gdt_1_1) %>%
  mutate(diff = gdt_1_1 - ctrl_1_1,
         rank = rank(diff))

# -- Output directory ----------------------------------------------------------

fig_dir <- file.path(script_dir, "../../output/figure1_S1")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

save_fig <- function(name, p, w = 7, h = 7) {
  ggsave(file.path(fig_dir, paste0(name, ".pdf")), p, width = w, height = h)
  ggsave(file.path(fig_dir, paste0(name, ".png")), p, width = w, height = h,
         dpi = 300)
  message("  Saved: ", file.path(fig_dir, name), " (.pdf + .png)")
}

# =============================================================================
# Figure 1C — beta-score rank plot (GDT201 1:1 vs NC 1:1)
# Blue = positively selected (beta > 0.5); Red = negatively selected (< -0.5)
# =============================================================================

message("\n--- Figure 1C: GDT201 1:1 vs NC 1:1 rank plot ---")

labels_1c <- df_1_1 %>%
  filter(abs(diff) > BETA_CUTOFF) %>%
  filter(!grepl("NonTargetingControlGuideForHuman", Gene))

keyvals <- ifelse(df_1_1$diff > BETA_CUTOFF, "blue",
                  ifelse(df_1_1$diff < -BETA_CUTOFF, "red", "grey"))

p_1c <- ggplot(df_1_1, aes(x = diff, y = rank)) +
  geom_point(colour = keyvals, size = 0.8) +
  geom_text_repel(data = labels_1c, aes(label = Gene),
                  size = 2.5, max.overlaps = 60) +
  xlab("Treatment - Control Beta Score") +
  ylab("Rank") +
  ylim(-500, 25000) +
  theme_bw()

save_fig("fig1C_rank_gdt201_1_1_vs_nc", p_1c)

# =============================================================================
# Figure 1D — Curated table of top positively selected genes
#
# Selected genes of biological interest with their functional annotation.
# Rank = position among all genes sorted by descending beta-score (diff).
# =============================================================================

message("\n--- Figure 1D: Curated gene table ---")

gene_annotations <- tribble(
  ~Gene,     ~Alias,            ~Function,
  "TM9SF3",  "TM9SF3",          "Golgi regulation",
  "MAN2A1",  "MAN2A1 (aManII)", "Glycosylation",
  "MGAT2",   "MGAT2 (GnT-II)",  "Glycosylation",
  "FUT4",    "FUT4",             "Fucosylation",
  "MGAT5",   "MGAT5 (GnT-V)",   "Glycosylation",
  "SLC35C1", "SLC35C1",          "Fucosylation",
  "PROCR",   "EPCR (PROCR)",     "TCR Ligand",
  "GMDS",    "GMDS",             "Fucosylation",
  "TSTA3",   "TSTA3",            "Fucosylation",
  "MAN1A2",  "MAN1A2",           "Glycosylation",
  "MGAT1",   "MGAT1 (GnT-I)",   "Glycosylation"
)

pos_ranked <- df_1_1 %>%
  arrange(desc(diff)) %>%
  mutate(Rank = row_number(),
         Beta_Score = round(diff, 3))

fig1d_table <- pos_ranked %>%
  inner_join(gene_annotations, by = "Gene") %>%
  dplyr::select(Rank, Beta_Score, Gene = Alias, Function) %>%
  arrange(Rank)

write_csv(fig1d_table, file.path(fig_dir, "fig1D_top_genes.csv"))
message("  Wrote: ", file.path(fig_dir, "fig1D_top_genes.csv"))
message("  Table:")
print(as.data.frame(fig1d_table), row.names = FALSE, right = FALSE)

message("\nDone. Output in: ", normalizePath(fig_dir, mustWork = FALSE))
