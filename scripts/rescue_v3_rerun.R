# ============================================================================
# rescue_v3_rerun.R  — Stage-2 (significant RTN) + kidney-RTN reconciliation
# Run on the server (mpoets2: DESeq2 1.46.0, ashr 2.2.63). One per tissue.
# Goal: reproduce Stage-1 RTN (kidney must hit 1,961 / muscle 491), then add
#       Stage 2 = Stage-1 ∩ padj(combi vs age) < 0.05. Nothing else changes.
# Conventions (CLAUDE.md): use the cleaned *_df tables, NEVER raw lfcShrink objects.
#   res_age        = age   vs ctrl   -> LFC_age,        padj_age
#   res_intv       = combi vs age    -> LFC_combi_age,  padj_combi_age  (primary rescue contrast)
#   res_combi_ctrl = combi vs ctrl   -> LFC_combi_ctrl                  (normalisation check)
# ============================================================================

library(DESeq2); library(dplyr)

tissue <- "kidney"   # <-- set to "kidney" or "muscle_v2" and run twice
outdir <- "results/rescue_v3"; dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# --- 1. Load the cleaned _df tables for this tissue --------------------------
# Point these at your saved DESeq2 outputs (deseq2_out/). Each has a gene id col + log2FoldChange + padj.
base <- if (tissue == "kidney") "original_data/bulk/kidney_bulk/deseq2_out" else
                                 "original_data/bulk/muscle_bulk/v2/deseq2_out"
res_age        <- read.csv(file.path(base, "res_age_df.csv"),        row.names = 1)
res_intv       <- read.csv(file.path(base, "res_intv_df.csv"),       row.names = 1)   # combi vs age
res_combi_ctrl <- read.csv(file.path(base, "res_combi_ctrl_df.csv"), row.names = 1)   # combi vs ctrl

m <- data.frame(
  gene           = rownames(res_age),
  LFC_age        = res_age$log2FoldChange,        padj_age       = res_age$padj,
  LFC_combi_age  = res_intv[rownames(res_age), "log2FoldChange"],
  padj_combi_age = res_intv[rownames(res_age), "padj"],
  LFC_combi_ctrl = res_combi_ctrl[rownames(res_age), "log2FoldChange"]
)

# --- 2. Stage 1 (RTN) + RECONCILE kidney to 1,961 ----------------------------
aging   <- m %>% filter(!is.na(padj_age), padj_age < 0.05)
reversal <- aging %>% filter(sign(LFC_age) != sign(LFC_combi_age))     # processing step (4,645 / 719)
RTN      <- reversal %>% filter(abs(LFC_combi_ctrl) < 0.5)             # Stage 1 (target: kidney 1,961 / muscle 491)

cat(sprintf("[%s] aging=%d  reversal=%d  RTN(Stage1)=%d\n",
            tissue, nrow(aging), nrow(reversal), nrow(RTN)))

# If kidney RTN != 1,961, the divergence is the |LFC_combi_ctrl|<0.5 step (shrunk vs MLE LFC).
# Diagnostic: how many pass with shrunk LFC (above) vs MLE/unshrunken LFC?  Needs the dds object:
#   dds <- readRDS(...);  mle <- as.data.frame(results(dds, contrast=c("condition","combi","ctrl")))
#   reversal$LFC_combi_ctrl_mle <- mle[reversal$gene, "log2FoldChange"]
#   cat("RTN with MLE filter:", sum(abs(reversal$LFC_combi_ctrl_mle) < 0.5, na.rm=TRUE), "\n")
#   hist(abs(reversal$LFC_combi_ctrl), breaks=80); abline(v=0.5)   # inspect 0.45–0.55 boundary
# Pin whichever recipe reproduces 1,961 and use it for RTN above. Do NOT proceed until kidney RTN==1,961.

# --- 3. Stage 2 (significant RTN) — the ONLY new number ----------------------
sigRTN <- RTN %>% filter(!is.na(padj_combi_age), padj_combi_age < 0.05)
cat(sprintf("[%s] Stage2 (significant RTN)=%d  (%.0f%% of Stage 1)\n",
            tissue, nrow(sigRTN), 100 * nrow(sigRTN) / nrow(RTN)))

# --- 4. Write outputs (flat names are what Fig 3.7D cell 146 auto-loads) ------
write.csv(RTN,    file.path(outdir, sprintf("%s_stage1_RTN.csv",    tissue)), row.names = FALSE)
write.csv(sigRTN, file.path(outdir, sprintf("%s_stage2_sigRTN.csv", tissue)), row.names = FALSE)
write.csv(
  data.frame(metric = c("aging_DEGs","reversal","RTN_stage1","sigRTN_stage2"),
             n = c(nrow(aging), nrow(reversal), nrow(RTN), nrow(sigRTN))),
  file.path(outdir, sprintf("%s_counts_summary.csv", tissue)), row.names = FALSE)

# After running for BOTH tissues: send counts_summary.csv + the two *_stage2_sigRTN.csv back.
