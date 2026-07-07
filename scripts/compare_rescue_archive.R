# Compare new contrasts vs archived kidney rescue tables.
arc_dir <- file.path(Sys.getenv("TEMP"), "kidney_archive")

comb <- read.csv(file.path(arc_dir, "kidney_combined_rescued_genes_DE_info.csv"), stringsAsFactors = FALSE)
r2n  <- read.csv(file.path(arc_dir, "kidney_rescued_to_normal.csv"), stringsAsFactors = FALSE)
nonr <- read.csv(file.path(arc_dir, "kidney_non_rescued_aged_DEGs.csv"), stringsAsFactors = FALSE)
colnames(r2n)[1]  <- "gene_id"
colnames(nonr)[1] <- "gene_id"

cat("================ ARCHIVE DEFINITION CHECKS ================\n")
cat(sprintf("combined table rows: %d\n", nrow(comb)))
cat("rescue_significant table:\n"); print(table(comb$rescue_significant, useNA="ifany"))
# Infer archive rule: how does rescue_significant relate to LFC signs?
opp <- sign(comb$log2FoldChange_aging) != sign(comb$log2FoldChange_intervention)
cat("\nrescue_significant vs (aging/intervention opposite sign):\n")
print(table(rescue_sig = comb$rescue_significant, opposite_sign = opp))
cat(sprintf("\nall combined padj_aging<0.05? %s (min %.3g, max %.3g)\n",
            all(comb$padj_aging < 0.05, na.rm=TRUE), min(comb$padj_aging,na.rm=TRUE), max(comb$padj_aging,na.rm=TRUE)))
cat(sprintf("rescued_to_normal rows: %d | non_rescued rows: %d\n", nrow(r2n), nrow(nonr)))

cat("\n================ RECONSTRUCT FROM NEW CONTRASTS ================\n")
age  <- read.csv("results/kidney_bulk/deseq2_out/res_age_df.csv", stringsAsFactors=FALSE)
intv <- read.csv("results/kidney_bulk/deseq2_out/res_intv_df.csv", stringsAsFactors=FALSE)        # combi vs age
cc   <- read.csv("results/kidney_bulk/deseq2_out/res_combi_ctrl_df.csv", stringsAsFactors=FALSE)  # combi vs ctrl

# aging DEGs
aged <- age[which(age$padj < 0.05), ]
cat(sprintf("new aging DEGs (padj<0.05): %d\n", nrow(aged)))

m <- merge(aged[,c("gene_id","log2FoldChange","padj")],
           intv[,c("gene_id","log2FoldChange","padj")], by="gene_id", suffixes=c("_aging","_intv"))
m <- merge(m, cc[,c("gene_id","log2FoldChange","padj")], by="gene_id")
colnames(m)[(ncol(m)-1):ncol(m)] <- c("log2FoldChange_cc","padj_cc")

# KEY: only classify genes with VALID intervention data (independent-filtering gate)
m <- m[!is.na(m$padj_intv), ]
cat(sprintf("aging DEGs with valid intervention data: %d\n", nrow(m)))

# rescue = sign change (intervention opposes aging) -- CLAUDE.md: sign-change only
m$rescued <- sign(m$log2FoldChange_aging) != sign(m$log2FoldChange_intv)
# rescued_to_normal = rescued AND combi vs ctrl NOT significant (restored)
m$rescued_to_normal <- m$rescued & (is.na(m$padj_cc) | m$padj_cc >= 0.05)

cat(sprintf("new rescued (sign-change): %d\n", sum(m$rescued, na.rm=TRUE)))
cat(sprintf("new non-rescued (same sign): %d\n", sum(!m$rescued, na.rm=TRUE)))
cat(sprintf("new rescued_to_normal: %d\n", sum(m$rescued_to_normal, na.rm=TRUE)))

cat("\n================ SET OVERLAP vs ARCHIVE ================\n")
ov <- function(a,b,labA,labB){
  i<-length(intersect(a,b)); cat(sprintf("%-22s n_new=%d n_arc=%d  overlap=%d  new-only=%d  arc-only=%d  Jaccard=%.3f\n",
      paste0(labA,"~",labB), length(a), length(b), i, length(setdiff(a,b)), length(setdiff(b,a)),
      i/length(union(a,b))))
}
ov(m$gene_id[m$rescued],           comb$GeneID,  "rescued",           "combined")
ov(m$gene_id[m$rescued_to_normal], r2n$gene_id,  "rescued_to_normal", "archive")
ov(m$gene_id[!m$rescued],          nonr$gene_id, "non_rescued",       "archive")
ov(aged$gene_id,                   comb$GeneID,  "aging_DEGs",        "combined")
