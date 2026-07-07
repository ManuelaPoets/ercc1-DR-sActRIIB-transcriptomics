# Compare new tximport+DESeq2 res_age vs archived kidney res_age table.
new <- read.csv("results/kidney_bulk/deseq2_out/res_age_df.csv", stringsAsFactors = FALSE)
arc <- read.csv(file.path(Sys.getenv("TEMP"), "kidney_archive", "kidney_res_age_archive.csv"),
                stringsAsFactors = FALSE)
colnames(arc)[1] <- "gene_id"

cat("== sizes ==\n")
cat(sprintf("new genes: %d | archive genes: %d\n", nrow(new), nrow(arc)))

cat("\n== padj<0.05 counts ==\n")
cat(sprintf("new: %d | archive: %d\n",
            sum(new$padj < 0.05, na.rm = TRUE),
            sum(arc$padj < 0.05, na.rm = TRUE)))

# significant gene-set overlap
ns <- new$gene_id[which(new$padj < 0.05)]
as <- arc$gene_id[which(arc$padj < 0.05)]
cat(sprintf("sig overlap: %d | new-only: %d | archive-only: %d | Jaccard: %.4f\n",
            length(intersect(ns, as)), length(setdiff(ns, as)),
            length(setdiff(as, ns)),
            length(intersect(ns, as)) / length(union(ns, as))))

# merge on gene_id
m <- merge(new[, c("gene_id","log2FoldChange","padj")],
           arc[, c("gene_id","log2FoldChange","padj")],
           by = "gene_id", suffixes = c(".new",".arc"))
cat(sprintf("\n== merged on gene_id: %d common genes ==\n", nrow(m)))

cat("\n== padj concordance (should be ~identical) ==\n")
cat(sprintf("Pearson r(padj): %.6f | max abs diff: %.3e\n",
            cor(m$padj.new, m$padj.arc, use = "complete.obs"),
            max(abs(m$padj.new - m$padj.arc), na.rm = TRUE)))

cat("\n== LFC: new(ashr-shrunk) vs archive(MLE) ==\n")
cat(sprintf("Pearson r(LFC): %.4f | Spearman: %.4f\n",
            cor(m$log2FoldChange.new, m$log2FoldChange.arc, use = "complete.obs"),
            cor(m$log2FoldChange.new, m$log2FoldChange.arc, method = "spearman", use = "complete.obs")))
sign_agree <- mean(sign(m$log2FoldChange.new) == sign(m$log2FoldChange.arc), na.rm = TRUE)
cat(sprintf("LFC sign agreement: %.4f\n", sign_agree))
