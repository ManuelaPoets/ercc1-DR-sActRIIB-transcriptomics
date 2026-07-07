# Pin down the archive's non_rescued definition.
arc_dir <- file.path(Sys.getenv("TEMP"), "kidney_archive")
nonr <- read.csv(file.path(arc_dir, "kidney_non_rescued_aged_DEGs.csv"), stringsAsFactors=FALSE)
colnames(nonr)[1] <- "gene_id"

age  <- read.csv("results/kidney_bulk/deseq2_out/res_age_df.csv", stringsAsFactors=FALSE)
intv <- read.csv("results/kidney_bulk/deseq2_out/res_intv_df.csv", stringsAsFactors=FALSE)
cc   <- read.csv("results/kidney_bulk/deseq2_out/res_combi_ctrl_df.csv", stringsAsFactors=FALSE)

aged <- age[which(age$padj < 0.05), ]
m <- merge(aged[,c("gene_id","log2FoldChange","padj")],
           intv[,c("gene_id","log2FoldChange","padj")], by="gene_id", suffixes=c("_aging","_intv"))
m <- merge(m, cc[,c("gene_id","log2FoldChange","padj")], by="gene_id")
colnames(m)[(ncol(m)-1):ncol(m)] <- c("log2FoldChange_cc","padj_cc")
m <- m[!is.na(m$padj_intv), ]
m$same_sign <- sign(m$log2FoldChange_aging) == sign(m$log2FoldChange_intv)

J <- function(a,b) length(intersect(a,b))/length(union(a,b))
test <- function(set, label){
  cat(sprintf("%-46s n=%5d  overlap=%4d  arc-only=%4d  Jaccard=%.3f\n",
      label, length(set), length(intersect(set, nonr$gene_id)),
      length(setdiff(nonr$gene_id, set)), J(set, nonr$gene_id)))
}
cat(sprintf("archive non_rescued n = %d\n\n", nrow(nonr)))
test(m$gene_id[m$same_sign], "same-sign (valid intv)")
test(m$gene_id[m$same_sign & m$padj_intv < 0.05], "same-sign & intv significant")
test(m$gene_id[m$same_sign & !is.na(m$padj_cc) & m$padj_cc < 0.05], "same-sign & combi-vs-ctrl significant")
test(m$gene_id[m$same_sign & m$padj_intv < 0.05 & !is.na(m$padj_cc) & m$padj_cc < 0.05], "same-sign & intv sig & cc sig")
# directionality: is archive non_rescued only up or only down in aging?
hit <- m[m$gene_id %in% nonr$gene_id, ]
cat(sprintf("\narchive non_rescued among my genes: %d; aging LFC up=%d down=%d\n",
            nrow(hit), sum(hit$log2FoldChange_aging>0), sum(hit$log2FoldChange_aging<0)))
cat(sprintf("of these, same_sign=%d  opposite_sign=%d\n", sum(hit$same_sign), sum(!hit$same_sign)))
