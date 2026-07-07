# Compare PAIRWISE outputs vs archived kidney tables (res_age + rescue tables).
arc_dir <- file.path(Sys.getenv("TEMP"), "kidney_archive")
od <- "results/kidney_bulk/deseq2_out_pairwise"

arc  <- read.csv(file.path(arc_dir, "kidney_res_age_archive.csv"), stringsAsFactors=FALSE); colnames(arc)[1]<-"gene_id"
comb <- read.csv(file.path(arc_dir, "kidney_combined_rescued_genes_DE_info.csv"), stringsAsFactors=FALSE)
r2n  <- read.csv(file.path(arc_dir, "kidney_rescued_to_normal.csv"), stringsAsFactors=FALSE); colnames(r2n)[1]<-"gene_id"
nonr <- read.csv(file.path(arc_dir, "kidney_non_rescued_aged_DEGs.csv"), stringsAsFactors=FALSE); colnames(nonr)[1]<-"gene_id"

age  <- read.csv(file.path(od,"res_age_df.csv"), stringsAsFactors=FALSE)
intv <- read.csv(file.path(od,"res_intv_df.csv"), stringsAsFactors=FALSE)
cc   <- read.csv(file.path(od,"res_combi_ctrl_df.csv"), stringsAsFactors=FALSE)

J <- function(a,b) length(intersect(a,b))/length(union(a,b))
ov <- function(a,b,lab){ i<-length(intersect(a,b)); cat(sprintf("%-22s n_new=%5d n_arc=%5d overlap=%5d new-only=%5d arc-only=%4d Jaccard=%.3f\n",
      lab,length(a),length(b),i,length(setdiff(a,b)),length(setdiff(b,a)),i/length(union(a,b)))) }

cat("== res_age: pairwise vs archive ==\n")
cat(sprintf("pairwise sig: %d | archive sig: %d\n", sum(age$padj<0.05,na.rm=TRUE), sum(arc$padj<0.05,na.rm=TRUE)))
ov(age$gene_id[which(age$padj<0.05)], arc$gene_id[which(arc$padj<0.05)], "res_age sig")

cat("\n== rescue reconstruction (pairwise) ==\n")
aged <- age[which(age$padj<0.05),]
m <- merge(aged[,c("gene_id","log2FoldChange","padj")], intv[,c("gene_id","log2FoldChange","padj")], by="gene_id", suffixes=c("_aging","_intv"))
m <- merge(m, cc[,c("gene_id","log2FoldChange","padj")], by="gene_id"); colnames(m)[(ncol(m)-1):ncol(m)]<-c("log2FoldChange_cc","padj_cc")
m <- m[!is.na(m$padj_intv),]
m$rescued <- sign(m$log2FoldChange_aging) != sign(m$log2FoldChange_intv)
m$rescued_to_normal <- m$rescued & (is.na(m$padj_cc) | m$padj_cc >= 0.05)
cat(sprintf("aging DEGs w/ valid intv: %d | rescued: %d | non-rescued: %d | rescued_to_normal: %d\n",
            nrow(m), sum(m$rescued), sum(!m$rescued), sum(m$rescued_to_normal)))
ov(m$gene_id[m$rescued],           comb$GeneID, "rescued~combined")
ov(m$gene_id[!m$rescued],          nonr$gene_id,"non_rescued~archive")
ov(m$gene_id[m$rescued_to_normal], r2n$gene_id, "rescued_to_normal~arc")
