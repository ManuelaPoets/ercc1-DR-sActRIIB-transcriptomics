# Faithful reproduction of the notebook's rescue classification, using the
# full-model ashr contrasts in deseq2_out/. Definitions taken verbatim from
# 1_DESeq2_kidney__multi-omic_rescue_pipeline.ipynb (cells ~1168-1200, 1180-1182).
arc_dir <- file.path(Sys.getenv("TEMP"), "kidney_archive")
od <- "results/kidney_bulk/deseq2_out"
rd <- function(f) { x <- read.csv(file.path(od,f), stringsAsFactors=FALSE); rownames(x)<-x$gene_id; x }

res_age        <- rd("res_age_df.csv")
res_intv       <- rd("res_intv_df.csv")        # combi vs age
res_intv1      <- rd("res_intv1_df.csv")       # DR vs age
res_intv2      <- rd("res_intv2_df.csv")       # sAct vs age
res_combi_ctrl <- rd("res_combi_ctrl_df.csv")  # combi vs ctrl

L <- function(df, ids) df[ids, "log2FoldChange"]

# deg_age: padj<0.05 & |LFC|>0
deg_age <- res_age[which(res_age$padj < 0.05 & abs(res_age$log2FoldChange) > 0), ]
ga <- rownames(deg_age)

# rescued by combi (res_intv): opposite sign  (which() => NA dropped)
rescued_combi_ids <- ga[which(sign(deg_age$log2FoldChange) != sign(L(res_intv, ga)))]
# multi-intervention sets (no which(): NA stays NA, harmless in setdiff/union)
rescued_combi <- ga[sign(deg_age$log2FoldChange) != sign(L(res_intv,  ga))]
rescued_DR    <- ga[sign(deg_age$log2FoldChange) != sign(L(res_intv1, ga))]
rescued_sAct  <- ga[sign(deg_age$log2FoldChange) != sign(L(res_intv2, ga))]
rescued_all   <- unique(c(rescued_combi, rescued_DR, rescued_sAct))
non_rescued_genes <- setdiff(ga, rescued_all)

# rescued_genes (by combi) + combi-vs-ctrl info
rescued_genes <- deg_age[rescued_combi_ids, ]
rescued_genes$log2FC_combi_ctrl <- L(res_combi_ctrl, rescued_combi_ids)
rescued_genes$padj_combi_ctrl   <- res_combi_ctrl[rescued_combi_ids, "padj"]
rescued_to_normal_ids <- rescued_combi_ids[which(abs(rescued_genes$log2FC_combi_ctrl) < 0.5)]

# variants for rescued_to_normal definition
r2n_lfc      <- rescued_combi_ids[which(abs(rescued_genes$log2FC_combi_ctrl) < 0.5)]
r2n_padj     <- rescued_combi_ids[which(rescued_genes$padj_combi_ctrl > 0.05)]
r2n_both     <- rescued_combi_ids[which(abs(rescued_genes$log2FC_combi_ctrl) < 0.5 & rescued_genes$padj_combi_ctrl > 0.05)]
r2n_padj_na  <- rescued_combi_ids[which(rescued_genes$padj_combi_ctrl > 0.05 | is.na(rescued_genes$padj_combi_ctrl))]
r2n_both_na  <- rescued_combi_ids[which(abs(rescued_genes$log2FC_combi_ctrl) < 0.5 & (rescued_genes$padj_combi_ctrl > 0.05 | is.na(rescued_genes$padj_combi_ctrl)))]

cat(sprintf("deg_age: %d | rescued(combi): %d | rescued_all(any): %d | non_rescued: %d | rescued_to_normal: %d\n",
            nrow(deg_age), length(rescued_combi_ids), length(na.omit(rescued_all)),
            length(non_rescued_genes), length(rescued_to_normal_ids)))

# ---- Compare to archive -----------------------------------------------------
r2n  <- read.csv(file.path(arc_dir,"kidney_rescued_to_normal.csv"), stringsAsFactors=FALSE); colnames(r2n)[1]<-"gene_id"
nonr <- read.csv(file.path(arc_dir,"kidney_non_rescued_aged_DEGs.csv"), stringsAsFactors=FALSE); colnames(nonr)[1]<-"gene_id"
ov <- function(a,b,lab){ i<-length(intersect(a,b)); cat(sprintf("%-20s n_new=%5d n_arc=%5d overlap=%5d new-only=%4d arc-only=%4d Jaccard=%.3f\n",
      lab,length(a),length(b),i,length(setdiff(a,b)),length(setdiff(b,a)),i/length(union(a,b)))) }
cat("\n== vs archive (non_rescued) ==\n")
ov(non_rescued_genes,     nonr$gene_id, "non_rescued")
cat("\n== rescued_to_normal: which gate matches archive (1961)? ==\n")
ov(r2n_lfc,     r2n$gene_id, "|LFC_cc|<0.5")
ov(r2n_padj,    r2n$gene_id, "padj_cc>0.05")
ov(r2n_both,    r2n$gene_id, "both")
ov(r2n_padj_na, r2n$gene_id, "padj>0.05|NA")
ov(r2n_both_na, r2n$gene_id, "both(padj|NA)")
