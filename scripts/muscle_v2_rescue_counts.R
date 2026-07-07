od <- "results/muscle_bulk/v2/deseq2_out"
rd <- function(f) read.csv(file.path(od, f), stringsAsFactors = FALSE)
age <- rd("res_age_df.csv"); iv <- rd("res_intv_df.csv")
d1 <- rd("res_intv1_df.csv"); d2 <- rd("res_intv2_df.csv")
L <- function(df, ids) df[match(ids, df$gene_id), "log2FoldChange"]
deg <- age[which(age$padj < 0.05 & abs(age$log2FoldChange) > 0), ]; ga <- deg$gene_id
opp <- function(x) ga[sign(deg$log2FoldChange) != sign(L(x, ga))]
rc <- na.omit(opp(iv)); rdd <- na.omit(opp(d1)); rs <- na.omit(opp(d2))
ra <- unique(c(rc, rdd, rs)); nr <- setdiff(ga, ra)
cat(sprintf("deg_age=%d rescued_combi=%d rescued_DR=%d rescued_sAct=%d rescued_all=%d non_rescued=%d (%.1f%%)\n",
            nrow(deg), length(rc), length(rdd), length(rs), length(ra), length(nr), 100*length(nr)/nrow(deg)))
