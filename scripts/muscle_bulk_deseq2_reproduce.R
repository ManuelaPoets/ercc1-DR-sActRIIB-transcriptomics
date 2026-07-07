# ============================================================================
# Muscle bulk RNA-seq: faithful reproduction of the archived muscle_v2 tables.
# Source = GEO featureCounts matrix (GSE268971), per notebook
# 2_DESeq2_muscle_multi-omic_rescue_pipeline.ipynb.
# Full model ~condition, ashr shrinkage, IDs version-stripped after results.
# ============================================================================
suppressMessages({ library(DESeq2); library(ashr) })

data_dir <- "results/muscle_bulk"
out_dir  <- file.path(data_dir, "deseq2_out"); dir.create(out_dir, showWarnings = FALSE)
arc_dir  <- file.path(Sys.getenv("TEMP"), "muscle_archive")

# ---- Counts + sample sheet (notebook cells 468-484) ------------------------
muscle_counts <- read.table(file.path(data_dir, "GSE268971_Muscle_counts_file.txt.gz"),
                            header = TRUE, sep = "\t", row.names = 1)
ss <- read.csv(file.path(data_dir, "muscle_sample_metadata_filtered.csv"),
               header = TRUE, sep = ",")
ss$sample_name_full <- paste0(ss$sample_name, "Aligned.out.bam")
rownames(ss) <- ss$sample_name_full
muscle_counts <- muscle_counts[, rownames(ss)]              # drop annotation cols, order to ss
ss$condition  <- factor(ss$condition, levels = c("ctrl","age","DR","sAct","combi"))
stopifnot(ncol(muscle_counts) == 15L, nrow(ss) == 15L)

# ---- DESeq2 (cells 524-535): filter rowSums>=10, full model --------------
dds <- DESeqDataSetFromMatrix(countData = muscle_counts, colData = ss, design = ~condition)
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)

shrink <- function(num, den) {
  r <- results(dds, contrast = c("condition", num, den))
  as.data.frame(lfcShrink(dds, contrast = c("condition", num, den), type = "ashr", res = r))
}
strip <- function(df) { rownames(df) <- sub("\\.[0-9]+$", "", rownames(df)); df }

res_age_df        <- strip(shrink("age",   "ctrl"))
res_intv_df       <- strip(shrink("combi", "age"))
res_intv1_df      <- strip(shrink("DR",    "age"))
res_intv2_df      <- strip(shrink("sAct",  "age"))
res_combi_ctrl_df <- strip(shrink("combi", "ctrl"))
for (nm in c("res_age_df","res_intv_df","res_intv1_df","res_intv2_df","res_combi_ctrl_df")) {
  d <- get(nm); d$gene_id <- rownames(d)
  write.csv(d[order(d$padj), ], file.path(out_dir, paste0(nm, ".csv")), row.names = FALSE)
}

# ---- Rescue classification (cells 677-771) ---------------------------------
L <- function(df, ids) df[ids, "log2FoldChange"]
deg_age <- res_age_df[which(res_age_df$padj < 0.05 & abs(res_age_df$log2FoldChange) > 0), ]
ga <- rownames(deg_age)
rescued_combi <- ga[sign(deg_age$log2FoldChange) != sign(L(res_intv_df,  ga))]
rescued_DR    <- ga[sign(deg_age$log2FoldChange) != sign(L(res_intv1_df, ga))]
rescued_sAct  <- ga[sign(deg_age$log2FoldChange) != sign(L(res_intv2_df, ga))]
rescued_all   <- unique(c(rescued_combi, rescued_DR, rescued_sAct))
non_rescued   <- setdiff(ga, rescued_all)
rescued_combi_ids <- ga[which(sign(deg_age$log2FoldChange) != sign(L(res_intv_df, ga)))]
lfc_cc <- L(res_combi_ctrl_df, rescued_combi_ids)
rescued_to_normal <- rescued_combi_ids[which(abs(lfc_cc) < 0.5)]

cat(sprintf("deg_age: %d | rescued(combi): %d | rescued_all(any): %d | non_rescued: %d | rescued_to_normal: %d\n",
            nrow(deg_age), length(rescued_combi_ids), length(na.omit(rescued_all)),
            length(non_rescued), length(rescued_to_normal)))

# ---- Compare to archive -----------------------------------------------------
rd <- function(f){ x<-read.csv(file.path(arc_dir,f),stringsAsFactors=FALSE); colnames(x)[1]<-"gene_id"; x }
a_age <- rd("muscle_v2_res_age_results_with_symbols.csv")
a_deg <- rd("muscle_v2_deg_age.csv")
a_nr  <- rd("muscle_v2_non_rescued_aged_DEGs.csv")
a_r2n <- rd("muscle_v2_rescued_to_normal.csv")
ov <- function(a,b,lab){ i<-length(intersect(a,b)); cat(sprintf("%-20s n_new=%5d n_arc=%5d overlap=%5d new-only=%4d arc-only=%4d Jaccard=%.3f\n",
      lab,length(a),length(b),i,length(setdiff(a,b)),length(setdiff(b,a)),i/length(union(a,b)))) }
cat("\n== vs archive ==\n")
cat(sprintf("res_age sig: new %d | archive %d\n", sum(res_age_df$padj<0.05,na.rm=TRUE), sum(a_age$padj<0.05,na.rm=TRUE)))
ov(rownames(deg_age),    a_deg$gene_id, "deg_age")
ov(non_rescued,          a_nr$gene_id,  "non_rescued")
ov(rescued_to_normal,    a_r2n$gene_id, "rescued_to_normal")
