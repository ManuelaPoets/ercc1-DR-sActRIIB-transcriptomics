# ============================================================================
# Muscle bulk v2: Salmon quant.sf -> tximport -> DESeq2, GRCm39.110.
# tx2gene built from the SAME GTF as the v2 notebook (2b_DESeq2_muscleV2_*).
# Full model ~condition, ashr shrinkage; rescue defs identical to kidney/muscle.
# Validates against archived muscle_v2_* tables.
# ============================================================================
suppressMessages({ library(tximport); library(DESeq2); library(ashr); library(txdbmaker) })

data_dir <- "results/muscle_bulk"
v2_dir   <- file.path(data_dir, "v2")
out_dir  <- file.path(v2_dir, "deseq2_out"); dir.create(out_dir, showWarnings = FALSE)
arc_dir  <- file.path(Sys.getenv("TEMP"), "muscle_archive")
gtf_path <- "genomes/Mus_musculus.GRCm39.110.gtf.gz"

# ---- tx2gene from GTF (verbatim from v2 notebook) ---------------------------
t2g_cache <- file.path(v2_dir, "tx2gene_GRCm39.110.tsv")
if (file.exists(t2g_cache)) {
  tx2gene <- read.delim(t2g_cache, stringsAsFactors = FALSE)
} else {
  txdb <- suppressWarnings(makeTxDbFromGFF(gtf_path))
  k <- keys(txdb, keytype = "TXNAME")
  tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")
  write.table(tx2gene, t2g_cache, sep = "\t", row.names = FALSE, quote = FALSE)
}
cat("tx2gene:", nrow(tx2gene), "transcripts |", length(unique(tx2gene$GENEID)), "genes\n")

# ---- Salmon files + coldata (condition from folder name) -------------------
qdirs   <- list.dirs(file.path(v2_dir, "salmon_out"), recursive = FALSE)
samples <- basename(qdirs)
files   <- file.path(qdirs, "quant.sf"); names(files) <- samples
stopifnot(all(file.exists(files)), length(files) == 15L)
coldata <- data.frame(row.names = samples,
                      condition = factor(sub("_[0-9]+$", "", samples),
                                         levels = c("ctrl","age","DR","sAct","combi")))
print(table(coldata$condition))

# ---- Import + DESeq2 --------------------------------------------------------
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
dds <- DESeqDataSetFromTximport(txi, colData = coldata, design = ~condition)
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)

shrink <- function(num, den) {
  r <- results(dds, contrast = c("condition", num, den))
  d <- as.data.frame(lfcShrink(dds, contrast = c("condition", num, den), type = "ashr", res = r))
  d$gene_id <- rownames(d); d
}
res_age_df        <- shrink("age",   "ctrl")
res_intv_df       <- shrink("combi", "age")
res_intv1_df      <- shrink("DR",    "age")
res_intv2_df      <- shrink("sAct",  "age")
res_combi_ctrl_df <- shrink("combi", "ctrl")
for (nm in c("res_age_df","res_intv_df","res_intv1_df","res_intv2_df","res_combi_ctrl_df")) {
  d <- get(nm); write.csv(d[order(d$padj), ], file.path(out_dir, paste0(nm, ".csv")), row.names = FALSE)
}

# ---- Rescue classification (same defs as kidney/muscle notebooks) ----------
L <- function(df, ids) df[match(ids, df$gene_id), "log2FoldChange"]
deg_age <- res_age_df[which(res_age_df$padj < 0.05 & abs(res_age_df$log2FoldChange) > 0), ]
ga <- deg_age$gene_id
rescued_combi <- ga[sign(deg_age$log2FoldChange) != sign(L(res_intv_df,  ga))]
rescued_DR    <- ga[sign(deg_age$log2FoldChange) != sign(L(res_intv1_df, ga))]
rescued_sAct  <- ga[sign(deg_age$log2FoldChange) != sign(L(res_intv2_df, ga))]
rescued_all   <- unique(c(rescued_combi, rescued_DR, rescued_sAct))
non_rescued   <- setdiff(ga, rescued_all)
rescued_combi_ids <- ga[which(sign(deg_age$log2FoldChange) != sign(L(res_intv_df, ga)))]
lfc_cc <- L(res_combi_ctrl_df, rescued_combi_ids)
rescued_to_normal <- rescued_combi_ids[which(abs(lfc_cc) < 0.5)]
cat(sprintf("\ndeg_age: %d | rescued(combi): %d | rescued_all: %d | non_rescued: %d | rescued_to_normal: %d\n",
            nrow(deg_age), length(rescued_combi_ids), length(na.omit(rescued_all)),
            length(non_rescued), length(rescued_to_normal)))

# ---- Compare to archive -----------------------------------------------------
rd <- function(f){ x<-read.csv(file.path(arc_dir,f),stringsAsFactors=FALSE); colnames(x)[1]<-"gene_id"; x }
a_age <- rd("muscle_v2_res_age_results_with_symbols.csv")
a_deg <- rd("muscle_v2_deg_age.csv"); a_nr <- rd("muscle_v2_non_rescued_aged_DEGs.csv"); a_r2n <- rd("muscle_v2_rescued_to_normal.csv")
ov <- function(a,b,lab){ i<-length(intersect(a,b)); cat(sprintf("%-20s n_new=%5d n_arc=%5d overlap=%5d new-only=%4d arc-only=%4d Jaccard=%.3f\n",
      lab,length(a),length(b),i,length(setdiff(a,b)),length(setdiff(b,a)),i/length(union(a,b)))) }
cat("\n== vs archive muscle_v2 ==\n")
cat(sprintf("res_age sig: new %d | archive %d\n", sum(res_age_df$padj<0.05,na.rm=TRUE), sum(a_age$padj<0.05,na.rm=TRUE)))
ov(deg_age$gene_id,   a_deg$gene_id, "deg_age")
ov(non_rescued,       a_nr$gene_id,  "non_rescued")
ov(rescued_to_normal, a_r2n$gene_id, "rescued_to_normal")
