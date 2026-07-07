# Independent re-run of muscle bulk DESeq2 from raw salmon quants.
# Goal: confirm (or refute) that sActRIIB has ~no standalone effect in muscle,
# without using any of Honey's derived tables.
suppressMessages({library(tximport); library(DESeq2)})

SALMON <- "original_data/bulk/muscle_bulk/v2/salmon_out"
GTF    <- "original_data/genomes/Mus_musculus.GRCm39.110.gtf.gz"

# --- tx2gene fresh from the GTF (manual parse of transcript lines) ---
cat("parsing GTF for tx2gene...\n")
con <- gzfile(GTF); L <- readLines(con); close(con)
tx <- L[grepl("\ttranscript\t", L, fixed = TRUE)]
tid <- sub('.*transcript_id "([^"]+)".*', "\\1", tx)
gid <- sub('.*gene_id "([^"]+)".*',       "\\1", tx)
gnm <- ifelse(grepl('gene_name "', tx), sub('.*gene_name "([^"]+)".*', "\\1", tx), gid)
tx2gene <- data.frame(TXNAME = tid, GENEID = gid)
g2sym   <- unique(data.frame(GENEID = gid, SYMBOL = gnm))
cat("  transcripts:", nrow(tx2gene), " genes:", length(unique(gid)), "\n")

# --- import all 15 samples ---
samples <- c(paste0("ctrl_", 1:3), paste0("age_", 1:3), paste0("DR_", 1:3),
             paste0("sAct_", 1:3), paste0("combi_", 1:3))
files <- file.path(SALMON, samples, "quant.sf")
stopifnot(all(file.exists(files)))
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
colnames(txi$counts) <- samples

condition <- factor(sub("_[0-9]+$", "", samples), levels = c("ctrl","age","DR","sAct","combi"))
dds <- DESeqDataSetFromTximport(txi, data.frame(condition), ~condition)
dds <- dds[rowSums(counts(dds)) >= 10, ]
cat("genes after >=10-count filter:", nrow(dds), "\n")
dds <- DESeq(dds, quiet = TRUE)

R <- function(a, b) as.data.frame(results(dds, contrast = c("condition", a, b)))  # MLE, Wald padj
res_age   <- R("age",   "ctrl")
res_dr    <- R("DR",    "age")
res_sact  <- R("sAct",  "age")
res_combi <- R("combi", "age")

sig <- function(d) sum(d$padj < 0.05, na.rm = TRUE)
mx  <- function(d) round(max(abs(d$log2FoldChange), na.rm = TRUE), 2)
cat("\n=== INDEPENDENT muscle contrasts (MLE LFC, Wald padj) ===\n")
cat(sprintf("  aging (age vs ctrl): %d DEGs\n", sig(res_age)))
cat(sprintf("  DR   vs age: %d sig | max|LFC| %.2f\n",  sig(res_dr),   mx(res_dr)))
cat(sprintf("  sAct vs age: %d sig | max|LFC| %.2f\n",  sig(res_sact), mx(res_sact)))
cat(sprintf("  combi vs age: %d sig | max|LFC| %.2f\n", sig(res_combi),mx(res_combi)))

# genes genuinely reversed by each intervention (aging DEG + intervention sig + opposite sign)
agdeg <- rownames(res_age)[which(res_age$padj < 0.05)]
rev_by <- function(rint) {
  g <- intersect(agdeg, rownames(rint)[which(rint$padj < 0.05)])
  g[sign(rint[g, "log2FoldChange"]) != sign(res_age[g, "log2FoldChange"])]
}
cat(sprintf("\n  aging DEGs genuinely REVERSED (sig + opposite sign): DR=%d  sAct=%d\n",
            length(rev_by(res_dr)), length(rev_by(res_sact))))
cat(sprintf("  |sAct LFC| of aging DEGs: median=%.4f  max=%.3f  (>0.25: %d / %d)\n",
            median(abs(res_sact[agdeg,"log2FoldChange"]), na.rm=TRUE),
            max(abs(res_sact[agdeg,"log2FoldChange"]), na.rm=TRUE),
            sum(abs(res_sact[agdeg,"log2FoldChange"])>0.25, na.rm=TRUE), length(agdeg)))

# cross-check Honey's 20 sActRIIB-exclusive genes in MY fresh run
her <- read.csv("results/Muscle_v2_analysis_results/3_Effect_comparison_tables/muscle_v2_intervention_impact_comparison_sAct.csv")
ids <- intersect(her$Gene, rownames(res_sact))
cat(sprintf("\n  Honey's 20 sActRIIB-exclusive genes -> in my fresh run: matched %d\n", length(ids)))
cat(sprintf("    their sAct LFC (my run): median=%.4f  max|LFC|=%.4f  any sig(padj<.05)? %d\n",
            median(res_sact[ids,"log2FoldChange"]), max(abs(res_sact[ids,"log2FoldChange"])),
            sum(res_sact[ids,"padj"]<0.05, na.rm=TRUE)))
