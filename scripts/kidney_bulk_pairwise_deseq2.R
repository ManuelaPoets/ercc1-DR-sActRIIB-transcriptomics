# ============================================================================
# Kidney bulk RNA-seq: per-pairwise DESeq2 (matches archive's modelling).
# Each contrast is fit on ONLY its two conditions (6 samples, n=3+3):
# separate dispersion estimation + independent filtering per comparison.
# Conventions: denominator = reference; ashr shrinkage; cleaned *_df outputs.
# ============================================================================

suppressMessages({ library(tximport); library(DESeq2); library(ashr) })

data_dir <- "results/kidney_bulk"
out_dir  <- file.path(data_dir, "deseq2_out_pairwise")
dir.create(out_dir, showWarnings = FALSE)

# ---- Sample sheet (5 aging conditions present) -----------------------------
coldata <- read.csv(file.path(data_dir, "Activin_kidney_info_selection.csv"),
                    stringsAsFactors = FALSE, check.names = FALSE)
coldata$condition <- trimws(coldata$condition)
keep_cond <- c("ctrl", "age", "DR", "sAct", "combi")
coldata   <- coldata[coldata$condition %in% keep_cond, ]
coldata$file <- file.path(data_dir, paste0(coldata$sample_name, "_quant.sf"))
coldata   <- coldata[file.exists(coldata$file), ]
stopifnot(nrow(coldata) == 15L)
rownames(coldata) <- coldata$sample_name

# ---- tx2gene ----------------------------------------------------------------
t2g <- read.delim(file.path(data_dir, "salmon_tx2gene.tsv"), header = FALSE,
                  stringsAsFactors = FALSE)
colnames(t2g) <- c("tx", "gene_id", "gene_symbol")
tx2gene  <- t2g[, c("tx", "gene_id")]
gene_sym <- unique(t2g[, c("gene_id", "gene_symbol")])
gene_sym <- gene_sym[!duplicated(gene_sym$gene_id), ]

# ---- Import once (all 15 samples) ------------------------------------------
files <- coldata$file; names(files) <- coldata$sample_name
txi   <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)

subset_txi <- function(txi, idx) list(
  abundance = txi$abundance[, idx, drop = FALSE],
  counts    = txi$counts[,    idx, drop = FALSE],
  length    = txi$length[,    idx, drop = FALSE],
  countsFromAbundance = txi$countsFromAbundance)

# ---- Per-pairwise contrast helper ------------------------------------------
run_pair <- function(num, den, name) {
  idx <- which(coldata$condition %in% c(num, den))
  cd  <- coldata[idx, ]
  cd$condition <- relevel(factor(cd$condition, levels = c(den, num)), ref = den)
  dds <- DESeqDataSetFromTximport(subset_txi(txi, idx), colData = cd, design = ~condition)
  dds <- dds[rowSums(counts(dds)) > 1, ]
  dds <- DESeq(dds)
  res <- results(dds, contrast = c("condition", num, den))
  res <- lfcShrink(dds, contrast = c("condition", num, den), type = "ashr", res = res)
  df  <- as.data.frame(res)
  df$gene_id     <- rownames(df)
  df$gene_symbol <- gene_sym$gene_symbol[match(df$gene_id, gene_sym$gene_id)]
  df  <- df[order(df$padj), c("gene_id","gene_symbol","baseMean",
                              "log2FoldChange","lfcSE","pvalue","padj")]
  write.csv(df, file.path(out_dir, paste0(name, ".csv")), row.names = FALSE)
  message(sprintf("%-16s %s vs %s | n=%d | tested(non-NA padj): %d | padj<0.05: %d",
                  name, num, den, nrow(cd), sum(!is.na(df$padj)),
                  sum(df$padj < 0.05, na.rm = TRUE)))
  df
}

res_age_df        <- run_pair("age",   "ctrl", "res_age_df")
res_intv1_df      <- run_pair("DR",    "age",  "res_intv1_df")
res_intv2_df      <- run_pair("sAct",  "age",  "res_intv2_df")
res_intv_df       <- run_pair("combi", "age",  "res_intv_df")
res_combi_ctrl_df <- run_pair("combi", "ctrl", "res_combi_ctrl_df")

message("Done. Outputs in ", normalizePath(out_dir))
