# ============================================================================
# Kidney bulk RNA-seq: Salmon quant.sf -> tximport -> DESeq2
# Honey / Thesis. Reads raw quant from results/kidney_bulk/.
# Conventions: ctrl reference; ashr shrinkage; cleaned *_df outputs.
# Standard 5-condition aging design (young intervention arms excluded).
# ============================================================================

suppressMessages({
  library(tximport)
  library(DESeq2)
  library(ashr)
})

# ---- Paths -----------------------------------------------------------------
data_dir <- "results/kidney_bulk"                 # run from thesis root
out_dir  <- file.path(data_dir, "deseq2_out")
dir.create(out_dir, showWarnings = FALSE)

# ---- Sample sheet ----------------------------------------------------------
coldata <- read.csv(file.path(data_dir, "Activin_kidney_info_selection.csv"),
                    stringsAsFactors = FALSE, check.names = FALSE)
coldata$condition <- trimws(coldata$condition)

# Keep only the 5 aging conditions we have quant files for
keep_cond <- c("ctrl", "age", "DR", "sAct", "combi")
coldata   <- coldata[coldata$condition %in% keep_cond, ]

# Match to available quant.sf files (subset + order together)
coldata$file <- file.path(data_dir, paste0(coldata$sample_name, "_quant.sf"))
coldata      <- coldata[file.exists(coldata$file), ]
stopifnot(nrow(coldata) == 15L)

rownames(coldata)   <- coldata$sample_name
coldata$condition   <- factor(coldata$condition, levels = keep_cond)  # ctrl = ref

files        <- coldata$file
names(files) <- coldata$sample_name

# ---- tx2gene (transcript -> gene_id; symbol kept for annotation) -----------
t2g <- read.delim(file.path(data_dir, "salmon_tx2gene.tsv"),
                  header = FALSE, stringsAsFactors = FALSE)
colnames(t2g) <- c("tx", "gene_id", "gene_symbol")
tx2gene <- t2g[, c("tx", "gene_id")]
gene_sym <- unique(t2g[, c("gene_id", "gene_symbol")])
gene_sym <- gene_sym[!duplicated(gene_sym$gene_id), ]

# ---- Import (gene-level) ----------------------------------------------------
txi <- tximport(files, type = "salmon", tx2gene = tx2gene,
                ignoreTxVersion = FALSE)

# ---- DESeq2 -----------------------------------------------------------------
dds <- DESeqDataSetFromTximport(txi, colData = coldata, design = ~condition)
dds <- estimateSizeFactors(dds)
dds <- dds[rowSums(counts(dds)) >= 10, ]   # original notebook pre-filter (>=10 reads)
dds <- DESeq(dds)

# ---- Contrast helper: ashr-shrunk, cleaned data.frame ----------------------
get_df <- function(dds, num, den, name) {
  res <- results(dds, contrast = c("condition", num, den))
  res <- lfcShrink(dds, contrast = c("condition", num, den),
                   type = "ashr", res = res)
  df  <- as.data.frame(res)
  df$gene_id     <- rownames(df)
  df$gene_symbol <- gene_sym$gene_symbol[match(df$gene_id, gene_sym$gene_id)]
  df  <- df[order(df$padj), ]
  df  <- df[, c("gene_id", "gene_symbol",
                "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj")]
  write.csv(df, file.path(out_dir, paste0(name, ".csv")), row.names = FALSE)
  message(sprintf("%-16s %s vs %s | padj<0.05: %d",
                  name, num, den, sum(df$padj < 0.05, na.rm = TRUE)))
  df
}

# ---- Standard five contrasts (CLAUDE.md) -----------------------------------
res_age_df        <- get_df(dds, "age",   "ctrl", "res_age_df")        # aging effect
res_intv1_df      <- get_df(dds, "DR",    "age",  "res_intv1_df")      # DR effect
res_intv2_df      <- get_df(dds, "sAct",  "age",  "res_intv2_df")      # sActRIIB effect
res_intv_df       <- get_df(dds, "combi", "age",  "res_intv_df")       # PRIMARY rescue
res_combi_ctrl_df <- get_df(dds, "combi", "ctrl", "res_combi_ctrl_df") # restoration check

saveRDS(dds, file.path(out_dir, "dds_kidney_bulk.rds"))
message("Done. Outputs in ", normalizePath(out_dir))
