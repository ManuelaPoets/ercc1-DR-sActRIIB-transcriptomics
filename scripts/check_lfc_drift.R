# Confirm: is rescued_to_normal residual due to ashr shrinkage drift near |LFC|<0.5?
arc_dir <- file.path(Sys.getenv("TEMP"), "kidney_archive")
r2n <- read.csv(file.path(arc_dir,"kidney_rescued_to_normal.csv"), stringsAsFactors=FALSE)
colnames(r2n)[1] <- "gene_id"
cc  <- read.csv("results/kidney_bulk/deseq2_out/res_combi_ctrl_df.csv", stringsAsFactors=FALSE)

m <- merge(r2n[,c("gene_id","log2FC_combi_ctrl")], cc[,c("gene_id","log2FoldChange")], by="gene_id")
colnames(m) <- c("gene_id","lfc_arc","lfc_new")
cat(sprintf("overlap genes (archive r2n): %d\n", nrow(m)))
cat(sprintf("Pearson r(LFC_cc): %.4f\n", cor(m$lfc_arc, m$lfc_new)))
cat(sprintf("mean |LFC| archive: %.4f | mean |LFC| new: %.4f  (new smaller => more shrunk)\n",
            mean(abs(m$lfc_arc)), mean(abs(m$lfc_new))))
cat(sprintf("median |LFC| archive: %.4f | new: %.4f\n", median(abs(m$lfc_arc)), median(abs(m$lfc_new))))
cat(sprintf("genes where |new| < |arc| (further shrunk): %d / %d (%.1f%%)\n",
            sum(abs(m$lfc_new) < abs(m$lfc_arc)), nrow(m),
            100*mean(abs(m$lfc_new) < abs(m$lfc_arc))))

# the 460 extras: my rescued-by-combi genes with |new LFC|<0.5 not in archive r2n.
# How close to the 0.5 boundary are their archive-equivalent magnitudes?
deg_dummy <- cc  # full combi_ctrl table for all genes
extra_band <- cc[which(abs(cc$log2FoldChange) >= 0.40 & abs(cc$log2FoldChange) < 0.50), ]
cat(sprintf("\nall genes with my |LFC_cc| in [0.40,0.50): %d (these are the boundary-sensitive band)\n",
            nrow(extra_band)))
