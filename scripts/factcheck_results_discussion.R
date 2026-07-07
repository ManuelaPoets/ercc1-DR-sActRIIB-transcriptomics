# Fact-check bulk numeric claims in Results/Discussion against archived res_age tables.
ka <- file.path(Sys.getenv("TEMP"), "kidney_archive", "kidney_res_age_archive.csv")
ma <- file.path(Sys.getenv("TEMP"), "muscle_archive", "muscle_v2_res_age_results_with_symbols.csv")
K <- read.csv(ka, stringsAsFactors = FALSE)
M <- read.csv(ma, stringsAsFactors = FALSE)
colnames(K)[1] <- "gene_id"; colnames(M)[1] <- "gene_id"
cat("kidney res_age cols:", paste(colnames(K), collapse=","), "\n")
cat("muscle res_age cols:", paste(colnames(M), collapse=","), "\n\n")

updown <- function(df, lab) {
  sig <- df[which(df$padj < 0.05), ]
  up <- sum(sig$log2FoldChange > 0); dn <- sum(sig$log2FoldChange < 0)
  cat(sprintf("%s: total sig=%d | up=%d | down=%d\n", lab, nrow(sig), up, dn))
}
cat("== up/down DEG counts (thesis: kidney 3,704/3,795 ; muscle 722/659) ==\n")
updown(K, "kidney")
updown(M, "muscle")

look <- function(df, sym) {
  i <- which(toupper(df$Symbol) == toupper(sym))
  if (!length(i)) return(c(NA, NA))
  i <- i[which.min(df$padj[i])]
  c(df$log2FoldChange[i], df$padj[i])
}
genes <- c("Tgfb1","Bax","Lmna","Stmn1","Hcls1","Maf","Pcbp2","Tnks2","Ctnna1","Ssr3","Eda2r","Gdf15","Ercc1")
cat("\n== cited gene LFCs (age vs ctrl) : kidney | muscle ==\n")
cat(sprintf("%-8s %10s %10s | %10s %10s\n","gene","K_LFC","K_padj","M_LFC","M_padj"))
for (g in genes) {
  k <- look(K, g); m <- look(M, g)
  cat(sprintf("%-8s %10.3f %10.2e | %10.3f %10.2e\n", g, k[1], k[2], m[1], m[2]))
}
cat("\nThesis-cited (shared-80, age vs ctrl ashr): Tgfb1 K+0.533/M+1.368 ; Bax K+0.579/M+0.657 ; Lmna K+0.504/M+0.672 ;\n")
cat(" Stmn1 K+0.356/M+1.452 ; Hcls1 K+0.769/M+1.026 ; Maf K-1.075/M-0.777 ; Pcbp2 K-0.557/M-0.366 ;\n")
cat(" Tnks2 K-0.514/M-0.325 ; Ctnna1 K+0.398/M+0.372 ; Ssr3 K-0.519/M-0.395 ; Eda2r K+4.04 ; Gdf15 M+1.97\n")
cat("(NB kidney archive LFC = MLE/un-shrunk, so magnitudes run larger than the ashr values cited; check SIGN + ballpark.)\n")
