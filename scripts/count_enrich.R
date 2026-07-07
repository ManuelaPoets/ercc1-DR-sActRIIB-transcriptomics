dir <- file.path(Sys.getenv("TEMP"), "enrich_check")
files <- c(k_aging="k_aging.csv", k_rescued="k_resc.csv", k_nonrescued="k_nonresc.csv", k_DR="k_dr.csv",
           m_aging="m_aging.csv", m_rescued="m_resc.csv", m_nonrescued="m_nonresc.csv", m_DR="m_dr.csv")
exp <- c(k_aging=2757, k_rescued=703, k_nonrescued=195, k_DR=276, m_aging=244, m_rescued=17, m_nonrescued=45, m_DR=105)
for (n in names(files)) {
  p <- file.path(dir, files[[n]])
  d <- tryCatch(read.csv(p, stringsAsFactors=FALSE), error=function(e) NULL)
  if (is.null(d)) { cat(sprintf("%-14s ERROR\n", n)); next }
  pc <- if ("p.adjust" %in% colnames(d)) "p.adjust" else if ("padj" %in% colnames(d)) "padj" else NA
  sig <- if (!is.na(pc)) sum(d[[pc]] < 0.05, na.rm=TRUE) else NA
  cat(sprintf("%-14s rows=%-5d sig=%-5d thesis=%-5d %s\n", n, nrow(d), sig, exp[[n]],
              ifelse(!is.na(sig) && sig==exp[[n]], "MATCH", "DIFF")))
}
