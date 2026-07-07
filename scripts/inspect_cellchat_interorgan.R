d <- "notebooks/aging_interorgan/cellchat_objects"
conds <- c("ctrl","age","DR","sAct","combi")
cat("=== Inter-organ 9-cell-type subset (_sub) per condition ===\n")
for (c in conds) {
  p <- file.path(d, paste0("cellchat_", c, "_object_sub.rds"))
  if (file.exists(p)) {
    s <- readRDS(p); n <- attr(s, "net")
    cat(sprintf("%-6s -> %d ct | sum weight = %.2f | count = %d\n",
                c, nrow(n$count), sum(n$weight), sum(n$count)))
  } else cat(c, "MISSING\n")
}
cat("\n=== Full merged (_object) per condition ===\n")
for (c in conds) {
  p <- file.path(d, paste0("cellchat_", c, "_object.rds"))
  if (file.exists(p)) {
    s <- readRDS(p); n <- attr(s, "net")
    cat(sprintf("%-6s -> %d ct | sum weight = %.2f | count = %d\n",
                c, nrow(n$count), sum(n$weight), sum(n$count)))
  } else cat(c, "MISSING\n")
}
