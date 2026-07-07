d <- "notebooks/aging_interorgan/cellchat_objects"
o <- readRDS(file.path(d, "cellchat_ctrl_object.rds"))
net <- attr(o, "net")
cat("full ctrl object:", nrow(net$count), "cell types | sum weight =", round(sum(net$weight),2), "| count =", sum(net$count), "\n")
cat(paste(rownames(net$count), collapse=" | "), "\n\n")
for (fn in c("cellchat_ctrl_object_sub.rds","cellchat_ctrl_object_sub_imm.rds")) {
  p <- file.path(d, fn)
  if (file.exists(p)) {
    s <- readRDS(p); n <- attr(s, "net")
    cat(fn, "->", nrow(n$count), "cell types | sum weight =", round(sum(n$weight),2), "| count =", sum(n$count), "\n")
    cat("   ", paste(rownames(n$count), collapse=", "), "\n")
  }
}
