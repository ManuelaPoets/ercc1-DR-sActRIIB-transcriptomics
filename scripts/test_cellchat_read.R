d <- "notebooks/aging_interorgan/cellchat_objects"
obj <- readRDS(file.path(d, "cellchat_ctrl_object.rds"))
cat("attr names:", paste(names(attributes(obj)), collapse=", "), "\n")
net <- attr(obj, "net")
cat("net is null:", is.null(net), "\n")
if (!is.null(net)) {
  cat("net names:", paste(names(net), collapse=", "), "\n")
  cat("sum count:", sum(net$count), " | sum weight:", round(sum(net$weight),3), "\n")
  cat("dim count:", paste(dim(net$count), collapse="x"), "\n")
}
