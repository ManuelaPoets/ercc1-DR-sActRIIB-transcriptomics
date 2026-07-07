# Inspect a kidney CellChat object WITHOUT loading CellChat (S4 slots = attributes)
k <- readRDS("D:/manuela/epyc_aging/snRNA_aging-mouse/cellchat_objects/cellchat_sg18_object.rds")
cat("attribute (slot) names:", paste(names(attributes(k)), collapse=", "), "\n")
net  <- attr(k, "net")
netP <- attr(k, "netP")
cat("net names:", paste(names(net), collapse=", "), "\n")
cat("sum(net$weight) =", round(sum(net$weight),3), "\n")
cat("sum(net$count)  =", sum(net$count), "\n")
cat("netP names:", paste(names(netP), collapse=", "), "\n")
pw <- dimnames(netP$prob)[[3]]
cat("n pathways:", length(pw), "\n")
cat("pathways:", paste(pw, collapse=", "), "\n")
for (p in c("VEGF","EGF","SEMA3","EPHA","EPHB","PDGF")) {
  if (p %in% pw) cat(sprintf("  %-7s total prob = %.4f\n", p, sum(netP$prob[,,p])))
}
