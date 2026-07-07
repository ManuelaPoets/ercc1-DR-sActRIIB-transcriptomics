d <- "notebooks/aging_interorgan/cellchat_objects"
conds <- c("ctrl","age","DR","sAct","combi")
want <- c("SEMA3","PECAM1","CADM","CDH5","THBS","COLLAGEN")

# collect per-pathway total interaction probability per condition
tab <- list()
npaths <- c()
for (cc in conds) {
  p <- file.path(d, paste0("cellchat_", cc, "_object_sub.rds"))
  o <- readRDS(p)
  netP <- attr(o, "netP")
  prob <- netP$prob                      # [sender, receiver, pathway]
  paths <- dimnames(prob)[[3]]
  tot <- apply(prob, 3, sum)             # total prob per pathway
  names(tot) <- paths
  tab[[cc]] <- tot
  npaths[cc] <- length(paths)
}

allpaths <- sort(unique(unlist(lapply(tab, names))))
cat("Pathways detected per condition:\n"); print(npaths)
cat("Union of pathways across conditions:", length(allpaths), "\n\n")

getv <- function(cc, pw) { v <- tab[[cc]][pw]; if (is.na(v)) 0 else v }

cat(sprintf("%-9s %8s %8s %8s %8s %8s\n","pathway","ctrl","age","DR","sAct","combi"))
for (pw in want) {
  cat(sprintf("%-9s %8.3f %8.3f %8.3f %8.3f %8.3f\n", pw,
      getv("ctrl",pw), getv("age",pw), getv("DR",pw), getv("sAct",pw), getv("combi",pw)))
}

# how many pathways restored by combi (combi closer to ctrl than age was, given aging decreased it)
cat("\n# pathways where combi > age AND aging decreased (ctrl>age):\n")
restored <- 0
for (pw in allpaths) {
  ctv <- getv("ctrl",pw); agv <- getv("age",pw); cbv <- getv("combi",pw)
  if (ctv > agv && cbv > agv) restored <- restored + 1
}
cat(restored, "of", length(allpaths), "\n")
