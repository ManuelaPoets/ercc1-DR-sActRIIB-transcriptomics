# Per-rgzj muscle CellChat full totals + cell counts (no CellChat pkg)
base <- "D:/manuela/epyc_aging/aging_sn_muscle/cellchat_objects"
cat("rgzj  cells  sum(net$weight)  sum(net$count)\n")
for (s in c("rgzj1","rgzj2","rgzj3","rgzj4","rgzj5")) {
  o <- readRDS(sprintf("%s/cellchat_%s_object.rds", base, s))
  net <- attr(o,"net")
  idents <- attr(o,"idents")
  nc <- if (is.factor(idents)) length(idents) else length(unlist(idents))
  cat(sprintf("%-6s %5d  %8.3f  %7d\n", s, nc, sum(net$weight), sum(net$count)))
}
cat("\nThesis CSV (named): ctrl 16.3/87, age 8.6/67, DR 9.0/68, sAct 13.4/100, combi 16.0/105\n")
cat("Physical: rgzj3=sAct, rgzj4=DR.  If rgzj3 weight > rgzj4 -> sAct>DR ranking, CSV 'sAct'=13.4 correctly = rgzj3.\n")
