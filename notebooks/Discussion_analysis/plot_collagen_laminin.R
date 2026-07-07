library(CellChat)

# ---- Adjust this path to your server ----
rds_dir <- "/data/bonn-epyc/projects/manuela/aging/aging_interorgan/cellchat_objects/"

cellchat.ctrl  <- readRDS(paste0(rds_dir, "cellchat_ctrl_centr_sub.rds"))
cellchat.age   <- readRDS(paste0(rds_dir, "cellchat_age_centr_sub.rds"))
cellchat.DR    <- readRDS(paste0(rds_dir, "cellchat_DR_centr_sub.rds"))
cellchat.sAct  <- readRDS(paste0(rds_dir, "cellchat_sAct_centr_sub.rds"))
cellchat.combi <- readRDS(paste0(rds_dir, "cellchat_combi_centr_sub.rds"))

object.list <- list(
  ctrl  = cellchat.ctrl,
  age   = cellchat.age,
  DR    = cellchat.DR,
  sAct  = cellchat.sAct,
  combi = cellchat.combi
)

out_dir <- "./interorgan_missing_plots/"
dir.create(out_dir, showWarnings = FALSE)

pathways <- c("COLLAGEN", "LAMININ")

for (pathway in pathways) {
  has_pw <- sapply(names(object.list), function(c) pathway %in% object.list[[c]]@netP$pathways)
  cat(pathway, "detected in:", paste(names(has_pw)[has_pw], collapse=", "), "\n")
  cat(pathway, "absent in:",   paste(names(has_pw)[!has_pw], collapse=", "), "\n")

  if (sum(has_pw) == 0) { cat("Skipping\n\n"); next }

  png(paste0(out_dir, pathway, "_circle_all_conditions.png"),
      width = 2400, height = 1600, bg = "white", res = 150)
  par(mfrow = c(2, 3), xpd = TRUE)
  for (i in 1:length(object.list)) {
    nm <- names(object.list)[i]
    if (has_pw[nm]) {
      netVisual_aggregate(object.list[[i]], signaling = pathway,
                          signaling.name = paste(pathway, nm))
    } else {
      plot.new(); title(main = paste(pathway, nm, "\n(not detected)"))
    }
  }
  dev.off()
  cat("Saved:", pathway, "\n\n")
}

cat("\n=== Probability summary ===\n")
for (pw in pathways) {
  cat("\n", pw, ":\n")
  for (nm in names(object.list)) {
    obj <- object.list[[nm]]
    if (pw %in% obj@netP$pathways) {
      cat("  ", nm, ":", round(sum(obj@netP$prob[,,pw]), 4), "\n")
    } else {
      cat("  ", nm, ": 0.000 (not detected)\n")
    }
  }
}
cat("\nDone.\n")
