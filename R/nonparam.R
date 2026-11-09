model = Surv(veteran$time, veteran$status)

# ctype = 1 gives Nelson-Aalen, conf.type = "log-log" keeps CI in [0, 1]
fitModel = survfit(model ~ celltype, data = veteran, conf.type = "log-log", ctype = 1)
karno    = survfit(model ~ karno1,   data = veteran, conf.type = "log-log", ctype = 1)
fitTrt   = survfit(model ~ trtFac,   data = veteran, conf.type = "log-log", ctype = 1)

lrCell = survdiff(model ~ celltype, data = veteran, rho = 0)
ppCell = survdiff(model ~ celltype, data = veteran, rho = 1)
lrKarn = survdiff(model ~ karno1,   data = veteran, rho = 0)
ppKarn = survdiff(model ~ karno1,   data = veteran, rho = 1)
lrTrt  = survdiff(model ~ trtFac,   data = veteran, rho = 0)
ppTrt  = survdiff(model ~ trtFac,   data = veteran, rho = 1)

extractSurvTable = function(fitObj, groupLabel, lrObj, ppObj, tau = 365) {
  tab = as.data.frame(summary(fitObj, rmean = tau)$table)
  pLr = pchisq(lrObj$chisq, df = length(lrObj$n) - 1, lower.tail = FALSE)
  pPp = pchisq(ppObj$chisq, df = length(ppObj$n) - 1, lower.tail = FALSE)

  rmeanCol = grep("rmean", colnames(tab), value = TRUE)[1]
  seCol    = grep("se\\(rmean\\)", colnames(tab), value = TRUE)[1]

  data.frame(
    variable   = groupLabel,
    stratum    = sub("^[a-zA-Z0-9_.]+=", "", rownames(tab)),
    n          = tab$records,
    events     = tab$events,
    median     = tab$median,
    medLow95   = tab[, "0.95LCL"],
    medUp95    = tab[, "0.95UCL"],
    rmst365    = round(tab[, rmeanCol], 2),
    rmstSe     = round(tab[, seCol], 2),
    rmstLow95  = round(tab[, rmeanCol] - 1.96 * tab[, seCol], 2),
    rmstUp95   = round(tab[, rmeanCol] + 1.96 * tab[, seCol], 2),
    logRankChi = round(lrObj$chisq, 2),
    logRankP   = signif(pLr, 3),
    petoChi    = round(ppObj$chisq, 2),
    petoP      = signif(pPp, 3),
    row.names  = NULL
  )
}

rmstSum = rbind(
  extractSurvTable(fitModel, "celltype", lrCell, ppCell),
  extractSurvTable(karno,    "karno1",   lrKarn, ppKarn),
  extractSurvTable(fitTrt,   "trt",      lrTrt,  ppTrt)
)

cat("\n--- veteran: non-parametric survival and RMST ---\n")
print(rmstSum)
write.csv(rmstSum, file.path(outDir, "rmstSum.csv"), row.names = FALSE)

cellCols = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a")
cellLabs = levels(veteran$celltype)

png(file.path(figDir, "kmCurve.png"), width = 1800, height = 820, res = 160)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3.0, 1.5), bg = "white")

plot(fitModel, col = cellCols, lwd = 2.5,
     xlab = "Time (days)",
     ylab = "Survival probability",
     main = "Kaplan-Meier Survival by Cell Type")
abline(h = 0.5, lty = 3, col = "grey45")
grid(col = "grey88")
legend("topright", legend = cellLabs, col = cellCols, lwd = 2.5, cex = 0.9, bty = "n")

plot(fitModel, fun = "cumhaz", col = cellCols, lwd = 2.5,
     xlab = "Time (days)",
     ylab = "Cumulative hazard",
     main = "Nelson-Aalen Cumulative Hazard")
grid(col = "grey88")
legend("topleft", legend = cellLabs, col = cellCols, lwd = 2.5, cex = 0.9, bty = "n")

dev.off()
