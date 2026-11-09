coxUnstrat = coxph(Surv(time, status) ~ trtFac + priorFac + celltype + karno1, data = veteran)
zphUnstrat = cox.zph(coxUnstrat)

modelStrat = coxph(Surv(time, status) ~ trtFac + priorFac + strata(celltype, karno1), data = veteran)
zphStrat   = cox.zph(modelStrat)

fitCox = coxph(Surv(time, status) ~ trt * prior + strata(celltype, karno1), data = veteran)

v    = vcov(fitCox)
bTrt = unname(coef(fitCox)["trt"])
bPri = unname(coef(fitCox)["prior"])
bInt = unname(coef(fitCox)["trt:prior"])

seTrt0 = sqrt(v["trt", "trt"])
hrTrt0 = exp(bTrt)
ciTrt0 = exp(bTrt + c(-1, 1) * qnorm(0.975) * seTrt0)
pTrt0  = 2 * pnorm(-abs(bTrt / seTrt0))

# prior is coded 0/10 in veteran, scale interaction contrast by 10
logHrTrt10 = bTrt + 10 * bInt
varTrt10   = v["trt", "trt"] + 100 * v["trt:prior", "trt:prior"] + 20 * v["trt", "trt:prior"]
seTrt10    = sqrt(varTrt10)
hrTrt10    = exp(logHrTrt10)
ciTrt10    = exp(logHrTrt10 + c(-1, 1) * qnorm(0.975) * seTrt10)
pTrt10     = 2 * pnorm(-abs(logHrTrt10 / seTrt10))

sumStrat = summary(modelStrat)
bMain    = coef(modelStrat)
seMain   = sumStrat$coefficients[, "se(coef)"]
ciMain   = exp(confint(modelStrat))
pMain    = sumStrat$coefficients[, "Pr(>|z|)"]

coxCont = data.frame(
  model     = c("Stratified Main Effects", "Stratified Main Effects",
                "Stratified Interaction", "Stratified Interaction", "Stratified Interaction"),
  contrast  = c("Test vs Standard (trt)",
                "Prior Yes vs No (prior = 10 vs 0)",
                "Test vs Standard | Prior = No (0)",
                "Test vs Standard | Prior = Yes (10)",
                "Interaction ratio (Prior Yes / Prior No)"),
  logHr     = round(c(bMain[1], bMain[2], bTrt, logHrTrt10, 10 * bInt), 4),
  seLogHr   = round(c(seMain[1], seMain[2], seTrt0, seTrt10, 10 * sqrt(v["trt:prior", "trt:prior"])), 4),
  hr        = round(c(exp(bMain[1]), exp(bMain[2]), hrTrt0, hrTrt10, exp(10 * bInt)), 4),
  low95     = round(c(ciMain[1, 1], ciMain[2, 1], ciTrt0[1], ciTrt10[1],
                      exp(10 * bInt - qnorm(0.975) * 10 * sqrt(v["trt:prior", "trt:prior"]))), 4),
  up95      = round(c(ciMain[1, 2], ciMain[2, 2], ciTrt0[2], ciTrt10[2],
                      exp(10 * bInt + qnorm(0.975) * 10 * sqrt(v["trt:prior", "trt:prior"]))), 4),
  pValue    = round(c(pMain[1], pMain[2], pTrt0, pTrt10,
                      summary(fitCox)$coefficients["trt:prior", "Pr(>|z|)"]), 4),
  zphP      = c(round(zphStrat$table["trtFac", "p"], 4),
                round(zphStrat$table["priorFac", "p"], 4),
                NA, NA, NA),
  row.names = NULL
)

cat("\n--- veteran: unstratified and stratified Schoenfeld tests ---\n")
print(zphUnstrat)
print(zphStrat)

cat("\n--- veteran: stratified Cox contrasts ---\n")
print(coxCont)
write.csv(coxCont, file.path(outDir, "coxCont.csv"), row.names = FALSE)

cellCols  = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a")
karnoCols = c("#d95f02", "#1b9e77")

png(file.path(figDir, "phDiag.png"), width = 1800, height = 1400, res = 160)
par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3.0, 1.5), bg = "white")

plot(fitModel, fun = "cloglog", col = cellCols, lwd = 2.5,
     xlab = "Time (days, log scale)", ylab = "log(-log S(t))",
     main = "Log-Log Survival by Cell Type")
grid(col = "grey88")
legend("topleft", legend = levels(veteran$celltype), col = cellCols, lwd = 2.5, bty = "n")

plot(karno, fun = "cloglog", col = karnoCols, lwd = 2.5,
     xlab = "Time (days, log scale)", ylab = "log(-log S(t))",
     main = "Log-Log Survival by Karnofsky Score")
grid(col = "grey88")
legend("topleft", legend = levels(veteran$karno1), col = karnoCols, lwd = 2.5, bty = "n")

plot(zphStrat[1], lwd = 2, col = "#1b9e77",
     main = "Schoenfeld Residuals: Treatment")
abline(h = 0, lty = 2, col = "grey40")
grid(col = "grey88")

plot(zphStrat[2], lwd = 2, col = "#7570b3",
     main = "Schoenfeld Residuals: Prior Therapy")
abline(h = 0, lty = 2, col = "grey40")
grid(col = "grey88")

dev.off()

coxNull  = coxph(Surv(time, status) ~ 1, data = veteran)
martNull = residuals(coxNull, type = "martingale")
dfbStrat = residuals(modelStrat, type = "dfbeta")

png(file.path(figDir, "coxResid.png"), width = 1800, height = 820, res = 160)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3.0, 1.5), bg = "white")

plot(veteran$karno, martNull, pch = 19, col = rgb(0.2, 0.4, 0.6, 0.65),
     xlab = "Karnofsky score",
     ylab = "Martingale residual",
     main = "Martingale Residuals vs Karnofsky Score")
lines(lowess(veteran$karno, martNull, f = 0.6), col = "#d95f02", lwd = 3)
abline(h = 0, lty = 2, col = "grey45")
abline(v = 50, lty = 3, col = "#e7298a", lwd = 2)
grid(col = "grey88")

plot(seq_len(nrow(veteran)), dfbStrat[, 1], type = "h", lwd = 1.8, col = "#7570b3",
     xlab = "Patient index",
     ylab = "dfbeta (trt)",
     main = "dfbeta Influence for Treatment")
points(seq_len(nrow(veteran)), dfbStrat[, 1], pch = 19, cex = 0.7, col = "#7570b3")
thresh = 2 / sqrt(nrow(veteran)) * sd(dfbStrat[, 1])
abline(h = c(-thresh, 0, thresh), lty = c(3, 2, 3), col = c("#d95f02", "grey45", "#d95f02"))
grid(col = "grey88")

dev.off()
