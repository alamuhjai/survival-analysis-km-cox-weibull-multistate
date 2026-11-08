weib2 = survreg(Surv(time, status) ~ age + sex + ph.ecog, data = lung, dist = "weibull")

if (requireNamespace("eha", quietly = TRUE)) {
  weib1 = eha::weibreg(Surv(time, status) ~ age + sex + ph.ecog, data = lung)
  bPh   = coef(weib1)[c("age", "sex", "ph.ecog")]
  sePh  = sqrt(diag(vcov(weib1)))[c("age", "sex", "ph.ecog")]
} else {
  sigmaAft = weib2$scale
  gammaPh  = 1 / sigmaAft
  bAft     = coef(weib2)[c("age", "sex", "ph.ecog")]
  bPh      = -bAft * gammaPh
  V        = vcov(weib2)
  # delta method SE for beta_PH = -beta_AFT / sigma
  sePh = sapply(seq_along(bAft), function(j) {
    grad = c(-gammaPh, -bPh[j])
    idx  = c(j + 1, nrow(V))
    sqrt(as.numeric(t(grad) %*% V[idx, idx] %*% grad))
  })
  names(sePh) = names(bPh)
  weib1 = list(coefficients = bPh, se = sePh, shape = gammaPh, scale = exp(coef(weib2)[1]),
               loglik = weib2$loglik)
}

sigmaAft = weib2$scale
gammaPh  = 1 / sigmaAft
bAft     = coef(weib2)[c("age", "sex", "ph.ecog")]
seAft    = sqrt(diag(vcov(weib2)))[c("age", "sex", "ph.ecog")]

deaths  = sum(lungData$event)
totTime = sum(lungData$time)
hz      = deaths / totTime
seLin   = sqrt(deaths) / totTime
ciLin   = (hz + c(-1, 1) * qnorm(0.975) * seLin) * 1000
ciLog   = (hz * exp(c(-1, 1) * qnorm(0.975) / sqrt(deaths))) * 1000

model1 = survreg(Surv(time, event) ~ age + sex + ph.ecog, data = lungData, dist = "exponential")
model2 = glm(event ~ age + sex + ph.ecog, offset = log(time), family = poisson, data = lungData)

bExpPh = -coef(model1)[c("age", "sex", "ph.ecog")]
bPois  =  coef(model2)[c("age", "sex", "ph.ecog")]

weibullDual = data.frame(
  covariate   = c("age (per year)", "sex (Female vs Male)", "ph.ecog (per 1-unit)"),
  betaPh      = round(bPh, 4),
  sePh        = round(sePh, 4),
  hrPh        = round(exp(bPh), 4),
  hrLow95     = round(exp(bPh - qnorm(0.975) * sePh), 4),
  hrUp95      = round(exp(bPh + qnorm(0.975) * sePh), 4),
  pPh         = signif(2 * pnorm(-abs(bPh / sePh)), 3),
  betaAft     = round(bAft, 4),
  seAft       = round(seAft, 4),
  afAft       = round(exp(bAft), 4),
  afLow95     = round(exp(bAft - qnorm(0.975) * seAft), 4),
  afUp95      = round(exp(bAft + qnorm(0.975) * seAft), 4),
  checkRatio  = round(-bPh / gammaPh, 4),
  poisGlmBeta = round(bPois, 4),
  poisGlmHr   = round(exp(bPois), 4),
  row.names   = NULL
)

cat("\n--- Study 2: lung Weibull PH vs AFT ---\n")
print(weibullDual)
write.csv(weibullDual, file.path(outDir, "weibullDual.csv"), row.names = FALSE)

png(file.path(figDir, "weibullFit.png"), width = 1800, height = 820, res = 160)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3.0, 1.5), bg = "white")

lungEcog1 = subset(lungData, ph.ecog == 1)
kmSex     = survfit(Surv(time, event) ~ sexFac, data = lungEcog1, conf.type = "log-log")
sexCols   = c("#d95f02", "#1b9e77")

plot(kmSex, col = sexCols, lwd = 2, lty = 2,
     xlab = "Time (days)", ylab = "Survival probability",
     main = "Weibull Fit vs Kaplan-Meier (lung)")
grid(col = "grey88")

tGrid   = seq(1, 1000, length.out = 300)
meanAge = mean(lungData$age)
for (sVal in c(1, 2)) {
  lpAft = coef(weib2)["(Intercept)"] + coef(weib2)["age"] * meanAge +
          coef(weib2)["sex"] * sVal + coef(weib2)["ph.ecog"] * 1
  sWeib = exp(-(tGrid / exp(lpAft))^gammaPh)
  lines(tGrid, sWeib, col = sexCols[sVal], lwd = 3)
}
legend("topright",
       legend = c("Male (KM)", "Female (KM)", "Male (Weibull)", "Female (Weibull)"),
       col = c(sexCols, sexCols), lty = c(2, 2, 1, 1), lwd = c(2, 2, 3, 3), cex = 0.88, bty = "n")

yPos   = c(3, 2, 1)
hrVals = weibullDual$hrPh
hrLow  = weibullDual$hrLow95
hrUp   = weibullDual$hrUp95
afVals = weibullDual$afAft
afLow  = weibullDual$afLow95
afUp   = weibullDual$afUp95

plot(NULL, xlim = c(0.35, 2.05), ylim = c(0.5, 3.5), yaxt = "n",
     xlab = "Effect size (HR vs AFT time factor)",
     ylab = "",
     main = "Weibull PH vs AFT Effects")
abline(v = 1, lty = 2, col = "grey40", lwd = 1.5)
grid(col = "grey88")
axis(2, at = yPos, labels = c("Age (+1 yr)", "Sex (Female)", "ECOG (+1 pt)"), las = 1, cex.axis = 0.9)

arrows(hrLow, yPos + 0.12, hrUp, yPos + 0.12, angle = 90, code = 3, length = 0.05, col = "#d95f02", lwd = 2.2)
points(hrVals, yPos + 0.12, pch = 19, col = "#d95f02", cex = 1.3)

arrows(afLow, yPos - 0.12, afUp, yPos - 0.12, angle = 90, code = 3, length = 0.05, col = "#1b9e77", lwd = 2.2)
points(afVals, yPos - 0.12, pch = 17, col = "#1b9e77", cex = 1.3)

legend("bottomright",
       legend = c("PH hazard ratio", "AFT time factor"),
       col = c("#d95f02", "#1b9e77"), pch = c(19, 17), lwd = 2.2, cex = 0.88, bty = "n")

dev.off()
