time = seq(0, 3, length = 300)

# direct (1 -> 3) vs illness-mediated (1 -> 2 -> 3) mortality
probFrom1  = (5 / 6) * (1 - exp(-1.2 * time))
probAlpha1 = (1 / 6) * (1 - exp(-1.2 * time)) + (exp(-1.2 * time) - exp(-1 * time))
probAlpha2 = (1 / 6) * (1 - exp(-1.2 * time)) - (1 / 4) * (exp(-1.2 * time) - exp(-2 * time))
probAlpha5 = (1 / 6) * (1 - exp(-1.2 * time)) - (1 / 19) * (exp(-1.2 * time) - exp(-5 * time))

transCounts = matrix(c(
  395, 400, 280,
    1,   1,   1,
 1040, 500, 320,
    2,   1,   3
), nrow = 4, byrow = TRUE,
dimnames = list(c("1->2", "1->3", "2->1", "2->3"),
                c("Yr1", "Yr2", "Yr3")))

timeAtRisk = matrix(c(
   980,  800,  700,
  2500, 2600, 2500
), nrow = 2, byrow = TRUE,
dimnames = list(c("State1", "State2"), c("Yr1", "Yr2", "Yr3")))

stateMap   = c(1, 1, 2, 2)
lambdaFull = transCounts / timeAtRisk[stateMap, ]
lambdaPool = rowSums(transCounts) / rowSums(timeAtRisk[stateMap, ])
seFull     = lambdaFull / sqrt(transCounts)

llKernel = function(n, S) ifelse(n == 0, 0, n * log(n / S)) - n
llFull   = sum(llKernel(transCounts, timeAtRisk[stateMap, ]))
llRed    = sum(llKernel(rowSums(transCounts), rowSums(timeAtRisk[stateMap, ])))
W        = 2 * (llFull - llRed)
df       = (ncol(transCounts) - 1) * nrow(transCounts)
pval     = pchisq(W, df, lower.tail = FALSE)

contrib = sapply(1:4, function(i) {
  full = sum(llKernel(transCounts[i, ], timeAtRisk[stateMap[i], ]))
  red  = llKernel(sum(transCounts[i, ]), sum(timeAtRisk[stateMap[i], ]))
  2 * (full - red)
})
names(contrib) = rownames(transCounts)

S1Tot   = sum(timeAtRisk[1, ])
S2Tot   = sum(timeAtRisk[2, ])
lam12   = lambdaPool["1->2"]; lam13 = lambdaPool["1->3"]
lam21   = lambdaPool["2->1"]; lam23 = lambdaPool["2->3"]
lamOut1 = lam12 + lam13
lamOut2 = lam21 + lam23
varSum1 = lam12 / S1Tot + lam13 / S1Tot
varSum2 = lam21 / S2Tot + lam23 / S2Tot

mu1 = 1 / lamOut1
mu2 = 1 / lamOut2
se1 = sqrt(varSum1 / lamOut1^4)
se2 = sqrt(varSum2 / lamOut2^4)
z   = qnorm(0.975)

# 3-year piecewise-constant restricted sojourn time mu_k(3) and delta-method SE
lamOut1Yr = lambdaFull["1->2", ] + lambdaFull["1->3", ]
lamOut2Yr = lambdaFull["2->1", ] + lambdaFull["2->3", ]
varOut1Yr = lambdaFull["1->2", ]^2 / transCounts["1->2", ] + lambdaFull["1->3", ]^2 / transCounts["1->3", ]
varOut2Yr = lambdaFull["2->1", ]^2 / transCounts["2->1", ] + lambdaFull["2->3", ]^2 / transCounts["2->3", ]

calcPiecewiseSoj = function(l, v) {
  l1 = l[1]; l2 = l[2]; l3 = l[3]
  t1 = (1 - exp(-l1)) / l1
  t2 = exp(-l1) * (1 - exp(-l2)) / l2
  t3 = exp(-l1 - l2) * (1 - exp(-l3)) / l3
  mu = t1 + t2 + t3

  d1 = (l1 * exp(-l1) - (1 - exp(-l1))) / l1^2 - t2 - t3
  d2 = exp(-l1) * (l2 * exp(-l2) - (1 - exp(-l2))) / l2^2 - t3
  d3 = exp(-l1 - l2) * (l3 * exp(-l3) - (1 - exp(-l3))) / l3^2
  se = sqrt(sum(c(d1, d2, d3)^2 * v))
  c(mu = unname(mu), se = unname(se))
}

pw1 = calcPiecewiseSoj(lamOut1Yr, varOut1Yr)
pw2 = calcPiecewiseSoj(lamOut2Yr, varOut2Yr)

matExp3 = function(Q, dt) {
  M    = diag(nrow(Q))
  term = diag(nrow(Q))
  for (m in 1:25) {
    term = (term %*% (Q * dt)) / m
    M    = M + term
  }
  M
}

makeQ = function(r) {
  l12 = lambdaFull["1->2", r]; l13 = lambdaFull["1->3", r]
  l21 = lambdaFull["2->1", r]; l23 = lambdaFull["2->3", r]
  matrix(c(-(l12 + l13), l12, l13,
           l21, -(l21 + l23), l23,
           0, 0, 0), nrow = 3, byrow = TRUE)
}

Q1  = makeQ(1); Q2 = makeQ(2); Q3 = makeQ(3)
P01 = matExp3(Q1, 1)
P02 = P01 %*% matExp3(Q2, 1)

uGrid = seq(0, 1, length.out = 201)
du    = uGrid[2] - uGrid[1]
wTrap = c(0.5, rep(1, length(uGrid) - 2), 0.5) * du

occ1 = 0; occ2 = 0
for (idx in seq_along(uGrid)) {
  u    = uGrid[idx]
  P_u1 = matExp3(Q1, u)
  P_u2 = P01 %*% matExp3(Q2, u)
  P_u3 = P02 %*% matExp3(Q3, u)
  occ1 = occ1 + wTrap[idx] * (P_u1[1, 1] + P_u2[1, 1] + P_u3[1, 1])
  occ2 = occ2 + wTrap[idx] * (P_u1[2, 2] + P_u2[2, 2] + P_u3[2, 2])
}

stateSoj = data.frame(
  quantity  = c(rownames(transCounts),
                "State 1 Sojourn (Piecewise 3-yr)", "State 2 Sojourn (Piecewise 3-yr)",
                "State 1 Sojourn (Pooled inf-horizon)", "State 2 Sojourn (Pooled inf-horizon)",
                "State 1 Total Occupation Psi_11(3)", "State 2 Total Occupation Psi_22(3)"),
  yr1       = c(round(lambdaFull[, 1], 4), rep(NA, 6)),
  yr2       = c(round(lambdaFull[, 2], 4), rep(NA, 6)),
  yr3       = c(round(lambdaFull[, 3], 4), rep(NA, 6)),
  estimate  = round(c(lambdaPool, pw1["mu"], pw2["mu"], mu1, mu2, occ1, occ2), 4),
  se        = round(c(lambdaPool / sqrt(rowSums(transCounts)), pw1["se"], pw2["se"], se1, se2, NA, NA), 4),
  low95     = round(c(lambdaPool - z * (lambdaPool / sqrt(rowSums(transCounts))),
                      pw1["mu"] - z * pw1["se"], pw2["mu"] - z * pw2["se"],
                      mu1 - z * se1, mu2 - z * se2, NA, NA), 4),
  up95      = round(c(lambdaPool + z * (lambdaPool / sqrt(rowSums(transCounts))),
                      pw1["mu"] + z * pw1["se"], pw2["mu"] + z * pw2["se"],
                      mu1 + z * se1, mu2 + z * se2, NA, NA), 4),
  lrtStatW  = c(round(contrib, 3), rep(round(W, 3), 6)),
  row.names = NULL
)

cat("\n--- Study 3: multistate transition summary ---\n")
print(stateSoj)
write.csv(stateSoj, file.path(outDir, "stateSoj.csv"), row.names = FALSE)

png(file.path(figDir, "stateOcc.png"), width = 1800, height = 820, res = 160)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3.0, 1.5), bg = "white")

plot(time, probFrom1, type = "l", col = "#1f78b4", ylim = c(0, 1), lwd = 3,
     xlab = "Time (years)", ylab = "Probability",
     main = "Mortality by Pathway (State 3)")
lines(time, probAlpha1, col = "#e31a1c", lwd = 2.5, lty = 2)
lines(time, probAlpha2, col = "#ff7f00", lwd = 2.5, lty = 3)
lines(time, probAlpha5, col = "#6a3d9a", lwd = 2.5, lty = 4)
grid(col = "grey88")
legend("topleft",
       legend = c("Direct (1 -> 3)",
                  "Via State 2 (alpha = 1)",
                  "Via State 2 (alpha = 2)",
                  "Via State 2 (alpha = 5)"),
       col = c("#1f78b4", "#e31a1c", "#ff7f00", "#6a3d9a"),
       lwd = c(3, 2.5, 2.5, 2.5), lty = c(1, 2, 3, 4), cex = 0.88, bty = "n")

yrs = 1:3
r12 = lambdaFull["1->2", ]; se12 = seFull["1->2", ]
r21 = lambdaFull["2->1", ]; se21 = seFull["2->1", ]

plot(yrs - 0.04, r12, type = "b", pch = 19, col = "#d95f02", lwd = 2.5,
     xlim = c(0.8, 3.2), ylim = c(0.05, 0.60), xaxt = "n",
     xlab = "Year",
     ylab = "Transition rate (per year)",
     main = "Annual Transition Rates by State")
axis(1, at = 1:3, labels = c("Year 1", "Year 2", "Year 3"))
grid(col = "grey88")

arrows(yrs - 0.04, r12 - z * se12, yrs - 0.04, r12 + z * se12,
       angle = 90, code = 3, length = 0.05, col = "#d95f02", lwd = 2)

lines(yrs + 0.04, r21, type = "b", pch = 17, col = "#1b9e77", lwd = 2.5, lty = 2)
arrows(yrs + 0.04, r21 - z * se21, yrs + 0.04, r21 + z * se21,
       angle = 90, code = 3, length = 0.05, col = "#1b9e77", lwd = 2)

abline(h = lambdaPool["1->2"], lty = 3, col = "#d95f02")
abline(h = lambdaPool["2->1"], lty = 3, col = "#1b9e77")

legend("topright",
       legend = c("Illness onset (1 -> 2)",
                  "Recovery (2 -> 1)",
                  "Pooled rate"),
       col = c("#d95f02", "#1b9e77", "grey40"),
       pch = c(19, 17, NA), lty = c(1, 2, 3), lwd = 2.2, cex = 0.88, bty = "n")

dev.off()
