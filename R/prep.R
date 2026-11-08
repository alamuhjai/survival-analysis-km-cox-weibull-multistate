library(survival)

figDir = "figures"
outDir = "results"
dir.create(figDir, showWarnings = FALSE, recursive = TRUE)
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)

data = veteran

# order karno1 levels so <= 50 is first; recode prior (0/10) and trt (1/2) as factors
veteran$karno1 = factor(ifelse(veteran$karno > 50, "karno > 50", "karno <= 50"),
                        levels = c("karno <= 50", "karno > 50"))
veteran$priorFac = factor(veteran$prior, levels = c(0, 10), labels = c("No", "Yes"))
veteran$trtFac   = factor(veteran$trt, levels = c(1, 2), labels = c("Standard", "Test"))

veteran$trtBin   = ifelse(veteran$trt == 2, 1, 0)
veteran$priorBin = ifelse(veteran$prior == 10, 1, 0)

lungData = subset(lung, !is.na(ph.ecog))
lungData$event  = as.integer(lungData$status == 2)
lungData$sexFac = factor(lungData$sex, levels = c(1, 2), labels = c("Male", "Female"))
