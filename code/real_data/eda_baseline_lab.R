library(reshape2); library(readxl)
real_dat <- readRDS("data/clean_data.RData")
lab_baseline <- read_excel("data/Lab_baseline.xlsx")
lab_baseline <- lab_baseline[, c("Unique Subject Identifier",
                                "Parameter", "Baseline Value")]
lab_baseline <- as.data.frame(lab_baseline)
lab_baseline_wide <- reshape(lab_baseline, idvar = "Unique Subject Identifier",
                             timevar = "Parameter", direction = "wide")

too_many_miss <- which(sapply(1:ncol(lab_baseline_wide), function(j){
  sum(is.na(lab_baseline_wide[, j]))
}) > 50)
lab_baseline_wide <- lab_baseline_wide[, -too_many_miss]
lab_baseline_wide <- lab_baseline_wide[complete.cases(lab_baseline_wide), ]

safety <- read_excel("data/Data_Safety.xlsx")
safety <- safety[safety$`Analysis Visit`=="Week 24", ]
safety <- safety[!is.na(safety$`Percent Change from Baseline`), ]
safety <- safety[, c(1, 11, 14)]
colnames(safety) <- c("Unique Subject Identifier", "baseline_lymph",
                      "pct_change_lymph")
safety_lab <- merge(safety, lab_baseline_wide, by = "Unique Subject Identifier", 
                    all = "left", all.x = TRUE, all.y = FALSE)
safety_lab <- safety_lab[complete.cases(safety_lab), ]
safety_lab <- safety_lab[, -c(1)]
safety_lab <- scale(safety_lab)

library(glmnet)

mod <- cv.glmnet(x = as.matrix(safety_lab[, -2]), y = safety_lab[, 2], 
                 family = "gaussian")
mod_coef <- coef(mod, s = mod$lambda.min)
print(mod_coef@Dimnames[[1]][mod_coef@i+1][which(abs(mod_coef@x) > 0.1)])
# This supports the inclusion of albumin and lymph_marker in safety model. 