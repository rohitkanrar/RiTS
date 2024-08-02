library(readxl)

region <- read_excel("data/Region.xlsx")
region <- region[, c(1, 3)]
colnames(region)[2] <- "region"

efficacy <- read_excel("data/Data_Efficacy.xlsx")
efficacy <- efficacy[efficacy$`Analysis Visit` == "Week 24", ]
efficacy <- efficacy[, c(1, 2, 11, 14)]
colnames(efficacy) <- c("Unique Subject Identifier", "arm",
                        "baseline_salt", "pct_change_salt")
efficacy <- efficacy[!is.na(efficacy$pct_change_salt), ]
efficacy <- as.data.frame(efficacy)

baseline <- read_excel("data/baseline.xlsx")
baseline <- baseline[, -c(2, 11, 15)]
colnames(baseline) <- c("Unique Subject Identifier", "height", "gender",
                        "weight", "race", "age", "bmi", "ethnicity", 
                        "duration" , "severity", "prior_trt", "smoke")
baseline <- as.data.frame(baseline)

retli_dat <- merge(efficacy, baseline, by = "Unique Subject Identifier", 
                   all = "left", all.x = TRUE, all.y = FALSE)
retli_dat <- merge(retli_dat, region, by = "Unique Subject Identifier", 
                   all = "left", all.x = TRUE, all.y = FALSE)

safety <- read_excel("data/Data_Safety.xlsx")
safety <- safety[safety$`Analysis Visit`=="Week 24", ]
safety <- safety[!is.na(safety$`Percent Change from Baseline`), ]
safety <- safety[, c(1, 11, 14)]
colnames(safety) <- c("Unique Subject Identifier", "baseline_lymph",
                      "pct_change_lymph")

lab_baseline <- read_excel("data/Lab_baseline.xlsx")
lab_baseline <- lab_baseline[lab_baseline$`Parameter Code` %in% 
                               c("L00025P", "L05024P"), ]
lab_baseline <- lab_baseline[, c("Unique Subject Identifier", "Parameter Code",
                                 "Baseline Value")]
lab_baseline <- as.data.frame(lab_baseline)
library(reshape2)
lab_baseline <- reshape(lab_baseline, idvar = "Unique Subject Identifier",
                        timevar = "Parameter Code", direction = "wide")
colnames(lab_baseline)[-1] <- c("albumin", "lymph_marker")


safety <- merge(safety, lab_baseline, by = "Unique Subject Identifier", 
                all = "left", all.x = TRUE, all.y = FALSE)

retli_dat <- merge(retli_dat, safety, by = "Unique Subject Identifier", 
                   all = "left", all.x = TRUE, all.y = FALSE)
retli_dat <- retli_dat[complete.cases(retli_dat), ]
lab <- c("Placebo" = "Arm 1", "Ritlecitinib 10 mg QD" = "Arm 2", 
         "Ritlecitinib 30 mg QD" = "Arm 3", "Ritlecitinib 50 mg QD" = "Arm 4",
         "Ritlecitinib 200/30 mg QD" = "Arm 5", 
         "Ritlecitinib 200/50 mg QD" = "Arm 6")
retli_dat[["arm"]] <- factor(retli_dat[["arm"]], levels = names(lab), 
                             labels = lab)
to_scale <- numeric(0)
for(j in 3:ncol(retli_dat)){
  if(colnames(retli_dat)[j] %in% c("arm", "pct_change_salt", "pct_change_lymph"))
    next
  if(is.character(retli_dat[, j])){
    retli_dat[, j] <- factor(retli_dat[, j])
  } else{
    to_scale <- c(to_scale, j)
  }
}
retli_dat[, to_scale] <- scale(retli_dat[, to_scale])
retli_dat[["pct_change_salt"]] <- -scale(retli_dat[["pct_change_salt"]]) 
retli_dat[["pct_change_lymph"]] <- scale(retli_dat[["pct_change_lymph"]]) 
retli_dat[["Unique Subject Identifier"]] <- NULL

saveRDS(retli_dat, "data/clean_data.RData")

rm(list = ls())
