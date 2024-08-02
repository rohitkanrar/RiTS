real_dat <- readRDS("data/clean_data.RData")
K <- length(unique(real_dat$arm))

mod_eff <- lm(pct_change_salt ~ ., subset(real_dat, select = -pct_change_lymph))
summary(mod_eff)

mod_safe <- lm(pct_change_lymph ~ ., subset(real_dat, select = -pct_change_salt))
summary(mod_safe)
