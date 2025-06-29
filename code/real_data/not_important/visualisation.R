ts_real <- readRDS("output/ts_real.RData")
rits_real <- readRDS("output/rits_real.RData")
rand_real <- readRDS("output/rand_real.RData")
real_dat <- readRDS("data/clean_data.RData")
source("code/real_data/application_function.R")

library(ggplot2)
library(patchwork)
library(reshape2)
N <- length(ts_real$trt)

# regret utility
df <- data.frame(patient = 1:N, ts = cumsum(ts_real$regret),
                 rand = cumsum(rand_real$regret), 
                 rits = cumsum(rits_real$regret))
df_long <- tidyr::gather(df, key = "strategy", value = "value", -patient)

# Plot using ggplot2
plot1 <- ggplot(df_long, aes(x = patient, y = value, color = strategy)) +
  geom_line() +
  geom_vline(xintercept = ts_real$tr_first, linetype = "dashed", color = "blue") +
  labs(title = "Utility Regret",
       x = "Patient", y = "Cumulative Regret") +
  theme(text = element_text(size = 8))

df <- data.frame(patient = 1:N, ts = cumsum(ts_real$subopt_trt),
                 rand = cumsum(rand_real$subopt_trt),
                 rits = cumsum(rits_real$subopt_trt))
df_long <- tidyr::gather(df, key = "strategy", value = "value", -patient)

# Plot using ggplot2
plot2 <- ggplot(df_long, aes(x = patient, y = value, color = strategy)) +
  geom_line() +
  geom_vline(xintercept = ts_real$tr_first, linetype = "dashed", color = "blue") +
  labs(title = "Utility Sub-optimal Treatment Count",
       x = "Patient", y = "Cumulative Sub-optimal Treatment Count") +
  theme(text = element_text(size = 8))

plot12 <- plot1 + plot2 + plot_layout(ncol = 2, guides = "collect")
ggsave("plot/regret_all_real.jpg", height = 3, width = 6, units = "in")


# regret benefit
df <- data.frame(patient = 1:N, ts = cumsum(ts_real$regret_benf),
                 rand = cumsum(rand_real$regret_benf), 
                 rits = cumsum(rits_real$regret_benf))
df_long <- tidyr::gather(df, key = "strategy", value = "value", -patient)

# Plot using ggplot2
plot1 <- ggplot(df_long, aes(x = patient, y = value, color = strategy)) +
  geom_line() +
  geom_vline(xintercept = ts_real$tr_first, linetype = "dashed", color = "blue") +
  labs(title = "Efficacy Regret",
       x = "Patient", y = "Cumulative Regret") +
  theme(text = element_text(size = 8))

df <- data.frame(patient = 1:N, ts = cumsum(ts_real$subopt_trt_benf),
                 rand = cumsum(rand_real$subopt_trt_benf),
                 rits = cumsum(rits_real$subopt_trt_benf))
df_long <- tidyr::gather(df, key = "strategy", value = "value", -patient)

# Plot using ggplot2
plot2 <- ggplot(df_long, aes(x = patient, y = value, color = strategy)) +
  geom_line() +
  geom_vline(xintercept = ts_real$tr_first, linetype = "dashed", color = "blue") +
  labs(title = "Efficacy Sub-optimal Treatment Count",
       x = "Patient", y = "Cumulative Sub-optimal Treatment Count") +
  theme(text = element_text(size = 8))

plot12 <- plot1 + plot2 + plot_layout(ncol = 2, guides = "collect")
ggsave("plot/regret_benf_real.jpg", height = 3, width = 6, units = "in")

# regret safety
df <- data.frame(patient = 1:N, ts = cumsum(ts_real$regret_safe),
                 rand = cumsum(rand_real$regret_safe), 
                 rits = cumsum(rits_real$regret_safe))
df_long <- tidyr::gather(df, key = "strategy", value = "value", -patient)

# Plot using ggplot2
plot1 <- ggplot(df_long, aes(x = patient, y = value, color = strategy)) +
  geom_line() +
  geom_vline(xintercept = ts_real$tr_first, linetype = "dashed", color = "blue") +
  labs(title = "Safety Regret",
       x = "Patient", y = "Cumulative Regret") +
  theme(text = element_text(size = 8))

df <- data.frame(patient = 1:N, ts = cumsum(ts_real$subopt_trt_safe),
                 rand = cumsum(rand_real$subopt_trt_safe),
                 rits = cumsum(rits_real$subopt_trt_safe))
df_long <- tidyr::gather(df, key = "strategy", value = "value", -patient)

# Plot using ggplot2
plot2 <- ggplot(df_long, aes(x = patient, y = value, color = strategy)) +
  geom_line() +
  geom_vline(xintercept = ts_real$tr_first, linetype = "dashed", color = "blue") +
  labs(title = "Safety Sub-optimal Treatment Count",
       x = "Patient", y = "Cumulative Sub-optimal Treatment Count") +
  theme(text = element_text(size = 8))

plot12 <- plot1 + plot2 + plot_layout(ncol = 2, guides = "collect")
ggsave("plot/regret_safe_real.jpg", height = 3, width = 6, units = "in")


### plotting real-time arm allocation 

my_colors <- c("#E69F00", "#56B4E9", "#009E73", 
               "#0072B2", "#D55E00", "#CC79A7")

# rand_real
df <- data.frame(patient_id = 1:length(rand_real$regret), 
                 dose = real_dat[["arm"]]
)
alloc_real_rand <- ggplot(df, aes(x = factor(patient_id), y = dose, 
                                  color = dose)) +
  geom_point(size = 0.5) +
  labs(x = "Patient", y = "Dose") +
  scale_color_manual(values = my_colors) +  # Specify manual color palette
  theme_minimal() + guides(color = FALSE) +
  theme(axis.text.x = element_blank(),  # Hide x-axis text
        axis.ticks.x = element_blank(),
        text = element_text(size = 12)) +
  labs(title = "Actual Randomized Allocation")


# ts_real
df <- data.frame(patient_id = 1:length(rand_real$regret), 
                 dose = factor(paste("Arm", ts_real$trt))
)
alloc_real_ts <- ggplot(df, aes(x = factor(patient_id), y = dose, 
                                color = dose)) +
  geom_point(size = 0.5) +
  labs(x = "Patient", y = "Dose") +
  scale_color_manual(values = my_colors) +  # Specify manual color palette
  theme_minimal() + guides(color = FALSE) +
  theme(axis.text.x = element_blank(),  # Hide x-axis text
        axis.ticks.x = element_blank(),
        text = element_text(size = 12))   +
  labs(title = "TS Allocation")

# rits_real
df <- data.frame(patient_id = 1:length(rand_real$regret), 
                 dose = factor(paste("Arm", rits_real$trt))
)
alloc_real_rits <- ggplot(df, aes(x = factor(patient_id), y = dose, 
                                  color = dose)) +
  geom_point(size = 0.5) +
  labs(x = "Patient", y = "Dose") +
  scale_color_manual(values = my_colors) +  # Specify manual color palette
  theme_minimal() + guides(color = FALSE) +
  theme(axis.text.x = element_blank(),  # Hide x-axis text
        axis.ticks.x = element_blank(),
        text = element_text(size = 12))   +
  labs(title = "RiTS Allocation")

plot_alloc <- alloc_real_rand + alloc_real_ts + alloc_real_rits + 
  plot_layout(ncol = 3)
ggsave("plot/arm_alloc_real.jpg", height = 4, width = 12, units = "in")


### plotting propensity
df <- as.data.frame(ts_real$log_dat$prpns_mat)
K <- length(unique(ts_real$trt))
colnames(df) <- paste("Arm", 1:K)
df$obs <- 1:N  # Adding a column for observation number

# Melt the data frame to long format for ggplot
df_long <- melt(df, id.vars = "obs", variable.name = "Arm", value.name = "Prpn")

# Plot using ggplot
plot_prop_ts <- ggplot(df_long, aes(x = obs, y = Prpn, color = Arm)) +
  geom_line(linewidth = 0.5) +
  labs(x = "Patient", y = "Propensity", 
       title = "Thompson Sampling") + ylim(0, 1) +
  theme(text = element_text(size = 12)) +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73", 
                                "#0072B2", "#D55E00", "#CC79A7"))

df <- as.data.frame(rits_real$log_dat$prpns_mat)
colnames(df) <- paste("Arm", 1:K)
df$obs <- 1:N  # Adding a column for observation number

# Melt the data frame to long format for ggplot
df_long <- melt(df, id.vars = "obs", variable.name = "Arm", value.name = "Prpn")

# Plot using ggplot
plot_prop_rits <- ggplot(df_long, aes(x = obs, y = Prpn, color = Arm)) +
  geom_line(linewidth = 0.5) +
  labs(x = "Patient", y = "Propensity", 
       title = "Risk-inclusive Thompson Sampling") + ylim(0, 1) +
  theme(text = element_text(size = 12)) +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73", 
                                "#0072B2", "#D55E00", "#CC79A7"))

trt_tab_real <- rbind("Rand" = table(real_dat[["arm"]]), 
                      "TS" = table(ts_real$trt), 
                      "RiTS" = table(rits_real$trt))
tab_grob <- gridExtra::tableGrob(trt_tab_real)
print(xtable::xtable(trt_tab_real))

plot12 <- plot_prop_ts + plot_prop_rits + tab_grob +
  plot_layout(ncol = 3, guides = "collect")
ggsave("plot/prop_all_real.jpg", height = 4, width = 20, units = "in")





# calculating sample averages
real_data <- readRDS("data/clean_data.RData")
mu_true <- sapply(1:K, function(k){
  mean(real_data[["pct_change_salt"]][real_data$arm == paste("Arm", k)])
}
  )

ate_start <- 30
ate <- ts_real$ate
# ATE
## TS
df_ts <- data.frame(
  Patient = ate_start:N,
  Arm = rep(1:K, each = N - (ate_start-1)),
  ATE = c(ate[, , 1]),
  ATEL = c(ate[, , 2]),
  ATEH = c(ate[, , 3]),
  TrueATE = rep(mu_true, each = N - (ate_start-1))
)

## RiTS
ate <- rits_real$ate
df_rits <- data.frame(
  Patient = ate_start:N,
  Arm = rep(1:K, each = N - (ate_start-1)),
  ATE = c(ate[, , 1]),
  ATEL = c(ate[, , 2]),
  ATEH = c(ate[, , 3]),
  TrueATE = rep(mu_true, each = N - (ate_start-1))
)

real_data <- readRDS("data/clean_data.RData")
placebo_arm <- 1
ate_rand <- get_ate_rand(real_data[["pct_change_salt"]], real_data[["arm"]], 
                         alpha = 0.05, placebo_arm = placebo_arm, 
                         ate_start = ate_start)
contr_rand <- ate_rand$contr
ate_rand <- ate_rand$ate
## rand
df_rand <- data.frame(
  Patient = ate_start:N,
  Arm = rep(1:K, each = N - (ate_start-1)),
  ATE = c(ate_rand[ate_start:N, , 1]),
  ATEL = c(ate_rand[ate_start:N, , 2]),
  ATEH = c(ate_rand[ate_start:N, , 3]),
  TrueATE = rep(mu_true, each = N - (ate_start-1))
)

min_ylim <- min(df_rand$ATEL, df_ts$ATEL, df_rits$ATEL, na.rm = TRUE)
max_ylim <- max(df_rand$ATEH, df_ts$ATEH, df_rits$ATEH, na.rm = TRUE)

plot3 <- ggplot(df_ts, aes(x = Patient, y = TrueATE, color = factor(Arm))) +
  geom_line() +
  geom_line(aes(y = ATEL), linetype = "dashed") +
  geom_line(aes(y = ATEH), linetype = "dashed") +
  geom_line(aes(y = ATE), linetype = "longdash") +
  # geom_hline(yintercept = mu[1], color = "red") +
  facet_wrap(~Arm) +
  labs(x = "Patient", y = "Average Arm Efficacy", 
       title = "AsympCS + TS") +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73", 
                                "#0072B2", "#D55E00", "#CC79A7"),
                     name = "Arm") +
  theme(text = element_text(size = 8)) +
  scale_y_continuous(limits = c(min_ylim, max_ylim))

plot4 <- ggplot(df_rits, aes(x = Patient, y = TrueATE, color = factor(Arm))) +
  geom_line() +
  geom_line(aes(y = ATEL), linetype = "dashed") +
  geom_line(aes(y = ATEH), linetype = "dashed") +
  geom_line(aes(y = ATE), linetype = "longdash") +
  # geom_hline(yintercept = mu[1], color = "red") +
  facet_wrap(~Arm) +
  labs(x = "Patient", y = "Average Arm Efficacy", 
       title = "AsympCS + RiTS") +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73", 
                                "#0072B2", "#D55E00", "#CC79A7"),
                     name = "Arm") +
  theme(text = element_text(size = 8)) +
  scale_y_continuous(limits = c(min_ylim, max_ylim))

# plot5 <- ggplot(df_rand, aes(x = Patient, y = TrueATE, color = factor(Arm))) +
#   geom_line() +
#   geom_line(aes(y = ATEL), linetype = "dashed") +
#   geom_line(aes(y = ATEH), linetype = "dashed") +
#   geom_line(aes(y = ATE), linetype = "longdash") +
#   # geom_hline(yintercept = mu[1], color = "red") +
#   facet_wrap(~Arm) +
#   labs(x = "Patient", y = "Average Arm Efficacy", 
#        title = "Simple Average + Rand") +
#   scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73", 
#                                 "#0072B2", "#D55E00", "#CC79A7"),
#                      name = "Arm") +
#   theme(text = element_text(size = 8)) +
#   scale_y_continuous(limits = c(min_ylim, max_ylim))

plot345 <- plot3 + plot4 + plot_layout(ncol = 2, guides = "collect")
ggsave("plot/ate_ts_rand_real.jpg", height = 3, width = 9, units = "in")


# contrast
## TS
true_contr <- mu_true - mu_true[placebo_arm]
true_contr <- true_contr[-placebo_arm]
contr_ts <- ts_real$contr
df_ts <- data.frame(
  Patient = ate_start:N,
  Arm = rep(2:K, each = N - (ate_start-1)),
  Contr = c(contr_ts[, , 1]),
  ContrL = c(contr_ts[, , 2]),
  ContrH = c(contr_ts[, , 3]),
  TrueContr = rep(true_contr, each = N - (ate_start-1))
)
contr_rits <- rits_real$contr
df_rits <- data.frame(
  Patient = ate_start:N,
  Arm = rep(2:K, each = N - (ate_start-1)),
  Contr = c(contr_rits[, , 1]),
  ContrL = c(contr_rits[, , 2]),
  ContrH = c(contr_rits[, , 3]),
  TrueContr = rep(true_contr, each = N - (ate_start-1))
)

df_rand <- data.frame(
  Patient = ate_start:N,
  Arm = rep(2:K, each = N - (ate_start-1)),
  Contr = c(contr_rand[ate_start:N, , 1]),
  ContrL = c(contr_rand[ate_start:N, , 2]),
  ContrH = c(contr_rand[ate_start:N, , 3]),
  TrueContr = rep(true_contr, each = N - (ate_start-1))
)

min_ylim <- min(df_rand$ContrL, df_ts$ContrL, df_rits$ContrL, na.rm = TRUE)
max_ylim <- max(df_rand$ContrH, df_ts$ContrH, df_rits$ContrH, na.rm = TRUE)

plot6 <- ggplot(df_ts, aes(x = Patient, y = TrueContr, color = factor(Arm))) +
  geom_line() +
  geom_line(aes(y = ContrL), linetype = "dashed") +
  geom_line(aes(y = ContrH), linetype = "dashed") +
  geom_line(aes(y = Contr), linetype = "longdash") +
  facet_wrap(~Arm) +
  labs(x = "Patient", y = "Effect Size", 
       title = "AsympCS + TS Contrast") +
  scale_color_manual(values = c("#56B4E9", "#009E73", 
                                "#0072B2", "#D55E00", "#CC79A7"),
                     name = "Arm") +
  theme(text = element_text(size = 8)) +
  scale_y_continuous(limits = c(min_ylim, max_ylim))

plot7 <- ggplot(df_rits, aes(x = Patient, y = TrueContr, color = factor(Arm))) +
  geom_line() +
  geom_line(aes(y = ContrL), linetype = "dashed") +
  geom_line(aes(y = ContrH), linetype = "dashed") +
  geom_line(aes(y = Contr), linetype = "longdash") +
  facet_wrap(~Arm) +
  labs(x = "Patient", y = "Effect Size", 
       title = "AsympCS + RiTS Contrast") +
  scale_color_manual(values = c("#56B4E9", "#009E73", 
                                "#0072B2", "#D55E00", "#CC79A7"),
                     name = "Arm") +
  theme(text = element_text(size = 8)) +
  scale_y_continuous(limits = c(min_ylim, max_ylim))

# plot8 <- ggplot(df_rand, aes(x = Patient, y = TrueContr, color = factor(Arm))) +
#   geom_line() +
#   geom_line(aes(y = ContrL), linetype = "dashed") +
#   geom_line(aes(y = ContrH), linetype = "dashed") +
#   geom_line(aes(y = Contr), linetype = "longdash") +
#   facet_wrap(~Arm) +
#   labs(x = "Patient", y = "Differene of Average Arm Efficacy over Placebo", 
#        title = "Simple Average + Rand Contrast") + 
#   scale_color_manual(values = c("#56B4E9", "#009E73", 
#                                 "#0072B2", "#D55E00", "#CC79A7"),
#                      name = "Arm") +
#   theme(text = element_text(size = 8)) +
#   scale_y_continuous(limits = c(min_ylim, max_ylim))

plot678 <- plot6 + plot7 + plot_layout(ncol = 2, guides = "collect")
ggsave("plot/contr_ts_rand_real.jpg", height = 3, width = 9, units = "in")
