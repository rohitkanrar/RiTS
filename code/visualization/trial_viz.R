library(ggplot2)
library(patchwork)
library(reshape2)

ts_out <- readRDS("output/ts_out.RData")
rits_out <- readRDS("output/rits_out.RData")
rand_out <- readRDS("output/rand_out.RData")
sim_dat <- readRDS("output/sim_dat.RData")

N <- length(ts_out$trt)
beta_true <- sim_dat$beta_true
mu_true <- sim_dat$mu_true
ate_start <- sim_dat$ate_start
tr_start <- sim_dat$tr_start
placebo_arm <- sim_dat$placebo_arm

# Trt Table
trt_tab <- rbind("Rand" = table(rand_out$trt), "TS" = table(ts_out$trt), 
                 "RiTS" = table(rits_out$trt))
colnames(trt_tab) <- paste("Arm", 1:4)
xtable::xtable(trt_tab)

# regret utility
df <- data.frame(patient = 1:N, ts = cumsum(ts_out$regret),
                 rand = cumsum(rand_out$regret), 
                 rits = cumsum(rits_out$regret))
df_long <- tidyr::gather(df, key = "strategy", value = "value", -patient)

# Plot using ggplot2
plot1 <- ggplot(df_long, aes(x = patient, y = value, color = strategy)) +
  geom_line() +
  geom_vline(xintercept = ts_out$tr_first, linetype = "dashed", color = "blue") +
  labs(title = "Utility Regret",
       x = "Patient", y = "Cumulative Regret") +
  theme(text = element_text(size = 8)) +
  scale_color_manual(
    values = c("rand" = "#CC79A7", "ts" = "#0072B2", "rits" = "#D55E00"),
    labels = c("ts" = "TS", "rand" = "Rand", "rits" = "RiTS")
  )

df <- data.frame(patient = 1:N, ts = cumsum(ts_out$subopt_trt),
                 rand = cumsum(rand_out$subopt_trt),
                 rits = cumsum(rits_out$subopt_trt))
df_long <- tidyr::gather(df, key = "strategy", value = "value", -patient)

# Plot using ggplot2
plot2 <- ggplot(df_long, aes(x = patient, y = value, color = strategy)) +
  geom_line() +
  geom_vline(xintercept = ts_out$tr_first, linetype = "dashed", color = "blue") +
  labs(title = "Utility Sub-optimal Treatment Count",
       x = "Patient", y = "Cumulative Sub-optimal Treatment Count") +
  theme(text = element_text(size = 8)) +
  scale_color_manual(
    values = c("rand" = "#CC79A7", "ts" = "#0072B2", "rits" = "#D55E00"),
    labels = c("ts" = "TS", "rand" = "Rand", "rits" = "RiTS")
  )

plot12 <- plot1 + plot2 + plot_layout(ncol = 2, guides = "collect")
ggsave("plot/regret_all.jpg", height = 3, width = 6, units = "in")


# regret benefit
df <- data.frame(patient = 1:N, ts = cumsum(ts_out$regret_benf),
                 rand = cumsum(rand_out$regret_benf), 
                 rits = cumsum(rits_out$regret_benf))
df_long <- tidyr::gather(df, key = "strategy", value = "value", -patient)

# Plot using ggplot2
plot1 <- ggplot(df_long, aes(x = patient, y = value, color = strategy)) +
  geom_line() +
  geom_vline(xintercept = ts_out$tr_first, linetype = "dashed", color = "blue") +
  labs(title = "Efficacy Regret",
       x = "Patient", y = "Cumulative Regret") +
  theme(text = element_text(size = 8)) +
  scale_color_manual(
    values = c("rand" = "#CC79A7", "ts" = "#0072B2", "rits" = "#D55E00"),
    labels = c("ts" = "TS", "rand" = "Rand", "rits" = "RiTS")
  )

df <- data.frame(patient = 1:N, ts = cumsum(ts_out$subopt_trt_benf),
                 rand = cumsum(rand_out$subopt_trt_benf),
                 rits = cumsum(rits_out$subopt_trt_benf))
df_long <- tidyr::gather(df, key = "strategy", value = "value", -patient)

# Plot using ggplot2
plot2 <- ggplot(df_long, aes(x = patient, y = value, color = strategy)) +
  geom_line() +
  geom_vline(xintercept = ts_out$tr_first, linetype = "dashed", color = "blue") +
  labs(title = "Efficacy Sub-optimal Treatment Count",
       x = "Patient", y = "Cumulative Sub-optimal Treatment Count") +
  theme(text = element_text(size = 8)) +
  scale_color_manual(
    values = c("rand" = "#CC79A7", "ts" = "#0072B2", "rits" = "#D55E00"),
    labels = c("ts" = "TS", "rand" = "Rand", "rits" = "RiTS")
  )

plot12 <- plot1 + plot2 + plot_layout(ncol = 2, guides = "collect")
ggsave("plot/regret_benf.jpg", height = 3, width = 6, units = "in")

# regret safety
df <- data.frame(patient = 1:N, ts = cumsum(ts_out$regret_safe),
                 rand = cumsum(rand_out$regret_safe), 
                 rits = cumsum(rits_out$regret_safe))
df_long <- tidyr::gather(df, key = "strategy", value = "value", -patient)

# Plot using ggplot2
plot1 <- ggplot(df_long, aes(x = patient, y = value, color = strategy)) +
  geom_line() +
  geom_vline(xintercept = ts_out$tr_first, linetype = "dashed", color = "blue") +
  labs(title = "Safety Regret",
       x = "Patient", y = "Cumulative Regret") +
  theme(text = element_text(size = 8)) +
  scale_color_manual(
    values = c("rand" = "#CC79A7", "ts" = "#0072B2", "rits" = "#D55E00"),
    labels = c("ts" = "TS", "rand" = "Rand", "rits" = "RiTS")
  )

df <- data.frame(patient = 1:N, ts = cumsum(ts_out$subopt_trt_safe),
                 rand = cumsum(rand_out$subopt_trt_safe),
                 rits = cumsum(rits_out$subopt_trt_safe))
df_long <- tidyr::gather(df, key = "strategy", value = "value", -patient)

# Plot using ggplot2
plot2 <- ggplot(df_long, aes(x = patient, y = value, color = strategy)) +
  geom_line() +
  geom_vline(xintercept = ts_out$tr_first, linetype = "dashed", color = "blue") +
  labs(title = "Safety Sub-optimal Treatment Count",
       x = "Patient", y = "Cumulative Sub-optimal Treatment Count") +
  theme(text = element_text(size = 8)) +
  scale_color_manual(
    values = c("rand" = "#CC79A7", "ts" = "#0072B2", "rits" = "#D55E00"),
    labels = c("ts" = "TS", "rand" = "Rand", "rits" = "RiTS")
  )

plot12 <- plot1 + plot2 + plot_layout(ncol = 2, guides = "collect")
ggsave("plot/regret_safe.jpg", height = 3, width = 6, units = "in")



# plotting propensity
K <- length(unique(rand_out$trt))
df <- as.data.frame(ts_out$log_dat$prpns_mat)
colnames(df) <- paste("Arm", 1:K)
df$obs <- 1:N  # Adding a column for observation number

# Melt the data frame to long format for ggplot
df_long <- melt(df, id.vars = "obs", variable.name = "Arm", value.name = "Prpn")

# Plot using ggplot
plot_prop_ts <- ggplot(df_long, aes(x = obs, y = Prpn, color = Arm)) +
  geom_line(linewidth = 0.5) +
  labs(x = "Patient", y = "Propensity", 
       title = "Thompson Sampling") +
  theme(text = element_text(size = 10)) +
  scale_color_manual(
    values = c("Arm 1" = "#E69F00", "Arm 2" = "#56B4E9", 
               "Arm 3" = "#009E73", "Arm 4" = "#CC79A7")
  )

df <- as.data.frame(rits_out$log_dat$prpns_mat)
colnames(df) <- paste("Arm", 1:K)
df$obs <- 1:N  # Adding a column for observation number

# Melt the data frame to long format for ggplot
df_long <- melt(df, id.vars = "obs", variable.name = "Arm", value.name = "Prpn")

# Plot using ggplot
plot_prop_bivts <- ggplot(df_long, aes(x = obs, y = Prpn, color = Arm)) +
  geom_line(linewidth = 0.5) +
  labs(x = "Patient", y = "Propensity", 
       title = "Risk-inclusive Thompson Sampling") +
  theme(text = element_text(size = 10)) +
  scale_color_manual(
    values = c("Arm 1" = "#E69F00", "Arm 2" = "#56B4E9", 
               "Arm 3" = "#009E73", "Arm 4" = "#CC79A7")
  )

plot12 <-plot_prop_ts + plot_prop_bivts + plot_layout(ncol = 2, guides = "collect")
ggsave("plot/prop_all.jpg", height = 4, width = 9, units = "in")


ate_start <- ts_out$tr_first
ate <- ts_out$ate
# ATE
## TS
df_ts <- data.frame(
  Patient = ate_start:N,
  Arm = rep(1:K, each = N - (ate_start-1)),
  ATE = c(ate[ate_start:N, , 1]),
  ATEL = c(ate[ate_start:N, , 2]),
  ATEH = c(ate[ate_start:N, , 3]),
  TrueATE = rep(mu_true, each = N - (ate_start-1))
)

## RiTS
ate <- rits_out$ate
df_rits <- data.frame(
  Patient = ate_start:N,
  Arm = rep(1:K, each = N - (ate_start-1)),
  ATE = c(ate[ate_start:N, , 1]),
  ATEL = c(ate[ate_start:N, , 2]),
  ATEH = c(ate[ate_start:N, , 3]),
  TrueATE = rep(mu_true, each = N - (ate_start-1))
)

ate_rand <- rand_out$ate
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
       title = "AW-AIPW + TS") +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#CC79A7"), 
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
       title = "AW-AIPW + RiTS") +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#CC79A7"), 
                     name = "Arm") +
  theme(text = element_text(size = 8)) +
  scale_y_continuous(limits = c(min_ylim, max_ylim))

plot5 <- ggplot(df_rand, aes(x = Patient, y = TrueATE, color = factor(Arm))) +
  geom_line() +
  geom_line(aes(y = ATEL), linetype = "dashed") +
  geom_line(aes(y = ATEH), linetype = "dashed") +
  geom_line(aes(y = ATE), linetype = "longdash") +
  # geom_hline(yintercept = mu[1], color = "red") +
  facet_wrap(~Arm) +
  labs(x = "Patient", y = "Average Arm Efficacy", 
       title = "Simple Average + Rand") +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#CC79A7"), 
                     name = "Arm") +
  theme(text = element_text(size = 8)) +
  scale_y_continuous(limits = c(min_ylim, max_ylim))

plot345 <- plot3 + plot4 + plot5 + plot_layout(ncol = 3, guides = "collect")
ggsave("plot/ate_ts_rand.jpg", height = 3, width = 9, units = "in")


# contrast
## TS
true_contr <- mu_true - mu_true[placebo_arm]
true_contr <- true_contr[-placebo_arm]
contr_ts <- ts_out$contr
df_ts <- data.frame(
  Patient = ate_start:N,
  Arm = rep(2:K, each = N - (ate_start-1)),
  Contr = c(contr_ts[ate_start:N, , 1]),
  ContrL = c(contr_ts[ate_start:N, , 2]),
  ContrH = c(contr_ts[ate_start:N, , 3]),
  TrueContr = rep(true_contr, each = N - (ate_start-1))
)
contr_rits <- rits_out$contr
df_rits <- data.frame(
  Patient = ate_start:N,
  Arm = rep(2:K, each = N - (ate_start-1)),
  Contr = c(contr_rits[ate_start:N, , 1]),
  ContrL = c(contr_rits[ate_start:N, , 2]),
  ContrH = c(contr_rits[ate_start:N, , 3]),
  TrueContr = rep(true_contr, each = N - (ate_start-1))
)
contr_rand <- rand_out$contr
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
  labs(x = "Patient", y = "Differene of Average Arm Efficacy over Placebo", 
       title = "AW-AIPW + TS Contrast") +
  scale_color_manual(values = c("#56B4E9", "#009E73", "#CC79A7"), 
                     name = "Arm") +
  theme(text = element_text(size = 8)) +
  scale_y_continuous(limits = c(min_ylim, max_ylim))

plot7 <- ggplot(df_rits, aes(x = Patient, y = TrueContr, color = factor(Arm))) +
  geom_line() +
  geom_line(aes(y = ContrL), linetype = "dashed") +
  geom_line(aes(y = ContrH), linetype = "dashed") +
  geom_line(aes(y = Contr), linetype = "longdash") +
  facet_wrap(~Arm) +
  labs(x = "Patient", y = "Differene of Average Arm Efficacy over Placebo", 
       title = "AW-AIPW + RiTS Contrast") +
  scale_color_manual(values = c("#56B4E9", "#009E73", "#CC79A7"), 
                     name = "Arm") +
  theme(text = element_text(size = 8)) +
  scale_y_continuous(limits = c(min_ylim, max_ylim))

plot8 <- ggplot(df_rand, aes(x = Patient, y = TrueContr, color = factor(Arm))) +
  geom_line() +
  geom_line(aes(y = ContrL), linetype = "dashed") +
  geom_line(aes(y = ContrH), linetype = "dashed") +
  geom_line(aes(y = Contr), linetype = "longdash") +
  facet_wrap(~Arm) +
  labs(x = "Patient", y = "Differene of Average Arm Efficacy over Placebo", 
       title = "Simple Average + Rand Contrast") +
  scale_color_manual(values = c("#56B4E9", "#009E73", "#CC79A7"), 
                     name = "Arm") +
  theme(text = element_text(size = 8)) +
  scale_y_continuous(limits = c(min_ylim, max_ylim))

plot678 <- plot6 + plot7 + plot8 + plot_layout(ncol = 3, guides = "collect")
ggsave("plot/contr_ts_rand.jpg", height = 3, width = 9, units = "in")