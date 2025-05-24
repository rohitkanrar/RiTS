file_choice <- 3
ts_sim <- readRDS(paste("output/ts_mis_sim_", file_choice, ".RData", sep = ""))
rand_sim <- readRDS(paste("output/rand_sim_", 1, ".RData", sep = ""))
rits_sim <- readRDS(paste("output/rits_mis_sim_", file_choice, ".RData", sep = ""))
n_iter <- length(ts_sim)
N <- length(ts_sim[[1]]$trt)
ind <- c(seq(50, 100, 10), 150, 200)
K <- length(unique(rand_sim[[1]]$trt))

library(ggplot2)
library(tidyr)
library(patchwork)

### Regret Combined
## Utility Regret
rand_reg <- t(sapply(1:n_iter, function(i){
  cumsum(rand_sim[[i]]$regret)[ind]
}))
ts_reg <- t(sapply(1:n_iter, function(i){
  cumsum(ts_sim[[i]]$regret)[ind]
}))
rits_reg <- t(sapply(1:n_iter, function(i){
  cumsum(rits_sim[[i]]$regret)[ind]
}))


# Combine matrices into a single data frame
df <- data.frame(
  Value = c(rand_reg, ts_reg, rits_reg),
  Method = rep(c("rand", "ts", "rits"), 
               each = nrow(rand_reg)*ncol(rand_reg)),
  Column = rep(rep(1:ncol(rand_reg), each = nrow(rand_reg), 
                   times = 3))
)
library(ggplot2)

# Plot grouped boxplot
sim_reg_plot1 <- ggplot(df, aes(x = factor(Column), 
                                y = Value, fill = Method)) +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.7) +
  labs(x = "Number of Participants", y = "Cumulative Regret", 
       fill = "Method") +
  scale_x_discrete(labels = ind) + ggtitle("Utility") +
  scale_fill_manual(
    values = c("rand" = "#CC79A7", "ts" = "#0072B2", "rits" = "#D55E00"),
    labels = c("ts" = "TS", "rand" = "Rand", "rits" = "RiTS")
  )

## Efficacy Regret
rand_reg <- t(sapply(1:n_iter, function(i){
  cumsum(rand_sim[[i]]$regret_benf)[ind]
}))
ts_reg <- t(sapply(1:n_iter, function(i){
  cumsum(ts_sim[[i]]$regret_benf)[ind]
}))
rits_reg <- t(sapply(1:n_iter, function(i){
  cumsum(rits_sim[[i]]$regret_benf)[ind]
}))


# Combine matrices into a single data frame
df <- data.frame(
  Value = c(rand_reg, ts_reg, rits_reg),
  Method = rep(c("rand", "ts", "rits"), 
               each = nrow(rand_reg)*ncol(rand_reg)),
  Column = rep(rep(1:ncol(rand_reg), each = nrow(rand_reg), 
                   times = 3))
)
library(ggplot2)

# Plot grouped boxplot
sim_reg_plot2 <- ggplot(df, aes(x = factor(Column), 
                                y = Value, fill = Method)) +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.7) +
  labs(x = "Number of Participants", y = "Cumulative Regret", 
       fill = "Method") +
  scale_x_discrete(labels = ind) + ggtitle("Efficacy") +
  scale_fill_manual(
    values = c("rand" = "#CC79A7", "ts" = "#0072B2", "rits" = "#D55E00"),
    labels = c("ts" = "TS", "rand" = "Rand", "rits" = "RiTS")
  ) 


## Safety Regret
rand_reg <- t(sapply(1:n_iter, function(i){
  cumsum(rand_sim[[i]]$regret_safe)[ind]
}))
ts_reg <- t(sapply(1:n_iter, function(i){
  cumsum(ts_sim[[i]]$regret_safe)[ind]
}))
rits_reg <- t(sapply(1:n_iter, function(i){
  cumsum(rits_sim[[i]]$regret_safe)[ind]
}))


# Combine matrices into a single data frame
df <- data.frame(
  Value = c(rand_reg, ts_reg, rits_reg),
  Method = rep(c("rand", "ts", "rits"), 
               each = nrow(rand_reg)*ncol(rand_reg)),
  Column = rep(rep(1:ncol(rand_reg), each = nrow(rand_reg), 
                   times = 3))
)
library(ggplot2)

# Plot grouped boxplot
sim_reg_plot3 <- ggplot(df, aes(x = factor(Column), 
                                y = Value, fill = Method)) +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.7) +
  labs(x = "Number of Participants", y = "Cumulative Regret", 
       fill = "Method") +
  scale_x_discrete(labels = ind)  + ggtitle("Safety") +
  scale_fill_manual(
    values = c("rand" = "#CC79A7", "ts" = "#0072B2", "rits" = "#D55E00"),
    labels = c("ts" = "TS", "rand" = "Rand", "rits" = "RiTS")
  ) 

sim_regret_plot <- sim_reg_plot1 + sim_reg_plot2 + sim_reg_plot3 +
  plot_layout(ncol = 3, guides = "collect")
ggsave(paste("plot/regret_mis_sim_bwplot_", file_choice, ".jpg", sep = ""), 
       height = 4, width = 12, units = "in")


# ------------------------------------------- #

## Arm Allocation Boxplot 
trt_rand <- t(sapply(1:n_iter, function(i) table(rand_sim[[i]]$trt)))
trt_ts <- t(sapply(1:n_iter, function(i) table(ts_sim[[i]]$trt)))
trt_rits <- t(sapply(1:n_iter, function(i) table(rits_sim[[i]]$trt)))

df <- data.frame(
  Frequency = c(trt_rand, trt_ts, trt_rits),
  Arm = rep(rep(1:4, each = nrow(trt_rand)), times = 3),
  Method = rep(c("rand", "ts", "rits"), each = nrow(trt_rand)*ncol(trt_rand))
)

ggplot(df, aes(x = factor(Arm), 
               y = Frequency, fill = Method)) +
  geom_boxplot() +
  labs(x = "Arm", y = "Frequency of Allocation", 
       fill = "Method") + 
  scale_fill_manual(
    values = c("rand" = "#CC79A7", "ts" = "#0072B2", "rits" = "#D55E00"),
    labels = c("ts" = "TS", "rand" = "Rand", "rits" = "RiTS")
  ) +
  ggtitle("Frequency of Arm Allocations by Different Methods Across 1000 Trials")
ggsave(paste("plot/freq_arm_alloc_mis_sim_", file_choice, ".jpg", sep = ""), 
       height = 4, width = 8, units = "in")

# -------------------------------------------- #

sim_choice <- readRDS("metadata/sim_choice.RData")
ate_start <- sim_choice$ate_start
ate_ind <- ind - ate_start
## Width
# Arm 2
contr_rand_width2 <- t(sapply(1:n_iter, function(iter){
  rand_sim[[iter]]$contr[ate_ind, 1, 3] - rand_sim[[iter]]$contr[ate_ind, 1, 2]
}))

contr_ts_width2 <- t(sapply(1:n_iter, function(iter){
  ts_sim[[iter]]$contr[ate_ind, 1, 3] - ts_sim[[iter]]$contr[ate_ind, 1, 2]
}))

contr_rits_width2 <- t(sapply(1:n_iter, function(iter){
  rits_sim[[iter]]$contr[ate_ind, 1, 3] - ts_sim[[iter]]$contr[ate_ind, 1, 2]
}))

df <- data.frame(
  Width = c(contr_rand_width2, contr_ts_width2, contr_rits_width2),
  Method = rep(c("rand", "ts", "rits"), 
               each = nrow(contr_rand_width2)*ncol(contr_rand_width2)),
  Column = rep(rep(1:ncol(contr_rand_width2), each = nrow(contr_rand_width2), 
                   times = 3))
)

sim_wid_plot1 <- ggplot(df, aes(x = factor(Column), 
                                y = Width, fill = Method)) +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.7) +
  labs(x = "Number of Participants", y = "Width", 
       fill = "Arm") + ylim(0, 15) +
  scale_fill_manual(
    values = c("rand" = "#CC79A7", "ts" = "#0072B2", "rits" = "#D55E00"),
    labels = c("ts" = "TS", "rand" = "Rand", "rits" = "RiTS")
  ) +
  scale_x_discrete(labels = ind) + ggtitle("Width for Arm 2 - Arm 1") 


# Arm 3
contr_rand_width3 <- t(sapply(1:n_iter, function(iter){
  rand_sim[[iter]]$contr[ate_ind, 2, 3] - rand_sim[[iter]]$contr[ate_ind, 2, 2]
}))

contr_ts_width3 <- t(sapply(1:n_iter, function(iter){
  ts_sim[[iter]]$contr[ate_ind, 2, 3] - ts_sim[[iter]]$contr[ate_ind, 2, 2]
}))

contr_rits_width3 <- t(sapply(1:n_iter, function(iter){
  rits_sim[[iter]]$contr[ate_ind, 2, 3] - ts_sim[[iter]]$contr[ate_ind, 2, 2]
}))

df <- data.frame(
  Width = c(contr_rand_width3, contr_ts_width3, contr_rits_width3),
  Method = rep(c("rand", "ts", "rits"), 
               each = nrow(contr_rand_width3)*ncol(contr_rand_width3)),
  Column = rep(rep(1:ncol(contr_rand_width3), each = nrow(contr_rand_width3), 
                   times = 3))
)

sim_wid_plot2 <- ggplot(df, aes(x = factor(Column), 
                                y = Width, fill = Method)) +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.7) +
  labs(x = "Number of Participants", y = "Width", 
       fill = "Arm") + ylim(0, 15) +
  scale_fill_manual(
    values = c("rand" = "#CC79A7", "ts" = "#0072B2", "rits" = "#D55E00"),
    labels = c("ts" = "TS", "rand" = "Rand", "rits" = "RiTS")
  ) +
  scale_x_discrete(labels = ind) + ggtitle("Width for Arm 3 - Arm 1") 


# Arm 4
contr_rand_width4 <- t(sapply(1:n_iter, function(iter){
  rand_sim[[iter]]$contr[ate_ind, 3, 3] - rand_sim[[iter]]$contr[ate_ind, 3, 2]
}))

contr_ts_width4 <- t(sapply(1:n_iter, function(iter){
  ts_sim[[iter]]$contr[ate_ind, 3, 3] - ts_sim[[iter]]$contr[ate_ind, 3, 2]
}))

contr_rits_width4 <- t(sapply(1:n_iter, function(iter){
  rits_sim[[iter]]$contr[ate_ind, 3, 3] - ts_sim[[iter]]$contr[ate_ind, 3, 2]
}))

df <- data.frame(
  Width = c(contr_rand_width4, contr_ts_width4, contr_rits_width4),
  Method = rep(c("rand", "ts", "rits"), 
               each = nrow(contr_rand_width4)*ncol(contr_rand_width4)),
  Column = rep(rep(1:ncol(contr_rand_width4), each = nrow(contr_rand_width4), 
                   times = 3))
)

sim_wid_plot3 <- ggplot(df, aes(x = factor(Column), 
                                y = Width, fill = Method)) +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.7) +
  labs(x = "Number of Participants", y = "Width", 
       fill = "Arm") + ylim(0, 15) +
  scale_fill_manual(
    values = c("rand" = "#CC79A7", "ts" = "#0072B2", "rits" = "#D55E00"),
    labels = c("ts" = "TS", "rand" = "Rand", "rits" = "RiTS")
  ) +
  scale_x_discrete(labels = ind) + ggtitle("Width for Arm 4 - Arm 1") 

sim_wid_plot <- sim_wid_plot1 + sim_wid_plot2 + sim_wid_plot3 +
  plot_layout(ncol = 3, guides = "collect")
ggsave(paste("plot/width_mis_sim_bwplot_", file_choice, ".jpg", sep = ""), 
       height = 4, width = 12, units = "in")

# ------------------------------------------------------ #

## bias
sim_dat <- readRDS("metadata/sim_dat.RData")
mu_true <- sim_dat$mu_true
contr_true <- mu_true - mu_true[1]
contr_true <- contr_true[setdiff(1:K, sim_dat$placebo_arm)]

# Arm 2
contr_rand_bias2 <- t(sapply(1:n_iter, function(iter){
  rand_sim[[iter]]$contr[ate_ind, 1, 1] - contr_true[1]
}))

contr_ts_bias2 <- t(sapply(1:n_iter, function(iter){
  ts_sim[[iter]]$contr[ate_ind, 1, 1] - contr_true[1]
}))

contr_rits_bias2 <- t(sapply(1:n_iter, function(iter){
  rits_sim[[iter]]$contr[ate_ind, 1, 1] - contr_true[1]
}))

df <- data.frame(
  Bias = c(contr_rand_bias2, contr_ts_bias2, contr_rits_bias2),
  Method = rep(c("rand", "ts", "rits"), 
               each = nrow(contr_rand_bias2)*ncol(contr_rand_bias2)),
  Column = rep(rep(1:ncol(contr_rand_bias2), each = nrow(contr_rand_bias2), 
                   times = 3))
)

sim_bias_plot1 <- ggplot(df, aes(x = factor(Column), 
                                 y = Bias, fill = Method)) +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.7) +
  labs(x = "Number of Participants", y = "Bias", 
       fill = "Arm") + ylim(-3, 3) +
  scale_fill_manual(
    values = c("rand" = "#CC79A7", "ts" = "#0072B2", "rits" = "#D55E00"),
    labels = c("ts" = "TS", "rand" = "Rand", "rits" = "RiTS")
  ) +
  scale_x_discrete(labels = ind) + ggtitle("Bias for Arm 2 - Arm 1") 


## bias
# Arm 3
contr_rand_bias3 <- t(sapply(1:n_iter, function(iter){
  rand_sim[[iter]]$contr[ate_ind, 2, 1] - contr_true[2]
}))

contr_ts_bias3 <- t(sapply(1:n_iter, function(iter){
  ts_sim[[iter]]$contr[ate_ind, 2, 1] - contr_true[2]
}))

contr_rits_bias3 <- t(sapply(1:n_iter, function(iter){
  rits_sim[[iter]]$contr[ate_ind, 2, 1] - contr_true[2]
}))

df <- data.frame(
  Bias = c(contr_rand_bias3, contr_ts_bias3, contr_rits_bias3),
  Method = rep(c("rand", "ts", "rits"), 
               each = nrow(contr_rand_bias3)*ncol(contr_rand_bias3)),
  Column = rep(rep(1:ncol(contr_rand_bias3), each = nrow(contr_rand_bias3), 
                   times = 3))
)

sim_bias_plot2 <- ggplot(df, aes(x = factor(Column), 
                                 y = Bias, fill = Method)) +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.7) +
  labs(x = "Number of Participants", y = "Bias", 
       fill = "Arm") + ylim(-3, 3) +
  scale_fill_manual(
    values = c("rand" = "#CC79A7", "ts" = "#0072B2", "rits" = "#D55E00"),
    labels = c("ts" = "TS", "rand" = "Rand", "rits" = "RiTS")
  ) +
  scale_x_discrete(labels = ind) + ggtitle("Bias for Arm 3 - Arm 1") 


## bias
# Arm 4
contr_rand_bias4 <- t(sapply(1:n_iter, function(iter){
  rand_sim[[iter]]$contr[ate_ind, 3, 1] - contr_true[3]
}))

contr_ts_bias4 <- t(sapply(1:n_iter, function(iter){
  ts_sim[[iter]]$contr[ate_ind, 3, 1] - contr_true[3]
}))

contr_rits_bias4 <- t(sapply(1:n_iter, function(iter){
  rits_sim[[iter]]$contr[ate_ind, 3, 1] - contr_true[3]
}))

df <- data.frame(
  Bias = c(contr_rand_bias4, contr_ts_bias4, contr_rits_bias4),
  Method = rep(c("rand", "ts", "rits"), 
               each = nrow(contr_rand_bias4)*ncol(contr_rand_bias4)),
  Column = rep(rep(1:ncol(contr_rand_bias4), each = nrow(contr_rand_bias4), 
                   times = 3))
)

sim_bias_plot3 <- ggplot(df, aes(x = factor(Column), 
                                 y = Bias, fill = Method)) +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.7) +
  labs(x = "Number of Participants", y = "Bias", 
       fill = "Arm") + ylim(-3, 3) +
  scale_fill_manual(
    values = c("rand" = "#CC79A7", "ts" = "#0072B2", "rits" = "#D55E00"),
    labels = c("ts" = "TS", "rand" = "Rand", "rits" = "RiTS")
  ) +
  scale_x_discrete(labels = ind) + ggtitle("Bias for Arm 4 - Arm 1") 

sim_bias_plot <- sim_bias_plot1 + sim_bias_plot2 + sim_bias_plot3 +
  plot_layout(ncol = 3, guides = "collect")
ggsave(paste("plot/bias_mis_sim_bwplot_", file_choice, ".jpg", sep = ""), 
       height = 4, width = 12, units = "in")

# ---------------------------------------------- #

# Cumulative Miscoverage Rates
source("code/function/misc.R")
library(reshape2)
K <- length(unique(rand_sim[[1]]$trt))
delays <- c(20, 70, 120)

for(delay in delays){
  # TS
  cum_miscov <- get_cum_mis_cov(ts_sim, mu_true = mu_true, 
                                contr_true = contr_true, delay = delay)
  
  df <- as.data.frame(cum_miscov[[2]] / length(rand_sim))
  colnames(df) <- paste("Arm", 2:K)
  df$obs <- (ate_start+delay):N  # Adding a column for observation number
  
  # Melt the data frame to long format for ggplot
  df_long <- melt(df, id.vars = "obs", variable.name = "Arm", value.name = "Miscov")
  
  # Plot using ggplot
  plot_cum_miscov_ts <- ggplot(df_long, aes(x = obs, y = Miscov, color = Arm)) +
    geom_line(linewidth = 0.5) + ylim(c(0, 0.2)) + 
    geom_hline(yintercept = sim_choice$alpha/(K-1), linetype = "dashed", color = "blue") +
    labs(x = "Patient", y = "Cumulative Miscoverage", 
         title = "TS Allocation") +
    theme(text = element_text(size = 10)) +
    scale_color_manual(
      values = c("Arm 1" = "#E69F00", "Arm 2" = "#56B4E9", 
                 "Arm 3" = "#009E73", "Arm 4" = "#CC79A7")
    )
  
  # RiTS
  cum_miscov <- get_cum_mis_cov(rits_sim, mu_true = mu_true, 
                                contr_true = contr_true, delay = delay)
  
  df <- as.data.frame(cum_miscov[[2]] / length(rand_sim))
  colnames(df) <- paste("Arm", 2:K)
  df$obs <- (ate_start+delay):N  # Adding a column for observation number
  
  
  # Melt the data frame to long format for ggplot
  df_long <- melt(df, id.vars = "obs", variable.name = "Arm", value.name = "Miscov")
  
  # Plot using ggplot
  plot_cum_miscov_rits <- ggplot(df_long, aes(x = obs, y = Miscov, color = Arm)) +
    geom_line(linewidth = 0.5) + ylim(c(0, 0.2)) + 
    geom_hline(yintercept = sim_choice$alpha/(K-1), linetype = "dashed", color = "blue") +
    labs(x = "Patient", y = "Cumulative Miscoverage", 
         title = "RiTS Allocation") +
    theme(text = element_text(size = 10)) +
    scale_color_manual(
      values = c("Arm 1" = "#E69F00", "Arm 2" = "#56B4E9", 
                 "Arm 3" = "#009E73", "Arm 4" = "#CC79A7")
    )
  
  plot_cum_miscov <- plot_cum_miscov_ts + plot_cum_miscov_rits +
    plot_layout(ncol = 2, guides = "collect")
  ggsave(paste("plot/cum_miscov_mis_contr_", file_choice, "_peek_", 
               delay+ate_start, ".jpg", sep = ""), 
         height = 4, width = 8, units = "in")
}
