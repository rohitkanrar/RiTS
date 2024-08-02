source("code/visualization/replicated_trial_tables.R")

colnames(tab_metric) <- paste(ind)
rownames(tab_metric) <- rep(c("Rand", "TS", "RiTS"), 2)
print(xtable::xtable(tab_metric))

covr_contr <- cbind(covr_contr_rand, covr_contr_ts, covr_contr_rits)
rownames(covr_contr) <- paste("Arm", 2:K, "- Arm 1")
colnames(covr_contr) <- rep(paste(ind), 3)
print(xtable::xtable(covr_contr))


width_contr <- cbind(width_contr_rand, width_contr_ts, width_contr_rits)
rownames(width_contr) <- paste("Arm", 2:K, "- Arm 1")
colnames(width_contr) <- rep(paste(ind), 3)
print(xtable::xtable(width_contr))


bias_contr <- cbind(bias_contr_rand, bias_contr_ts, bias_contr_rits)
rownames(bias_contr) <- paste("Arm", 2:K, "- Arm 1")
colnames(bias_contr) <- rep(paste(ind), 3)
print(xtable::xtable(bias_contr))



# stop_n <- data.frame(Rand = stop_n_rand, TS = stop_n_ts,
#                      RiTS = stop_n_rits)

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
ggsave("plot/regret_sim_bwplot.jpg", height = 4, width = 12, units = "in")

### Sub-opt plot Combined
## Utility Sub-optimal Treatment Count
rand_sub <- t(sapply(1:n_iter, function(i){
  cumsum(rand_sim[[i]]$subopt_trt)[ind]
}))
ts_sub <- t(sapply(1:n_iter, function(i){
  cumsum(ts_sim[[i]]$subopt_trt)[ind]
}))
rits_sub <- t(sapply(1:n_iter, function(i){
  cumsum(rits_sim[[i]]$subopt_trt)[ind]
}))


# Combine matrices into a single data frame
df <- data.frame(
  Value = c(rand_sub, ts_sub, rits_sub),
  Method = rep(c("rand", "ts", "rits"), 
               each = nrow(rand_sub)*ncol(rand_sub)),
  Column = rep(rep(1:ncol(rand_sub), each = nrow(rand_sub), 
                   times = 3))
)
library(ggplot2)

# Plot grouped boxplot
sim_sub_plot1 <- ggplot(df, aes(x = factor(Column), 
                                y = Value, fill = Method)) +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.7) +
  labs(x = "Number of Participants", y = "Cumulative Sub-optimal Treatment Count", 
       fill = "Method") +
  scale_x_discrete(labels = ind) + ggtitle("Utility") +
  scale_fill_manual(
    values = c("rand" = "#CC79A7", "ts" = "#0072B2", "rits" = "#D55E00"),
    labels = c("ts" = "TS", "rand" = "Rand", "rits" = "RiTS")
  ) 

## Efficacy Sub-optimal Treatment Count
rand_sub <- t(sapply(1:n_iter, function(i){
  cumsum(rand_sim[[i]]$subopt_trt_benf)[ind]
}))
ts_sub <- t(sapply(1:n_iter, function(i){
  cumsum(ts_sim[[i]]$subopt_trt_benf)[ind]
}))
rits_sub <- t(sapply(1:n_iter, function(i){
  cumsum(rits_sim[[i]]$subopt_trt_benf)[ind]
}))


# Combine matrices into a single data frame
df <- data.frame(
  Value = c(rand_sub, ts_sub, rits_sub),
  Method = rep(c("rand", "ts", "rits"), 
               each = nrow(rand_sub)*ncol(rand_sub)),
  Column = rep(rep(1:ncol(rand_sub), each = nrow(rand_sub), 
                   times = 3))
)
library(ggplot2)

# Plot grouped boxplot
sim_sub_plot2 <- ggplot(df, aes(x = factor(Column), 
                                y = Value, fill = Method)) +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.7) +
  labs(x = "Number of Participants", y = "Cumulative Sub-optimal Treatment Count", 
       fill = "Method") +
  scale_x_discrete(labels = ind) + ggtitle("Efficacy") +
  scale_fill_manual(
    values = c("rand" = "#CC79A7", "ts" = "#0072B2", "rits" = "#D55E00"),
    labels = c("ts" = "TS", "rand" = "Rand", "rits" = "RiTS")
  ) 


## Safety Sub-optimal Treatment Count
rand_sub <- t(sapply(1:n_iter, function(i){
  cumsum(rand_sim[[i]]$subopt_trt_safe)[ind]
}))
ts_sub <- t(sapply(1:n_iter, function(i){
  cumsum(ts_sim[[i]]$subopt_trt_safe)[ind]
}))
rits_sub <- t(sapply(1:n_iter, function(i){
  cumsum(rits_sim[[i]]$subopt_trt_safe)[ind]
}))


# Combine matrices into a single data frame
df <- data.frame(
  Value = c(rand_sub, ts_sub, rits_sub),
  Method = rep(c("rand", "ts", "rits"), 
               each = nrow(rand_sub)*ncol(rand_sub)),
  Column = rep(rep(1:ncol(rand_sub), each = nrow(rand_sub), 
                   times = 3))
)
library(ggplot2)

# Plot grouped boxplot
sim_sub_plot3 <- ggplot(df, aes(x = factor(Column), 
                                y = Value, fill = Method)) +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.7) +
  labs(x = "Number of Participants", y = "Cumulative Sub-optimal Treatment Count", 
       fill = "Method") +
  scale_x_discrete(labels = ind)  + ggtitle("Safety") +
  scale_fill_manual(
    values = c("rand" = "#CC79A7", "ts" = "#0072B2", "rits" = "#D55E00"),
    labels = c("ts" = "TS", "rand" = "Rand", "rits" = "RiTS")
  ) 

sim_subopt_plot <- sim_sub_plot1 + sim_sub_plot2 + sim_sub_plot3 +
  plot_layout(ncol = 3, guides = "collect")
ggsave("plot/subopt_sim_bwplot.jpg", height = 4, width = 12, units = "in")


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
ggsave("plot/freq_Arm_alloc_sim.jpg", width = 8, height = 4)



## Width
# Arm 2
contr_rand_width2 <- t(sapply(1:n_iter, function(iter){
  rand_sim[[iter]]$contr[ind, 1, 3] - rand_sim[[iter]]$contr[ind, 1, 2]
}))

contr_ts_width2 <- t(sapply(1:n_iter, function(iter){
  ts_sim[[iter]]$contr[ind, 1, 3] - ts_sim[[iter]]$contr[ind, 1, 2]
}))

contr_rits_width2 <- t(sapply(1:n_iter, function(iter){
  rits_sim[[iter]]$contr[ind, 1, 3] - ts_sim[[iter]]$contr[ind, 1, 2]
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
       fill = "Arm") + #ylim(0, 15) +
  scale_fill_manual(
    values = c("rand" = "#CC79A7", "ts" = "#0072B2", "rits" = "#D55E00"),
    labels = c("ts" = "TS", "rand" = "Rand", "rits" = "RiTS")
  ) +
  scale_x_discrete(labels = ind) + ggtitle("Width for Arm 2 - Arm 1") 


# Arm 3
contr_rand_width3 <- t(sapply(1:n_iter, function(iter){
  rand_sim[[iter]]$contr[ind, 2, 3] - rand_sim[[iter]]$contr[ind, 2, 2]
}))

contr_ts_width3 <- t(sapply(1:n_iter, function(iter){
  ts_sim[[iter]]$contr[ind, 2, 3] - ts_sim[[iter]]$contr[ind, 2, 2]
}))

contr_rits_width3 <- t(sapply(1:n_iter, function(iter){
  rits_sim[[iter]]$contr[ind, 2, 3] - ts_sim[[iter]]$contr[ind, 2, 2]
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
       fill = "Arm") + #ylim(0, 15) +
  scale_fill_manual(
    values = c("rand" = "#CC79A7", "ts" = "#0072B2", "rits" = "#D55E00"),
    labels = c("ts" = "TS", "rand" = "Rand", "rits" = "RiTS")
  ) +
  scale_x_discrete(labels = ind) + ggtitle("Width for Arm 3 - Arm 1") 


# Arm 4
contr_rand_width4 <- t(sapply(1:n_iter, function(iter){
  rand_sim[[iter]]$contr[ind, 3, 3] - rand_sim[[iter]]$contr[ind, 3, 2]
}))

contr_ts_width4 <- t(sapply(1:n_iter, function(iter){
  ts_sim[[iter]]$contr[ind, 3, 3] - ts_sim[[iter]]$contr[ind, 3, 2]
}))

contr_rits_width4 <- t(sapply(1:n_iter, function(iter){
  rits_sim[[iter]]$contr[ind, 3, 3] - ts_sim[[iter]]$contr[ind, 3, 2]
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
       fill = "Arm") + #ylim(0, 15) +
  scale_fill_manual(
    values = c("rand" = "#CC79A7", "ts" = "#0072B2", "rits" = "#D55E00"),
    labels = c("ts" = "TS", "rand" = "Rand", "rits" = "RiTS")
  ) +
  scale_x_discrete(labels = ind) + ggtitle("Width for Arm 4 - Arm 1") 

sim_wid_plot <- sim_wid_plot1 + sim_wid_plot2 + sim_wid_plot3 +
  plot_layout(ncol = 3, guides = "collect")
ggsave("plot/width_sim_bwplot.jpg", height = 4, width = 12, units = "in")


## bias
# Arm 2
contr_rand_bias2 <- t(sapply(1:n_iter, function(iter){
  rand_sim[[iter]]$contr[ind, 1, 1] - contr_true[1]
}))

contr_ts_bias2 <- t(sapply(1:n_iter, function(iter){
  ts_sim[[iter]]$contr[ind, 1, 1] - contr_true[1]
}))

contr_rits_bias2 <- t(sapply(1:n_iter, function(iter){
  rits_sim[[iter]]$contr[ind, 1, 1] - contr_true[1]
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
  rand_sim[[iter]]$contr[ind, 2, 1] - contr_true[2]
}))

contr_ts_bias3 <- t(sapply(1:n_iter, function(iter){
  ts_sim[[iter]]$contr[ind, 2, 1] - contr_true[2]
}))

contr_rits_bias3 <- t(sapply(1:n_iter, function(iter){
  rits_sim[[iter]]$contr[ind, 2, 1] - contr_true[2]
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
  rand_sim[[iter]]$contr[ind, 3, 1] - contr_true[3]
}))

contr_ts_bias4 <- t(sapply(1:n_iter, function(iter){
  ts_sim[[iter]]$contr[ind, 3, 1] - contr_true[3]
}))

contr_rits_bias4 <- t(sapply(1:n_iter, function(iter){
  rits_sim[[iter]]$contr[ind, 3, 1] - contr_true[3]
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
ggsave("plot/bias_sim_bwplot.jpg", height = 4, width = 12, units = "in")



# stop_n_long <- pivot_longer(stop_n, cols = everything(), names_to = "Variable", values_to = "Value")
# 
# # Plot using ggplot
# ggplot(stop_n_long, aes(x = Variable, y = Value, fill = Variable)) +
#   geom_boxplot() +
#   theme(text = element_text(size = 8)) +
#   labs(x = "Method", y = "Sample Size", 
#        title = "Minimum sample sizes required by different methods to stop the trial")
# ggsave("plot/stop_trial.jpg", height = 4, width = 6, units = "in")
