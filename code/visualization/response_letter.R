ts_sim_low <- readRDS("output/ts_sim_dgp_low_min_prpn_0.05_tr_start_24.RData")
ts_sim_low <- ts_sim_low[1:100]
rand_sim_low <- readRDS("output/rand_sim_dgp_low_min_prpn_0.005_tr_start_24.RData")
rand_sim_low <- rand_sim_low[1:100]
rits_sim_low <- readRDS("output/rits_sim_dgp_low_min_prpn_0.05_tr_start_24.RData")
rits_sim_low <- rits_sim_low[1:100]

plot(ts_sim_low[[1]]$contr_standard[, 1, 1], col = "red", type = "l", 
     ylim = c(min(ts_sim_low[[1]]$contr_standard), max(ts_sim_low[[1]]$contr_standard)))
lines(ts_sim_low[[1]]$contr_standard[, 2, 1], col = "blue")
lines(ts_sim_low[[1]]$contr_standard[, 3, 1], col = "green")

lines(ts_sim_low[[1]]$contr_standard[, 1, 2], col = "red", lty=2)
lines(ts_sim_low[[1]]$contr_standard[, 1, 3], col = "red", lty=2)
lines(ts_sim_low[[1]]$contr_standard[, 2, 2], col = "blue", lty=2)
lines(ts_sim_low[[1]]$contr_standard[, 2, 3], col = "blue", lty=2)
lines(ts_sim_low[[1]]$contr_standard[, 3, 2], col = "green", lty=2)
lines(ts_sim_low[[1]]$contr_standard[, 3, 3], col = "green", lty=2)
abline(h = 0.1)

library(ggplot2)
custom_colors <- c("Arm 2 - Arm 1" = "steelblue", 
                   "Arm 3 - Arm 1" = "darkgreen", 
                   "Arm 4 - Arm 1" = "purple")
custom_fill_colors <- c("Arm 2 - Arm 1" = "skyblue", 
                        "Arm 3 - Arm 1" = "lightgreen", 
                        "Arm 4 - Arm 1" = "plum")
theme_set(theme_gray(base_size = 14))

i <- 5
df2 <- data.frame(x = as.numeric(dimnames(ts_sim_low[[i]]$contr)[[1]]),
                  center = ts_sim_low[[i]]$contr_standard[, 1, 1],
                  lower = ts_sim_low[[i]]$contr_standard[, 1, 2],
                  upper = ts_sim_low[[i]]$contr_standard[, 1, 3],
                  lineid = "Arm 2 - Arm 1"
                  )
df3 <- data.frame(x = as.numeric(dimnames(ts_sim_low[[i]]$contr)[[1]]),
                  center = ts_sim_low[[i]]$contr_standard[, 2, 1],
                  lower = ts_sim_low[[i]]$contr_standard[, 2, 2],
                  upper = ts_sim_low[[i]]$contr_standard[, 2, 3],
                  lineid = "Arm 3 - Arm 1"
                  )
df4 <- data.frame(x = as.numeric(dimnames(ts_sim_low[[i]]$contr)[[1]]),
                  center = ts_sim_low[[i]]$contr_standard[, 3, 1],
                  lower = ts_sim_low[[i]]$contr_standard[, 3, 2],
                  upper = ts_sim_low[[i]]$contr_standard[, 3, 3],
                  lineid = "Arm 4 - Arm 1"
                  )
df <- rbind(df2, df3, df4)
df$lineid <- factor(df$lineid)

ggplot(df, aes(x = x, group = lineid)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = lineid), alpha = 0.3) +
  geom_line(aes(y = center, colour = lineid), linewidth = 1) +
  scale_color_manual(name = "AsympCS", values = custom_colors) +
  scale_fill_manual(name = "AsympCS", values = custom_fill_colors) +
  # Customize the plot appearance
  labs(
    title = "Estimated Effect Sizes and AsympCSs at Different Stages of a Trial",
    x = "Number of Participant",
    y = "Effect Size"
  ) + geom_hline(yintercept = 0.1, color = "black", linetype = "dashed", 
                 linewidth = 0.8)
ggsave("plot/example_asympcs.jpg", height = 6, width = 12, 
       units = "in")
