source("code/visualization/viz_function.R")

out_dir <- "output_git/"

if (!dir.exists("plot")) {
  dir.create("plot", recursive = TRUE)
}

required_files <- c(
  paste0(out_dir, "ts_sim_dgp_high_min_prpn_0.05_tr_start_24.RData"),
  paste0(out_dir, "rand_sim_dgp_high_min_prpn_0.005_tr_start_24.RData"),
  paste0(out_dir, "rits_sim_dgp_high_min_prpn_0.05_tr_start_24.RData"),
  paste0(out_dir, "ts_sim_dgp_low_min_prpn_0.05_tr_start_24.RData"),
  paste0(out_dir, "rand_sim_dgp_low_min_prpn_0.005_tr_start_24.RData"),
  paste0(out_dir, "rits_sim_dgp_low_min_prpn_0.05_tr_start_24.RData")
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    "The following files are missing. Please rerun the Trial Optimization and Trial Evaluation chunks first:\n",
    paste(missing_files, collapse = "\n")
  )
}

# High SNR
ts_sim_high <- readRDS(
  paste0(out_dir, "ts_sim_dgp_high_min_prpn_0.05_tr_start_24.RData")
)

rand_sim_high <- readRDS(
  paste0(out_dir, "rand_sim_dgp_high_min_prpn_0.005_tr_start_24.RData")
)

rits_sim_high <- readRDS(
  paste0(out_dir, "rits_sim_dgp_high_min_prpn_0.05_tr_start_24.RData")
)

# Low SNR
ts_sim_low <- readRDS(
  paste0(out_dir, "ts_sim_dgp_low_min_prpn_0.05_tr_start_24.RData")
)

rand_sim_low <- readRDS(
  paste0(out_dir, "rand_sim_dgp_low_min_prpn_0.005_tr_start_24.RData")
)

rits_sim_low <- readRDS(
  paste0(out_dir, "rits_sim_dgp_low_min_prpn_0.05_tr_start_24.RData")
)

sim_choice <- readRDS("metadata/sim_choice.RData")

n_iter <- length(ts_sim_high)
N <- length(ts_sim_high[[1]]$trt)
K <- length(unique(rand_sim_high[[1]]$trt))

# Frequency of Arm Allocation Plot
df_high <- gen_freq_arm_alloc_df(
  out_rand = rand_sim_high,
  out_ts = ts_sim_high,
  out_rits = rits_sim_high
)

df_low <- gen_freq_arm_alloc_df(
  out_rand = rand_sim_low,
  out_ts = ts_sim_low,
  out_rits = rits_sim_low
)

df_high[["dgp"]] <- "High-SNR"
df_low[["dgp"]] <- "Low-SNR"

# The original gen_freq_arm_alloc() function requires a df_null argument.
# For this small example, we only generated High-SNR and Low-SNR cases,
# so we pass df_null = NULL explicitly.
alloc_plot <- gen_freq_arm_alloc(
  df_high = df_high,
  df_low = df_low,
  df_null = NULL,
  ylims = c(0, 150)
)

# Proportion of trials where Arm 4 is the winner
winner_high <- gen_winner_curve_df(
  rand_out = rand_sim_high,
  ts_out = ts_sim_high,
  rits_out = rits_sim_high,
  true_best_arm = 4,
  include_ipw = FALSE
)

winner_low <- gen_winner_curve_df(
  rand_out = rand_sim_low,
  ts_out = ts_sim_low,
  rits_out = rits_sim_low,
  true_best_arm = 4,
  include_ipw = FALSE
)

winner_high[["dgp"]] <- "High-SNR"
winner_low[["dgp"]] <- "Low-SNR"

winner <- rbind(winner_high, winner_low)
winner[["type"]] <- "Winner (Arm 4)"

# Proportion of trials where stopping criteria is met
power_high <- gen_power_curve_df(
  rand_out = rand_sim_high,
  ts_out = ts_sim_high,
  rits_out = rits_sim_high,
  min_thresh = 0.1,
  include_ipw = FALSE
)

power_low <- gen_power_curve_df(
  rand_out = rand_sim_low,
  ts_out = ts_sim_low,
  rits_out = rits_sim_low,
  min_thresh = 0.1,
  include_ipw = FALSE
)

power_high[["dgp"]] <- "High-SNR"
power_low[["dgp"]] <- "Low-SNR"

power_df <- rbind(power_high, power_low)
power_df[["type"]] <- "Stopping Criteria"

metric_plots <- gen_metrics_plot(
  df_winner = winner,
  df_power = power_df
)

metric_alloc_plot_example <- alloc_plot + metric_plots + patchwork::plot_layout(ncol = 2)

ggsave(
  filename = "plot/metric_alloc_plot_example.jpg",
  plot = metric_alloc_plot_example,
  height = 4,
  width = 15,
  units = "in"
)

metric_alloc_plot_example
