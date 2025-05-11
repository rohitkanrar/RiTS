library(ggplot2); library(tidyr); library(patchwork)
# Function to generate cumulative regret plot as in Appendix A.1.
gen_cum_reg_bwplot <- function(out_rand, out_ts, out_rits, ind, criteria){
  n_iter <- length(out_rand)
  if(criteria == "Utility"){
    rand_reg <- t(sapply(1:n_iter, function(i){
      cumsum(out_rand[[i]]$regret)[ind]
    }))
    ts_reg <- t(sapply(1:n_iter, function(i){
      cumsum(out_ts[[i]]$regret)[ind]
    }))
    rits_reg <- t(sapply(1:n_iter, function(i){
      cumsum(out_rits[[i]]$regret)[ind]
    }))
  } else if(criteria == "Efficacy"){
    rand_reg <- t(sapply(1:n_iter, function(i){
      cumsum(out_rand[[i]]$regret_benf)[ind]
    }))
    ts_reg <- t(sapply(1:n_iter, function(i){
      cumsum(out_ts[[i]]$regret_benf)[ind]
    }))
    rits_reg <- t(sapply(1:n_iter, function(i){
      cumsum(out_rits[[i]]$regret_benf)[ind]
    }))
  } else if(criteria == "Safety"){
    rand_reg <- t(sapply(1:n_iter, function(i){
      cumsum(out_rand[[i]]$regret_safe)[ind]
    }))
    ts_reg <- t(sapply(1:n_iter, function(i){
      cumsum(out_ts[[i]]$regret_safe)[ind]
    }))
    rits_reg <- t(sapply(1:n_iter, function(i){
      cumsum(out_rits[[i]]$regret_safe)[ind]
    }))
  } else{
    stop("Criteria should be one of Utility, Efficacy or Safety.")
  }
  
  df <- data.frame(
    Value = c(rand_reg, ts_reg, rits_reg),
    Method = rep(c("rand", "ts", "rits"), 
                 each = nrow(rand_reg)*ncol(rand_reg)),
    Column = rep(rep(1:ncol(rand_reg), each = nrow(rand_reg), 
                     times = 3))
  )
  sim_reg_plot <- ggplot(df, aes(x = factor(Column), 
                                  y = Value, fill = Method)) +
    geom_boxplot(position = position_dodge(width = 0.8), width = 0.7) +
    labs(x = "Number of Participants", y = "Cumulative Regret", 
         fill = "Method") +
    scale_x_discrete(labels = ind) + ggtitle(criteria) +
    scale_fill_manual(
      values = c("rand" = "#CC79A7", "ts" = "#0072B2", "rits" = "#D55E00"),
      labels = c("ts" = "TS", "rand" = "Rand", "rits" = "RiTS")
    )
  return(sim_reg_plot)
}

# Function to generate freq of arm alloc plot as in Appendix A.1.
gen_freq_arm_alloc <- function(out_rand, out_ts, out_rits, ylims, titl){
  n_iter <- length(out_rand)
  trt_rand <- t(sapply(1:n_iter, function(i) table(out_rand[[i]]$trt)))
  trt_ts <- t(sapply(1:n_iter, function(i) table(out_ts[[i]]$trt)))
  trt_rits <- t(sapply(1:n_iter, function(i) table(out_rits[[i]]$trt)))
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
    ) + ylim(ylims) + ggtitle(titl)
}

# Function to generate width box plot as in Appendix A.2.
gen_width_bwplot <- function(out_rand, out_ts, out_rits, ate_ind, arm, ylims){
  n_iter <- length(out_rand)
  titl <- paste("Arm ",  arm, " - Arm 1", sep = "")
  contr_std_width <- t(sapply(1:n_iter, function(iter){
    out_rand[[iter]]$contr_standard[ate_ind, (arm - 1), 3] - 
      out_rand[[iter]]$contr_standard[ate_ind, (arm - 1), 2]
  }))
  contr_rand_width <- t(sapply(1:n_iter, function(iter){
    out_rand[[iter]]$contr[ate_ind, (arm - 1), 3] - 
      out_rand[[iter]]$contr[ate_ind, (arm - 1), 2]
  }))
  contr_ts_width <- t(sapply(1:n_iter, function(iter){
    out_ts[[iter]]$contr[ate_ind, (arm - 1), 3] - 
      out_ts[[iter]]$contr[ate_ind, (arm - 1), 2]
  }))
  contr_rits_width <- t(sapply(1:n_iter, function(iter){
    out_rits[[iter]]$contr[ate_ind, (arm - 1), 3] - 
      out_rits[[iter]]$contr[ate_ind, (arm - 1), 2]
  }))
  
  df <- data.frame(
    Width = c(contr_std_width, contr_rand_width, contr_ts_width, 
              contr_rits_width),
    Method = rep(c("std", "rand", "ts", "rits"), 
                 each = nrow(contr_rand_width)*ncol(contr_rand_width)),
    Column = rep(rep(1:ncol(contr_rand_width), each = nrow(contr_rand_width), 
                     times = 4))
  )
  wid_plot <- ggplot(df, aes(x = factor(Column), 
                                  y = Width, fill = Method)) +
    geom_boxplot(position = position_dodge(width = 0.8), width = 0.7) +
    labs(x = "Number of Participants", y = "Width", 
         fill = "Arm") + ylim(ylims) +
    scale_fill_manual(
      values = c("std" = "#009E73", "rand" = "#CC79A7", "ts" = "#0072B2", 
                 "rits" = "#D55E00"),
      labels = c("std" = "Std", "ts" = "TS", "rand" = "Rand", "rits" = "RiTS")
    ) +
    scale_x_discrete(labels = ind) + ggtitle(titl) 
  return(wid_plot)
}

# Function to generate bias box plot as in Appendix A.2.
gen_bias_bwplot <- function(out_rand, out_ts, out_rits, ate_ind, arm, 
                            contr_true, ylims){
  n_iter <- length(out_rand)
  titl <- paste("Arm ",  arm, " - Arm 1", sep = "")
  contr_std_bias <- t(sapply(1:n_iter, function(iter){
    out_rand[[iter]]$contr_standard[ate_ind, arm-1, 1] - contr_true[arm-1]
  }))
  contr_rand_bias <- t(sapply(1:n_iter, function(iter){
    out_rand[[iter]]$contr[ate_ind, arm-1, 1] - contr_true[arm-1]
  }))
  contr_ts_bias <- t(sapply(1:n_iter, function(iter){
    out_ts[[iter]]$contr[ate_ind, arm-1, 1] - contr_true[arm-1]
  }))
  contr_rits_bias <- t(sapply(1:n_iter, function(iter){
    out_rits[[iter]]$contr[ate_ind, arm-1, 1] - contr_true[arm-1]
  }))
  
  df <- data.frame(
    Bias = c(contr_std_bias, contr_rand_bias, contr_ts_bias, contr_rits_bias),
    Method = rep(c("std", "rand", "ts", "rits"), 
                 each = nrow(contr_rand_bias)*ncol(contr_rand_bias)),
    Column = rep(rep(1:ncol(contr_rand_bias), each = nrow(contr_rand_bias), 
                     times = 4))
  )
  
  sim_bias_plot1 <- ggplot(df, aes(x = factor(Column), 
                                   y = Bias, fill = Method)) +
    geom_boxplot(position = position_dodge(width = 0.8), width = 0.7) +
    labs(x = "Number of Participants", y = "Bias", 
         fill = "Arm") + ylim(ylims) +
    scale_fill_manual(
      values = c("std" = "#009E73", "rand" = "#CC79A7", "ts" = "#0072B2", 
                 "rits" = "#D55E00"),
      labels = c("std" = "Std", "ts" = "TS", "rand" = "Rand", "rits" = "RiTS")
    ) +
    scale_x_discrete(labels = ind) + ggtitle(titl) 
}

# Function to generate cumulative miscoverage plot as in Appendix A.3.

gen_cum_miscov_plot <- function(out, ate_true, contr_true, alpha, 
                                ate_start, titl){
  n_iter <- length(out); K <- length(ate_true); N <- length(out[[1]]$trt)
  cum_miscov <- get_cum_mis_cov(out, mu_true = mu_true, 
                                contr_true = contr_true, delay = 0)
  
  df <- as.data.frame(cum_miscov[[2]] / n_iter)
  colnames(df) <- paste("Arm", 2:K)
  df$obs <- as.numeric(dimnames(out[[1]]$ate)[[1]])  
  df_long <- melt(df, id.vars = "obs", variable.name = "Arm", 
                  value.name = "Miscov")
  
  ggplot(df_long, aes(x = obs, y = Miscov, color = Arm)) +
    geom_line(linewidth = 0.5) + ylim(c(0, 0.1)) + 
    geom_hline(yintercept = sim_choice$alpha/(K-1), linetype = "dashed", 
               color = "blue") +
    labs(x = "Patient", y = "Cumulative Miscoverage", 
         title = titl) +
    theme(text = element_text(size = 10)) +
    scale_color_manual(
      values = c("Arm 1" = "#E69F00", "Arm 2" = "#56B4E9", 
                 "Arm 3" = "#009E73", "Arm 4" = "#CC79A7")
    )
}

gen_winner_curve <- function(rand_out, ts_out, rits_out, titl = "High SNR", 
                             true_best_arm = 4){
  n_iter <- length(rand_out)
  true_best_arm <- true_best_arm - 1
  n_comp <- dim(rand_sim_high[[1]]$contr)[1]
  power_std <- numeric(n_comp)
  power_rand <- numeric(n_comp)
  power_ts <- numeric(n_comp)
  power_rits <- numeric(n_comp)
  for(i in 1:n_iter){
    power_std <- power_std + as.numeric(
      sapply(1:n_comp, function(t) which.max(rand_out[[i]]$contr_standard[t, , 1]) == true_best_arm)
    )
    power_rand <- power_rand + as.numeric(
      sapply(1:n_comp, function(t) which.max(rand_out[[i]]$contr[t, , 1]) == true_best_arm)
    )
    power_ts <- power_ts + as.numeric(
      sapply(1:n_comp, function(t) which.max(ts_out[[i]]$contr[t, , 1]) == true_best_arm)
    )
    power_rits <- power_rits + as.numeric(
      sapply(1:n_comp, function(t) which.max(rits_out[[i]]$contr[t, , 1]) == true_best_arm)
    )
  }
  power_std <- power_std / n_iter
  power_rand <- power_rand / n_iter
  power_ts <- power_ts / n_iter
  power_rits <- power_rits / n_iter
  
  df <- data.frame(std = power_std, rand = power_rand, ts = power_ts,
                   rits = power_rits, 
                   participant = as.numeric(dimnames(rand_out[[1]]$contr)[[1]]))
  df_long <- tidyr::gather(df, key = "Methods", value = "value", -participant)
  out_plot <- ggplot(df_long, aes(x = participant, y = value, color = Methods)) +
    geom_line() + labs(x = "Participant", y = "Proportion") +
    scale_color_manual(
      values = c("rand" = "#CC79A7", "ts" = "#0072B2", "rits" = "#D55E00", 
                 "std" = "#009E73"),
      labels = c("ts" = "TS", "rand" = "Rand", "rits" = "RiTS", "std" = "Std")
    ) + ggtitle(titl)
  return(out_plot)
}

gen_power_curve <- function(rand_out, ts_out, rits_out, titl = "High SNR", 
                            min_thresh = 0.1){
  n_iter <- length(rand_out)
  n_comp <- dim(rand_sim_high[[1]]$contr)[1]
  
  power_std <- numeric(n_comp)
  power_rand <- numeric(n_comp)
  power_ts <- numeric(n_comp)
  power_rits <- numeric(n_comp)
  
  for(i in 1:n_iter){
    power_std <- power_std + sapply(1:n_comp, function(t){
      rej <- sapply(1:3, function(k){
        zero_in_intv(rand_out[[i]]$contr_standard[t, k, 2:3], zero = min_thresh)
      })
      as.numeric(FALSE %in% rej)
    })
    power_rand <- power_rand + sapply(1:n_comp, function(t){
      rej <- sapply(1:3, function(k){
        zero_in_intv(rand_out[[i]]$contr[t, k, 2:3], zero = min_thresh)
      })
      as.numeric(FALSE %in% rej)
    })
    power_ts <- power_ts + sapply(1:n_comp, function(t){
      rej <- sapply(1:3, function(k){
        zero_in_intv(ts_out[[i]]$contr[t, k, 2:3], zero = min_thresh)
      })
      as.numeric(FALSE %in% rej)
    })
    power_rits <- power_rits + sapply(1:n_comp, function(t){
      rej <- sapply(1:3, function(k){
        zero_in_intv(rits_out[[i]]$contr[t, k, 2:3], zero = min_thresh)
      })
      as.numeric(FALSE %in% rej)
    })
  }
  
  power_std <- power_std / n_iter
  power_rand <- power_rand / n_iter
  power_ts <- power_ts / n_iter
  power_rits <- power_rits / n_iter
  
  df <- data.frame(std = power_std, rand = power_rand, ts = power_ts,
                   rits = power_rits, 
                   participant = as.numeric(dimnames(rand_out[[1]]$contr)[[1]]))
  df_long <- tidyr::gather(df, key = "Methods", value = "value", -participant)
  out_plot <- ggplot(df_long, aes(x = participant, y = value, color = Methods)) +
    geom_line() + labs(x = "Participant", y = "Power") +
    scale_color_manual(
      values = c("rand" = "#CC79A7", "ts" = "#0072B2", "rits" = "#D55E00", 
                 "std" = "#009E73"),
      labels = c("ts" = "TS", "rand" = "Rand", "rits" = "RiTS", "std" = "Std")
    ) + ggtitle(titl)
  return(out_plot)
}
