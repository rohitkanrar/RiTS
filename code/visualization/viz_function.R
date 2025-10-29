library(ggplot2); library(tidyr); library(patchwork)
update_geom_defaults("boxplot", list(
  outlier.size = 0.25,
  color = "black",
  linewidth = 0.25
))
# Function to generate cumulative regret plot as in Appendix A.1.
gen_cum_reg_df <- function(out_rand, out_ts, out_rits, ind, criteria){
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
  df["criteria"] <- criteria
  return(df)
}
gen_cum_reg_bwplot <- function(df_high, df_low, df_null, ind){
  df <- rbind(df_high, df_low, df_null)
  sim_reg_plot <- ggplot(df, aes(x = factor(Column), 
                                  y = Value, fill = Method)) +
    geom_boxplot(position = position_dodge(width = 0.8), width = 0.7, 
                 outlier.size = 0.5, outlier.shape = NA) +
    facet_grid(dgp ~ criteria, scales = "free_y") + 
    labs(x = "Number of participants", y = "Cumulative regret", 
         fill = "Method") +
    scale_x_discrete(labels = ind) +
    scale_fill_manual(
      values = c("rand" = "#CC79A7", "ts" = "#0072B2", "rits" = "#D55E00"),
      labels = c("ts" = "TS", "rand" = "Rand", "rits" = "RiTS")
    ) + theme(legend.position = "top")
  return(sim_reg_plot)
}

# Function to generate freq of arm alloc plot as in Appendix A.1.
gen_freq_arm_alloc_df <- function(out_rand, out_ts, out_rits){
  # browser()
  n_iter <- length(out_rand); K <- length(unique(out_rand[[1]]$trt))
  trt_rand <- t(sapply(1:n_iter, function(i) table(out_rand[[i]]$trt)))
  trt_ts <- t(sapply(1:n_iter, function(i) table(out_ts[[i]]$trt)))
  trt_rits <- t(sapply(1:n_iter, function(i) table(out_rits[[i]]$trt)))
  df <- data.frame(
    Frequency = c(trt_rand, trt_ts, trt_rits),
    Arm = rep(rep(1:K, each = nrow(trt_rand)), times = 3),
    Method = rep(c("rand", "ts", "rits"), each = nrow(trt_rand)*ncol(trt_rand))
  )
  return(df)
}
gen_freq_arm_alloc <- function(df_high, df_low, df_null, ylims){
  if(!is.null(df_low)){
    df <- rbind(df_high, df_low, df_null)
  } else{
    df <- df_high
  }
  out_plot <- ggplot(df, aes(x = factor(Arm), 
                 y = Frequency, fill = Method)) +
    geom_boxplot(outlier.size = 0.5, outlier.shape = NA) +
    labs(x = "Arm", y = "Number of participants", 
         fill = "Method") + 
    scale_fill_manual(
      values = c("rand" = "#CC79A7", "ts" = "#0072B2", "rits" = "#D55E00"),
      labels = c("ts" = "TS", "rand" = "Rand", "rits" = "RiTS")
    ) + ylim(ylims) + theme(legend.position = "top")
  if(!is.null(df_low)){
    out_plot <- out_plot + facet_wrap(~dgp)
  }
  return(out_plot)
}

gen_width_df <- function(out_rand, out_ts, out_rits, ate_ind){
  n_iter <- length(out_rand); K <- length(unique(out_rand[[1]]$trt))
  ind <- dimnames(out_rand[[1]]$contr)[[1]][ate_ind]
  df <- data.frame()
  for(arm in 2:K){
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
    df_tmp <- data.frame(
      Width = c(contr_std_width, contr_rand_width, contr_ts_width, 
                contr_rits_width),
      Method = rep(c("ttest", "rand", "ts", "rits"), 
                   each = nrow(contr_rand_width)*ncol(contr_rand_width)),
      Column = rep(rep(1:ncol(contr_rand_width), each = nrow(contr_rand_width), 
                       times = 4)),
      Arm = paste("Arm", arm, "- Arm 1")
    )
    df <- rbind(df, df_tmp)
  }
  df
}
require(ggh4x)
# Function to generate width box plot as in Appendix A.2.
gen_width_bwplot <- function(df_high, df_low, ind, ylims){
  df_high["dgp"] <- "High-SNR"; df_low["dgp"] <- "Low-SNR"
  df <- rbind(df_high, df_low)
  y_limits <- list(
    scale_y_continuous(limits = ylims[[1]]),
    scale_y_continuous(limits = ylims[[2]])
  )
  wid_plot <- ggplot(df, aes(x = factor(Column), 
                                  y = Width, fill = Method)) + 
    facet_grid(dgp ~ Arm, scales = "free_y") +
    geom_boxplot(position = position_dodge(width = 0.8), width = 0.7,
                 outlier.size = 0.5, outlier.shape = NA) +
    labs(x = "Number of participants", y = "Width", 
         fill = "Method") +
    scale_fill_manual(
      values = c("ttest" = "#009E73", "rand" = "#CC79A7", "ts" = "#0072B2", 
                 "rits" = "#D55E00"),
      labels = c("ttest" = "Rand-OF", "ts" = "TS-AIPW", 
                 "rand" = "Rand-AIPW", "rits" = "RiTS-AIPW")
    ) +
    scale_x_discrete(labels = ind) + theme(legend.position = "top")
  if(!is.null(ylims)){
    wid_plot <- wid_plot + facetted_pos_scales(y = y_limits)
  }
  return(wid_plot)
}

# Function to generate bias box plot as in Appendix A.2.
gen_bias_df <- function(out_rand, out_ts, out_rits, ate_ind, contr_true){
  n_iter <- length(out_rand); K <- length(unique(out_rand[[1]]$trt))
  ind <- dimnames(out_rand[[1]]$contr)[[1]][ate_ind]
  df <- data.frame()
  for(arm in 2:K){
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
    
    df_tmp <- data.frame(
      Bias = c(contr_std_bias, contr_rand_bias, contr_ts_bias, contr_rits_bias),
      Method = rep(c("ttest", "rand", "ts", "rits"), 
                   each = nrow(contr_rand_bias)*ncol(contr_rand_bias)),
      Column = rep(rep(1:ncol(contr_rand_bias), each = nrow(contr_rand_bias), 
                       times = 4)),
      Arm = paste("Arm", arm, "- Arm 1")
    )
    df <- rbind(df, df_tmp)
  }
  df
}
gen_bias_bwplot <- function(df_high, df_low, ind, ylims){
  df_high["dgp"] <- "High-SNR"; df_low["dgp"] <- "Low-SNR"
  df <- rbind(df_high, df_low)
  y_limits <- list(
    scale_y_continuous(limits = ylims[[1]]),
    scale_y_continuous(limits = ylims[[2]])
  )
  bias_plot <- ggplot(df, aes(x = factor(Column), 
                                   y = Bias, fill = Method)) +
    geom_boxplot(position = position_dodge(width = 0.8), width = 0.7,
                 outlier.size = 0.5, outlier.shape = NA) +
    facet_grid(dgp ~ Arm, scales = "free_y") +
    labs(x = "Number of participants", y = "Bias", 
         fill = "Method") +
    scale_fill_manual(
      values = c("ttest" = "#009E73", "rand" = "#CC79A7", "ts" = "#0072B2", 
                 "rits" = "#D55E00"),
      labels = c("ttest" = "Rand-OF", "ts" = "TS-AIPW", 
                 "rand" = "Rand-AIPW", "rits" = "RiTS-AIPW")
    ) +
    scale_x_discrete(labels = ind) + theme(legend.position = "top")
  if(!is.null(ylims)){
    bias_plot <- bias_plot + facetted_pos_scales(y = y_limits)
  }
  bias_plot
}

# Function to generate cumulative miscoverage plot as in Appendix A.3.
library(reshape2)
gen_cum_miscov_df <- function(out_rand, out_ts, out_rits, mu_true, contr_true,
                              delay_aipw = 77){
  n_iter <- length(out_rand); K <- length(mu_true); N <- length(out_rand[[1]]$trt)
  cum_miscov <- get_cum_mis_cov(out_rand, mu_true = mu_true, 
                                contr_true = contr_true, delay_ipw = 0,
                                delay_aipw = delay_aipw, 
                                need_ipw = FALSE, need_std = TRUE)
  df_rand <- as.data.frame(cum_miscov[["contr"]] / n_iter)
  colnames(df_rand) <- paste("Arm", 2:K, "- Arm 1")
  obs <- dimnames(out_rand[[1]]$contr)[[1]] 
  obs <- obs[(delay_aipw+1):length(obs)]
  df_rand$obs <- as.numeric(obs)
  df_rand <- melt(df_rand, id.vars = "obs", variable.name = "Arm", 
                  value.name = "Miscov")
  df_rand["Method"] <- "Rand-AIPW"
  
  df_std <- as.data.frame(cum_miscov[["contr_std"]] / n_iter)
  colnames(df_std) <- paste("Arm", 2:K, "- Arm 1")
  df_std$obs <- as.numeric(obs)
  df_std <- melt(df_std, id.vars = "obs", variable.name = "Arm", 
                  value.name = "Miscov")
  df_std["Method"] <- "Rand-OF"
  
  cum_miscov <- get_cum_mis_cov(out_ts, mu_true = mu_true, 
                                contr_true = contr_true, delay_ipw = 0,
                                delay_aipw = delay_aipw, need_ipw = FALSE, 
                                need_std = FALSE)
  df_ts <- as.data.frame(cum_miscov[["contr"]] / n_iter)
  colnames(df_ts) <- paste("Arm", 2:K, "- Arm 1")
  df_ts$obs <- as.numeric(obs)
  df_ts <- melt(df_ts, id.vars = "obs", variable.name = "Arm", 
                  value.name = "Miscov")
  df_ts["Method"] <- "TS-AIPW"
  
  cum_miscov <- get_cum_mis_cov(out_rits, mu_true = mu_true, 
                                contr_true = contr_true, delay_ipw = 0,
                                delay_aipw = delay_aipw, need_ipw = FALSE, 
                                need_std = FALSE)
  df_rits <- as.data.frame(cum_miscov[["contr"]] / n_iter)
  colnames(df_rits) <- paste("Arm", 2:K, "- Arm 1")
  df_rits$obs <- as.numeric(obs)
  df_rits <- melt(df_rits, id.vars = "obs", variable.name = "Arm", 
                  value.name = "Miscov")
  df_rits["Method"] <- "RiTS-AIPW"
  
  rbind(df_std, df_rand, df_ts, df_rits)
}

gen_cum_miscov_plot <- function(df_high, df_low, alpha, ate_start){
  df_high["dgp"] <- "High-SNR"; df_low["dgp"] <- "Low-SNR"
  df <- rbind(df_high, df_low)
  df$Method <- factor(df$Method, levels = c("Rand-OF", "Rand-AIPW", 
                                            "TS-AIPW", "RiTS-AIPW"))
  
  ggplot(df, aes(x = obs, y = Miscov, color = Arm)) +
    geom_line(linewidth = 0.5) + ylim(c(0, 0.1)) + 
    geom_hline(yintercept = sim_choice$alpha/(K-1), linetype = "dashed", 
               color = "blue") +
    facet_grid(dgp ~ Method) +
    labs(x = "Number of participants", y = "Cumulative miscoverage") +
    theme(text = element_text(size = 10)) +
    scale_color_manual(
      name = "Estimand",
      values = c("Arm 1" = "#E69F00", "Arm 2 - Arm 1" = "#56B4E9", 
                 "Arm 3 - Arm 1" = "#009E73", "Arm 4 - Arm 1" = "#CC79A7"),
      labels = c("Arm 2 - Arm 1" = expression(paste(Delta, " (2)")),
                 "Arm 3 - Arm 1" = expression(paste(Delta, " (3)")),
                 "Arm 4 - Arm 1" = expression(paste(Delta, " (4)"))
                 )
    ) + theme(legend.position = "top")
}

expand_standard_array <- function(target_array, standard_array){
  target_times_chr <- dimnames(target_array)[[1]]
  standard_times_chr <- dimnames(standard_array)[[1]]
  target_times_num <- as.numeric(target_times_chr)
  standard_times_num <- as.numeric(standard_times_chr)
  source_indices <- findInterval(target_times_num, standard_times_num)
  source_indices[source_indices == 0] <- 1
  new_array <- standard_array[source_indices, , , drop = FALSE]
  dimnames(new_array) <- dimnames(target_array)
  return(new_array)
}
expand_standard_vector <- function(target_vector, standard_vector){
  target_times_chr <- names(target_vector)
  standard_times_chr <- names(standard_vector)
  target_times_num <- as.numeric(target_times_chr)
  standard_times_num <- as.numeric(standard_times_chr)
  source_indices <- findInterval(target_times_num, standard_times_num)
  source_indices[source_indices == 0] <- 1
  new_vector <- standard_vector[source_indices, drop = FALSE]
  names(new_vector) <- names(target_vector)
  return(new_vector)
}

gen_winner_curve_df <- function(rand_out, ts_out, rits_out, 
                             true_best_arm = 4, include_std = TRUE, 
                             include_ipw = TRUE){
  n_iter <- length(rand_out)
  true_best_arm <- true_best_arm - 1
  n_comp <- dim(rand_out[[1]]$contr)[1]
  n_comp_obrien <- dim(rand_out[[1]]$contr_standard)[1]
  power_std <- numeric(n_comp_obrien)
  power_rand <- numeric(n_comp)
  power_ts <- numeric(n_comp)
  power_rits <- numeric(n_comp)
  if(include_ipw){
    power_rand_ipw <- numeric(n_comp)
    power_ts_ipw <- numeric(n_comp)
    power_rits_ipw <- numeric(n_comp)
  }
  for(i in 1:n_iter){
    if(include_std){
      power_std <- power_std + as.numeric(
        sapply(1:n_comp_obrien, 
               function(t) which.max(rand_out[[i]]$contr_standard[t, , 1]) == true_best_arm)
      )
    }
    power_rand <- power_rand + as.numeric(
      sapply(1:n_comp, function(t) which.max(rand_out[[i]]$contr[t, , 1]) == true_best_arm)
    )
    power_ts <- power_ts + as.numeric(
      sapply(1:n_comp, function(t) which.max(ts_out[[i]]$contr[t, , 1]) == true_best_arm)
    )
    power_rits <- power_rits + as.numeric(
      sapply(1:n_comp, function(t) which.max(rits_out[[i]]$contr[t, , 1]) == true_best_arm)
    )
    if(include_ipw){
      power_rand_ipw <- power_rand_ipw + as.numeric(
        sapply(1:n_comp, function(t) which.max(rand_out[[i]]$contr_ipw[t, , 1]) == true_best_arm)
      )
      power_ts_ipw <- power_ts_ipw + as.numeric(
        sapply(1:n_comp, function(t) which.max(ts_out[[i]]$contr_ipw[t, , 1]) == true_best_arm)
      )
      power_rits_ipw <- power_rits_ipw + as.numeric(
        sapply(1:n_comp, function(t) which.max(rits_out[[i]]$contr_ipw[t, , 1]) == true_best_arm)
      )
    }
    names(power_std) <- dimnames(rand_out[[i]]$contr_standard)[[1]]
    names(power_rand) <- dimnames(rand_out[[i]]$contr)[[1]]
    power_std_step <- expand_standard_vector(target_vector = power_rand, 
                                            standard_vector = power_std)
  }
  names(power_std) <- NULL; names(power_rand) <- NULL; names(power_std_step) <- NULL
  power_std <- power_std / n_iter
  power_std_step <- power_std_step / n_iter
  power_rand <- power_rand / n_iter
  power_ts <- power_ts / n_iter
  power_rits <- power_rits / n_iter
  if(include_ipw){
    power_rand_ipw <- power_rand_ipw / n_iter
    power_ts_ipw <- power_ts_ipw / n_iter
    power_rits_ipw <- power_rits_ipw / n_iter
  }
  
  df <- data.frame(rand = power_rand, ts = power_ts,
                   rits = power_rits, 
                   participant = as.numeric(dimnames(rand_out[[1]]$contr)[[1]]))
  if(include_std){
    df[["ttest"]] <- power_std_step
  }
  if(include_ipw){
    df[["rand_ipw"]] <- power_rand_ipw
    df[["ts_ipw"]] <- power_ts_ipw
    df[["rits_ipw"]] <- power_rits_ipw
  }
  df_long <- tidyr::gather(df, key = "Methods", value = "value", -participant)
  df_long
}
require(colorspace)
gen_winner_curve <- function(df_high, df_low, df_null){
  if(!is.null(df_low)){
    df <- rbind(df_high, df_low, df_null)
  } else{
    df <- df_high
  }
  out_plot <- ggplot(df, aes(x = participant, y = value, color = Methods)) +
    geom_line() + labs(x = "Number of participants", y = "Proportion of trials") + 
    scale_color_manual(
      name = "Method",
      values = c("rand" = "#CC79A7", "ts" = "#0072B2", "rits" = "#D55E00", 
                 "ttest" = "#009E73", "rand_ipw" = "#000000",
                 "ts_ipw" = "#56B4E9", "rits_ipw" = "#E69F00"),
      labels = c("ts" = "TS-AIPW", "rand" = "Rand-AIPW", "rits" = "RiTS-AIPW",
                 "ts_ipw" = "TS-IPW", "rand_ipw" = "Rand-IPW", 
                 "rits_ipw" = "RiTS-IPW", "ttest" = "Rand-OF")
    ) + theme(legend.position = "top")
  if(!is.null(df_low)){
    out_plot <- out_plot + facet_wrap(~dgp)
  }
  return(out_plot)
}

source("code/function/misc.R")
gen_power_curve_df <- function(rand_out, ts_out, rits_out, min_thresh = 0.1,
                               include_std = TRUE, include_ipw = TRUE){
  n_iter <- length(rand_out); K <- length(unique(rand_out[[1]]$trt))
  n_comp <- dim(rand_out[[1]]$contr)[1]
  n_comp_obrien <- dim(rand_out[[1]]$contr_standard)[1]
  power_std <- numeric(n_comp_obrien)
  power_rand <- numeric(n_comp)
  power_ts <- numeric(n_comp)
  power_rits <- numeric(n_comp)
  if(include_ipw){
    power_rand_ipw <- numeric(n_comp)
    power_ts_ipw <- numeric(n_comp)
    power_rits_ipw <- numeric(n_comp)
  }
  
  for(i in 1:n_iter){
    if(include_std){
      power_std <- power_std + sapply(1:n_comp_obrien, function(t){
        rej <- sapply(1:(K-1), function(k){
          zero_in_intv(rand_out[[i]]$contr_standard[t, k, 2:3], zero = min_thresh)
        })
        as.numeric(FALSE %in% rej)
      })
    }
    power_rand <- power_rand + sapply(1:n_comp, function(t){
      rej <- sapply(1:(K-1), function(k){
        zero_in_intv(rand_out[[i]]$contr[t, k, 2:3], zero = min_thresh)
      })
      as.numeric(FALSE %in% rej)
    })
    power_ts <- power_ts + sapply(1:n_comp, function(t){
      rej <- sapply(1:(K-1), function(k){
        zero_in_intv(ts_out[[i]]$contr[t, k, 2:3], zero = min_thresh)
      })
      as.numeric(FALSE %in% rej)
    })
    power_rits <- power_rits + sapply(1:n_comp, function(t){
      rej <- sapply(1:(K-1), function(k){
        zero_in_intv(rits_out[[i]]$contr[t, k, 2:3], zero = min_thresh)
      })
      as.numeric(FALSE %in% rej)
    })
    if(include_ipw){
      power_rand_ipw <- power_rand_ipw + sapply(1:n_comp, function(t){
        rej <- sapply(1:(K-1), function(k){
          zero_in_intv(rand_out[[i]]$contr_ipw[t, k, 2:3], zero = min_thresh)
        })
        as.numeric(FALSE %in% rej)
      })
      power_ts_ipw <- power_ts_ipw + sapply(1:n_comp, function(t){
        rej <- sapply(1:(K-1), function(k){
          zero_in_intv(ts_out[[i]]$contr_ipw[t, k, 2:3], zero = min_thresh)
        })
        as.numeric(FALSE %in% rej)
      })
      power_rits_ipw <- power_rits_ipw + sapply(1:n_comp, function(t){
        rej <- sapply(1:(K-1), function(k){
          zero_in_intv(rits_out[[i]]$contr_ipw[t, k, 2:3], zero = min_thresh)
        })
        as.numeric(FALSE %in% rej)
      })
    }
    names(power_std) <- dimnames(rand_out[[i]]$contr_standard)[[1]]
    names(power_rand) <- dimnames(rand_out[[i]]$contr)[[1]]
    power_std_step <- expand_standard_vector(target_vector = power_rand, 
                                             standard_vector = power_std)
  }
  power_std <- power_std / n_iter
  power_std_step <- power_std_step / n_iter
  power_rand <- power_rand / n_iter
  power_ts <- power_ts / n_iter
  power_rits <- power_rits / n_iter
  if(include_ipw){
    power_rand_ipw <- power_rand_ipw / n_iter
    power_ts_ipw <- power_ts_ipw / n_iter
    power_rits_ipw <- power_rits_ipw / n_iter
  }
  
  df <- data.frame(rand = power_rand, ts = power_ts,
                   rits = power_rits,  
                   participant = as.numeric(dimnames(rand_out[[1]]$contr)[[1]]))
  if(include_std){
    df[["ttest"]] <- power_std_step
  }
  if(include_ipw){
    df[["rand_ipw"]] <- power_rand_ipw
    df[["ts_ipw"]] <- power_ts_ipw
    df[["rits_ipw"]] <- power_rits_ipw
  }
  df_long <- tidyr::gather(df, key = "Methods", value = "value", -participant)
  return(df_long)
}

gen_power_curve <- function(df_high, df_low, df_null){
  if(!is.null(df_low)){
    df <- rbind(df_high, df_low, df_null)
  } else{
    df <- df_high
  }
  out_plot <- ggplot(df, aes(x = participant, y = value, color = Methods)) +
    geom_line() + labs(x = "Number of participants", y = "Power") +
    scale_color_manual(
      name = "Method",
      values = c("rand" = "#CC79A7", "ts" = "#0072B2", "rits" = "#D55E00", 
                 "ttest" = "#009E73", "rand_ipw" = "#000000",
                 "ts_ipw" = "#56B4E9", "rits_ipw" = "#E69F00"),
      labels = c("ts" = "TS-AIPW", "rand" = "Rand-AIPW", "rits" = "RiTS-AIPW",
                 "ts_ipw" = "TS-IPW", "rand_ipw" = "Rand-IPW", 
                 "rits_ipw" = "RiTS-IPW", "ttest" = "Rand-OF")
    ) + theme(legend.position = "top")
  if(!is.null(df_low)){
    out_plot <- out_plot + facet_wrap(~dgp)
  }
  return(out_plot)
}

gen_metrics_plot <- function(df_winner, df_power, dgp_exists = TRUE){
  df <- rbind(df_winner, df_power)
  out_plot <- ggplot(df, aes(x = participant, y = value, color = Methods)) +
    geom_line() + labs(x = "Number of participants", y = "Proportion of trials") + 
    scale_color_manual(
      name = "Method",
      values = c("rand" = "#CC79A7", "ts" = "#0072B2", "rits" = "#D55E00", 
                 "ttest" = "#009E73", "rand_ipw" = "#000000",
                 "ts_ipw" = "#56B4E9", "rits_ipw" = "#E69F00"),
      labels = c("ts" = "TS-AIPW", "rand" = "Rand-AIPW", "rits" = "RiTS-AIPW",
                 "ts_ipw" = "TS-IPW", "rand_ipw" = "Rand-IPW", 
                 "rits_ipw" = "RiTS-IPW", "ttest" = "Rand-OF")
    ) + theme(legend.position = "top")
  if(dgp_exists){
    out_plot <- out_plot + facet_grid(type~dgp)
  } else{
    out_plot <- out_plot + facet_wrap(~type)
  }
  return(out_plot)
}

gen_summary_for_table <- function(sim, K, ate_ind, contr_true, need_std = FALSE,
                                  need_ipw = FALSE){
  n_iter <- length(sim)
  width_contr <- matrix(0, K-1, length(ate_ind))
  bias_contr <- matrix(0, K-1, length(ate_ind))
  rmse_contr <- matrix(0, K-1, length(ate_ind))
  if(need_std){
    width_contr_std <- matrix(0, K-1, length(ate_ind))
    bias_contr_std <- matrix(0, K-1, length(ate_ind))
    rmse_contr_std <- matrix(0, K-1, length(ate_ind))
  }
  if(need_ipw){
    width_contr_ipw <- matrix(0, K-1, length(ate_ind))
    bias_contr_ipw <- matrix(0, K-1, length(ate_ind))
    rmse_contr_ipw <- matrix(0, K-1, length(ate_ind))
  }
  for(k in 2:K){
    for(i in 1:length(ate_ind)){
      for(iter in 1:n_iter){
        width_contr[k-1, i] <- width_contr[k-1, i] + 
          abs(sim[[iter]]$contr[ate_ind[i], k-1, 3] - 
                sim[[iter]]$contr[ate_ind[i], k-1, 2])
        bias_contr[k-1, i] <- bias_contr[k-1, i] + 
          sim[[iter]]$contr[ate_ind[i], k-1, 1] - contr_true[k-1]
        rmse_contr[k-1, i] <- rmse_contr[k-1, i] + 
          (sim[[iter]]$contr[ate_ind[i], k-1, 1] - contr_true[k-1])^2
        if(need_std){
          width_contr_std[k-1, i] <- width_contr_std[k-1, i] + 
            abs(sim[[iter]]$contr_standard[ate_ind[i], k-1, 3] - 
                  sim[[iter]]$contr_standard[ate_ind[i], k-1, 2])
          bias_contr_std[k-1, i] <- bias_contr_std[k-1, i] + 
            sim[[iter]]$contr_standard[ate_ind[i], k-1, 1] - contr_true[k-1]
          rmse_contr_std[k-1, i] <- rmse_contr_std[k-1, i] + 
            (sim[[iter]]$contr_standard[ate_ind[i], k-1, 1] - contr_true[k-1])^2
        }
        if(need_ipw){
          width_contr_ipw[k-1, i] <- width_contr_ipw[k-1, i] + 
            abs(sim[[iter]]$contr_ipw[ate_ind[i], k-1, 3] - 
                  sim[[iter]]$contr_ipw[ate_ind[i], k-1, 2])
          bias_contr_ipw[k-1, i] <- bias_contr_ipw[k-1, i] + 
            sim[[iter]]$contr_ipw[ate_ind[i], k-1, 1] - contr_true[k-1]
          rmse_contr_ipw[k-1, i] <- rmse_contr_ipw[k-1, i] + 
            (sim[[iter]]$contr_ipw[ate_ind[i], k-1, 1] - contr_true[k-1])^2
        }
      }
      width_contr[k-1, i] <- width_contr[k-1, i] / n_iter
      bias_contr[k-1, i] <- bias_contr[k-1, i] / n_iter
      rmse_contr[k-1, i] <- sqrt(rmse_contr[k-1, i] / n_iter)
      if(need_std){
        width_contr_std[k-1, i] <- width_contr_std[k-1, i] / n_iter
        bias_contr_std[k-1, i] <- bias_contr_std[k-1, i] / n_iter
        rmse_contr_std[k-1, i] <- sqrt(rmse_contr_std[k-1, i] / n_iter)
      }
      if(need_ipw){
        width_contr_ipw[k-1, i] <- width_contr_ipw[k-1, i] / n_iter
        bias_contr_ipw[k-1, i] <- bias_contr_ipw[k-1, i] / n_iter
        rmse_contr_ipw[k-1, i] <- sqrt(rmse_contr_ipw[k-1, i] / n_iter)
      }
    }
  }
  out_list <- list(width = width_contr, bias = bias_contr, 
       rmse = rmse_contr)
  if(need_std){
    out_list[["width_std"]] <- width_contr_std
    out_list[["bias_std"]] <- bias_contr_std
    out_list[["rmse_std"]] <- rmse_contr_std
  }
  if(need_ipw){
    out_list[["width_ipw"]] <- width_contr_ipw
    out_list[["bias_ipw"]] <- bias_contr_ipw
    out_list[["rmse_ipw"]] <- rmse_contr_ipw
  }
  out_list
}

gen_bias_rmse_tab <- function(summ_rand, summ_ts, summ_rits, ate_ind, ind){
  bias_tab <- matrix(NA, 7*(K-1), length(ate_ind))
  for(k in 1:(K-1)){
    bias_tab[7*(k-1)+1, ] <- summ_rand$bias_std[k, ]
    bias_tab[7*(k-1)+2, ] <- summ_rand$bias[k, ]
    bias_tab[7*(k-1)+3, ] <- summ_ts$bias[k, ]
    bias_tab[7*(k-1)+4, ] <- summ_rits$bias[k, ]
    bias_tab[7*(k-1)+5, ] <- summ_rand$bias_ipw[k, ]
    bias_tab[7*(k-1)+6, ] <- summ_ts$bias_ipw[k, ]
    bias_tab[7*(k-1)+7, ] <- summ_rits$bias_ipw[k, ]
  }
  rmse_tab <- matrix(NA, 7*(K-1), length(ate_ind))
  for(k in 1:(K-1)){
    rmse_tab[7*(k-1)+1, ] <- summ_rand$rmse_std[k, ]
    rmse_tab[7*(k-1)+2, ] <- summ_rand$rmse[k, ]
    rmse_tab[7*(k-1)+3, ] <- summ_ts$rmse[k, ]
    rmse_tab[7*(k-1)+4, ] <- summ_rits$rmse[k, ]
    rmse_tab[7*(k-1)+5, ] <- summ_rand$rmse_ipw[k, ]
    rmse_tab[7*(k-1)+6, ] <- summ_ts$rmse_ipw[k, ]
    rmse_tab[7*(k-1)+7, ] <- summ_rits$rmse_ipw[k, ]
  }
  est_err_tab <- matrix(paste(round(bias_tab, 2), "(", round(rmse_tab, 2), ")", sep = ""), 
                        nrow = nrow(bias_tab))
  colnames(est_err_tab) <- paste(ind)
  rownames(est_err_tab) <- paste(c("Rand-OF", "Rand-AIPW", "TS-AIPW", "RiTS-AIPW",
                                   "Rand-IPW", "TS-IPW", "RiTS-IPW"), 
                                 "(Arm", rep(2:4, each = 7), ")", sep = "")
  est_err_tab
}

gen_cum_miscov_for_tab <- function(out_rand, out_ts, out_rits, 
                                   mu_true, contr_true, delay_aipw = 77,
                                   delay_ipw = 0){
  n_iter <- length(out_rand); K <- length(mu_true); N <- length(out_rand[[1]]$trt)
  end_cum_miscov <- numeric((K-1)*7)
  names(end_cum_miscov) <- paste("Arm ", rep(2:K, each = 7), ":", 
                                 c("Rand+T-test", "Rand+AIPW", 
                                   "TS+AIPW", "RiTS+AIPW", "Rand+IPW", 
                                   "TS+IPW", "RiTS+IPW"))
  cum_miscov_rand <- get_cum_mis_cov(out_rand, mu_true = mu_true, 
                                contr_true = contr_true, delay_aipw = delay_aipw,
                                delay_ipw = delay_ipw, need_ipw = TRUE,
                                need_std = TRUE)
  end_ind <- nrow(cum_miscov_rand[[2]])
  cum_miscov_ts <- get_cum_mis_cov(out_ts, mu_true = mu_true,
                                   contr_true = contr_true, delay_aipw = delay_aipw,
                                   delay_ipw = delay_ipw, need_ipw = TRUE,
                                   need_std = FALSE)
  cum_miscov_rits <- get_cum_mis_cov(out_rits, mu_true = mu_true,
                                     contr_true = contr_true, delay_aipw = delay_aipw,
                                     delay_ipw = delay_ipw, need_ipw = TRUE,
                                     need_std = FALSE)
  for(k in 1:(K-1)){
    end_cum_miscov[seq(1, (K-1)*7, 7)] <- cum_miscov_rand[["contr_std"]][end_ind, ]
    end_cum_miscov[seq(2, (K-1)*7, 7)] <- cum_miscov_rand[["contr"]][end_ind, ]
    end_cum_miscov[seq(3, (K-1)*7, 7)] <- cum_miscov_ts[["contr"]][end_ind, ]
    end_cum_miscov[seq(4, (K-1)*7, 7)] <- cum_miscov_rits[["contr"]][end_ind, ]
    end_cum_miscov[seq(5, (K-1)*7, 7)] <- cum_miscov_rand[["contr_ipw"]][end_ind, ]
    end_cum_miscov[seq(6, (K-1)*7, 7)] <- cum_miscov_ts[["contr_ipw"]][end_ind, ]
    end_cum_miscov[seq(7, (K-1)*7, 7)] <- cum_miscov_rits[["contr_ipw"]][end_ind, ]
  }
  end_cum_miscov/n_iter
}
