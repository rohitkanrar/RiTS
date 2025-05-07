require(ggplot2); require(reshape2)
n_comp <- 151
power_std <- numeric(n_comp)
power_rand <- numeric(n_comp)
power_ts <- numeric(n_comp)
power_rits <- numeric(n_comp)
for(i in 1:n_iter){
  power_std <- power_std + as.numeric(
    sapply(1:n_comp, function(t) which.max(rand_sim_high[[i]]$contr_standard[t, , 1]) == 3)
  )
  power_rand <- power_rand + as.numeric(
    sapply(1:n_comp, function(t) which.max(rand_sim_high[[i]]$contr[t, , 1]) == 3)
  )
  power_ts <- power_ts + as.numeric(
    sapply(1:n_comp, function(t) which.max(ts_sim_high[[i]]$contr[t, , 1]) == 3)
  )
  power_rits <- power_rits + as.numeric(
    sapply(1:n_comp, function(t) which.max(rits_sim_high[[i]]$contr[t, , 1]) == 3)
  )
}
power_std <- power_std / n_iter
power_rand <- power_rand / n_iter
power_ts <- power_ts / n_iter
power_rits <- power_rits / n_iter

df <- data.frame(std = power_std, rand = power_rand, ts = power_ts,
                 rits = power_rits, participant = 50:200)
df_long <- tidyr::gather(df, key = "Methods", value = "value", -participant)
ggplot(df_long, aes(x = participant, y = value, color = Methods)) +
  geom_line() + labs(x = "Participant", y = "Proportion") +
  scale_color_manual(
    values = c("rand" = "#CC79A7", "ts" = "#0072B2", "rits" = "#D55E00", 
               "std" = "#009E73"),
    labels = c("ts" = "TS", "rand" = "Rand", "rits" = "RiTS", "std" = "Std")
  )
ggsave("plot/prop_select.jpg", height = 4, width = 6, units = "in")


zero_in_intv <- function(intv, zero = 0.1){
  a <- intv[1]; b <- intv[2]
  if(a < zero && b > zero){
    return(TRUE)
  } else{
    return(FALSE)
  }
}

power_std <- numeric(n_comp)
power_rand <- numeric(n_comp)
power_ts <- numeric(n_comp)
power_rits <- numeric(n_comp)

for(i in 1:n_iter){
  power_std <- power_std + sapply(1:n_comp, function(t){
    rej <- sapply(1:3, function(k){
      zero_in_intv(rand_sim_high[[i]]$contr_standard[t, k, 2:3])
    })
    as.numeric(FALSE %in% rej)
  })
  power_rand <- power_rand + sapply(1:n_comp, function(t){
    rej <- sapply(1:3, function(k){
      zero_in_intv(rand_sim_high[[i]]$contr[t, k, 2:3])
    })
    as.numeric(FALSE %in% rej)
  })
  power_ts <- power_ts + sapply(1:n_comp, function(t){
    rej <- sapply(1:3, function(k){
      zero_in_intv(ts_sim_high[[i]]$contr[t, k, 2:3])
    })
    as.numeric(FALSE %in% rej)
  })
  power_rits <- power_rits + sapply(1:n_comp, function(t){
    rej <- sapply(1:3, function(k){
      zero_in_intv(rits_sim_high[[i]]$contr[t, k, 2:3])
    })
    as.numeric(FALSE %in% rej)
  })
}

power_std <- power_std / n_iter
power_rand <- power_rand / n_iter
power_ts <- power_ts / n_iter
power_rits <- power_rits / n_iter

df <- data.frame(std = power_std, rand = power_rand, ts = power_ts,
                 rits = power_rits, participant = 50:200)
df_long <- tidyr::gather(df, key = "Methods", value = "value", -participant)
ggplot(df_long, aes(x = participant, y = value, color = Methods)) +
  geom_line() + labs(x = "Participant", y = "Power") +
  scale_color_manual(
    values = c("rand" = "#CC79A7", "ts" = "#0072B2", "rits" = "#D55E00", 
               "std" = "#009E73"),
    labels = c("ts" = "TS", "rand" = "Rand", "rits" = "RiTS", "std" = "Std")
  )
ggsave("plot/power.jpg", height = 4, width = 6, units = "in")
