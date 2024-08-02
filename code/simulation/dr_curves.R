library(ggplot2)
library(patchwork)
library(reshape2)

mu1 <- function(x){
  5*(1 - 0.005 * (x+0.5)^2 - 0.005 * (x-0.5)^2)
}
mu2 <- function(x){
  5*(1.35 - 0.1 * (x+0.5)^2 - 0.1 * (x-0.5)^2)
}
mu3 <- function(x){
  5*(1.35 - 0.05 * (x+0.5)^2 - 0.05 * (x-0.5)^2)
}
mu4 <- function(x){
  5*(1.6 - 0.1 * (x+0.5)^2 - 0.1 * (x-0.5)^2)
}

x_vals <- seq(-3, 3, by = 0.01)

df_mu <- data.frame(x = x_vals, mu1 = mu1(x_vals), mu2 = mu2(x_vals),
                    mu3 = mu3(x_vals), mu4 = mu4(x_vals))


eff <- ggplot(df_mu, aes(x = x)) +
  geom_line(aes(y = mu1, color = "Arm 1"), linetype = "solid") +  
  geom_line(aes(y = mu2, color = "Arm 2"), linetype = "solid") +  
  geom_line(aes(y = mu3, color = "Arm 3"), linetype = "solid") +
  geom_line(aes(y = mu4, color = "Arm 4"), linetype = "solid") +
  labs(x = "z", y = "Average Efficacy") +
  theme(text = element_text(size = 8)) +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#CC79A7"))

nu1 <- function(x){
  5*1
}
nu2 <- function(x){
  5*(1 - 0.05 * x^2)
}
nu3 <- function(x){
  5*(1 - 0.005 * x^2)
}
nu4 <- function(x){
  5*(1 - 0.3 * x^2)
}

df_nu <- data.frame(x = x_vals, nu1 = nu1(x_vals), nu2 = nu2(x_vals),
                    nu3 = nu3(x_vals), nu4 = nu4(x_vals))

safe <- ggplot(df_nu, aes(x = x)) +
  geom_line(aes(y = nu1, color = "Arm 1"), linetype = "solid") +  
  geom_line(aes(y = nu2, color = "Arm 2"), linetype = "solid") +  
  geom_line(aes(y = nu3, color = "Arm 3"), linetype = "solid") +
  geom_line(aes(y = nu4, color = "Arm 4"), linetype = "solid") +
  labs(x = "z", y = "Average Safety") +
  theme(text = element_text(size = 8)) +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#CC79A7"))

w <- 1
df_util <- data.frame(x = x_vals, util1 = (df_mu$mu1 + df_nu$nu1 * w)/(w+1),
                      util2 = (df_mu$mu2 + df_nu$nu2 * w)/(w+1),
                      util3 = (df_mu$mu3 + df_nu$nu3 * w)/(w+1),
                      util4 = (df_mu$mu4 + df_nu$nu4 * w)/(w+1))

util <- ggplot(df_util, aes(x = x)) +
  geom_line(aes(y = util1, color = "Arm 1"), linetype = "solid") +  
  geom_line(aes(y = util2, color = "Arm 2"), linetype = "solid") +  
  geom_line(aes(y = util3, color = "Arm 3"), linetype = "solid") +
  geom_line(aes(y = util4, color = "Arm 4"), linetype = "solid") +
  labs(x = "z", y = "Average Utility") +
  theme(text = element_text(size = 8)) +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#CC79A7"))

dr_curves <- eff + safe + util + plot_layout(ncol = 3, guides = "collect")
ggsave("plot/dr_curves.jpg", height = 2, width = 8, units = "in")


df_avg <- data.frame(Arm = 1:4, Efficacy = 5*c(0.9875, 1.1, 1.225, 1.35),
                     Safety = 5*c(1, 0.995, 0.95, 0.7),
                     Utility = 5*c(0.99375, 1.0475, 1.0875, 1.025))
ate <- ggplot(df_avg, aes(x = Arm)) +
  geom_line(aes(y = Efficacy, color = "Efficacy"), linetype = "solid") +  
  geom_line(aes(y = Safety, color = "Safety"), linetype = "solid") +  
  geom_line(aes(y = Utility, color = "Utility"), linetype = "solid") +
  labs(x = "Arm", y = "Average Effect") +
  theme(text = element_text(size = 8)) +
  scale_color_manual(
    values = c("Utility" = "#000000", "Efficacy" = "#D55E00", 
               "Safety" = "#009E73")
  )
ggsave("plot/dr_avg.jpg", height = 2, width = 3, units = "in")
