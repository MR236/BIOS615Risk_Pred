library(RISKmis)
library(survival)
library(doRNG)
library(glmnet)
library(ggplot2)
library(data.table)
library(dplyr)
n <- 500
m <- 500
beta0 <- c(0.7,0.7,0.7,-0.5,-0.5,-0.5,0.3,0.3,0.3,0)
sigmaZ <- 1
misT <- noise.loglinear.S
nsim <- 1000
estimates_new <- NULL
t1 <- Sys.time()
for (i in 1:nsim) {
if (i %% 20 == 0){print(i)}
sim_data <- sim.gen.rev(n, beta0, sigmaZ,m,misT)
labelled_data <- sim_data$data[1:m,]
betadelta <- init.beta.new(delta = labelled_data$delta,C=labelled_data$C,Z = labelled_data$Z)$beta
estimates_new <- rbind(estimates_new, betadelta)
}
t2 <- Sys.time()
time_new <- as.numeric(t2-t1)/nsim
# time_new (m = 250, nsim=1000) = 0.09s
# time_new (m = 500, nsim=1000) = 0.17s
# time_new (m = 2000, nsim=1000) = 0.91s


estimates_old <- NULL
t1 <- Sys.time()
for (i in 1:nsim) {
  if (i %% 2 == 0){print(i)}
  sim_data <- sim.gen.rev(n, beta0, sigmaZ,m,misT)
  labelled_data <- sim_data$data[1:m,]
  h = sd(labelled_data$C)/(sum(labelled_data$delta))^0.26
  KC = dnorm(as.matrix(dist(labelled_data$C/h,diag=T,upper=T)))/h
  betadelta <- init.beta(delta = labelled_data$delta,Z = labelled_data$Z,KC=KC, link=expit, dlink=dexpit)
  estimates_old <- rbind(estimates_old, betadelta)
}
t2 <- Sys.time()
time_old <- as.numeric(t2-t1)/nsim
# time_old (m = 250, nsim=1000)  = 0.2233s
# time_old (m = 500, nsim=1000)  = 1.1589s
# time_old (m = 2000, nsim=1000) = 16.2469s


old_results <- data.frame(estimates_old)
old_results_long <- melt(setDT(old_results), value.name="estimate")
old_results_long$method <- "Old"


new_results <- data.frame(estimates_new)
new_results_long <- melt(setDT(new_results), value.name="estimate")
new_results_long$method <- "New"

full_data <- rbind(old_results_long, new_results_long)
true_vals <- data.frame(variable = c("X1", "X2","X3","X4","X5",
                                     "X6", "X7", "X8","X9","X10"),
                        true_value = beta0)
full_data <- left_join(full_data, true_vals)
full_data <- mutate(full_data, bias = estimate-true_value)

write.csv(full_data, "SimResultsBetaDelta.csv")

full_data$variable <- factor(full_data$variable, levels = c("X1", "X2","X3","X4","X5",
                                                            "X6", "X7", "X8","X9","X10"))
ggplot(aes(x = variable, fill = method, y = bias), data=full_data) + 
  geom_boxplot() + 
  theme_bw() + 
  xlab("Coefficient") + ylab("Bias") + labs(fill="Method")

ggsave("InitialBiasPlot.jpg", height=4, width=7, dpi=400)

stored_data <- read.csv("SimResultsBetaDelta.csv")
MSE <- stored_data %>%
  mutate(sqr_error = bias^2) %>%
  group_by(variable, method) %>%
  summarise(Bias=mean(bias), MSE = mean(sqr_error)) %>% 
  ungroup()

Overall_MSE <- MSE %>% 
  group_by(method) %>%
  summarize(mean_MSE = mean(MSE), mean_bias = mean(Bias))
