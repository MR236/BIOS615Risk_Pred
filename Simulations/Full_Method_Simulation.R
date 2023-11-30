library(RISKmis)
library(survival)
library(doRNG)
library(glmnet)
library(ggplot2)
library(data.table)
library(dplyr)
library(MASS)
New_Method_Estimates <- NULL
Old_Method_Estimates <- NULL
New_Times <- NULL
Old_Times <- NULL
for (j in c(1)) {
n <- 500
m <- 500
beta0 <- c(0.7,0.7,0.7,-0.5,-0.5,-0.5,0.3,0.3,0.3,0)
sigmaZ <- 1
misT <- noise.loglinear.S
estimates_new <- NULL
sim_data <- sim.gen.rev(n, beta0, sigmaZ,m,misT)

data_labelled <- sim_data$data
data_unlabelled <- sim_data$dataS

t1 <- Sys.time()
result <- Full_Estimator_New(data_labelled, data_unlabelled, n=n, m=m, Nperturb=50)
t2 <- Sys.time()
time_new <- as.numeric(t2-t1, units="secs")
lowerCI <- NULL
upperCI <- NULL
for (i in c(1:10)) {
  lowerCI <- c(lowerCI, quantile(result[[2]][i,], 0.025))
  upperCI <- c(upperCI, quantile(result[[2]][i,], 0.975))
}
new_data <- data.frame(Est = result[[1]], lower=lowerCI, upper=upperCI, Coefficient = c(1:10), Sim=j)

t1 <- Sys.time()
result_old <- Full_Estimator_Original(data_labelled, data_unlabelled, n=n,m=m, Nperturb=50)
t2 <- Sys.time()
time_old <- as.numeric(t2-t1, units="secs")


lowerCI_old <- NULL
upperCI_old <- NULL
for (i in c(1:10)) {
  lowerCI_old <- c(lowerCI_old, quantile(result_old[[2]][i,], 0.025))
  upperCI_old <- c(upperCI_old, quantile(result_old[[2]][i,], 0.975))
}
old_data <- data.frame(Est = result_old[[1]], lower=lowerCI_old, upper=upperCI_old, Coefficient = c(1:10), Sim=j)
New_Times <- c(New_Times, time_new)
Old_Times <- c(Old_Times, time_old)
New_Method_Estimates <- rbind(New_Method_Estimates, new_data)
Old_Method_Estimates <- rbind(Old_Method_Estimates, old_data)
}
write.csv(New_Method_Estimates, "New_Method_Estimates1.csv")
write.csv(Old_Method_Estimates, "Old_Method_Estimates1.csv")
Time_Data <- rbind(data.frame(Type = "New", Time = New_Times), data.frame(Type = "Old", Time = Old_Times))
