source("Code/Author_Original_Functions.R")
source("Code/Author_Utility_Functions.R")
source("Code/Project_New_Functions.R")
source("Code/Full_Implementations.R")

library(survival)
library(doRNG)
library(glmnet)
library(ggplot2)
library(data.table)
library(dplyr)
library(MASS)
repeat_sims <- FALSE
if (repeat_sims) {
  New_Method_Estimates <- read.csv("Results/New_Method_Estimates1.csv")
  Old_Method_Estimates <- read.csv("Results/Old_Method_Estimates1.csv")
  Times <- read.csv("Results/Time_Data1.csv")
  New_Times <- filter(Times, Type == "New")$Time
  Old_Times <- filter(Times, Type == "Old")$Time
} else {New_Method_Estimates <- NULL
Old_Method_Estimates <- NULL
New_Times <- NULL
Old_Times <- NULL}
for (j in c(1:100)) {
  n <- 1000
  m <- 500
  beta0 <- c(0.7,0.7,0.7,-0.5,-0.5,-0.5,0.3,0.3,0.3,0)
  sigmaZ <- 1
  misT <- noise.loglinear.S
  estimates_new <- NULL
  sim_data <- sim.gen.rev(n, beta0, sigmaZ,m,misT)
  
  data_labelled <- sim_data$data
  data_unlabelled <- sim_data$dataS
  
  t1 <- Sys.time()
  result <- Full_Estimator_New(data_labelled, data_unlabelled, n=n, m=m, Nperturb=500)
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
  result_old <- Full_Estimator_Original(data_labelled, data_unlabelled, n=n,m=m, Nperturb=500)
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
write.csv(New_Method_Estimates, "Results/New_Method_Estimates2.csv", row.names=F)
write.csv(Old_Method_Estimates, "Results/Old_Method_Estimates2.csv",row.names=F)
Time_Data <- rbind(data.frame(Type = "New", Time = New_Times), data.frame(Type = "Old", Time = Old_Times))
write.csv(Time_Data, "Results/Time_Data2.csv",row.names=F)