library(RISKmis)
library(survival)
library(doRNG)
library(glmnet)
library(ggplot2)
library(data.table)
library(dplyr)
n <- 1000
m <- 250
beta0 <- c(0.7,0.7,0.7,-0.5,-0.5,-0.5,0.3,0.3,0.3,0)
p <- length(beta0)
sigmaZ <- 1
misT <- noise.loglinear.S
nsim <- 50
#------------------------------------------
#    initial estimator beta_delta
#------------------------------------------
sim_data <- sim.gen.rev(n, beta0, sigmaZ,m,misT)
labelled_data <- sim_data$data[1:m,]
betadelta <- init.beta.new(delta = labelled_data$delta,C=labelled_data$C,Z = labelled_data$Z)$beta


#------------------------------------------
#   score for beta_delta
#------------------------------------------
data <- sim_data$data
dataS <- sim_data$dataS
beta.std = betadelta/sqrt(sum(betadelta^2))
lp = drop(data$Z %*% beta.std)
sdlp = sd(lp)

h1 = sdlp/(sum(dataS$delta1))^0.3
h2 = sdlp/(sum(dataS$delta2))^0.3
t1 <- Sys.time()
for (i in 1:nsim) {
  if (i %% 20 == 0){print(i)}
  Sk = rep(0,p*2)
Sk[1:p] = Sk_sym_new(lp,data$Z,
                         dataS$X1,dataS$delta1,
                         data$C,h1)
Sk[p+1:p] = Sk_sym_new(lp,data$Z,
                           dataS$X2,dataS$delta2,
                           data$C,h2)}
t2 <- Sys.time()
time_new <- as.numeric(t2-t1) /nsim
# N = 10000, 1.83s


t1 <- Sys.time()
for (i in 1:nsim) {
  if (i %% 2 == 0){print(i)}
  Sk = rep(0,p*2)
  Sk[1:p] = Sk_sym(lp,data$Z,
                       dataS$X1,dataS$delta1,
                       data$C,dnorm,h1)
  Sk[p+1:p] = Sk_sym(lp,data$Z,
                         dataS$X2,dataS$delta2,
                         data$C,dnorm,h2)}
t2 <- Sys.time()
time_old <- as.numeric(t2-t1)/nsim
# N = 10000, 7.14s
