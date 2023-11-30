library(RISKmis)
library(survival)
library(doRNG)
library(glmnet)
library(ggplot2)
library(data.table)
library(dplyr)
n <- 20000
m <- 500
beta0 <- c(0.7,0.7,0.7,-0.5,-0.5,-0.5,0.3,0.3,0.3,0)
p <- length(beta0)
sigmaZ <- 1
misT <- noise.loglinear.S
nsim <- 40
errors <- NULL
#------------------------------------------
#    initial estimator beta_delta
#------------------------------------------
for (i in 1:nsim) {
print(i)
sim_data <- sim.gen.rev(n, beta0, sigmaZ,m,misT)
labelled_data <- sim_data$data[1:m,]
betadelta <- init.beta.new(delta = labelled_data$delta,C=labelled_data$C,Z = labelled_data$Z)


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

Sk_new = rep(0,p*2)
Sk_new[1:p] = Sk_sym_new(lp,data$Z,
                       dataS$X1,dataS$delta1,
                       data$C,h1)
Sk_new[p+1:p] = Sk_sym_new(lp,data$Z,
                         dataS$X2,dataS$delta2,
                         data$C,h2)
Sk = rep(0,p*2)
Sk[1:p] = Sk_sym(lp,data$Z,
                   dataS$X1,dataS$delta1,
                   data$C,dnorm,h1)
Sk[p+1:p] = Sk_sym(lp,data$Z,
                     dataS$X2,dataS$delta2,
                     data$C,dnorm,h2)
errors <- c(errors, mean(abs(Sk_new/Sk - 1)))}
