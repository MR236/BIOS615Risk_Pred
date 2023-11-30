#' Author's Original Implementation of Full Method for Two Surrogate Outcomes
#' @param data labelled data with columns delta (event indicator), C (censoring indicator), Z (covariates)
#' @param dataS unlabelled data containing event time and censoring indicator for surrogate outcomes
#' @param link link function for score equation to be solved
#' @param dlink derivative of link function
#' @param inv.link inverse of link function
#' @param Nperturb number of perturbations in perturbation analysis
#' @param min.lam minimum value of lambda for adaptive ridge regression
#' @param maxsoft maximum value of lambda for adaptive ridge regression
#' @param Nfold number of fold in cross-validation analysis
#' @param m size of labelled data subset
#' @param n size of full data
#' @param p number of covariates
#' @param Nlam number of lambdas to test in adaptive ridge regression
#' @return A list containing an estimate betaSSL, and set of perturbed estimates betaSSLk for inference
#' @export
Full_Estimator_Original <- function(data, 
                                    dataS,
                                    link = expit, 
                                    dlink=dexpit, 
                                    inv.link=logit, 
                                    Nperturb=1000,
                                    min.lam = 1e-4,
                                    maxsoft = 0.5,
                                    Nfold = 5,
                                    m = nrow(data),
                                    n = nrow(dataS),
                                    p = ncol(data$Z),
                                    Nlam = 100) {
  foldid = c(rep_len(1:Nfold,m)[sample(m)],
             rep_len(Nfold:1,n-m)[sample(n-m)])
  
  #------------------------------------------
  #    beta_delta and its perturbation
  #------------------------------------------
  
  h = sd(data$C[1:m])/(sum(data$delta[1:m]))^0.26
  KC = dnorm(as.matrix(dist(data$C[1:m]/h,diag=T,upper=T)))/h
  
  betadelta = init.beta(data$delta[1:m],data$Z[1:m,],KC,
                        link = link, dlink = dlink)
  
  init.cv = matrix(0,p,Nfold)
  lp.init.cv = matrix(0,n,Nfold)
  for(ifold in 1:Nfold)
  {
    i.label = which(foldid[1:m] != ifold)
    i.train = which(foldid != ifold)
    init.cv[,ifold] = init.beta(data$delta[i.label],data$Z[i.label,],KC[i.label,i.label],
                                link = link, dlink = dlink)
    lp.init.cv[,ifold] = drop(data$Z %*% init.cv[,ifold])
  }
  
  
  
  
  #------------------------------------------
  #    Sk
  #------------------------------------------
  
  beta.std = betadelta/sqrt(sum(betadelta^2))
  lp = drop(data$Z %*% beta.std)
  sdlp = sd(lp)
  
  # for (isetup in 1:Nsetup)
  # {
  
  h1 = sdlp/(sum(dataS$delta1))^0.3
  h2 = sdlp/(sum(dataS$delta2))^0.3
  Sk = rep(0,p*2)
  Sk[1:p] = Sk_sym(lp,data$Z,
                   dataS$X1,dataS$delta1,
                   data$C,dnorm,h1)
  Sk[p+1:p] = Sk_sym(lp,data$Z,
                     dataS$X2,dataS$delta2,
                     data$C,dnorm,h2)
  
  #------------------------------------------
  #    Sk perturbation
  #------------------------------------------
  
  
  Skb =  matrix(0, Nperturb,p*2)
  
  h1 = sdlp/(sum(dataS$delta1))^0.3
  h2 = sdlp/(sum(dataS$delta2))^0.3
  
  betak =  matrix(0,p, Nperturb)
  for (iperturb in 1:Nperturb) 
  {
    V = rbeta(n,0.5,1.5)*4
    betak[,iperturb] = 
      init.beta.perturb(data$delta[1:m],data$Z[1:m,],KC,
                        V[1:m], 
                        init = betadelta,
                        link = link, dlink = dlink)
    betak.std = betak[,iperturb]/sqrt(sum(betak[,iperturb]^2))
    lp = data$Z %*% betak.std
    Skb[iperturb,1:p] = c(Sk_sym_perturb(lp,data$Z,
                                         dataS$X1,dataS$delta1,
                                         data$C,dnorm,h1,
                                         matrix(V)))
    Skb[iperturb,p+1:p] = c(Sk_sym_perturb(lp,data$Z,
                                           dataS$X2,dataS$delta2,
                                           data$C,dnorm,h2,
                                           matrix(V)))
  }
  
  #------------------------------------------
  # CV for adaptive soft-thresholding initial
  #------------------------------------------ 
  ada.factor = 1/(abs(betadelta)*apply(data$Z, 2, sd))
  lambda = c(0,exp(seq(log(min.lam),log(max(betadelta^2)),length.out = Nlam))[-Nlam])
  
  #beta.cv = matrix(0,p,Nfold)
  cv.concord = matrix(0,Nfold,Nlam)
  for(ifold in 1:Nfold)
  {
    i.train = which(foldid != ifold)
    i.test = which(foldid == ifold)
    beta.cv = init.cv[,ifold] 
    beta.thres.cv = sign(beta.cv)*matrix(pmax(0,abs(beta.cv)-outer(ada.factor,lambda,"*")),p)
    lp.test = data$Z[i.test,] %*% beta.thres.cv
    tmp.concord = (cv.concordance(lp.test, 
                                  dataS$X1[i.test],dataS$delta1[i.test],
                                  data$C[i.test]) + 
                     cv.concordance(lp.test, 
                                    dataS$X2[i.test],dataS$delta2[i.test],
                                    data$C[i.test]))
    cv.concord[ifold,] = tmp.concord[,1]/tmp.concord[,2]
  }
  
  best.ilam = which.max(apply(cv.concord,2,mean))
  lambda.best = lambda[best.ilam]
  thres.best = lambda.best * ada.factor
  betadelta.thres = sign(betadelta) * pmax(0,abs(betadelta) - thres.best)
  
  #------------------------------------------
  #    Estimating W matrix
  #------------------------------------------
  
  W.hat.adaPCA.ridge = W_hat_adaPCA_ridge(betak, Skb)
  
  #------------------------------------------
  # CV for adaptive soft-thresholding SSL1
  #------------------------------------------
  
  W.hat = W.hat.adaPCA.ridge$W.hat
  
  betaSSL = betadelta - drop(W.hat %*% Sk)
  ada.factor = 1/(abs(betaSSL)*apply(data$Z, 2, sd))
  lambda = c(0,exp(seq(log(min.lam),log(min(max(betaSSL^2),maxsoft)),length.out = Nlam))[-Nlam])
  
  #beta.cv = matrix(0,p,Nfold)
  cv.concord = matrix(0,Nfold,Nlam)
  for(ifold in 1:Nfold)
  {
    i.train = which(foldid != ifold)
    i.test = which(foldid == ifold)
    Sk.cv = rep(0,p*2)
    Sk.cv[1:p] = Sk_sym(lp.init.cv[i.train,ifold],data$Z[i.train,],
                        dataS$X1[i.train],dataS$delta1[i.train],
                        data$C[i.train],dnorm,h1)
    Sk.cv[p+1:p] = Sk_sym(lp.init.cv[i.train,ifold],data$Z[i.train,],
                          dataS$X2[i.train],dataS$delta2[i.train],
                          data$C[i.train],dnorm,h2)
    beta.cv = init.cv[,ifold] - drop (W.hat %*% Sk.cv)
    beta.thres.cv = sign(beta.cv)*matrix(pmax(0,abs(beta.cv)-outer(ada.factor,lambda,"*")),p)
    lp.test = data$Z[i.test,] %*% beta.thres.cv
    tmp.concord = (cv.concordance(lp.test, 
                                  dataS$X1[i.test],dataS$delta1[i.test],
                                  data$C[i.test]) + 
                     cv.concordance(lp.test, 
                                    dataS$X2[i.test],dataS$delta2[i.test],
                                    data$C[i.test]))
    cv.concord[ifold,] = tmp.concord[,1]/tmp.concord[,2]
  }
  
  best.ilam = which.max(apply(cv.concord,2,mean))
  lambda.best = lambda[best.ilam]
  thres.best = lambda.best * ada.factor
  
  betaSSLk = betak - W.hat %*% t(Skb)
  
  
  return(list(betaSSL, betaSSLk))
}

#' New Implementation of Full Method for Two Surrogate Outcomes
#' @param data labelled data with columns delta (event indicator), C (censoring indicator), Z (covariates)
#' @param dataS unlabelled data containing event time and censoring indicator for surrogate outcomes
#' @param link link function for score equation to be solved
#' @param dlink derivative of link function
#' @param inv.link inverse of link function
#' @param Nperturb number of perturbations in perturbation analysis
#' @param min.lam minimum value of lambda for adaptive ridge regression
#' @param maxsoft maximum value of lambda for adaptive ridge regression
#' @param Nfold number of fold in cross-validation analysis
#' @param m size of labelled data subset
#' @param n size of full data
#' @param p number of covariates
#' @param Nlam number of lambdas to test in adaptive ridge regression
#' @param knots number of knots for spline in first-step estimator
#' @return A list containing an estimate betaSSL, and set of perturbed estimates betaSSLk for inference
#' @export
Full_Estimator_New <- function(data, 
                               dataS,
                               link = expit, 
                               dlink=dexpit, 
                               inv.link=logit, 
                               Nperturb=1000,
                               min.lam = 1e-4,
                               maxsoft = 0.5,
                               Nfold = 5,
                               knots = min(9,floor(m/500)),
                               m = nrow(data),
                               n = nrow(dataS),
                               p = ncol(data$Z),
                               Nlam=100) {
  foldid = c(rep_len(1:Nfold,m)[sample(m)],
             rep_len(Nfold:1,n-m)[sample(n-m)])
  
  #------------------------------------------
  #    beta_delta and its perturbation
  #------------------------------------------
  beta_res = init.beta.new(delta = data$delta[1:m],C=data$C[1:m],Z = data$Z[1:m,],knots=knots)
  betadelta = beta_res$beta
  
  init.cv = matrix(0,p,Nfold)
  lp.init.cv = matrix(0,n,Nfold)
  init.cv_m = matrix(0,p,Nfold)
  lp.init.cv_m = matrix(0,n,Nfold)
  for(ifold in 1:Nfold)
  {
    i.label = which(foldid[1:m] != ifold)
    i.train = which(foldid != ifold)
    init.cv[,ifold] = init.beta.new(data$delta[i.label],C=data$C[i.label],Z = data$Z[i.label,],knots=knots)$beta
    lp.init.cv[,ifold] = drop(data$Z %*% init.cv[,ifold])
    
  }
  
  #------------------------------------------
  #    Sk
  #------------------------------------------
  
  beta.std = betadelta/sqrt(sum(betadelta^2))
  lp = drop(data$Z %*% beta.std)
  sdlp = sd(lp)
  
  h1 = sdlp/(sum(dataS$delta1))^0.3
  h2 = sdlp/(sum(dataS$delta2))^0.3
  Sk = rep(0,p*2)
  Sk[1:p] = Sk_sym_new(lp,data$Z,
                       dataS$X1,dataS$delta1,
                       data$C,h1)
  Sk[p+1:p] = Sk_sym_new(lp,data$Z,
                         dataS$X2,dataS$delta2,
                         data$C,h2)
  
  
  #------------------------------------------
  #    Sk perturbation
  #------------------------------------------
  Skb =  matrix(0, Nperturb,p*2)
  
  betak =  matrix(0,p, Nperturb)
  for (iperturb in 1:Nperturb) 
  {
    # replaced original perturbation analysis with bootstrap
    i_boot = sample(c(1:m), m, replace=TRUE)
    lambdas = c(0.01,0.1,0.5,1,2,10,100,200,500,1000)*beta_res$lambda
    betak[,iperturb] = init.beta.new(delta = data$delta[i_boot],C=data$C[i_boot],Z = data$Z[i_boot,],knots=knots, lambda=lambdas)$beta
    
    betak.std = betak[,iperturb]/sqrt(sum(betak[,iperturb]^2))
    
    lp = data$Z %*% betak.std
    
    i_boot_full = c(i_boot, sample(c((m+1):n), n-m, replace=TRUE))
    lp_boot = lp[i_boot_full]
    Z_boot = data$Z[i_boot_full,]
    X1_boot = dataS$X1[i_boot_full]
    delta1_boot = dataS$delta1[i_boot_full]
    C_boot = data$C[i_boot_full]
    Skb[iperturb,1:p] = c(Sk_sym_new(lp_boot,Z_boot,
                                     X1_boot,delta1_boot,
                                     C_boot,h1))
    X2_boot = dataS$X2[i_boot_full]
    delta2_boot = dataS$delta2[i_boot_full]
    Skb[iperturb,p+1:p] = c(Sk_sym_new(lp_boot,Z_boot,
                                       X2_boot,delta2_boot,
                                       C_boot,h2))
  }
  #------------------------------------------
  # CV for adaptive soft-thresholding initial
  #------------------------------------------
  ada.factor = 1/(abs(betadelta)*apply(data$Z, 2, sd))
  lambda = c(0,exp(seq(log(min.lam),log(max(betadelta^2)),length.out = Nlam))[-Nlam])
  
  beta.cv = matrix(0,p,Nfold)
  cv.concord = matrix(0,Nfold,Nlam)
  for(ifold in 1:Nfold)
  {
    i.train = which(foldid != ifold)
    i.test = which(foldid == ifold)
    beta.cv = init.cv[,ifold]
    beta.thres.cv = sign(beta.cv)*matrix(pmax(0,abs(beta.cv)-outer(ada.factor,lambda,"*")),p)
    lp.test = data$Z[i.test,] %*% beta.thres.cv
    tmp.concord = (cv.concordance(lp.test,
                                  dataS$X1[i.test],dataS$delta1[i.test],
                                  data$C[i.test]) +
                     cv.concordance(lp.test,
                                    dataS$X2[i.test],dataS$delta2[i.test],
                                    data$C[i.test]))
    cv.concord[ifold,] = tmp.concord[,1]/tmp.concord[,2]
  }
  
  best.ilam = which.max(apply(cv.concord,2,mean))
  lambda.best = lambda[best.ilam]
  thres.best = lambda.best * ada.factor
  betadelta.thres = sign(betadelta) * pmax(0,abs(betadelta) - thres.best)
  
  #------------------------------------------
  #    Estimating W matrix
  #------------------------------------------
  
  W.hat.adaPCA.ridge = W_hat_adaPCA_ridge(betak, Skb)
  
  #------------------------------------------
  # CV for adaptive soft-thresholding SSL1
  #------------------------------------------
  
  W.hat = W.hat.adaPCA.ridge$W.hat
  
  betaSSL = betadelta - drop(W.hat %*% Sk)
  ada.factor = 1/(abs(betaSSL)*apply(data$Z, 2, sd))
  lambda = c(0,exp(seq(log(min.lam),log(min(max(betaSSL^2),maxsoft)),length.out = Nlam))[-Nlam])
  
  #beta.cv = matrix(0,p,Nfold)
  cv.concord = matrix(0,Nfold,Nlam)
  for(ifold in 1:Nfold)
  {
    i.train = which(foldid != ifold)
    i.test = which(foldid == ifold)
    Sk.cv = rep(0,p*2)
    Sk.cv[1:p] = Sk_sym_new(lp.init.cv[i.train,ifold],data$Z[i.train,],
                            dataS$X1[i.train],dataS$delta1[i.train],
                            data$C[i.train],h1)
    Sk.cv[p+1:p] = Sk_sym_new(lp.init.cv[i.train,ifold],data$Z[i.train,],
                              dataS$X2[i.train],dataS$delta2[i.train],
                              data$C[i.train],h2)
    beta.cv = init.cv[,ifold] - drop (W.hat %*% Sk.cv)
    beta.thres.cv = sign(beta.cv)*matrix(pmax(0,abs(beta.cv)-outer(ada.factor,lambda,"*")),p)
    lp.test = data$Z[i.test,] %*% beta.thres.cv
    tmp.concord = (cv.concordance(lp.test,
                                  dataS$X1[i.test],dataS$delta1[i.test],
                                  data$C[i.test]) +
                     cv.concordance(lp.test,
                                    dataS$X2[i.test],dataS$delta2[i.test],
                                    data$C[i.test]))
    cv.concord[ifold,] = tmp.concord[,1]/tmp.concord[,2]
  }
  
  best.ilam = which.max(apply(cv.concord,2,mean))
  lambda.best = lambda[best.ilam]
  thres.best = lambda.best * ada.factor
  
  betaSSLk = betak - W.hat %*% t(Skb)
  return(list(betaSSL, betaSSLk))
}