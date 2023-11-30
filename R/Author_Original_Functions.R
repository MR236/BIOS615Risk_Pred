#' Author's Original Implementation of First Step Estimator (only labelled data)
#' @param delta event indicator
#' @param Z covariates (design matrix)
#' @param KC kernel matrix evaluated based on censoring times
#' @param init initial values for coefficients beta, 0 by default
#' @param tol convergence tolerance for beta
#' @param maxit maximum number of iterations
#' @param min.factor minimum multiplicative decrease at a given iteration to not enter line search
#' @param ls.factor decrease in step size when conducting line search
#' @param max.move maximum step size within line search as proportion of newton step
#' @param link link function for score equation to be solved
#' @param dlink derivative of link function
#' @return A list of coefficients beta corresponding to each covariate column in Z
#' @export
init.beta = function(delta,Z, KC, init = rep(0,ncol(Z)), tol=1e-7,
                     maxit = 100, min.factor = 0.75,
                     ls.factor = 0.75, max.move = 1,
                     link = expit, dlink = dexpit)
{
  n = nrow(Z)
  KCd = drop(KC%*%delta)
  hC = rep(0,n)
  oldscore = NULL
  
  for(k in 1:maxit) 
  {
    lp = drop(Z%*%init)
    hC.flag = rep(TRUE,n)
    gij = wZbar = matrix(0,n,n)
    hHess = rep(0,n)
    for(kk in 1:maxit)
    {
      # print(hC[hC.flag])
      lij = outer(hC[hC.flag],lp,"+")
      gij[hC.flag,] = link(lij)
      tmp = KC[hC.flag,]*gij[hC.flag,]
      wZbar[hC.flag,] = KC[hC.flag,]*dlink(lij)
      if(sum(hC.flag)>=2)
      {
        hscore = apply(tmp,1,sum)-KCd[hC.flag]
        hHess[hC.flag] = apply(wZbar[hC.flag,],1,sum)
      }else
      {
        hscore = sum(tmp)-KCd[hC.flag]
        hHess[hC.flag] = sum(wZbar[hC.flag,])
      }
      
      dhC = hscore/hHess[hC.flag]
      dhC = sign(dhC)*pmin(abs(dhC),max.move)
      kk.flag = abs(hscore) > tol
      if(!any(kk.flag))
        break
      hC[hC.flag][kk.flag] = hC[hC.flag][kk.flag] - dhC[kk.flag]
      hC.flag[hC.flag] = kk.flag
    }
    if(kk >= maxit)
      stop("Numerical error when computing h0(Ci)")
    Zbar =  (wZbar%*%Z) / hHess 
    
    gi = link(hC+lp)
    bscore = drop(t(Z)%*% (delta - gi))
    if(!is.null(oldscore))
      if(((sum(oldscore^2)*min.factor) <= sum(bscore^2)))
      {
        init = init+dinit
        dinit = dinit*ls.factor
        if(max(abs(dinit))<tol)
        {
          if(max(abs(oldscore)) > 1e-6)
            warning(paste("Algorithm stops in line-search. Target tol: ",
                          tol, ". Current tol: ", max(abs(oldscore)),
                          ". ", sep = ''))
          break
        }
        init = init - dinit
        next
      }
    oldscore = bscore
    bHess = t(dlink(hC+lp)*Z) %*% (Zbar-Z)
    dinit = solve(bHess,bscore)
    if(all(abs(bscore)<tol))
      break
    # print(rbind(init,bscore,dinit))
    init = init - dinit
  }
  if(k >=maxit)
    stop("Numerical error when computing beta_delta")
  
  return(init)
}

#' Author's Original Implementation of Score Function for Second Step Estimator
#' Corresponds to a single surrogate outcome (run one time per surrogate)
#' @param lp linear predictor for each patient (X transpose beta)
#' @param Z covariates (design matrix)
#' @param Xk surrogate event times
#' @param Dk surrogate event indicators
#' @param Ct censoring times for surrogate event (requires to compute inverse probability of censoring weights)
#' @param K kernel function (e.g. dnorm)
#' @param h kernel bandwith
#' @return Score function value for each first step estimator beta (later to be combined to correct initial estimator)
#' @export
Sk_sym = function(lp, Z, Xk, Dk, Ct, K, h)
{
  n = nrow(Z)
  
  # ECDF of G
  X.order = order(Xk)
  C.sort = sort(Ct)
  Ctail = 0
  wZ = rep(0,n)
  w = 0
  next.X = n
  for(i.X in (n-1):1)
  {
    i = X.order[i.X]
    if(Xk[i]<Xk[X.order[next.X]])
      next.X = i.X
    while (Ctail < n) 
    {
      if(C.sort[n-Ctail] < Xk[i])
        break
      Ctail = Ctail + 1
    }
    
    if( (next.X == n) | (Dk[i]== 0))
      next
    GXi2 = (Ctail/n)^2
    w = w + (n-next.X)/GXi2
    
    js = X.order[(next.X+1):n]
    Kbz = K((lp[i]-lp[js])/h)/h
    wZ[i] = sum(Kbz)/GXi2
    wZ[js] = wZ[js] - Kbz/GXi2
  }
  
  return(drop(wZ%*%Z)/w)
}


#' Author's Original Implementation of Score Function for Second Step Estimator
#' Corresponds to a single surrogate outcome (run one time per surrogate)
#' @param lp linear predictor for each patient (X transpose beta)
#' @param Z covariates (design matrix)
#' @param Xk surrogate event times
#' @param Dk surrogate event indicators
#' @param Ct censoring times for surrogate event (requires to compute inverse probability of censoring weights)
#' @param K kernel function (e.g. dnorm)
#' @param h kernel bandwith
#' @param V vector of perturbation values
#' @return Score function value for each first step estimator beta, for each set of perturbation values
#' @export
Sk_sym_perturb = function(lp, Z, Xk, Dk, Ct, K, h, V)
{
  n = nrow(Z)
  n.perturb = ncol(V)
  
  # ECDF of G
  X.order = order(Xk)
  C.order = order(Ct)
  GX = rep(0,n)
  Ctail.count = 0
  Ctail.sum = rep(0,n.perturb)
  wZ = matrix(0,n.perturb,n)
  w = rep(0,n.perturb)
  next.X = n
  for(i.X in (n-1):1)
  {
    # print(i.X)
    i = X.order[i.X]
    if(Xk[i]<Xk[X.order[next.X]])
      next.X = i.X
    while (Ctail.count < n) 
    {
      if(Ct[C.order[n-Ctail.count]] < Xk[i])
        break
      Ctail.sum = Ctail.sum + V[C.order[n-Ctail.count],]
      Ctail.count = Ctail.count + 1
    }
    GXi2 = Ctail.sum^2
    
    if( (next.X == n) | (Dk[i]== 0))
      next    
    js = X.order[(next.X+1):n]
    njs = n-next.X
    
    Vij.GXi2 = V[i,]*matrix(t(V[js,]),n.perturb)/GXi2
    w = w + apply(Vij.GXi2,1,sum)
    
    Kbz.GXi2 = K((lp[i,]-
                    matrix(t(lp[js,]),n.perturb))/h)/h * Vij.GXi2
    wZ[,i] = apply(Kbz.GXi2,1,sum)
    wZ[,js] = wZ[,js] - Kbz.GXi2
  }
  
  return(drop(wZ%*%Z)/w)
}

#' Author's Perturbation Analysis for First Step Estimator (only labelled data)
#' @param delta event indicator
#' @param Z covariates (design matrix)
#' @param KC kernel matrix evaluated based on censoring times
#' @param V vector of perturbation values
#' @param init initial values for coefficients beta, 0 by default
#' @param tol convergence tolerance for beta
#' @param maxit maximum number of iterations
#' @param min.factor minimum multiplicative decrease at a given iteration to not enter line search
#' @param ls.factor decrease in step size when conducting line search
#' @param max.move maximum step size within line search as proportion of newton step
#' @param link link function for score equation to be solves
#' @param dlink derivative of link function
#' @return A list of coefficients beta corresponding to each covariate column in Z, for each set of perturbations
#' @export
init.beta.perturb = function(delta,Z, KC, 
                             V, init = rep(0,ncol(Z)), tol=1e-7,
                             maxit = 100,min.factor = 0.75,
                             ls.factor = 0.75,max.move = 1,
                             link = expit, dlink = dexpit)
{
  n = nrow(Z)
  KC = t(V*KC)
  KCd = drop(KC%*%delta)
  hC = rep(0,n)
  oldscore = NULL
  max.dbeta = max.move
  
  for(k in 1:maxit) 
  {
    # print(k)
    # print(init)
    # print(oldscore)
    lp = drop(Z%*%init)
    hC.flag = rep(TRUE,n)
    gij = wZbar = matrix(0,n,n)
    hHess = rep(0,n)
    for(kk in 1:maxit)
    {
      lij = outer(hC[hC.flag],lp,"+")
      gij[hC.flag,] = link(lij)
      tmp = KC[hC.flag,]*gij[hC.flag,]
      wZbar[hC.flag,] = KC[hC.flag,]*dlink(lij)
      if(sum(hC.flag)>=2)
      {
        hscore = apply(tmp,1,sum)-KCd[hC.flag]
        hHess[hC.flag] = apply(wZbar[hC.flag,],1,sum)
      }else
      {
        hscore = sum(tmp)-KCd[hC.flag]
        hHess[hC.flag] = sum(wZbar[hC.flag,])
      }
      
      dhC = hscore/hHess[hC.flag]
      dhC = sign(dhC)*pmin(abs(dhC),max.move)
      kk.flag = abs(hscore) > tol
      if(!any(kk.flag))
        break
      hC[hC.flag][kk.flag] = hC[hC.flag][kk.flag] - dhC[kk.flag]
      hC.flag[hC.flag] = kk.flag
    }
    if(kk >= maxit)
      stop("Numerical error when computing h0(Ci)")
    Zbar =  (wZbar%*%Z) / hHess 
    
    gi = link(hC+lp)
    bscore = drop(t(V*Z)%*% (delta - gi))
    if(!is.null(oldscore))
      if((sum(oldscore^2)*min.factor) <= sum(bscore^2))
      {
        init = init+dinit
        dinit = dinit*ls.factor
        if(max(abs(dinit))<tol)
        {
          if(max(abs(oldscore)) > 1e-1)
            stop(paste("Algorithm stops in line-search. Target tol: ",
                       tol, ". Current tol: ", max(abs(oldscore)),
                       ". ", sep = ''))
          
          if(max(abs(oldscore)) > 1e-6)
            warning(paste("Algorithm stops in line-search. Target tol: ",
                          tol, ". Current tol: ", max(abs(oldscore)),
                          ". ", sep = ''))
          break
        }
        init = init - dinit
        next
      }
    oldscore = bscore
    bHess = t(V*dlink(hC+lp)*Z) %*% (Zbar-Z)
    dinit = solve(bHess,bscore)
    dsize = sqrt(sum(dinit^2))
    if(dsize > max.dbeta)
    {
      dinit = dinit/dsize*max.dbeta
    }else
    {
      max.dbeta = dsize
    }
    # print(max(abs(bscore)))
    if(all(abs(bscore)<tol))
      break
    init = init - dinit
  }
  if(k >=maxit)
    stop("Numerical error when computing beta_delta")
  
  return(init)
}


#' Author's Implementation of Aggregated Rank Correlation to Conduct Cross-Validation for Second-Step Estimator
#' @param lp linear predictor for each patient (X transpose beta)
#' @param Xk surrogate event times
#' @param Dk surrogate event indicators
#' @param Ct censoring times for surrogate event (requires to compute inverse probability of censoring weights)
#' @return Aggregated rank correlation
#' @export
cv.concordance = function(lp,Xk, Dk, Ct)
{
  n = nrow(lp)
  nlam = ncol(lp)
  
  # ECDF of G
  X.order = order(Xk)
  C.sort = sort(Ct)
  Ctail = 0
  wc = rep(0,nlam)
  w = 0
  next.X = n
  for(i.X in (n-1):1)
  {
    i = X.order[i.X]
    if(Xk[i]<Xk[X.order[next.X]])
      next.X = i.X
    while (Ctail < n) 
    {
      if(C.sort[n-Ctail] < Xk[i])
        break
      Ctail = Ctail + 1
    }
    
    if( (next.X == n) | (Dk[i]== 0))
      next
    GXi2 = (Ctail/n)^2
    w = w + (n-next.X)/GXi2
    
    js = X.order[(next.X+1):n]
    wc = wc + apply(lp[i,] > t(lp[js,]),1,sum)/GXi2
  }
  
  return(cbind(wc,w))
}

#' Author's Implementation of Adaptive Ridge Regression for Second-Step Estimator
#' @param betak set of first step estimators obtained from perturbation analysis
#' @param Skb set of score functions for first step estimators obtained from perturbation
#' @param a number of surrogate outcomes
#' @return A list containing W.hat (optimal loading matrix), optimal lambda, penalty factor
#' @export
W_hat_adaPCA_ridge = function(betak,Skb, a = 2)
{
  p = nrow(betak)
  B = ncol(betak)
  Kp = ncol(Skb)
  
  Skb.sd = apply(Skb, 2, sd)
  Skb.std = scale(Skb)
  eig.var = eigen(var(Skb.std), symmetric = T)
  PSkb =  Skb.std %*% eig.var$vectors
  ada.pen = 1/eig.var$values^a
  ada.pen = ada.pen/mean(ada.pen)
  
  W.hat = matrix(0,p,Kp)
  lambda = rep(0, p)
  
  for(j in 1:p)
  {
    tmp.fit = cv.glmnet(PSkb,betak[j,],alpha = 0, penalty.factor = ada.pen)
    W.hat[j,]=drop(diag(1/Skb.sd) %*% eig.var$vectors %*% coef(tmp.fit,s="lambda.min")[-1])
    lambda[j] = tmp.fit$lambda.min
  }
  
  return(list(W.hat = W.hat, lambda =lambda, 
              ada.pen = ada.pen))
}


#' Author's Implementation of Adaptive Lasso Regression for Second-Step Estimator
#' @param betak set of first step estimators obtained from perturbation analysis
#' @param Skb set of score functions for first step estimators obtained from perturbation
#' @param a number of surrogate outcomes
#' @return A list containing W.hat (optimal loading matrix), optimal lambda, penalty factor
#' @export
W_hat_adaPCA_lasso = function(betak,Skb, a=1)
{
  p = nrow(betak)
  B = ncol(betak)
  Kp = ncol(Skb)
  
  Skb.sd = apply(Skb, 2, sd)
  Skb.std = scale(Skb)
  eig.var = eigen(var(Skb.std), symmetric = T)
  PSkb =  Skb.std %*% eig.var$vectors
  ada.pen = 1/eig.var$values^a
  ada.pen = ada.pen/mean(ada.pen)
  
  W.hat = matrix(0,p,Kp)
  lambda = rep(0, p)
  
  for(j in 1:p)
  {
    tmp.fit = cv.glmnet(PSkb,betak[j,], penalty.factor = ada.pen)
    W.hat[j,]=drop(diag(1/Skb.sd) %*% eig.var$vectors %*% coef(tmp.fit,s="lambda.min")[-1])
    lambda[j] = tmp.fit$lambda.min
  }
  
  return(list(W.hat = W.hat, lambda =lambda, 
              ada.pen = ada.pen))
}


#' Author's Original Function for Generating Simulated Data
#' @param n total number of observations
#' @param beta0 true effect of each covariate
#' @param sigmaZ variance-covariance matrix for covariates
#' @param m number of labelled observations
#' @param misT noise function relating surrogate outcome times to true outcome times
#' @param Nperturb additional perturbations fo be generated (default none)
#' @param inv.link inverse of link function relating censoring time to probability of event
#' @param inv.h inverse of link function relating linear predictor to survival time
#' @return A list containing two datasets: data (covariates and true outcomes) and dataS (surrogate outcomes)
#' @export
sim.gen.rev = function(n,beta0,sigmaZ,m,
                       misT, Nperturb=0,
                       inv.link = logit, 
                       inv.h = function(x) exp(x/3)*4)
{
  p = length(beta0)
  Z = mvrnorm(n,rep(0,p),sigmaZ^2*(0.2+0.8*diag(1,p)))
  u = runif(n)
  Tevent = inv.h((inv.link(u) - Z%*%beta0))
  C = runif(n,0,12)
  delta = (Tevent<=C)
  if(m < n)
    delta[(m+1):n] = NA
  
  Tstar = misT(Tevent)
  
  X1 = pmin(Tstar[,1],C)
  delta1 = (Tstar[,1]<=C)
  
  X2 = pmin(Tstar[,2],C)
  delta2 = (Tstar[,2]<=C)
  
  dataS = data.frame(X1=X1,delta1=delta1,
                     X2=X2,delta2=delta2)
  
  data = data.frame(delta = delta, C=C,
                    Tevent = Tevent, Tstar = Tstar)
  data$Z = Z
  if(Nperturb > 0)
    data$V = matrix(rbeta(n*Nperturb,0.5,1.5)*4,n)
  
  return(list(data=data,dataS=dataS))
}


#' Author's Original Function for applying low log-linear noise to event times
#' @param t event times
#' @return matrix with two mismeasured surrogate times
#' @export
noise.loglinear.S = function(t)
{
  mu = c(0,0.5,-0.25,0)
  sigma = c(0.5,0.15,0.35,0.45)
  
  n = length(t)
  
  Di = rbinom(n,1,0.5)  
  Di2 = rbinom(n,1,0.5)
  
  error = matrix(rnorm(n*4),n,4)
  
  epsilon = (Di*(error[,1]*sigma[1]+mu[1]) +
               (1-Di)*(error[,2]*sigma[2]+mu[2]))
  Tstar = exp(epsilon)*t
  
  epsilon2 = (Di2*(error[,3]*sigma[3]+mu[3]) +
                (1-Di2)*(error[,4]*sigma[4]+mu[4]))
  Tstar2 = exp(epsilon2)*t
  
  return(cbind(Tstar,Tstar2))
}

#' Author's Original Function for applying high log-linear noise to event times
#' @param t event times
#' @return matrix with two mismeasured surrogate times
#' @export
noise.loglinear.L = function(t)
{
  mu = c(1,-0.5,0,1.5)
  sigma = c(1.5,0.5,1,0.5)
  
  n = length(t)
  
  Di = rbinom(n,1,0.5)  
  Di2 = rbinom(n,1,0.5)
  
  error = matrix(rnorm(n*4),n,4)
  
  epsilon = (Di*(error[,1]*sigma[1]+mu[1]) +
               (1-Di)*(error[,2]*sigma[2]+mu[2]))
  Tstar = exp(epsilon)*t
  
  epsilon2 = (Di2*(error[,3]*sigma[3]+mu[3]) +
                (1-Di2)*(error[,4]*sigma[4]+mu[4]))
  Tstar2 = exp(epsilon2)*t
  
  return(cbind(Tstar,Tstar2))
}

#' Author's Original Function for applying low exponential noise to event times
#' @param t event times
#' @return matrix with two mismeasured surrogate times
#' @export
noise.exp.S = function(t)
{
  n = length(t)
  r1 = 0.9
  r2 = 0.05
  
  Tstar = pmax(t+ r2*(rexp(n)-0.5),0)
  Tstar2 = pmax(t+ r2*(rchisq(n,df=1)-0.5),0)
  
  return(cbind(Tstar,Tstar2))
}

#' Author's Original Function for applying high exponential noise to event times
#' @param t event times
#' @return matrix with two mismeasured surrogate times
#' @export
noise.exp.L = function(t)
{
  n = length(t)
  r1 = 0.45
  r2 = 0.5
  
  Tstar = pmax(t+ r2*(rexp(n)-0.5),0)
  Tstar2 = pmax(t+ r2*(rchisq(n,df=1)-0.5),0)
  
  return(cbind(Tstar,Tstar2))
}

