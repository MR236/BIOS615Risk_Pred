#' New Implementation of First Step Estimator based on Penalized Splines
#' @param delta event indicator
#' @param C censoring times
#' @param Z covariates (design matrix)
#' @param df degrees of freedom for spline
#' @param knots number of knots for spline
#' @param lambda penalty terms lambda to consider in cross-validation
#' @return A list of coefficients beta corresponding to each covariate column in Z
#' @export
init.beta.new  = function(delta,C,Z, df=3, knots = min(9,floor(length(delta)/500)), lambda=NULL)
{
  Y <- delta
  if (knots == 0) {knot_locations <- NULL}
  else {knot_locations <- as.numeric(quantile(C,seq(0,1,1/(knots+1))))
  knot_locations <- knot_locations[2:length(knot_locations)]}
  # construct spline using censoring times
  spline_X <-  splines::bs(C, knots=knot_locations, degree=df)
  X <- cbind(spline_X,Z)
  # find optimal penalty term using cross-validation
  cv.lasso <- cv.glmnet(X, Y, alpha = 1, family = "binomial", lambda = lambda)
  coef_lasso = coef(cv.lasso, cv.lasso$lambda.min)
  # extract coefficients corresponding to covariates rather than spline
  Z_coef = coef_lasso[c((length(coef_lasso)-ncol(Z)+1):length(coef_lasso))]
  return(list(beta = Z_coef, lambda = cv.lasso$lambda.min))
}


#' New Implementation of Score Function for Second Step Estimator
#' Corresponds to a single surrogate outcome (run one time per surrogate)
#' @param lp linear predictor for each patient (X transpose beta)
#' @param Z covariates (design matrix)
#' @param Xk surrogate event times
#' @param Dk surrogate event indicators
#' @param Ct censoring times for surrogate event (requires to compute inverse probability of censoring weights)
#' @param h kernel bandwith
#' @return Score function value for each first step estimator beta
#' @export
Sk_sym_new = function(lp, Z, Xk, Dk, Ct, h)
{
  n = nrow(Z)
  # selection of a second bandwith b, outside of which kernel evaluations can be
  # truncated without affecting estimation
  bandwith <- function(h, tol=10^-6) {
    val <- 2*h^2*log(sqrt(2*pi)*h*tol)
    return(sqrt(abs(val)))
  }
  b <- bandwith(h)
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
    # Toss out kernel evaluations that are beyond the bandwith (sufficiently far apart)
    js = js[lp[i]-lp[js] < b]
    # Replace dnorm with direct kernel evaluation (much more computationally efficient)
    Kbz = exp(-((lp[i]-lp[js])/h)^2/2)
    wZ[i] = sum(Kbz)/GXi2
    wZ[js] = wZ[js] - Kbz/GXi2
  }
  return(drop(wZ%*%Z)/(w*sqrt(2*pi)*h))
}
