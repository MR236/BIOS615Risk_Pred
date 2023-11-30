
#' Exponential Link Function (Inverse Logit)
#' @param x value
#' @return expit(x)
#' @export
expit <- function(x){
  1/(1+exp(-x))
}

#' Derivative of Exponential Link Function
#' @param x value
#' @return dexpit(x)
#' @export
dexpit <- function(x){
  expit(x)*(1-expit(x))
}

#' Logit Link Function
#' @param x value
#' @return logit(x)
#' @export
logit <- function(x){
  log(x/(1-x))
}

#' Derivative of Normal PDF
#' @param x value
#' @return ddnorm(x)
#' @export
ddnorm = function(x)
{
  -x*exp(-x^2/2)/sqrt(2*pi)
}

#' Second Derivative of Normal PDF
#' @param x value
#' @return dddnorm(x)
#' @export
dddnorm = function(x)
{
  (x^2-1)*exp(-x^2/2)/sqrt(2*pi)
}

#' Loglog Function
#' @param x value
#' @return -log(-log(x))
#' @export
loglog = function(x)
{
  -log(-log(x))
}

#' Expexp Function
#' @param x value
#' @return exp(-exp(-x))
#' @export
expexp = function(x)
{
  exp(-exp(-x))
}

#' Derivative of Expexp Function
#' @param x value
#' @return d/dx exp(-exp(-x))
#' @export
dexpexp = function(x)
{
  exp(-exp(-x))*exp(-x)
}

#' Weighted Mean Function
#' @param x value
#' @param wgt weight
#' @return weighted mean
#' @export
wtd.mean = function(x,wgt)
{
  drop(x%*%wgt)/sum(wgt)
}

#' Weighted Standard Deviation Function
#' @param x value
#' @param wgt weight
#' @return weighted sd
#' @export
wtd.sd = function(x, wgt)
{
  sqrt( drop(((x - wtd.mean(x,wgt))^2)%*% wgt) / sum(wgt))
}