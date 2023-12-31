% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Author_Original_Functions.R
\name{init.beta}
\alias{init.beta}
\title{Author's Original Implementation of First Step Estimator (only labelled data)}
\usage{
init.beta(
  delta,
  Z,
  KC,
  init = rep(0, ncol(Z)),
  tol = 1e-07,
  maxit = 100,
  min.factor = 0.75,
  ls.factor = 0.75,
  max.move = 1,
  link = expit,
  dlink = dexpit
)
}
\arguments{
\item{delta}{event indicator}

\item{Z}{covariates (design matrix)}

\item{KC}{kernel matrix evaluated based on censoring times}

\item{init}{initial values for coefficients beta, 0 by default}

\item{tol}{convergence tolerance for beta}

\item{maxit}{maximum number of iterations}

\item{min.factor}{minimum multiplicative decrease at a given iteration to not enter line search}

\item{ls.factor}{decrease in step size when conducting line search}

\item{max.move}{maximum step size within line search as proportion of newton step}

\item{link}{link function for score equation to be solved}

\item{dlink}{derivative of link function}
}
\value{
A list of coefficients beta corresponding to each covariate column in Z
}
\description{
Author's Original Implementation of First Step Estimator (only labelled data)
}
