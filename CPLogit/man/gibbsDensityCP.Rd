% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{gibbsDensityCP}
\alias{gibbsDensityCP}
\title{Grouped Bayesian linear regression - only regressors.}
\usage{
gibbsDensityCP(
  X_list0,
  Y_list0,
  Z_prior,
  psi_par,
  H,
  H_upper,
  alpha,
  b_prior,
  Q_prior,
  n_iter0
)
}
\arguments{
\item{X_list0}{list of n_i x p matrices of covariates (including intercept)}

\item{Y_list0}{list of n_i x 1 vectors of response variable}

\item{Z_prior}{prior information on clustering}

\item{psi_par}{value for the penalizations parameter}

\item{H}{number or clusters}

\item{H_upper}{upper limit for number of clusters}

\item{alpha}{value for the Dirichlet Process concentration parameter}

\item{b_prior}{mean vector of the normal multivariate prior for the coefficients}

\item{Q_prior}{covariance matrix of the normal multivariate prior for the coefficients}

\item{n_iter0}{total number of iterations}
}
\description{
This function does a gibbs sampling for a linear regression model with only regressors coef
note: this is a trick to clustering means of a normal distributed data
}
