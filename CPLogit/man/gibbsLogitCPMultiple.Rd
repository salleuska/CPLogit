% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{gibbsLogitCPMultiple}
\alias{gibbsLogitCPMultiple}
\title{Grouped logistic regression using Centered Partition process - multiple prior partitions}
\usage{
gibbsLogitCPMultiple(
  X_list0,
  Y_list0,
  Z_prior,
  psi_par,
  H,
  H_upper,
  alpha,
  a_prior,
  tau_prior,
  b_prior,
  Q_prior,
  n_iter0
)
}
\arguments{
\item{X_list0}{list of n_i x p matrices of covariates (including intercept)}

\item{Y_list0}{list of n_i x 1 vectors of response variable}

\item{Z_prior}{matrix of possible partitions for clustering allocation}

\item{psi_par}{value for the penalizations parameter}

\item{H}{number or clusters}

\item{H_upper}{upper limit for number of clusters}

\item{alpha}{value for the Dirichlet Process concentration parameter}

\item{a_prior}{mean of the normal prior on intercept}

\item{tau_prior}{inverse of variance of the normal prior on intercept}

\item{b_prior}{mean vector of the normal multivariate prior for the coefficients}

\item{Q_prior}{covariance matrix of the normal multivariate prior for the coefficients}

\item{n_iter0}{total number of iterations}
}
\description{
This function does a gibbs sampling for a logistic regression model using a DP prior for 
clustering regression coefficients. Sampling for the logistic regressions uses Polya-Gamma 
data augmentation.
}
