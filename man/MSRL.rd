\name{MSRL.cv}
\alias{MSRL.cv}
\title{Fit the multivariate square-root lasso.}
\usage{
MSRL.cv(X, Y, nlambda, lambda.vec = NULL,
    weighted = FALSE, standardize = FALSE, ADMM = FALSE, nfolds = NULL, 
    delta = .01, tol = 1e-8, quiet = TRUE, inner.quiet=TRUE)
}
\arguments{

 \item{X}{An \eqn{n \times p} matrix of predictors (unstandardized).}
  \item{Y}{An \eqn{n \times q} matrix of responses (unstandardized).}
  \item{nlambda}{The number of candidate tuning parameters to be considered.}
  \item{lambda.vec}{A vector of positive tuning parameters to be considered (overriding internal candidate set). Should be in descending order.}
  \item{weighted}{Should the penalty be weighted by response variables' standard deviations?}
  \item{standardize}{Should predictors be standardized for model fitting?}
  \item{ADMM}{Determines whether the algorithm starts using ADMM or APG. By default, if \eqn{n > q}, will start with APG and move to ADMM when appropriate.}
  \item{nfolds}{Number of folds for cross-validation.}
  \item{delta}{Ratio of maximum to minimum candidate tuning parameters. }
  \item{tol}{Convergence tolerance parameter}
  \item{quiet}{Print progress?  }
  \item{quiet}{Print iteration-by-iteration progress? Should be used for diagnostics.}
}
\value{
 	\item{beta}{A sparse matrix with estimated regression coefficients. Use MSRL.coef or MSRL.predictor for coefficient extraction or prediction.}
	\item{sparsity.mat}{A matrix indicating the degree of sparsity at each estimated \eqn{\beta}}
	\item{err.pred}{A matrix of sum-of-squared prediction errors from cross-validation.}
	\item{err.wpred}{A matrix of weighted sum-of-squared prediction errors from cross-validation.}
	\item{err.spec}{A matrix of spectral norm prediction errors from cross-validation.}
	\item{err.nuc}{A matrix of nuclear norm prediction errors from cross-validation.}
	\item{Y.mean}{Means for the response variables; used by MSRL.coef and MSRL.predictor.}
	\item{X.mean}{Means for the predictor variables; used by MSRL.coef and MSRL.predictor.}
	\item{Y.sd}{Standard deviations for the response variables; used by MSRL.coef and MSRL.predictor.}
	\item{X.sd}{Standard deviations for the predictor variables; used by MSRL.coef and MSRL.predictor.}
	\item{lambda.vec}{A vector of the candidate tuning parameters}
	\item{lam.min}{The tuning parameter value which had the minimum cross-validation sum-of-squared prediction errors.}
	\item{standardize}{Were predictors standardized?}
         }


\description{
  A function for fitting the multivariate square-root lasso and performing cross-validation.
}
