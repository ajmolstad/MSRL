\name{MSRL.cv}
\alias{MSRL.cv}
\title{Fit the multivariate square-root lasso}
\description{A function for performing cross-validation to select tuning parameters for the multivariate square-root lasso}
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
  \item{delta}{Ratio of minimum to maximum candidate tuning parameters. }
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
\details{This function is used to tune the multivariate square-root lasso using the \eqn{L_1}-vector norm penalty. Specifically, this function performs cross-validation for the estimator 
\deqn{\text{arg min}_{\beta \in \mathbb{R}^{p \times q}} \frac{1}{\sqrt{n}}\|Y - X\beta\|_* + \lambda \sum_{j,k}w_{j,k}|\beta_{j,k}|.}
When \code{weighted = FALSE}, the \eqn{w_{j,k} = 1}. When \code{weighted = TRUE}, \eqn{w_{j,k} = 1/\hat\sigma_{k}}, where \eqn{\hat\sigma_k} is the standard deviation of the $k$th response in the training set. 
For cross-validation, we measure sum-of-squared prediction error, weighted (by inverse variances of the responses) sum of squared prediction error, spectral norm prediction error, and nuclear norm prediction error. 
}

\example{
set.seed(1)
p <- 50
n <- 100
q <- 20
ntest <- 100

# --------------------------------
# Generate predictors
# --------------------------------
SigmaX <- matrix(0, nrow=p, ncol=p)
for(k in 1:p){
  for(j in 1:p){
    SigmaX[j,k] <- .5^abs(j-k)
  }
}
diag(SigmaX) <- 1
eo <- eigen(SigmaX)
SigmaXsqrt <- eo$vec\%*\%diag(eo$val^.5)\%*\%t(eo$vec)
X <- matrix(rnorm(n*p), nrow=n)\%*\%SigmaXsqrt
Xtest <- matrix(rnorm(ntest*p), nrow=ntest)\%*\%SigmaXsqrt
# -----------------------------------
# Generate regression coefficients
# -----------------------------------
beta <- matrix(0, nrow=p, ncol=q)
beta[sample(1:(p*q), floor(p*q)*.04, replace=TRUE)] <- 1

# -------------------------------------
# Generate responses
# -------------------------------------
Sigma <- matrix(0, nrow=q, ncol=q)
for(k in 1:q){
  for(j in 1:q){
    Sigma[j,k] <- .9^abs(j-k)
  }
}
diag(Sigma) <- 1
eo <- eigen(Sigma)
Sigmasqrt <- eo$vec\%*\%diag(eo$val^.5)\%*\%t(eo$vec)
Y <- X\%*\%beta + matrix(rnorm(n*q), nrow=n)\%*\%Sigmasqrt
Ytest <- Xtest\%*\%beta + matrix(rnorm(ntest*q), nrow=ntest)\%*\%Sigmasqrt

mod1 <- MSRL.cv(X = X, Y = Y, nlambda = 10, standardize = FALSE, 
                ADMM = FALSE, nfolds = NULL,  weighted = FALSE,  delta = .25, 
                tol = 1e-8, quiet = FALSE, inner.quiet = TRUE)
fitted.preds <- MSRL.predict(Xnew = Xtest, fit = mod1, lambda = mod1$lam.min)
str(fitted.preds)                
                
}