\name{MSRL.predict}
\alias{MSRL.predict}
\title{Make predictions using a fitted model}
\description{A function for making predictions using the multivariate square-root lasso model fit along the solution path.}
\usage{
MSRL.predict(Xnew, fit, lambda = NULL)
}
\arguments{
	\item{Xnew}{An \eqn{n_{\rm new} \times p} matrix of predictors.}
	\item{fit}{An object of class MSRL obtained from MSRL.cv}
	\item{lambda}{The value of the tuning parameter at which to make predictions. Must belong to the set fit\eqref{\$}lambda.vec}.}
}
\value{
	\item{pred}{Predicted values; each row corresponding to a row of Xnew.}
	\item{beta0}{The estimated intercept at the tuning parameter value \eqn{\lambda}, the third input argument.}
	\item{beta}{The estimated regression coefficient matrix at the tuning parameter value \eqn{\lambda}, the third input argument.}
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

