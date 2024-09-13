\name{MSRL.coef}
\alias{MSRL.coef}
\title{Extract estimated regression coefficients from the multivariate square-root lasso.}
\description{A function for extracting regression coefficients along the solution path for the multivariate square-root lasso.}
\usage{
MSRL.coef(fit, lambda = NULL)
}
\arguments{
	\item{fit}{An object of class MSRL obtained from MSRL.cv}
	\item{lambda}{The value of the tuning parameter at which to make predictions. Must belong to the set fit\eqref{\$}lambda.vec}.}
}
\value{
	\item{beta0}{The estimated intercept at the tuning parameter value \eqn{\lambda}, the third input argument.}
	\item{beta}{The estimated regression coefficient matrix at the tuning parameter value \eqn{\lambda}, the third input argument.}
}

\description{
  Extract estimated regression coefficients from the multivariate square-root lasso.
}
