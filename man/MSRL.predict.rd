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

