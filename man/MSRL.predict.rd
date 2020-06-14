\name{MSRL.predict}
\alias{MSRL.predict}
\title{Make predictions using the multivariate square-root lasso.}
\usage{
MSRL.predict(Xnew, fit, lambda = NULL)
}
\arguments{
	\item{Xnew}{An \eqn{n_{\rm new} \times p} matrix of predictors.}
	\item{fit}{An object of class MSRL obtained from MSRL.cv}
	\item{lambda}{The value of the tuning parameter at which to make predictions. Must belong to the set fit\eqref{\$}lambda.vec}.}
 % \item{time}{An \eqn{n}-variate or containing the failure/censoring times on the original scale -- not log-transformed.}
 %  \item{status}{An \eqn{n}-variate binary vector of same length as time: 0 indicates censored, 1 indicates failure. }
 %  \item{Z}{An \eqn{n \times q} design matrix for the linear mean function. Note that the first column should contain all ones to account for the intercept. We recommend construction using \code{model.matrix}. }
 %  \item{K}{Candidate kernel matrices in the form of an array of dimension \eqn{n \times n \times M}. The algorithm will work best if these kernels have diagonal entries on similar scales. }
 %  \item{tol}{The convergence tolerance. Default is \code{1e-7}.}
 %  \item{max.iter.MM}{The maximum number of iterations for the inner M-step algorithm.}
 %  \item{max.iter}{The maximum number of total EM-iterations.}
 %  \item{kern.type}{A character argument that can be set to either \code{K+I} or \code{multi-K}, indicating which algorithm to use. We highly recommend using \code{multi-K} in all cases, although if \eqn{M=1}, \code{K+I} can be faster, but less stable.}
 %  \item{quiet}{\code{TRUE/FALSE} -- print algorithm progress?}
 %  \item{initializer}{Only used when \code{kern.type = "multi-K"}: either \eqn{0,1} indicating the type of initialization used for the Monte-Carlo EM. The default is 0,  sets all variance components equal and gets an initial estimate of \eqn{\hat{\beta}} using the updating equation from Algorithm 2. The code 1 which fits a variance components model to the IPW-mean-imputed dataset.}
 %  \item{max.samples}{An upper bound on \eqn{s_k}, the Monte-Carlo sample size for the \eqn{k}th iteration. Note that the final imputed values of log-survival for censored subjects will be the average of \code{max.samples} Monte-Carlo draws.}
}
\value{
	\item{pred}{Predicted values; each row corresponding to a row of Xnew.}
	\item{beta0}{The estimated intercept at the tuning parameter value \eqn{\lambda}, the third input argument.}
	\item{beta}{The estimated regression coefficient matrix at the tuning parameter value \eqn{\lambda}, the third input argument.}
  %  \item{beta}{\eqn{\hat{\beta}}: The estimated regression coefficient vector corresponding to the columns of \code{Z}. }
  % \item{sigma2}{\eqn{\hat{\sigma}^2}: The estimated variance components: a vector of length \eqn{M+1}, with the final element corresponding to the variance of \eqn{\epsilon}. }
  % \item{Tout}{The log-failure and imputed log-failure times obtained from our MCEM algorithm. These are primarily to be used in the prediction function.}
  % \item{Yimpute}{The mean-imputed values of the training time-to-failures based on the method of Datta et al (2005). }
}

\description{
  Make predictions using the multivariate square-root lasso.
}
