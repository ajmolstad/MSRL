# ---------------------------------------
# Set model parameters
# --------------------------------------
uu <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
nreps <- 100
t0 <- expand.grid(rhos = rep(c(0.3, 0.5, 0.7, 0.9, 0.95), each=nreps), Model = c(1,2))
t1 <- expand.grid(rhos = rep(c(3, 5, 10, 25, 50), each=nreps), Model=c(3, 4))
t2 <- expand.grid(rhos = rep(c(2, 3, 4, 5, 6), each=nreps), Model=c(5,6))
params <- rbind(t0, t1, t2)
rho <- params[uu,1]
Model <- params[uu,2]
quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
} 


if(rho > 1){
  savename <- paste("~/blue/MSRL_R1/Simulations/L1/Results/Model", Model, "_rho", rho, "_", uu %% nreps + 1, ".RDS", sep="")
} else {
  savename <- paste("~/blue/MSRL_R1/Simulations/L1/Results/Model", Model, "_rho", rho*100, "_", uu %% nreps + 1, ".RDS", sep="")
}

# -----------------------------
# Preliminaries 
# -----------------------------
p <- 500
q <- 50
n <- 200
nval <- 200
ntest <- 1000

# ------------------------------
# Process results 
# ------------------------------
TPR <- function(beta, betahat){
  return(sum(beta!=0 & betahat !=0)/sum(beta!=0))
}

FPR <- function(beta, betahat){
  return(sum(beta==0 & betahat !=0)/sum(beta==0))
}

named.list = function(...) { 
  l = list(...)
  names(l) = as.character(match.call()[-1] )
  l
}


Results <- NULL
# ------------------------------
# Readin packages 
# ------------------------------
library(Rcpp)
library(Matrix)
library(flare)
library(MRCE)
library(glasso)
library(CVXR)

source("~/blue/MSRL_R1/Functions/MSRL_Functions.R")
source("~/blue/MSRL_R1/Functions/getResults.R")
source("~/blue/MSRL_R1/Functions/SUR.R")
sourceCpp("~/blue/MSRL_R1/Functions/APG_MRCE.cpp")

# url <- "http://cran.r-project.org/src/contrib/Archive/camel/camel_0.2.0.tar.gz"
# pkgFile <- "camel_0.2.0.tar.gz"
# download.file(url = url, destfile = pkgFile)
# install.packages(pkgs=pkgFile, type="source", repos=NULL)
# library(camel)

set.seed(uu)
# -------------------------------
# generate predictors 
# -------------------------------
SigmaX <- matrix(0, nrow=p, ncol=p)
for(k in 1:p){
  for(j in 1:p){
    SigmaX[j,k] <- .5^abs(j-k)
  }
}
diag(SigmaX) <- 1
eo <- eigen(SigmaX)
X <- matrix(rnorm(n*p), nrow=n)%*%eo$vec%*%diag(eo$val^.5)%*%t(eo$vec)
Xtest <-  matrix(rnorm((ntest)*p), nrow=ntest)%*%eo$vec%*%diag(eo$val^.5)%*%t(eo$vec)
Xval <-  matrix(rnorm((nval)*p), nrow=nval)%*%eo$vec%*%diag(eo$val^.5)%*%t(eo$vec)

# -------------------------------
# generate coefficient matrix 
# -------------------------------
beta <- matrix(0, nrow=p, ncol=q)
for(j in 1:q){
  nums <- 5
  temp <- sample(1:p, nums)
  beta[temp,j] <- rnorm(5)
}
source("~/blue/MSRL_R1/Functions/GenerateData.R")
# ---------------------------------------------
# Our method using AccPG
# ---------------------------------------------

OursCV.comptimeADMM <- system.time(
temp <- MSRL.cv(Y = Y, X = X, ADMM = TRUE, nlambda = 100, delta = .1, tol = 1e-10,
  quiet = FALSE, nfolds = NULL, inner.quiet=TRUE))

OursCV.comptimePGD <- system.time(
temp <- MSRL.cv(Y = Y, X = X, ADMM = FALSE, nlambda = 100, delta = .1, tol = 1e-10, 
  quiet = FALSE, nfolds = NULL, inner.quiet=TRUE))

Results <- append(Results, list("OursCV.comptimeADMM" = OursCV.comptimeADMM,
  "OursCV.comptimePGD" = OursCV.comptimePGD))

prederr.val <- rep(0, length(temp$lambda.vec))
for(jj in 1:length(temp$lambda.vec)){
  beta.mat <- MSRL.coef(temp, lam = temp$lambda.vec[jj])$beta
  Ypred.val <- MSRL.predict(Xval, temp, lam = temp$lambda.vec[jj])$pred
  prederr.val[jj] <- sum((Ypred.val - Yval)^2)/(nval*q)
}

sparsity.cutoff <- sum(beta.mat!=0)
get.min.prederr <- which(prederr.val == min(prederr.val))
Results <- append(Results, list("lambda.CV" = temp$lambda.vec[get.min.prederr]))
beta.mat <- MSRL.coef(temp, lam = temp$lambda.vec[get.min.prederr])$beta

Results <- append(Results, getResults(beta = beta, Ytest = Ytest, Xtest = Xtest, 
  Ytrain = Y, Xtrain = X, SigmaX = SigmaX, SigmaYinv = SigmaYinv, beta.est = beta.mat, listname = "oursCV"))

cat("# -------------------------- ", "\n")
cat("# Ours-CV", "\n")
cat("# -------------------------- ", "\n")

# ----------------------------------------------
# Compute single TP version of our estimator
# ----------------------------------------------
lambda.selected <- temp$lambda.vec[get.min.prederr]
OursCVSingle.comptime <- system.time(
temp2 <- MSRL.cv(Y = Y, X = X, ADMM = FALSE, lambda.vec = lambda.selected, delta = .15, tol = 1e-10,
  quiet = FALSE, nfolds = NULL, inner.quiet=TRUE))

Xtemp <- (X - rep(1,n)%*%t(colMeans(X)))
Ytemp <- (Y - rep(1,n)%*%t(colMeans(Y)))
beta.cvx <- Variable(p, q)
obj <- cvxr_norm(Ytemp - Xtemp%*%beta.cvx, "nuc") + sqrt(n)*lambda.selected*p_norm(beta.cvx, 1)
prob <- Problem(Minimize(obj))
CVXSingle.comptime <- system.time(result <- solve(prob))

OursCVSingleADMM.comptime <- system.time(
temp2 <- MSRL.cv(Y = Y, X = X, ADMM = TRUE, lambda.vec = lambda.selected, delta = .15, tol = 1e-10,
  quiet = FALSE, nfolds = NULL, inner.quiet=TRUE))

Results <- append(Results, list("OursCVSingle.comptime" = OursCVSingle.comptime,
  "CVXSingle.comptime" = CVXSingle.comptime, "OursCVSingleADMM.comptime" = OursCVSingleADMM.comptime))

cat("# ------------------------------ ", "\n")
cat("OursCVSingle.comptime = ",  OursCVSingle.comptime, "\n")
cat("CVXSingle.comptime = ",  CVXSingle.comptime, "\n")
cat("# ------------------------------ ", "\n")

# -----------------------------------------------
# Tuning parameter diagnostics
# -----------------------------------------------
Xtemp <- (X - rep(1,n)%*%t(colMeans(X)))
Xtemp <- Xtemp/tcrossprod(rep(1, n), apply(Xtemp, 2, sd))
epsilon <- svd((Y - X%*%beta))
epsilontemp <- tcrossprod(epsilon$u, epsilon$v)
lambda.opt1 <- n^(-.5)*max(abs(crossprod(Xtemp, epsilontemp)))
lambda.opt <- 1.01*lambda.opt1

# --- lambda empirical -----
temp.set <- rep(0, 10000)
for(k in 1:10000){
  tempmat <- svd(matrix(rnorm(n*q), nrow=n, ncol=q))
  temp.set[k] <- n^(-.5)*max(abs(crossprod(Xtemp, tcrossprod(tempmat$u, tempmat$v))))
}
lambda.q95 <- 1.01*quantile(temp.set, .95)
lambda.q85 <- 1.01*quantile(temp.set, .85)
lambda.q75 <- 1.01*quantile(temp.set, .75)
lambda.q50 <- 1.01*quantile(temp.set, .5)

Results <- append(Results, list("lambda.q95" = lambda.q95, 
                                "lambda.q85" = lambda.q85, 
                                "lambda.q75" = lambda.q75, 
                                "lambda.q50" = lambda.q50, 
                                "lambda.opt" = lambda.opt))


# ---------------------------------------------
# q95 tuning
# ---------------------------------------------
Oursq95Single.comptime <- system.time(temp <- MSRL.cv(Y = Y, X = X, ADMM = TRUE, lambda = lambda.q95, nlambda = NULL, delta = .1, tol = 1e-12, quiet = FALSE, nfolds = NULL, inner.quiet=TRUE, standardize=TRUE)
)

beta.mat <- tcrossprod(1/temp$X.sd, rep(1, q))*matrix(temp$beta,byrow=FALSE, nrow=p, ncol=q)
Results <- append(Results, list("Oursq95Single.comptime" = Oursq95Single.comptime))
Results <- append(Results, getResults(beta = beta, Ytest = Ytest, Xtest = Xtest,
              Ytrain = Y, Xtrain = X, SigmaX = SigmaX, SigmaYinv = SigmaYinv, beta.est = beta.mat, listname = "ourslam95"))

cat("# -------------------------- ", "\n")
cat("# Ours-q95", "\n")
cat("# -------------------------- ", "\n")

# ---------------------------------------------
# q95 tuning
# ---------------------------------------------
Oursq85Single.comptime <- system.time(temp <- MSRL.cv(Y = Y, X = X, ADMM = TRUE, lambda = lambda.q85, nlambda = NULL, delta = .1, tol = 1e-12, quiet = FALSE, nfolds = NULL, inner.quiet=TRUE, standardize=TRUE)
)

beta.mat <- tcrossprod(1/temp$X.sd, rep(1, q))*matrix(temp$beta,byrow=FALSE, nrow=p, ncol=q)
Results <- append(Results, list("Oursq85Single.comptime" = Oursq85Single.comptime))
Results <- append(Results, getResults(beta = beta, Ytest = Ytest, Xtest = Xtest,
                                      Ytrain = Y, Xtrain = X, SigmaX = SigmaX, SigmaYinv = SigmaYinv, beta.est = beta.mat, listname = "ourslam85"))

cat("# -------------------------- ", "\n")
cat("# Ours-q85", "\n")
cat("# -------------------------- ", "\n")



# ---------------------------------------------
# q95 tuning
# ---------------------------------------------
Oursq75Single.comptime <- system.time(temp <- MSRL.cv(Y = Y, X = X, ADMM = TRUE, lambda = lambda.q75, nlambda = NULL, delta = .1, tol = 1e-12, quiet = FALSE, nfolds = NULL, inner.quiet=TRUE, standardize=TRUE)
)

beta.mat <- tcrossprod(1/temp$X.sd, rep(1, q))*matrix(temp$beta,byrow=FALSE, nrow=p, ncol=q)
Results <- append(Results, list("Oursq75Single.comptime" = Oursq75Single.comptime))
Results <- append(Results, getResults(beta = beta, Ytest = Ytest, Xtest = Xtest,
                                      Ytrain = Y, Xtrain = X, SigmaX = SigmaX, SigmaYinv = SigmaYinv, beta.est = beta.mat, listname = "ourslam75"))

cat("# -------------------------- ", "\n")
cat("# Ours-q75", "\n")
cat("# -------------------------- ", "\n")



# ---------------------------------------------
# q95 tuning
# ---------------------------------------------
Oursq50Single.comptime <- system.time(temp <- MSRL.cv(Y = Y, X = X, ADMM = TRUE, lambda = lambda.q50, nlambda = NULL, delta = .1, tol = 1e-12, quiet = FALSE, nfolds = NULL, inner.quiet=TRUE, standardize=TRUE)
)

beta.mat <- tcrossprod(1/temp$X.sd, rep(1, q))*matrix(temp$beta,byrow=FALSE, nrow=p, ncol=q)
Results <- append(Results, list("Oursq50Single.comptime" = Oursq50Single.comptime))
Results <- append(Results, getResults(beta = beta, Ytest = Ytest, Xtest = Xtest,
                                      Ytrain = Y, Xtrain = X, SigmaX = SigmaX, SigmaYinv = SigmaYinv, beta.est = beta.mat, listname = "ourslam50"))

cat("# -------------------------- ", "\n")
cat("# Ours-q50", "\n")
cat("# -------------------------- ", "\n")


# ---------------------------------------------
# q95-half tuning
# ---------------------------------------------
Oursq95_2Single.comptime <- system.time(temp <- MSRL.cv(Y = Y, X = X, ADMM = FALSE, lambda = lambda.q95/2, nlambda = NULL, delta = .1, tol = 1e-12, quiet = FALSE, nfolds = NULL, inner.quiet=TRUE, standardize=TRUE)
)
beta.mat <- tcrossprod(1/temp$X.sd, rep(1, q))*matrix(temp$beta, byrow=FALSE, nrow=p, ncol=q)
Results <- append(Results, list("Oursq95_2Single.comptime" = Oursq95_2Single.comptime))
Results <- append(Results, getResults(beta = beta, Ytest = Ytest, Xtest = Xtest,
  Ytrain = Y, Xtrain = X, SigmaX = SigmaX, SigmaYinv = SigmaYinv, beta.est = beta.mat, listname = "ourslam95.2"))

cat("# -------------------------- ", "\n")
cat("# Ours-q95-2", "\n")
cat("# -------------------------- ", "\n")



# ---------------------------------------------
# theoretical tuning
# ---------------------------------------------
OursOptSingle.comptime <- system.time(temp <- MSRL.cv(Y = Y, X = X, ADMM = FALSE, lambda = lambda.opt, nlambda = NULL, delta = .1, tol = 1e-12, quiet = FALSE, nfolds = NULL, inner.quiet=TRUE, standardize=TRUE)
)

beta.mat <- tcrossprod(1/temp$X.sd, rep(1, q))*matrix(temp$beta,byrow=FALSE, nrow=p, ncol=q)
Results <- append(Results, list("OursOptSingle.comptime" = OursOptSingle.comptime))
Results <- append(Results, getResults(beta = beta, Ytest = Ytest, Xtest = Xtest,
  Ytrain = Y, Xtrain = X, SigmaX = SigmaX, SigmaYinv = SigmaYinv, beta.est = beta.mat, listname = "ourslamtheory"))

cat("# -------------------------- ", "\n")
cat("# Ours-theory", "\n")
cat("# -------------------------- ", "\n")


# ---------------------------------------------------
# One-tuning parameter Lasso 
# ---------------------------------------------------
library(glmnet)
max.lam <- 0
min.lam <- Inf
for(kk in 1:q){
  store <- glmnet(y = Y[,kk], x = X)
  max.lam <- max(max.lam, max(store$lambda))
  min.lam <- min(min.lam, min(store$lambda))
  # cat(kk, "\n")
}

lasso.u.comptime <- system.time({

lambda <- seq(max.lam, min.lam, length=100)
beta.temp <- array(0, dim=c(p,q, length(lambda)))
pred.test.lasso <- rep(0, length(lambda))
pred.val.lasso <- rep(0, length(lambda))
weightedpred.lasso <- rep(0, length(lambda))

for(kk in 1:length(lambda)){
  
  pred.val.temp <- matrix(0, nrow=nval, ncol=q)
  pred.test.temp <- matrix(0, nrow=ntest, ncol=q)
  
  for(jj in 1:q){
    
    store <- glmnet(y = Y[,jj], x = X, lambda = lambda[kk], standardize=FALSE)
    beta.temp[,jj,kk] <- coef(store)[-1]
    pred.val.temp[,jj] <- predict(store, newx = Xval)
    pred.test.temp[,jj] <- predict(store, newx = Xtest)
    
  }
  
  pred.val.lasso[kk] <- sum((pred.val.temp - Yval)^2)/(nval*q)
  pred.test.lasso[kk] <- sum((pred.test.temp - Ytest)^2)/(ntest*q)
  weightedpred.lasso[kk] <- sum(diag((pred.test.temp - Ytest)%*%SigmaYinv%*%t(pred.test.temp - Ytest)))/(ntest*q)
  
}

keep <- max(which(pred.val.lasso == min(pred.val.lasso)))
beta.mat <- beta.temp[,,keep]
})

Results <- append(Results, list("lasso.u.comptime" = lasso.u.comptime))

S <- n^(-1)*crossprod(Y - rep(1, n)%*%t(colMeans(Y)) - (X - rep(1, n)%*%t(colMeans(X)))%*%beta.mat)

Results <- append(Results, getResults(beta = beta, Ytest = Ytest, Xtest = Xtest, 
  Ytrain = Y, Xtrain = X, SigmaX = SigmaX, SigmaYinv = SigmaYinv, beta.est = beta.mat, listname = "lasso.u"))

cat("# -------------------------- ", "\n")
cat("# Lasso-m", "\n")
cat("# -------------------------- ", "\n")

# -----------------------------------
# Separate Sqrt Lasso estimators
# -----------------------------------

Ypred <- matrix(0, nrow=ntest, ncol=q)
keep.inds <- rep(0, q)
lam.max <- 0
lam.min <- Inf
for(jj in 1:q){

  store <- quiet(slim(X, Y[,jj], lambda = NULL, nlambda = 20,
                method="lq", q = 2, res.sd = FALSE, lambda.min.value = .05,
                prec = 1e-4, max.ite = 1e5, verbose = FALSE))
  lam.max <- max(lam.max, max(store$lambda))
  lam.min <- min(lam.min, min(store$lambda))
  pred.val.temp <- rep(0, length(store$lambda))
  for(kk in 1:length(store$lambda)){
    Ypred.temp <- c(predict(store, newdata = Xval, lambda.idx = kk, verbose=FALSE))[[1]]
    pred.val.temp[kk] <- sum((Yval[,jj] - Ypred.temp)^2)/nval
  }

  keep.ind <- which(pred.val.temp == min(pred.val.temp))
  keep.inds[jj] <- keep.ind
  if(keep.ind == 20){
    cat("TP boundary for sqrtLasso", "\n")
  }
  beta.mat[,jj] <- store$beta[,keep.ind]
  Ypred[,jj] <- c(predict(store, newdata = Xtest,  lambda.idx = keep.ind))[[1]]
  cat(jj,"\n")

}

Results <- append(Results, getResults(beta = beta, Ytest = Ytest, Xtest = Xtest, 
  Ytrain = Y, Xtrain = X, SigmaX = SigmaX, SigmaYinv = SigmaYinv, beta.est = beta.mat, listname = "sqrtlasso.u"))

cat("# -------------------------- ", "\n")
cat("# sqrt-Lasso-m", "\n")
cat("# -------------------------- ", "\n")

# -----------------------------------
# Separate Sqrt Lasso estimators
# -----------------------------------
Calibrated.comptime <- system.time({
  Ypred <- matrix(0, nrow=ntest, ncol=q)
  lambda <- 10^seq(log(lam.max, base=10), log(lam.min, base=10), length=100)
  pred.val.temp <- matrix(0, nrow=q, ncol=100)
  Ypred.array <- array(0, dim=c(ntest, q, 100))
  beta.array <- array(0, dim=c(p, q, 100))
  for(jj in 1:q){
    
    store <- quiet(slim(X, Y[,jj], lambda = lambda,
                  method="lq", q = 2, res.sd = FALSE, 
                  prec = 1e-4, max.ite = 1e5, verbose = FALSE))
    for(kk in 1:length(store$lambda)){
      Ypred.temp <- c(predict(store, newdata = Xval, lambda.idx = kk, verbose=FALSE))[[1]]
      pred.val.temp[jj,kk] <- sum((Yval[,jj] - Ypred.temp)^2)/nval[1]
      beta.array[,jj,kk] <- store$beta[,kk]
      Ypred.array[,jj,kk] <- c(predict(store, newdata = Xtest, lambda.idx = kk, verbose=FALSE))[[1]]
    }
    
    cat(jj,"\n")
    
  }
})

Results <- append(Results, list("Calibrated.comptime" = Calibrated.comptime))


keep <- which(colSums(pred.val.temp) == min(colSums(pred.val.temp)))
beta.mat <- beta.array[,,keep]
Results <- append(Results, getResults(beta = beta, Ytest = Ytest, Xtest = Xtest, 
  Ytrain = Y, Xtrain = X, SigmaX = SigmaX, SigmaYinv = SigmaYinv, beta.est = beta.mat, listname = "sqrtlasso.m"))

cat("# -------------------------- ", "\n")
cat("# sqrt-Lasso-u", "\n")
cat("# -------------------------- ", "\n")

# ---------------------------------------------
# MRCE
# ---------------------------------------------
maxup <- max(abs(n^(-1)*crossprod(X - rep(1, n)%*%t(colMeans(X)), Y - rep(1, n)%*%t(colMeans(Y)))%*%SigmaYinv))
minup <- maxup*.005
lam2.vec <- 10^seq(log(maxup,10), log(minup,10), length=100)
beta.array <- array(0, dim=c(p,q,100))
pred.err.val <- rep(Inf, 100)
for(kk in 1:100){
  temp <- mrce(X = X, Y = Y, lam2=lam2.vec[kk],
               method=c("fixed.omega"),
               cov.tol=1e-4, cov.maxit=1e3, omega=SigmaYinv,
               maxit.out=1e3, maxit.in=1e3, tol.out=1e-8,
               tol.in=1e-10, kfold=5, silent=TRUE, eps=1e-4)
  pred.err.val[kk] <- sum(((Yval - rep(1, dim(Yval)[1])%*%t(colMeans(Y))) - (Xval - rep(1, dim(Xval)[1])%*%t(colMeans(X)))%*%temp$Bhat)^2)/(nval*q)
  if(kk > 5){
    if(pred.err.val[kk] > pred.err.val[kk-1] && pred.err.val[kk-1] > pred.err.val[kk-2] && pred.err.val[kk-2] > pred.err.val[kk-3] && kk > 30){
      break
    }
  }
  beta.array[,,kk] <- temp$Bhat
  cat(kk, "\n")
  
}

keep <- which(pred.err.val == min(pred.err.val))
beta.mat <- beta.array[,,keep]
Results <- append(Results, getResults(beta = beta, Ytest = Ytest, Xtest = Xtest, 
  Ytrain = Y, Xtrain = X, SigmaX = SigmaX, SigmaYinv = SigmaYinv, beta.est = beta.mat, listname = "MRCE.opt"))


# -------------------------------------------
# Approximate-MRCE 
# -------------------------------------------
MRCEapCV.comptime <- system.time({
  
  lambda.max <- max(abs(S - diag(diag(S))))
  lambda.min <- .005*lambda.max
  rholist <- seq(lambda.max, lambda.min, length=25)
  penMat <- matrix(1, nrow=q,ncol=q)
  diag(penMat) <- 0
  store <- QUIC(S, rho = penMat, path = rholist)
  Xtrain <- X - rep(1, n)%*%t(colMeans(X))
  XtXtrain <- crossprod(Xtrain)
  XtXeigen <- max(eigen(XtXtrain)$val)
  Ytrain <- Y - rep(1, n)%*%t(colMeans(Y))
  cat("Through glassopath", "\n")
  beta.array <- array(0, dim=c(p,q, length(rholist), 100))
  pred.err.val <- matrix(Inf, nrow=length(rholist), ncol=100)

  for(jj in 1:length(rholist)){
    maxup <- max(abs(n^(-1)*crossprod(X - rep(1, n)%*%t(colMeans(X)), Y - rep(1, n)%*%t(colMeans(Y)))%*%as.matrix(store$X[,,jj])))
    minup <- maxup*(.025/jj)
    lam2.vec <- 10^seq(log(maxup + 1e-8,10), log(minup,10), length=100)
    Qeig <- max(eigen(as.matrix(store$X[,,jj]))$val)
    beta.init <- matrix(0, nrow=p, ncol=q)
    for(kk in 1:100){
      temp <- AccProxGrad(Y = Ytrain, X = Xtrain,  Omega=as.matrix(store$X[,,jj]), weights = matrix(1, nrow=p, ncol=q),
                          betainit = beta.init, lambda = 2*lam2.vec[kk], alpha = 1, 
                          maxiter = 1e3, tol = 1e-10, XtX = XtXtrain, XtXeig = XtXeigen, Qeig = Qeig)
      
      pred.err.val[jj,kk] <- sum(((Yval - rep(1, dim(Yval)[1])%*%t(colMeans(Y))) - (Xval - rep(1, dim(Xval)[1])%*%t(colMeans(X)))%*%temp)^2)/(nval*q)
      if(kk > 30){
        if(sum(sort(pred.err.val[jj, kk:(kk-4)], decreasing = TRUE, index=TRUE)$ix == 1:5) == 5){
          break
        }
      }
      beta.init <- temp
      beta.array[,,jj,kk] <- beta.init
      cat(kk, ":", sum(temp!=0), "\n")
    }
    cat(jj, "\n")
  }
  
})

Results <- append(Results, list("MRCEapCV.comptime" = MRCEapCV.comptime))

keep <- which(pred.err.val == min(pred.err.val), arr.ind=TRUE)
beta.mat <- beta.array[,,keep[1,1], keep[1,2]]
Results <- append(Results, getResults(beta = beta, Ytest = Ytest, Xtest = Xtest, 
  Ytrain = Y, Xtrain = X, SigmaX = SigmaX, SigmaYinv = SigmaYinv, beta.est = beta.mat, listname = "MRCE.ap"))


# ---------------------------------------------
# save results
# ---------------------------------------------
saveRDS(Results,  file=savename)


q("no")