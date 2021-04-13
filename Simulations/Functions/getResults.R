
getResults <- function(beta, Ytest, Xtest, Ytrain, Xtrain, SigmaX, SigmaYinv, beta.est, listname){
  
  tpr <- TPR(beta, beta.est)
  fpr <- FPR(beta, beta.est)
  mse <- sum((beta - beta.est)^2)
  modelerr <- sum(diag(t(beta - beta.est)%*%SigmaX%*%(beta - beta.est)))
  Ypred <- (Xtest - tcrossprod(rep(1, dim(Xtest)[1]), colMeans(Xtrain)))%*%beta.est
  residuals <- (Ytest - tcrossprod(rep(1, dim(Ytest)[1]), colMeans(Ytrain))) - Ypred
  prederr <- mean(residuals^2)
  weightedprederr <- sum(diag(residuals%*%SigmaYinv%*%t(residuals)))/(length(Ytest))
  
  # -------------------------------------
  # Refit OLS
  # -------------------------------------
  beta.refit <- matrix(0, nrow=p, ncol=q)
  for(k in 1:q){
    if(sum(beta.est[,k]!=0)>0){
      temp <- lm(Ytrain[,k] ~ Xtrain[,which(beta.est[,k]!=0)])
      beta.refit[which(beta.est[,k]!=0),k] <- coefficients(temp)[-1]
    }
  }
  mse.RFOLS <- sum((beta - beta.refit)^2)
  modelerr.RFOLS <- sum(diag(t(beta - beta.refit)%*%SigmaX%*%(beta - beta.refit)))
  Ypred <- (Xtest - tcrossprod(rep(1, dim(Xtest)[1]), colMeans(Xtrain)))%*%beta.refit
  residuals <- (Ytest - tcrossprod(rep(1, dim(Ytest)[1]), colMeans(Ytrain))) - Ypred
  prederr.RFOLS <- mean(residuals^2)
  weightedprederr.RFOLS <- sum(diag(residuals%*%SigmaYinv%*%t(residuals)))/(length(Ytest))
  
  # -------------------------------------
  # Refit SUR
  # -------------------------------------
  sparsity.mat <- matrix(1, nrow=p, ncol=q)
  sparsity.mat[which(beta.est==0)] <- 0
  beta.SUR <- SURfit(Y = Ytrain - rep(1, n)%*%t(colMeans(Ytrain)), X = Xtrain - rep(1, n)%*%t(colMeans(Xtrain)), sparsity.mat = sparsity.mat)
  mse.RFSUR <- sum((beta - beta.SUR)^2)
  modelerr.RFSUR <- sum(diag(t(beta - beta.SUR)%*%SigmaX%*%(beta - beta.SUR)))
  Ypred <- (Xtest - tcrossprod(rep(1, dim(Xtest)[1]), colMeans(Xtrain)))%*%beta.SUR
  residuals <- (Ytest - tcrossprod(rep(1, dim(Ytest)[1]), colMeans(Ytrain))) - Ypred
  prederr.RFSUR <- mean(residuals^2)
  weightedprederr.RFSUR <- sum(diag(residuals%*%SigmaYinv%*%t(residuals)))/(length(Ytest))
  
  temp <- named.list(tpr, fpr, mse, modelerr, prederr, weightedprederr,
                     mse.RFOLS, modelerr.RFOLS, prederr.RFOLS, weightedprederr.RFOLS,
                     mse.RFSUR, modelerr.RFSUR, prederr.RFSUR, weightedprederr.RFSUR
  )
  
  names(temp) <- paste(listname, ".", names(temp), sep="")
  return(temp)
  
}