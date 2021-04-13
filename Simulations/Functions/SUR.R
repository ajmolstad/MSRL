# ---------------------------------------------
# Write a SUR-style blockwise descent
# ---------------------------------------------

ProjGradSUR <- function(X, Y, beta, Omega, sparsity.mat, tol = 1e-8){
  
  n <- dim(X)[1]
  gradFunc <- function(X, Y, Omega, beta, n){
    temp <- 2*n^(-1)*crossprod(X, Y - X%*%beta)%*%Omega
    return(temp)
  }
  ObjFunc <- function(X, Y, Omega, beta, n){
    XB <- crossprod(t(X), beta)
    return((1/n)*sum(diag(crossprod(t(Y - XB), tcrossprod(Omega, (Y - XB))))))
  }
  L <- (4/n)*max(eigen(Omega)$val)*max(eigen(crossprod(X))$val)
  orig.obj <- abs(ObjFunc(X, Y, Omega, beta, n))
  old.obj <- ObjFunc(X, Y, Omega, beta, n)
    
    
  for(k in 1:100){
      update <- beta + (1/L)*gradFunc(X, Y, Omega, beta, n)
      update[which(sparsity.mat == 0)] <- 0
      new.obj <- ObjFunc(X, Y, Omega, update, n)
      beta <- update
      if(abs(old.obj - new.obj) < tol*orig.obj){
        break
      }
      old.obj <- new.obj
      #cat(new.obj, "\n")
  }
  
  return(beta)
}

SURfit <- function(X, Y, sparsity.mat){
  
  n <- dim(X)[1]
  p <- dim(X)[2]
  q <- dim(Y)[2]
  Omega <- diag(1, q)
  beta <- matrix(0, nrow=p, ncol=q)
  
  out.obj <- function(Y, X, beta, Omega){
    n <- dim(X)[1]
    XB <- crossprod(t(X), beta)
    return((1/n)*sum(diag(crossprod(t(Y - XB), tcrossprod(Omega, (Y - XB))))) - .5*determinant(Omega, logarithm=TRUE)$modulus[1]) + (1e-4/2)*sum(Omega^2)
  }
  
  old.obj <- out.obj(Y, X, beta, Omega)
  orig.obj <- abs(old.obj)
  
  for(k in 1:200){
    betaUp <- ProjGradSUR(X, Y, beta, Omega, sparsity.mat, tol = 1e-8)
    S <- (n)^(-1)*crossprod(Y - X%*%betaUp)
    eo <- eigen(S)
    OmegaUp <- (1/(2*1e-4))*(eo$vec%*%(diag(-eo$val + sqrt(eo$val^2 + 4*1e-4)))%*%t(eo$vec))
    new.obj <- out.obj(Y, X, betaUp, OmegaUp)
    beta <- betaUp
    Omega <- OmegaUp
    cat(new.obj, "\n")
    if(abs(new.obj - old.obj) < 1e-6*orig.obj){
      break
    }
    old.obj <- new.obj
  }
  
  return(beta)
}


