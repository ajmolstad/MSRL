
# if(Model == 1){
  
#   Sigma <- matrix(rho, nrow=q, ncol=q)
#   diag(Sigma) <- 1
#   Sigma <- diag(sqrt(seq(.5, 3, length=q)))%*%Sigma%*%diag(sqrt(seq(.5, 3, length=q)))
#   eo <- eigen(Sigma)
  
#   SigmaYinv <- eo$vec%*%diag(eo$val^-1)%*%t(eo$vec)
#   SigmaYinv[which(abs(SigmaYinv) < 1e-10)] <- 0
#   Y <- X%*%beta + matrix(rnorm(n*q), nrow=n)%*%eo$vec%*%diag(eo$val^.5)%*%t(eo$vec)
#   Ytest <- Xtest%*%beta + matrix(rnorm(ntest*q), nrow=ntest)%*%eo$vec%*%diag(eo$val^.5)%*%t(eo$vec)
#   Yval <- Xval%*%beta + matrix(rnorm(nval*q), nrow=nval)%*%eo$vec%*%diag(eo$val^.5)%*%t(eo$vec)

# }  


# if(Model == 2){
  
#   Sigma <- matrix(rho, nrow=q, ncol=q)
#   diag(Sigma) <- 1
#   Sigma <- diag(sqrt(seq(.5, 3, length=q)))%*%Sigma%*%diag(sqrt(seq(.5, 3, length=q)))
#   eo <- eigen(Sigma)
  
#   SigmaYinv <- eo$vec%*%diag(eo$val^-1)%*%t(eo$vec)
#   SigmaYinv[which(abs(SigmaYinv) < 1e-10)] <- 0
#   Y <- X%*%beta + matrix(rt(n*q, df = 3), nrow=n)%*%eo$vec%*%diag(eo$val^.5)%*%t(eo$vec)
#   Ytest <- Xtest%*%beta + matrix(rt(ntest*q, df = 3), nrow=ntest)%*%eo$vec%*%diag(eo$val^.5)%*%t(eo$vec)
#   Yval <- Xval%*%beta + matrix(rt(nval*q, df = 3), nrow=nval)%*%eo$vec%*%diag(eo$val^.5)%*%t(eo$vec)
  
# }




if(Model == 1){
  
  Sigma <- matrix(rho, nrow=q, ncol=q)
  diag(Sigma) <- 1
  Sigma <- 3*Sigma
  eo <- eigen(Sigma)
  
  SigmaYinv <- eo$vec%*%diag(eo$val^-1)%*%t(eo$vec)
  SigmaYinv[which(abs(SigmaYinv) < 1e-10)] <- 0
  Y <- X%*%beta + matrix(rnorm(n*q), nrow=n)%*%eo$vec%*%diag(eo$val^.5)%*%t(eo$vec)
  Ytest <- Xtest%*%beta + matrix(rnorm(ntest*q), nrow=ntest)%*%eo$vec%*%diag(eo$val^.5)%*%t(eo$vec)
  Yval <- Xval%*%beta + matrix(rnorm(nval*q), nrow=nval)%*%eo$vec%*%diag(eo$val^.5)%*%t(eo$vec)

}  


if(Model == 2){
  
  Sigma <- matrix(rho, nrow=q, ncol=q)
  diag(Sigma) <- 1
  Sigma <- 3*Sigma
  eo <- eigen(Sigma)
  
  SigmaYinv <- eo$vec%*%diag(eo$val^-1)%*%t(eo$vec)
  SigmaYinv[which(abs(SigmaYinv) < 1e-10)] <- 0
  Y <- X%*%beta + matrix(rt(n*q, df = 3), nrow=n)%*%eo$vec%*%diag(eo$val^.5)%*%t(eo$vec)
  Ytest <- Xtest%*%beta + matrix(rt(ntest*q, df = 3), nrow=ntest)%*%eo$vec%*%diag(eo$val^.5)%*%t(eo$vec)
  Yval <- Xval%*%beta + matrix(rt(nval*q, df = 3), nrow=nval)%*%eo$vec%*%diag(eo$val^.5)%*%t(eo$vec)
  
}



# if(Model == 3){
  

#    V <- matrix(rnorm(q*rho), nrow=q)
#   V <- V*(apply(V, 1, function(x){1/sqrt(sum(x^2))})%*%t(rep(1, rho)))*(sqrt(seq(.05, 2.45, length=q))%*%t(rep(1, rho)))
#   Sigma <- tcrossprod(V) + diag(.45, q)
#   eo <- eigen(Sigma)


#   V <- matrix(rnorm(q*rho), nrow=q)
#   V <- V*(apply(V, 1, function(x){1/sqrt(sum(x^2))})%*%t(rep(1, rho)))*sqrt(2.5)
#   Sigma <- tcrossprod(V) + diag(0.5, q)
#   eo <- eigen(Sigma)
  
#   SigmaYinv <- eo$vec%*%diag(eo$val^-1)%*%t(eo$vec)
#   SigmaYinv[which(abs(SigmaYinv) < 1e-10)] <- 0
#   Y <- X%*%beta + matrix(rnorm(n*q), nrow=n)%*%eo$vec%*%diag(eo$val^.5)%*%t(eo$vec)
#   Ytest <- Xtest%*%beta + matrix(rnorm(ntest*q), nrow=ntest)%*%eo$vec%*%diag(eo$val^.5)%*%t(eo$vec)
#   Yval <- Xval%*%beta + matrix(rnorm(nval*q), nrow=nval)%*%eo$vec%*%diag(eo$val^.5)%*%t(eo$vec)
  
# }


# if(Model == 4){
  
#   V <- matrix(rnorm(q*rho), nrow=q)
#   V <- V*(apply(V, 1, function(x){1/sqrt(sum(x^2))})%*%t(rep(1, rho)))*(sqrt(seq(.05, 2.45, length=q))%*%t(rep(1, rho)))
#   Sigma <- tcrossprod(V) + diag(.45, q)
#   eo <- eigen(Sigma)
  
#   SigmaYinv <- eo$vec%*%diag(eo$val^-1)%*%t(eo$vec)
#   SigmaYinv[which(abs(SigmaYinv) < 1e-10)] <- 0
#   Y <- X%*%beta + matrix(rt(n*q, df = 3), nrow=n)%*%eo$vec%*%diag(eo$val^.5)%*%t(eo$vec)
#   Ytest <- Xtest%*%beta + matrix(rt(ntest*q, df = 3), nrow=ntest)%*%eo$vec%*%diag(eo$val^.5)%*%t(eo$vec)
#   Yval <- Xval%*%beta + matrix(rt(nval*q, df = 3), nrow=nval)%*%eo$vec%*%diag(eo$val^.5)%*%t(eo$vec)
  
# }

if(Model == 3){
  
  V <- matrix(rnorm(q*rho), nrow=q)
  V <- V*(apply(V, 1, function(x){1/sqrt(sum(x^2))})%*%t(rep(1, rho)))*sqrt(1.45)
  Sigma <- tcrossprod(V) + diag(0.05, q)
  eo <- eigen(Sigma)
  
  SigmaYinv <- eo$vec%*%diag(eo$val^-1)%*%t(eo$vec)
  SigmaYinv[which(abs(SigmaYinv) < 1e-10)] <- 0
  Y <- X%*%beta + matrix(rnorm(n*q), nrow=n)%*%eo$vec%*%diag(eo$val^.5)%*%t(eo$vec)
  Ytest <- Xtest%*%beta + matrix(rnorm(ntest*q), nrow=ntest)%*%eo$vec%*%diag(eo$val^.5)%*%t(eo$vec)
  Yval <- Xval%*%beta + matrix(rnorm(nval*q), nrow=nval)%*%eo$vec%*%diag(eo$val^.5)%*%t(eo$vec)
  
}

if(Model == 4){
  
  V <- matrix(rnorm(q*rho), nrow=q)
  V <- V*(apply(V, 1, function(x){1/sqrt(sum(x^2))})%*%t(rep(1, rho)))*sqrt(1.45)
  Sigma <- tcrossprod(V) + diag(0.05, q)
  eo <- eigen(Sigma)
  
  SigmaYinv <- eo$vec%*%diag(eo$val^-1)%*%t(eo$vec)
  SigmaYinv[which(abs(SigmaYinv) < 1e-10)] <- 0
  Y <- X%*%beta + matrix(rt(n*q, df = 3), nrow=n)%*%eo$vec%*%diag(eo$val^.5)%*%t(eo$vec)
  Ytest <- Xtest%*%beta + matrix(rt(ntest*q, df = 3), nrow=ntest)%*%eo$vec%*%diag(eo$val^.5)%*%t(eo$vec)
  Yval <- Xval%*%beta + matrix(rt(nval*q, df = 3), nrow=nval)%*%eo$vec%*%diag(eo$val^.5)%*%t(eo$vec)
  
}  

# if(Model == 5){
  
#   l <- c(2, 5, 10, 25, 50, 100)
#   eigVec <- svd(matrix(rnorm(q^2), nrow=q))$u
#   eigVals <- seq(1/l[rho], 1, length=q)
#   Cor <- eigVec%*%diag(eigVals)%*%t(eigVec)
#   Sigma <- diag(sqrt(seq(.5, 3, length=q)))%*%Cor%*%diag(sqrt(seq(.5, 3, length=q)))
#   eo <- eigen(Sigma)
  
#   SigmaYinv <- eo$vec%*%diag(eo$val^-1)%*%t(eo$vec)
#   SigmaYinv[which(abs(SigmaYinv) < 1e-10)] <- 0
#   Y <- X%*%beta + matrix(rnorm(n*q), nrow=n)%*%eo$vec%*%diag(eo$val^.5)%*%t(eo$vec)
#   Ytest <- Xtest%*%beta + matrix(rnorm(ntest*q), nrow=ntest)%*%eo$vec%*%diag(eo$val^.5)%*%t(eo$vec)
#   Yval <- Xval%*%beta + matrix(rnorm(nval*q), nrow=nval)%*%eo$vec%*%diag(eo$val^.5)%*%t(eo$vec)
  
# }

# if(Model == 6){
  
  # l <- c(2, 5, 10, 25, 50, 100)
  # eigVec <- svd(matrix(rnorm(q^2), nrow=q))$u
  # eigVals <- seq(1/l[rho], 1, length=q)
  # Cor <- eigVec%*%diag(eigVals)%*%t(eigVec)
  # Sigma <- diag(sqrt(seq(.5, 3, length=q)))%*%Cor%*%diag(sqrt(seq(.5, 3, length=q)))
  # eo <- eigen(Sigma)
  
#   SigmaYinv <- eo$vec%*%diag(eo$val^-1)%*%t(eo$vec)
#   SigmaYinv[which(abs(SigmaYinv) < 1e-10)] <- 0
#   Y <- X%*%beta + matrix(rt(n*q, df = 3), nrow=n)%*%eo$vec%*%diag(eo$val^.5)%*%t(eo$vec)
#   Ytest <- Xtest%*%beta + matrix(rt(ntest*q, df = 3), nrow=ntest)%*%eo$vec%*%diag(eo$val^.5)%*%t(eo$vec)
#   Yval <- Xval%*%beta + matrix(rt(nval*q, df = 3), nrow=nval)%*%eo$vec%*%diag(eo$val^.5)%*%t(eo$vec)
  
# }

if(Model == 5){
  
  l <- c(2, 5, 10, 25, 50, 100)
  eigVec <- svd(matrix(rnorm(q^2), nrow=q))$u
  eigVals <- seq(1/l[rho], 1, length=q)
  Cor <- eigVec%*%diag(eigVals)%*%t(eigVec)
  Sigma <- 2*Cor
  eo <- eigen(Sigma)
  
  SigmaYinv <- eo$vec%*%diag(eo$val^-1)%*%t(eo$vec)
  SigmaYinv[which(abs(SigmaYinv) < 1e-10)] <- 0
  Y <- X%*%beta + matrix(rnorm(n*q), nrow=n)%*%eo$vec%*%diag(eo$val^.5)%*%t(eo$vec)
  Ytest <- Xtest%*%beta + matrix(rnorm(ntest*q), nrow=ntest)%*%eo$vec%*%diag(eo$val^.5)%*%t(eo$vec)
  Yval <- Xval%*%beta + matrix(rnorm(nval*q), nrow=nval)%*%eo$vec%*%diag(eo$val^.5)%*%t(eo$vec)
  
}

if(Model == 6){
  
  l <- c(2, 5, 10, 25, 50, 100)
  eigVec <- svd(matrix(rnorm(q^2), nrow=q))$u
  eigVals <- seq(1/l[rho], 1, length=q)
  Cor <- eigVec%*%diag(eigVals)%*%t(eigVec)
  Sigma <- 2*Cor
  eo <- eigen(Sigma)
  
  SigmaYinv <- eo$vec%*%diag(eo$val^-1)%*%t(eo$vec)
  SigmaYinv[which(abs(SigmaYinv) < 1e-10)] <- 0
  Y <- X%*%beta + matrix(rt(n*q, df = 3), nrow=n)%*%eo$vec%*%diag(eo$val^.5)%*%t(eo$vec)
  Ytest <- Xtest%*%beta + matrix(rt(ntest*q, df = 3), nrow=ntest)%*%eo$vec%*%diag(eo$val^.5)%*%t(eo$vec)
  Yval <- Xval%*%beta + matrix(rt(nval*q, df = 3), nrow=nval)%*%eo$vec%*%diag(eo$val^.5)%*%t(eo$vec)
  
}

