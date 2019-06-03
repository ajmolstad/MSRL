# -----------------------------------------------------
# Complete set of functions to fit NN estimator
# Please contact amolstad@fredhutch.org with issues  
# --------------------------------------------------------

AccPG <- function(y, x, beta, lam1, weight, tol, maxiter = 1e4, quiet = FALSE){
	    
	    betakm1 <- matrix(0, nrow=p, ncol=q)
	    betak <- matrix(0, nrow=p, ncol=q)
	    alphak <- 1
	    alphakm1 <- 1
	    k.iter <- 1
	    t <- 10
	    max.iter <- maxiter
	    obj <- rep(0,max.iter)
	    inner.temp <- svd(y - crossprod(t(x), betak))
	    h.old <-(1/sqrt(n))*sum(inner.temp$d) + lam1*sum(abs(weight*betak))
	    h.orig <- h.old
	    outeriter <- TRUE
	    wlam <- weight*lam1

	    while(outeriter){
	      
		    # -----------------------------
		    # Step 1
		    # -----------------------------
		    theta <- betak + ((alphak - 1)/alphakm1)*(betak - betakm1)

			# -----------------------------
		    # Step 2
		    # -----------------------------
		    temp <-  svd(y - crossprod(t(x), theta))
		    inneriter <- TRUE
		    t <- 1
	      
	        while(inneriter){
	        
	        	keep <- crossprod(x, tcrossprod(temp$u, temp$v))
	        	keep1 <- (t/sqrt(n))*keep + theta
	        	betakp1 <- pmax(abs(keep1) - weight*t*lam1, 0)*sign(keep1)
	        	inner.temp <- svd(y - crossprod(t(x), betakp1), nu = min(n,q), nv = min(n,q))
	        
	        	h <- (1/(sqrt(n)))*sum(inner.temp$d) + lam1*sum(abs(weight*betakp1))
	        	g <- (1/(sqrt(n)))*(sum(temp$d) - sum(keep*(betakp1 - theta))) + (1/(2*t))*sum((betakp1 - theta)^2) + lam1*sum(abs(weight*betakp1))
	        
	        	if(h < g){
	          		inneriter <- FALSE
	        	} else {
	          		t <- t/2
	        	}

	        	if(t < 1e-10){
	          		inneriter <- FALSE
	          		outeriter <- FALSE
	        	}
	        
	      	}
	      
		    if(h < h.old){
		    	h.old <- h
		    } else {
		    	betakp1 <- betak
		    }
	      
	      	inner.temp <- svd(y - crossprod(t(x), betakp1))
	      	betakm1 <- betak
	      	betak <- betakp1
	      	if(sum(inner.temp$d < .001) > 0){
	        	stop("Error!!!!")
	      	}
	      
	      	alphakm1 <- alphak 
	      	alphak <- (1 + sqrt(1 + 4*alphakm1^2))/2
	      
	      	obj[k.iter] <- h.old
	      
	      	if(k.iter > 5){
	        
	        	# --------------------------------
	        	# check convergence 
	        	# -------------------------------
	        	keep <- (1/sqrt(n))*crossprod(x, tcrossprod(inner.temp$u, inner.temp$v))
	        	t1 <- sum((abs(keep) > wlam)[which(betak == 0)])/sum(betak == 0)

	        	if(any(betak!=0)){
	          		t2 <- max(abs(keep - wlam*sign(betak))[which(betak != 0)])
	        	} else {
	          		t2 <- 0
	        	}

	        	# --------------------------------
	        	# print results if not quiet 
	        	# -------------------------------
	        	if(!quiet){
			        if(k.iter %% 50 == 0){
			          cat(k.iter, ": t1 = ", t1, "t2 = ", t2, "\n")
			        }
			    }

		        if(t1 == 0 && t2 < tol){
		          outeriter <- FALSE
		        }

		        if(obj[k.iter - 3] - obj[k.iter] < tol/abs(h.orig)){
		          outeriter <- FALSE
		        }

	      	}
	      
	      	k.iter <- k.iter + 1
	      
	      	if(k.iter > max.iter){
	        	outeriter <- FALSE
	      	}

	    }
	    
	    return(list("beta" = betakp1))

	}

 NN_ADMM <- function(Y, X, Gamma, Omega, beta, lambda, weight, tol, maxiter, rho, eta, quiet = FALSE){
      
      	PenProxP1 <- function(input, tau){
        	pmax(abs(input) - tau, 0)*sign(input)
      	}
      
	    n <- dim(Y)[1]
		p <- dim(X)[2]
		q <- dim(Y)[2]
		Gammakm1 <- Gamma
		Gammak <- Gamma
		Omegakm1 <- Omega
		Omegak <- Omega
		betakm1 <- beta
		Xbetakm1 <- crossprod(t(X), betakm1)
		iterating <- TRUE
		lam <- lambda*sqrt(n)
		k.iter <- 1
		wlam <- weight*lam
      
		while(iterating){

			# ---- Omega update 
			temp <- svd(Y + rho^(-1)*Gammakm1 - Xbetakm1, nu = min(n,q), nv = min(n,q))
			Omegak <- tcrossprod(temp$u*tcrossprod(rep(1, dim(temp$u)[1]), pmax(temp$d - (1/(rho)), 0)), temp$v)

			# -- beta update 
			H <- crossprod(X, Y + rho^(-1)*Gammakm1 - Omegak - Xbetakm1)
			betak <- PenProxP1(betakm1 + H/eta, wlam/(rho*eta))
			Xbetakm1 <- crossprod(t(X),  betak)

			# --- Gamma udpate
			Gammak <- Gammakm1 + ((1 + sqrt(5))/2)*rho*(Y - Xbetakm1 - Omegak)

			r1 <- sum((Y - Xbetakm1 - Omegak)^2); 
			s1 <- sum(crossprod(rho*X, Gammak - Gammakm1)^2)

			if(k.iter%%500 == 0){
				cat("r1 = ", r1, "; s1 = ", s1, "\n")
			}

			if(r1 < sqrt(n)*tol & s1 < tol){
				iterating = FALSE;
			}
		# update step size
		if(k.iter%%10 == 0){
			if(r1 > 10*s1){
				rho = rho*2;
			}

			if(s1 > 10*r1){
				rho = rho/2;
			}
		}

		Gammakm1 <- Gammak 
		Omegakm1 <- Omegak
		betakm1 <- betak
		k.iter <- k.iter + 1

		}

      	return(list("beta" = betak, "Gamma" = Gammak, "Omega" = Omegak))
    
    }
    


MSRL.cv <- function(X, Y, nlambda, lambda.vec = NULL,
		weighted = FALSE, standardize = FALSE, ADMM = FALSE, nfolds = NULL, 
		delta = .01, tol = 1e-8, quiet = TRUE, inner.quiet=TRUE){

	    
    # ---------------------------------
	# preliminaries
	# ---------------------------------
	n <- dim(X)[1]
	p <- dim(X)[2]
	q <- dim(Y)[2]
    ADMM.temp <- ADMM

	# --------------------------------
	# construct weight matrices
	# --------------------------------
	if(weighted){
		weight.mat <- rep(1, dim(X)[2])%*%t(apply(Y, 2, sd))
		weight <- weight.mat/max(weight.mat)
	} else {
		weight <- matrix(1, nrow=p, ncol=q)
	}

	# -----------------------------------------
	# standardize if necessary 
	# -----------------------------------------
	if(!standardize){

		x <- X - rep(1, n)%*%t(apply(X, 2, mean))
		xtx <- crossprod(x)
		y <- Y - rep(1, n)%*%t(apply(Y, 2, mean))
		xty <- crossprod(x, y)
		
	} else {
	  
		x <- (X - rep(1, n)%*%t(apply(X, 2, mean)))/(rep(1, n)%*%t(apply(X, 2, sd)))
		xtx <- crossprod(x)
		y <- (Y - rep(1, n)%*%t(apply(Y, 2, mean)))
		xty <- crossprod(x, y)

	}

	# -------------------------------------------
	# get candidate tuning parameters
	# -------------------------------------------
	if(is.null(lambda.vec)){

		lambda.vec <- rep(0, length=nlambda)
		temp <- svd(y)
		if(sum(weight==0) == 0){
			lambda.max <- 1.01*(1/sqrt(n))*max(abs(crossprod(x, temp$u%*%t(temp$v)))*weight^(-1))
		} else {
			lambda.max <- 1.01*(1/sqrt(n))*max(abs(crossprod(x, temp$u%*%t(temp$v))))
		}

		lambda.min <- delta*lambda.max
		for(kk in 1:nlambda){
			lambda.vec[kk] <- lambda.max^((nlambda-kk)/(nlambda-1))*lambda.min^((kk-1)/(nlambda-1))
		}
	}

	# ------------------------------------------------------------------------	
	# compute solution path 
	# ------------------------------------------------------------------------
	beta.full <- Matrix(0, nrow = p*q, ncol = length(lambda.vec), sparse=TRUE)
	sparsity.mat <- rep(0, length(lambda.vec))
	beta.old <- matrix(0, nrow=p, ncol=q)
	Gamma <- y
	temp <- svd(y)
	Omega <- tcrossprod(temp$u, temp$v)
  	xtxeig <- max(eigen(xtx)$val)

	for(kk in 1:length(lambda.vec)){

		# -------------------------------------------------
		# Compute using PGD if appropriate
		# -------------------------------------------------
		if(!ADMM.temp){
			
			temp <- try(AccPG(y = y, x = x, beta = beta.old, lam1 = lambda.vec[kk], weight = weight, 
  			tol = tol, maxiter = 1e4, quiet=inner.quiet))

  			if(class(temp) == "try-error"){

				ADMM.temp <- TRUE
				temp <- NN_ADMM(Y = y, X = x, Gamma = Gamma, Omega = Omega, 
				              beta = beta.old, lambda = lambda.vec[kk], weight = weight, 
				              tol = tol, maxiter = 1e4, rho = .0001, eta = 1.00001*xtxeig, quiet=inner.quiet)

				Gamma <- temp$Gamma
				Omega <- temp$Omega
				beta.old <- temp$beta
  			  	beta.full[,kk] <- c(beta.old)
  			  	sparsity.mat[kk] <- sum(beta.old!=0)
  			  
  			} else { 
  			  
  			  	beta.old <- temp$beta
  			  	beta.full[,kk] <- c(beta.old)
  			  	sparsity.mat[kk] <- sum(beta.old!=0)
  			  	if(!quiet){
  			  		cat(kk, ": non-zero = ", sum(beta.old!=0), "\n")
  			  		cat("# ------------------------------ ", "\n")
  			  	}
  			  
  			}
			
		} else {
		  
		  	temp <- NN_ADMM(Y = y, X = x, Gamma = Gamma, Omega = Omega, 
	          	beta = beta.old, lambda = lambda.vec[kk], weight = weight, 
	        	tol = tol, maxiter = 1e4, rho =  .0001, eta = 1.00001*xtxeig, quiet=inner.quiet)

		  	Gamma <- temp$Gamma
		  	Omega <- temp$Omega
          	beta.old <- temp$beta
    	  	beta.full[,kk] <- c(beta.old)
        	sparsity.mat[kk] <- sum(beta.old!=0)
        	if(!quiet){
       			cat(kk, ": non-zero = ", sum(beta.old!=0), "\n")
       			cat("# ------------------------------ ", "\n")
			}
		}
	}
		
	beta.full <- beta.full[,1:kk]
	sparsity.mat <- sparsity.mat[1:kk]
	lambda.vec <- lambda.vec[1:kk]

	# -----------------------------------------
	# Perform cross-validation 
	# -----------------------------------------
	if(!is.null(nfolds)){

		# -----------------------------------------
		# save metrics for cross-validation 
		# -----------------------------------------
		fold <- sample(rep(1:nfolds, length=n))
		cv.index <- split(1:n, fold)
		errs_wpred <- matrix(0, nrow=length(lambda.vec), ncol=nfolds)
		errs_pred <- matrix(0, nrow=length(lambda.vec), ncol=nfolds)
		errs_spec <- matrix(0, nrow=length(lambda.vec), ncol=nfolds)
		errs_nuc <- matrix(0, nrow=length(lambda.vec), ncol=nfolds)

		for(k in 1:nfolds){

			if(!standardize){
				
				# --- center X and Y
				ntrain <- dim(X[-cv.index[[k]],])[1]
				x.inner <- X[-cv.index[[k]], ] - rep(1, ntrain)%*%t(apply(X[-cv.index[[k]],], 2, mean))
				xtx.inner <- crossprod(x.inner)
				y.inner <- Y[-cv.index[[k]], ] - rep(1, ntrain)%*%t(apply(Y[-cv.index[[k]],], 2, mean))
				var.y <- apply(y.inner, 2, var)
				xty.inner <- crossprod(x.inner, y.inner)

				n.test <- length(cv.index[[k]])
				x.test <- X[cv.index[[k]], ] - rep(1, n.test)%*%t(apply(X[-cv.index[[k]],], 2, mean))
				y.test <- Y[cv.index[[k]], ] - rep(1, n.test)%*%t(apply(Y[-cv.index[[k]],], 2, mean))
				Gamma <- matrix(0, nrow=ntrain, ncol=q)
				Omega <- matrix(0, nrow=ntrain, ncol=q)
				beta.old <- matrix(0, nrow=p, ncol=q)
				xtxeig <- max(eigen(xtx.inner)$val)

				if(weighted){
			    	weight.mat <- rep(1, dim(X)[2])%*%t(apply(y.inner, 2, mad))
			    	weight <- weight.mat/max(weight.mat)
			    } else {
			    	weight <- matrix(1, nrow=p, ncol=q)
			    }
			
			} else {	

				# --- center X and Y --------------------
				ntrain <- dim(X[-cv.index[[k]],])[1]
				x.inner <- (X[-cv.index[[k]], ] - rep(1, ntrain)%*%t(apply(X[-cv.index[[k]],], 2, mean)))/(rep(1, ntrain)%*%t(apply(X[-cv.index[[k]], ], 2, sd)))
				xtx.inner <- crossprod(x.inner)
				y.inner <- (Y[-cv.index[[k]], ] - rep(1, ntrain)%*%t(apply(Y[-cv.index[[k]],], 2, mean)))
				var.y <- apply(y.inner, 2, var)
				xty.inner <- crossprod(x.inner, y.inner)

				# ------------------------------------------------
				n.test <- length(cv.index[[k]])
				x.test <- (X[cv.index[[k]], ] - rep(1, n.test)%*%t(apply(X[-cv.index[[k]],], 2, mean)))/(rep(1, n.test)%*%t(apply(X[-cv.index[[k]], ], 2, sd)))
				y.test <- (Y[cv.index[[k]], ] - rep(1, n.test)%*%t(apply(Y[-cv.index[[k]],], 2, mean)))
				Gamma <- matrix(0, nrow=ntrain, ncol=q)
				Omega <- matrix(0, nrow=ntrain, ncol=q)
				beta.old <- matrix(0, nrow=p, ncol=q)
				xtxeig <- max(eigen(xtx.inner)$val)

				if(weighted){
			    	weight.mat <- rep(1, dim(X)[2])%*%t(apply(y.inner, 2, sd))
			    	weight <- weight.mat/max(weight.mat)
			    } else {
			    	weight <- matrix(1, nrow=p, ncol=q)
			    }


			}

			Gamma <- y.inner
			temp <- svd(y.inner)
			Omega <- tcrossprod(temp$u, temp$v)

			# ------------------------------------------------
			for(kk in 1:length(lambda.vec)){


			    # -------------------------------------------------
				# Compute using PGD if appropriate
				# -------------------------------------------------
				if(!ADMM){
					
					temp <- try(AccPG(y = y.inner, x = x.inner, beta = beta.old, lam1 = lambda.vec[kk], weight = weight, 
		  			tol = tol, maxiter = 1e4, quiet=inner.quiet))

		  			if(class(temp) == "try-error"){

						ADMM <- TRUE
						temp <- NN_ADMM(Y = y.inner, X = x.inner, Gamma = Gamma, Omega = Omega, 
						              beta = beta.old, lambda = lambda.vec[kk], weight = weight, 
						              tol = tol, maxiter = 1e4, rho = .0001, eta = 1.00001*xtxeig, quiet=inner.quiet)

						Gamma <- temp$Gamma
						Omega <- temp$Omega
						beta.old <- temp$beta
		  			  
		  			} else { 
		  			  
		  			  	beta.old <- temp$beta
		  			  
		  			}
					
				} else {
				  
				  	temp <- NN_ADMM(Y = y.inner, X = x.inner, Gamma = Gamma, Omega = Omega, 
			          	beta = beta.old, lambda = lambda.vec[kk], weight = weight, 
			        	tol = tol, maxiter = 1e4, rho =  .0001, eta = 1.00001*xtxeig, quiet=inner.quiet)

				  	Gamma <- temp$Gamma
				  	Omega <- temp$Omega
		          	beta.old <- temp$beta

				}

				residual <- (y.test - x.test%*%temp$beta)
				errs_pred[kk, k] <- sum((residual)^2)
	   			errs_wpred[kk, k] <- sum(diag((residual)%*%diag(1/apply(Y[-cv.index[[k]], ], 2, var))%*%t(residual)))
	   			inner.temp <- svd(residual)$d
	   			errs_spec[kk, k] <- max(inner.temp)
	   			errs_nuc[kk, k] <- sum(inner.temp)
	   			# cat(kk, "\n")
			}
		  
		 	if(!quiet){
		 		cat("# ------------------------------ ", "\n")
		 		cat("Through CV fold", k, "\n")
		 		cat("# ------------------------------ ", "\n")

		 	}
		}

	} else {
		errs_pred <- NULL
		errs_wpred <- NULL
		errs_spec <- NULL
		errs_nuc <- NULL
	}

	if(!is.null(errs_wpred)){

		inds <- which(rowSums(errs_wpred) == min(rowSums(errs_wpred)), arr.ind = TRUE)
		lam.min <- lambda.vec[inds]

	} else {

		lam.min <- NULL

	}

		fit <-  list("beta" = beta.full, 
					"sparsity.mat" = sparsity.mat,
		     		"err.pred" = errs_pred, 
		     		"err.wpred" = errs_wpred, 
		     		"err.spec" = errs_spec, 
		     		"err.nuc" = errs_nuc, 
		     		"Y.offset" = apply(Y, 2, mean), 
		     		"X.offset" = apply(X, 2, mean),
		     		"Y.sd" = apply(Y, 2, sd), 
		     		"X.sd" = apply(X, 2, sd),
					"lambda.vec" = lambda.vec, 
					"lam.min" = lam.min,
					"standardize" = standardize)

		class(fit) <- "MSRL"
		return(fit)

	}
	

	MSRL.predict <- function(Xnew, fit, lambda = NULL){

		if(is.null(lambda)){
			lambda <- fit$lam.min
			if(is.null(fit$lam.min)){
				stop('No tuning parameters selected by CV')
			}

		} 
		
		lam.ind <- which(fit$lambda.vec == lambda)
		p <- length(fit$X.offset)
		q <- length(fit$Y.offset)
		if(!fit$standardize){
			beta.vec <- fit$beta[,lam.ind]
			beta.mat <- matrix(beta.vec, byrow=FALSE, nrow=p, ncol=q)
		} else {
			beta.vec <- fit$beta[,lam.ind]
			beta.mat <- tcrossprod(1/fit$X.sd, rep(1, q))*matrix(beta.vec, byrow=FALSE, nrow=p, ncol=q)
		}
		
		# --- get intercept 	
		B0 <- fit$Y.offset - crossprod(beta.mat, fit$X.offset)
		if(dim(Xnew)[1] > 1){
			preds <- tcrossprod(rep(1, dim(Xnew)[1]), B0) + tcrossprod(Xnew, t(beta.mat))
		} else {
			preds <- B0 + crossprod(beta.mat, Xnew)
		}

		return(list("pred" = preds, "beta0" = B0, "beta" = beta.mat))

	}


	MSRL.coef <- function(fit, lambda = NULL){

		if(is.null(lambda)){
			lambda <- fit$lam.min
			if(is.null(fit$lam.min)){
				stop('No tuning parameters selected by CV')
			}
		} 

		lam.ind <- which(fit$lambda.vec == lambda)
		p <- length(fit$X.offset)
		q <- length(fit$Y.offset)
		if(!fit$standardize){
			beta.vec <- fit$beta[,lam.ind]
			beta.mat <- matrix(beta.vec, byrow=FALSE, nrow=p, ncol=q)
		} else {
			beta.vec <- fit$beta[,lam.ind]
			beta.mat <- tcrossprod(1/fit$X.sd, rep(1, q))*matrix(beta.vec, byrow=FALSE, nrow=p, ncol=q)
		}
		# --- get intercept 	
		B0 <- fit$Y.offset - crossprod(beta.mat, fit$X.offset)
		return(list("beta0" = B0, "beta" = beta.mat))
	}
