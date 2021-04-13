#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

//[[Rcpp::depends("RcppArmadillo")]]
//[[Rcpp::export]]
double EvalObjFuncCpp(arma::mat Y, arma::mat X, arma::mat Omega, arma::mat weights,
  arma::mat beta, double lambda, double alpha){

  double n = X.n_rows;
  double p = X.n_cols; 
  double tau = (1-alpha)*lambda; 
  mat temp = Y - X*beta; 
  double keep = 0.0; 
  for(int i=0; i<n; ++i){
    keep += pow(n, -1)*(temp(span(i,i), span::all)*Omega*temp(span(i,i),span::all).t()).eval()(0,0); 
  }
  keep += lambda*alpha*accu(weights%abs(beta)); 

  for(int i=0; i<p; ++i){
    keep += tau*sqrt(accu(pow(beta(span(i,i),span::all), 2))); 
  }
  return keep; 

}

//[[Rcpp::depends("RcppArmadillo")]]
//[[Rcpp::export]]
arma::mat SFmat( arma::mat beta, arma::mat lam){
	int p = beta.n_rows;
	int q = beta.n_cols;
	mat A(p,q);  A.zeros();
	for(int j=0; j < p; j++){
		for(int k=0; k < q; k++){
			
			if(beta.at(j,k) > 0 && fabs(beta.at(j,k)) > lam.at(j,k)){
				A.at(j,k) =  beta.at(j,k) - lam.at(j,k);
			}

			if(beta.at(j,k) < 0 && fabs(beta.at(j,k)) > lam.at(j,k)){
				A.at(j,k) = beta.at(j,k) + lam.at(j,k);
			}
		}
	}
	return A;
}


//[[Rcpp::depends("RcppArmadillo")]]
//[[Rcpp::export]]
arma::mat rowSFmat( arma::mat beta, double lam){
	int p = beta.n_rows;
	int q = beta.n_cols;
	mat A(p,q);  A.zeros();
	double temp = 0.0;
	for(int j=0; j < p; j++){
		  temp = sqrt(accu(pow(beta(span(j,j), span::all), 2))); 
			if(temp > lam){
				A(span(j,j),span::all) = beta(span(j,j), span::all)*(1 - lam/temp); 
			}
	}
	return A; 
}


//[[Rcpp::depends("RcppArmadillo")]]
//[[Rcpp::export]]
arma::mat AccProxGrad(arma::mat Y, arma::mat X, arma::mat Omega, arma::mat weights, 
	arma::mat betainit, double lambda, double alpha, 
	int maxiter, double tol, arma::mat XtX, double XtXeig, double Qeig){

	int n = X.n_rows;
	mat XtY = X.t()*Y; 
  	double L = n/(4*XtXeig*Qeig);
	mat weightmat = L*alpha*lambda*weights; 
	double weightscale = L*(1-alpha)*lambda; 
	mat betakm1 = betainit;
	mat betak = betainit;
	mat betakp1 = betainit;
	double upobjfuncval = 0; 
	double resid = 0; 
	double wk = 0.0;
	double kiter = 0.0;
	bool iterating = TRUE;
	double origobjfuncval = EvalObjFuncCpp(Y, X, Omega, weights, betainit, lambda, alpha);
	double oldobjfuncval = origobjfuncval; 
	mat temp = betainit; 

   	while(iterating){

  	 	temp = betak + (wk)*(betak - betakm1);
  	 	temp = temp + 2*(L/n)*((XtY - XtX*temp)*Omega); 
  	 	temp = SFmat(temp, weightmat); 
		betakp1 = rowSFmat(temp, weightscale);
	    upobjfuncval = EvalObjFuncCpp(Y, X, Omega, weights, betakp1, lambda, alpha);
	  	resid = fabs(upobjfuncval - oldobjfuncval)/fabs(origobjfuncval);
      
      	// if(kiter> 100){
        if(resid < tol){
  				iterating = FALSE;
  			}
		// }
		wk = (kiter + 1)/(kiter + 4);
		kiter = kiter + 1;
		oldobjfuncval = upobjfuncval;
		betakm1 = betak;
		betak = betakp1;

		if(kiter > maxiter){
		  iterating = FALSE;
		}

 	}


  	return betakp1;
 }





