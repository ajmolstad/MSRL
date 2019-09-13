import numpy as np
from rpy2 import robjects
import rpy2.robjects.packages as rpackages
from rpy2.robjects import numpy2ri, pandas2ri

numpy2ri.activate()
pandas2ri.activate()
msrl = rpackages.importr("MSRL")


def MSRL_cv(X, Y, nlambda, tol, quiet=False):
    msrl_fit = msrl.MSRL_cv(X=X, Y=Y, nlambda=nlambda, tol=tol, quiet=quiet)
    n, p = X.shape
    q = Y.shape[1]
    fit_dict = dict(zip(msrl_fit.names, list(msrl_fit)))

    as_matrix = robjects.r['as']
    betas = np.array(as_matrix(fit_dict["beta"], "matrix"))
    betas = np.array([beta.reshape((p, q), order='F') for beta in betas.T])
    sparsity = np.array(fit_dict["sparsity.mat"])
    lambdas = np.array(fit_dict["lambda.vec"])

    return betas, lambdas, sparsity
