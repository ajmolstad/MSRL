import numpy as np
from msrl import MSRL_cv

n, p, q = 20, 30, 10

np.random.seed(18)
X = np.random.randn(n, p)
Y = np.random.randn(n, q)

betas, lambdas, sparsity = MSRL_cv(X, Y, nlambda=10, tol=1e-12)
print(betas.shape)
