# from rpy2 import robjects
from setuptools import setup
import rpy2.robjects.packages as rpackages


utils = rpackages.importr('utils')
utils.chooseCRANmirror(ind=1)
# installing devtools may take a while the first time (10 miutes)
if not rpackages.isinstalled("MSRL"):
    if not rpackages.isinstalled('devtools'):
        utils.install_packages("devtools")
    devtools = rpackages.importr('devtools')
    devtools.install_github("ajmolstad/MSRL")

desc = "Perform cross-validation, prediction, and extract coefficients from" \
    "fitted models for the multivariate square-root lasso."

setup(
    name="msrl",
    version='0.1',
    description=desc,
    author_email="amolstad@ufl.edu"
)
