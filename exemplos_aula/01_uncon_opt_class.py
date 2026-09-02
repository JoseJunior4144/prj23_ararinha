'''
INSTITUTO TECNOLÓGICO DE AERONÁUTICA
PROGRAMA DE ESPECIALIZAÇÃO EM ENGENHARIA AERONÁUTICA
OTIMIZAÇÃO MULTIDISCIPLINAR

Código com aplicações de ferramentas do Python para otimização
não restringida

Maj. Ney Sêcco 2025
'''

# IMPORTS
import numpy as np
from scipy.optimize import minimize, differential_evolution

# EXECUTION

#==============================================

# Define number of input variables
nvar = 100

# Define objective function
def objfun(x):

    fun = np.sum(100*(x[1:]-x[:-1]**2)**2 + (1-x[:-1])**2)

    return fun

def objfun_grad(x):

    grad = np.zeros(len(x))
    grad[0] = -400*(x[1]-x[0]**2)*x[0] - 2*(1-x[0])
    grad[1:-1] = -400*(x[2:]-x[1:-1]**2)*x[1:-1] - 2*(1-x[1:-1]) + 200*(x[1:-1] - x[:-2]**2)
    grad[-1] = 200*(x[-1]-x[-2]**2)

    return grad

#==============================================

# Define initial guess
x0 = np.zeros(nvar)

###! BFGS - With gradient

# Set optimization options
options = {'maxiter': 200_000}

# Run the optimization algorithm
# result = minimize(fun = objfun, x0 = x0, jac = objfun_grad,
#          method = 'BFGS', tol = 1e-6, options = options)

# # Print results
# print(result)

###! BFGS - Without gradient
# result2 = minimize(fun = objfun, x0 = x0, 
#                    method = 'BFGS', tol = 1e-6, options = options)

# print(result2)

###! NELDER-MEAD
options = {'maxiter': 200_000,
           'fatol': 1e-12,
           'adaptative': True}

# result3 = minimize(fun = objfun, x0 = x0, 
#                    method = 'Nelder-Mead', options = options)

# print(result3)

###! DIFFERENTIAL EVOLUTION

# Define bounds to determine initial population
bounds = [[-5, 5]] * nvar

# Solve the optimization problem
result4 = differential_evolution(func = objfun, bounds = bounds, seed = 1,
                                 maxiter = 5_000, polish = False, atol = 1e-12,
                                 popsize = 300)

# Print results
print(result4)