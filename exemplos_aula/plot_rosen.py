'''
Plot 3D da função de Rosenbrock para 2 variáveis (x1, x2)
'''

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# Objective function (Rosenbrock, generalized)
def objfun(x):
    return np.sum(100*(x[1:]-x[:-1]**2)**2 + (1-x[:-1])**2)

# Grid
x1 = np.linspace(-2.0, 2.0, 200)
x2 = np.linspace(-1.0, 3.0, 200)
X1, X2 = np.meshgrid(x1, x2)
Y = 100*(X2 - X1**2)**2 + (1 - X1)**2

# Plot
fig = plt.figure(figsize=(9, 7))
ax = fig.add_subplot(projection='3d')
ax.plot_surface(X1, X2, Y, cmap='viridis', norm=LogNorm(), alpha=0.9,
                linewidth=0, antialiased=True, rstride=2, cstride=2)
ax.contour(X1, X2, Y, levels=np.logspace(-1, 3.5, 15), norm=LogNorm(),
           offset=0, cmap='viridis', linewidths=0.6)

# Minimum at (1, 1)
ax.scatter([1], [1], [objfun(np.array([1.0, 1.0]))], color='red', s=60,
           label='Mínimo (1, 1)')

ax.set_xlabel('$x_1$')
ax.set_ylabel('$x_2$')
ax.set_zlabel('$f(x_1, x_2)$')
ax.set_title('Função de Rosenbrock (2 variáveis)')
ax.view_init(elev=35, azim=-125)
ax.legend()

plt.tight_layout()
plt.show()
