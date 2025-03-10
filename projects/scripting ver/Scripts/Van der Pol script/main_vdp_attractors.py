from VdP_Attractors import *
from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt


# First solving VdP (Van der Pol)

# Coordinates of initial condition
x0 = 0.5
y0 = 0.5
# Integration time
# Results references : 5.0 s -> half heart cycle, 8.3 s -> full cycle, system convergence with odeint ~ 12.7 s
tmax = 12.70
# dt : Time step ; N number of time points ; t : 1D array of time points
dt = 0.1
N = int(tmax/dt)
t = np.linspace(0, tmax, N)
# Another syntax to use odeint. Using lambda to avoid TypeError
P = odeint(lambda v, t: f(v), [x0, y0], t) # renvoit tableau 2 dim, de taille 2*N, P[0][j] = x(j) et P[1][j] = y(j)


# Plotting the vectorial phase portrait

xmin, xmax, dx = -4, 4, 0.1
ymin, ymax, dy = -4, 4, 0.1
# Create a network of vectors
X, Y = np.meshgrid(np.arange(xmin, xmax, dx), np.arange(ymin, ymax, dy) )
# Use the solution f to plot the phase portrait
f1, f2 = f([X, Y])
plt.streamplot(X, Y, f1, f2)
plt.title('Vectorial field (x, y)')
plt.xlabel('x')
plt.ylabel('y')
# Format Van der Pol parameters using '$\' or $, and by adding '$' after printing format
plt.title("Phase portrait : $\epsilon = %2.f$, $w_0 = %2.f$"%(epsilon, w_0))
plt.plot(0, 0, 'og')
plt.show()


# Then progressively plotting the basic phase portrait

for j in range(N):
    plt.clf()
    plt.plot(P[:j+1, 0], P[:j+1, 1], 'b-', label='Trajectory')
    # Plot in red the current point of the VdP solution
    plt.plot(P[j, 0], P[j, 1], 'ro')
    # Format Van der Pol parameters using '$\' or $, and by adding '$' after printing format
    plt.title( "Solution y(t) = x'(t), y'(t) = VdP(x(t), y(t))\n $\epsilon = %2.f$, pulsation = %2.f rad/s, t = %3.1f s"
               %(epsilon, w_0, t[j]) )
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    # For real-time plotting
    plt.pause(0.1)
plt.show()