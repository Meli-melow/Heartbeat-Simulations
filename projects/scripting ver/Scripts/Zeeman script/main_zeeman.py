from zeeman_cor import *
import numpy as np
from matplotlib import pyplot as plt
from scipy.integrate import odeint

#%%

# Zeeman system parameters
eps = 0.1
T = 1.0
# Initial point, for instance try x0 = ((1/3)**0.5 and y0 = 2*((1/3)**0.5)**3)
x0, y0 = zeeman_initial_point()

# Prepare the grid where to plot the system behavior
x_axis = [-2, 2]
y_axis = [-1, 1]
X, Y = np.meshgrid(np.arange(x_axis[0], x_axis[1] + 0.01, 0.1), np.arange(y_axis[0], y_axis[1] + 0.01, 0.1))
# 1st equation evaluated for 2500 points time points
x_points = np.linspace(x_axis[0], x_axis[1], 2501)
y_points = T * x_points - x_points ** 3
# Use 2500 time points
t_points = np.linspace(0, 20, 2501)


# Plot 3 graphs together, first associated to the systole for several values of x_d
# x-axis : length of myocardial cells
# y-axis : global cardiac electrochemical activity
print("Plotting Zeeman system during systole in process...")
# Example of average length values of myocardial cells during systole
Xs = np.array([-0.5, 0, 0.5]) * T ** 0.5
for k in range(3):
    x_s = Xs[k]
    plt.subplot(1, 3, k+1)
    plt.xlim(x_axis)
    plt.ylim(y_axis)
    # Plot the first equation
    plt.plot(x_points, y_points, 'k', lw=2)
    dXdt, dYdt = zeeman_no_peacemaker([X, Y], eps, T, x_s)
    plt.streamplot(X, Y, dXdt, dYdt)
    # Plot the second equation, compute the derivative of x with odeint
    sol_0 = odeint(lambda v, t: zeeman_no_peacemaker(v, eps, T, x_s), [x0, y0], t_points)
    x_sol, y_sol = sol_0.T
    plt.plot(x_sol, y_sol, 'r', lw=1.5)
    plt.grid()
    a = x_extrema(T)
    plt.plot(x0, y0, 'og')
    # Plot equilibrium points
    plt.plot(a, T*a - a**3, 'or')
    plt.plot(-a, a**3 - T*a, 'or')
    # Format Zeeman parameters using '$\' or $, and by adding '$' after printing format
    plt.title("Full cardiac cycle without peacemaker\n$\epsilon = %.2g$, $x_s = %.2f$, $T = %.2f$\n"%(eps, x_s, T)
              + "initial point = (%.3f, %.3f)"%(x0, y0))
plt.show()

# Plot 3 graphs together depicting the diastole
# Example of average length values of myocardial cells during diastole for several values of x_s
Xd = np.array([-2.5, 1, 3.5]) * T ** 0.5
print("Plotting Zeeman system during diastole in process...")
for k in range(3):
    x_s = Xs[k]
    x_d = Xd[k]
    plt.subplot(1, 3, k+1)
    plt.xlim(x_axis)
    plt.ylim(y_axis)
    # Plot the first equation
    plt.plot(x_points, y_points, 'k', lw=2)
    dXdt, dYdt = zeeman_constant_peacemaker([X, Y], eps, T, x_s, x_d)
    plt.streamplot(X, Y, dXdt, dYdt)
    # Plot the second equation, compute the derivative of x with odeint
    sol_1 = odeint(lambda y, t: zeeman_constant_peacemaker(y, eps, T, x_s, x_d), [x0, y0], t_points)
    x_sol, y_sol = sol_1.T
    plt.plot(x_sol, y_sol, 'r', lw=1.5)
    plt.grid()
    a = x_extrema(T)
    plt.plot(x0, y0, 'og')
    # Plot equilibrium points
    plt.plot(a, T*a - a**3, 'or')
    plt.plot(-a, a**3 - T*a, 'or')
    # Format Zeeman parameters using '$\' or $, and by adding '$' after printing format
    plt.title("Full cardiac cycle with constant peacemaker\n(value = 1)\n"
              "$\epsilon = %.2g$, $x_s = %.2f$, $x_d = %.2f$, $T = %.2f$\n"%(eps, x_s, x_d, T)
              + "initial point = (%.3f, %.3f)"%(x0, y0))
plt.show()