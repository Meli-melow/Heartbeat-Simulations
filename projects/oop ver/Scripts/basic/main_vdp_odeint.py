from VdP_Odeint import *
from scipy.integrate import odeint
from matplotlib.pyplot import *
import matplotlib.pyplot as plt
import numpy as np
# Might be required in a Spyder IDE to use odeint : from scipy import *


# Simulating VdP system with 10000 time points up to 50 s
vdp_system = VDPSystem(10.0, 0.99)
t = VdPOdeint.time_points
# Generate initial conditions
Ci = np.array([[0.5*np.cos(k*np.pi/5), -0.1*np.pi*np.sin(k*np.pi/5)] for k in range(10)])
#Another example : Ci2 = np.array([[2.5*np.cos(k*np.pi/5), -0.5*np.pi*np.sin(k*np.pi/5)] for k in range(10)])
# Identify each Van der Pol system with one color and a label
colors = {1: 'b', 2: 'r', 3: 'g', 4: 'c', 5: 'm', 6: 'y', 7: 'k', 8: 'orange', 9: 'purple', 10: 'gray'}
labels = ['n°1', 'n°2', 'n°3', 'n°4', 'n°5', 'n°6', 'n°7', 'n°8', 'n°9', 'n°10']

#plt.clf()

# Plot the solution and the phase portrait of each Van der Pol system
for i in range(10):
    solve_odeint_i = VdPOdeint(vdp_system, Ci[i])
    #Information on current plotting : print(f"iteration {i}")
    #solve_odeint1.toStringSolver(vdp_system)
    Ys = odeint(solve_odeint_i.vdp_odeint, solve_odeint_i.init_point, t)
    #print("\n")
    plt.subplot(221)
    plt.plot(t, Ys[:, 0], colors[i+1], label=labels[i])
    plt.plot(0, Ci[i][0], 'o', color=colors[i+1])
    # Format Van der Pol parameters using '$\' or $, and by adding '$' after printing format
    plt.title("x(t)\n$\epsilon = %2.f$, $w_0 = %2.f$"%(vdp_system.epsilon, vdp_system.w_0))
    plt.ylabel("x"); plt.xlabel("t"); grid()
    plt.subplot(222)
    plt.plot(Ys[:, 1], Ys[:, 0], colors[i+1])
    # Plot initial conditions
    plt.plot(Ci[i][1], Ci[i][0], 'o', color=colors[i+1])
    # Format Van der Pol parameters using '$\' or $, and by adding '$' after printing format
    plt.title("Phase portrait of Van der Pol systems\n$\epsilon = %2.f$, $w_0 = %2.f$"%(vdp_system.epsilon, vdp_system.w_0))
    plt.ylabel("x"); plt.xlabel("x_dot"); grid()

# Note that for the time representations, some points are on top of others
#for instance gray initial point is on top the red one
plt.legend(labels, loc='best')
plt.show()