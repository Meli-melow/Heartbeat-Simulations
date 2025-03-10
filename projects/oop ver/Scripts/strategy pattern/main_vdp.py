from VdP_Solver import *
from VdPEulerStrategy import VdPEulerSolver
from VdPOdeintStrategy import VdPOdeintSolver
from scipy.integrate import odeint
# Might be required in a Spyder IDE to use odeint : from scipy import *
import matplotlib.pyplot as plt
# Might be required in a Spyder IDE for plotting several figures at once : import subprocess


vdp_system = VDPSystem(2.0, 1.0)
t = VdPSolver.time_points

# Example of Euler resolution
solve_euler = VdPEulerSolver(vdp_system)
solve_euler.to_string_solver()
S1 = solve_euler.vdp_evaluate(vdp_system)
x1, v1 = S1[:, 0], S1[:, 1]

# Example of Odeint resolution
solve_odeint = VdPOdeintSolver(vdp_system)
solve_odeint.to_string_solver()
S2 = odeint(solve_odeint.vdp_evaluate, VdPSolver.x_0, VdPEulerSolver.time_points)
x2, v2 = S2[:, 0], S2[:, 1]

eps = vdp_system.epsilon; w_0 = vdp_system.w_0; init_sol = VdPSolver.x_0;
t_max = VdPSolver.integration_time; N = VdPSolver.nb_points

# Plotting the Euler solution
#plt.clf()
plt.subplot(221)
plt.plot(t, x1)
# Format Van der Pol parameters using '$\' or $, and by adding '$' after printing format
plt.title("x(t)\n" + "$\epsilon = %2.f$, $w_0 = %2.f$ rad/s, $x_0 = (%s, %s)$, $t_{total} = %3.1f$ s, N = %d points"
          % (eps, w_0, str(init_sol[0]), str(init_sol[1]), t_max, N))
plt.xlabel("t")
plt.ylabel("x(t) with Euler")
plt.subplot(222)
plt.plot(x1, v1)
plt.title("Euler phase portrait")
plt.xlabel("x")
plt.ylabel("x_dot")
plt.plot(init_sol[0], init_sol[1], 'og')

# Plotting the Odeint solution
plt.subplot(223)
plt.plot(t, x2)
# Format Zeeman parameters using '$\' or $, and by adding '$' after printing format
plt.title("x(t)\n" + "$\epsilon = %2.f$, $w_0 = %2.f$ rad/s, $x_0 = (%s, %s)$, $t_{total} = %3.1f$ s, N = %d points"
          % (eps, w_0, str(init_sol[0]), str(init_sol[1]), t_max, N))
plt.xlabel("t")
plt.ylabel("x(t) with odeint")
plt.subplot(224)
plt.plot(x2, v2)
plt.title("Odeint phase portrait")
plt.xlabel("x")
plt.ylabel("x_dot")
plt.plot(init_sol[0], init_sol[1], 'og')

plt.show()