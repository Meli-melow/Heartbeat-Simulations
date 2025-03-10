from VdP_Euler import *
import matplotlib.pyplot as plt
# Might be required in a Spyder IDE for plotting several figures at once : import subprocess


vdp_system = VDPSystem(2.0, 1.0)
solve_euler = VdPEuler(np.array([2, 1]), 50, 10000)
solve_euler.create_time_points()
solve_euler.to_string_solver(vdp_system)
# Example of Euler solving
S = solve_euler.euler(vdp_system)
t1 = solve_euler.time_points
x1, v1 = S[:, 0], S[:, 1]

# Plotting the solution and its phase portrait
#plt.clf()
plt.subplot(221)
plt.plot(t1, x1)
eps = vdp_system.epsilon; w_0 = vdp_system.w_0;
t_max = solve_euler.integration_time; N = solve_euler.nb_points
# Format Van der Pol parameters using '$\' or $, and by adding '$' after printing format
plt.title("Van der Pol solution x(t)\n" + "$\epsilon = %2.f$, $w_0 = %2.f$ rad/s, $x_0 = (%s, %s)$, $t_{total} = %3.1f$ s, N = %d points"
          % (eps, w_0, str(x1[0]), str(v1[0]), t_max, N))
plt.xlabel("t (s)")
plt.ylabel("x(t)")

plt.subplot(222)
plt.plot(x1, v1)
plt.title("Corresponding phase portrait")
plt.xlabel("x")
plt.ylabel("x_dot")
# Plot the initial condition
plt.plot(x1[0], v1[0], 'og')

plt.show()