from VdP_Solver import *


class VdPEulerSolver(VdPSolver):
    def to_string_solver(self):
        """
        Prints in the terminal information about the Euler Van der Pol solver.
        """

        print("Solving method : euler")
        super().to_string_solver()

    def vdp_evaluate(self, x=None, t=None):
        """
        This method computes the Euler's method' strategy.
        It directly builds up the solution, a two-dimensional time vector that has a total of 2*VDPSolver.nb_points.
        Which is why the ``x`` parameter is empty (``x=None``).
        Unlike odeint, it does not directly rely on the time vector (``t=None``)

        Returns
        -------
        array_like[float]
            The solution computed with Euler's method.
        """

        # Euler time step
        h = self.integration_time / (self.nb_points - 1)
        euler_sol = [[VdPSolver.x_0[0], VdPEulerSolver.x_0[1]]]
        # Compute the system for each time point
        # Upper bound excluded since the time vector is build with np.linespace
        for i in range(VdPSolver.nb_points-1):
            sol_i = euler_sol[-1]
            # Use a numpy array to allow scalar multiplication with h
            euler_sol += [sol_i + h * np.array([sol_i[1], super().vdp_system_equation(sol_i)])]
        return np.array(euler_sol)