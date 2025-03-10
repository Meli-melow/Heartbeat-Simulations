from VdP_Solver import *


class VdPOdeintSolver(VdPSolver):
    def to_string_solver(self):
        """
        Prints in the terminal information about the odeint Van der Pol solver.
        """

        print("Solving method : odeint")
        super().to_string_solver()

    def vdp_evaluate(self, x, t):
        """
        This method computes the odeint strategy.

        Parameters
        ----------
        x : array_like[float]
            The two-dimensional array associated to the Van der Pol system.
        t: array_like[float]
            The time point array required to use the scipy.odeint method.

        Returns
        -------
        array_like[float]
            The solution computed with Euler's method.
        """

        return np.array([x[1], super().vdp_system_equation(x)])