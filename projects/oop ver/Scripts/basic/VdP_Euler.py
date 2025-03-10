from TIPE import *
import numpy as np


class VdPEuler:
    """
    Given an instance of the Van der Pol model, this class solves it using Euler's method.

    Fields
    ------------
    x: array-like[float]
        The solution, first coordinate of zny Van der Pol system.
        A one-dimensional containing ``self.nb_points`` points.
    x_dot: array-like[float]
        The derivative of the solution, second coordinate of zny Van der Pol system.
        A one-dimensional containing ``self.nb_points`` points.
    integration_time: float
            Total integration time in seconds.
    nb_points: int
        Total number of time_points.
    time_points: array-like[float]
        A one-dimensional vectors containing ``self.nb_points`` time points.
    """

    def __init__(self, x_0, integration_time: float, nb_points: int):
        """
        Parameters
        ----------
        x_0: array_like[float]
            Two-dimensional vector that is the Van der Pol initial condition.
            First component is the first value of the Van der Pol solution.
            Second component is the derivative of the latter.
        integration_time: float
            Total integration time in seconds.
        nb_points: int
            Total number of time_points.
        """

        # Add the first value of the solution and its derivative
        self.x = [x_0[0]]
        self.x_dot = [x_0[1]]
        self.nb_points = nb_points
        self.integration_time = integration_time
        self.time_points = []

    def create_time_points(self):
        """
        Create a time points vector containing ``self.nb_points-1`` points.
        """

        self.time_points = np.linspace(0, self.integration_time, self.nb_points)

    def to_string_solver(self, system: VDPSystem):
        """
        Print in the terminal information about the chosen Van der Pol system and the solving parameters.

        Parameters
        ----------
        system: VDPSystem
            Instance of Van der Pol parameters.
        """

        print("Van der Pol parameters : {epsilon = " + f"{system.epsilon}, pulsation = {system.w_0}" + "}\n"
              + "Solving parameters : method : euler ; initial point = (" + f"{self.x[0]}, {self.x_dot[0]}) ; integration time = " + f"{self.integration_time}"
              + f"; number of points : {self.nb_points} ; time step : {self.integration_time / self.nb_points}")

    def vdp_evaluate(self, system: VDPSystem, x):
        """
        Compute the derivative of the Van der Pol system.

        Parameters
        ----------
        system: VDPSystem
            Instance of Van der Pol parameters.
        x: array_like[float]
            The two-dimensional system being solved.

        Returns
        -------
        array_like[float]
            Computation of the derivative of the Van der Pol system.
        """

        x1, x2 = x                                               # changement de variable matriciel
        return np.array([x2, system.epsilon * (1 - x1 ** 2) * x2 - system.w_0 * x1])       # Allow products with a scalar

    def euler(self, system: VDPSystem):
        """
        Compute the derivative of the Van der Pol system for given point in time.
        The solution is built up directly.

        Parameters
        ----------
        system: VDPSystem
            Instance of Van der Pol parameters.

        Returns
        -------
        array_like[float]
            Computation of the derivative of the Van der Pol system.
        """

        # Euler's step
        h = self.integration_time / (self.nb_points - 1)
        euler_sol = [[self.x[0], self.x_dot[0]]]
        # Upper bound is excluded because the time vector is built with np.linespace
        for i in range(self.nb_points-1):
            # Use latest values of the Van der Pol solution
            sol_i = euler_sol[-1]
            # Add the next value of the solution, image of the next time point
            euler_sol += [sol_i + h * self.vdp_evaluate(system, sol_i)]
        return np.array(euler_sol)