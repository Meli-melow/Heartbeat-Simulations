from VdP_System import *
from numpy import floating, ndarray, dtype
from typing import Any
from abc import abstractmethod
import numpy as np


class VdPSolver:
    """
    This class solves the Van der Pol equation using the strategy pattern.
    The solving relies on both Euler's method and odeint.

    Fields : class attributes in order to compare both Euler's method and odeint solvers' accuracy.
    ------
    integration_time : float
        Total integration time.
    nb_points : int
        Total number of points.
    time_points : array_like[float]
        Array of time points.
    x_0 : array_like[float]
        Initial condition.
    """

    integration_time = 50.0
    nb_points = 10000
    time_points = np.linspace(0, integration_time, nb_points)
    x_0 = np.array([2.0, 1.0])

    def __init__(self, system: VDPSystem):
        """
        Parameters
        ----------
        system : VdPSystem
            Instance of the Van der Pol model.
        """

        self.system = system

    def to_string_solver(self):
        """
        Print in the terminal information about the used the Van der Pol parameters and solver.
        """

        self.system.to_string_system()
        print(f"Initial point = {VdPSolver.x_0}) ; integration time = " + f"{VdPSolver.integration_time}"
              + f"; number of points : {VdPSolver.nb_points} ; time step : {VdPSolver.integration_time / VdPSolver.nb_points}")

    def vdp_system_equation(cls, v):
        """
        Compute the equation of the Van der Pol system equation.
        This is a class method as both solvers rely on the Van der Pol equation.

        Parameters
        ----------
        v: array_like[float]
            Two-dimensional array associated to the Van der Pol system.

        Return
        ------
        float
            Evaluation of the second derivative of the Van der Pol system.
        """

        x, xp = v[0], v[1]
        return cls.system.epsilon * (1 - x ** 2) * xp - cls.system.w_0 * x

    @abstractmethod
    def vdp_evaluate(self, x, t):
        """
        Evaluate the Van der Pol system for given time points.
        This method is an abstract strategy used to solve the Van der Pol equation using either Euler's method or odeint.

        Parameters
        ----------
        x : array_like[float]
            The array associated to the Van der Pol system.
        t : array_like[float]
            Array of time points.
        """
        pass