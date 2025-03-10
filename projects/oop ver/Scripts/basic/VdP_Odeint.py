from TIPE import *
from numpy import floating, ndarray, dtype
from typing import Any
import numpy as np


class VdPOdeint:
    """
    Given an instance of the Van der Pol model, this class solves it using the odeint from scipy.
    The definition is adapted for computing Van der Pol systems, with multiples initial conditions.

    Class Variables
    ---------------
    integration_time : float
        Value of the total integration time.
    nb_points: int
        Total number of time points.
    time_points: array_like[float]
        Vector containing the time points.
    """

    integration_time = 50.0
    nb_points = 10000
    time_points = np.linspace(0, integration_time, nb_points)

    def __init__(self, system: VDPSystem, x_0: ndarray[Any, dtype]):
        """
        Parameters
        ----------
        system: VDPSystem
            Instance of Van der Pol parameters.
        x_0: ndarray[Any, dtype]
            Initial condition of the Van der Pol system.
        """

        self.system = system
        self.init_point = x_0

    def to_string_solver(self, system: VDPSystem):
        """
        Print in the terminal information about the chosen Van der Pol system and the solving parameters.

        Parameters
        ----------
        system: VDPSystem
            Instance of Van der Pol parameters.
        """

        print("Van der Pol parameters : {epsilon = " + f"{system.epsilon}, pulsation = {system.w_0}" + "}\n"
              + "Solving parameters : method : odeint ; initial point = " + f"{self.init_point} ; integration time = " + f"{VdPOdeint.integration_time}"
              + f"; number of points : {VdPOdeint.nb_points} ; time step : {VdPOdeint.integration_time / VdPOdeint.nb_points}")

    def vdp_odeint(self, v: ndarray[Any, floating], t):
        """
        Compute the derivative of the Van der Pol system at a given point in time.
        The time vector is a necessary parameter for scipy.odeint().

        Parameters
        ----------
        v: ndarray[Any, floating]
            The two-dimensional system being solved.
        t: array_like[floating]
            Vector containing the time points.

        Returns
        -------
        array_like[float]
            Computation of the derivative of the Van der Pol system.
        """

        x, x_dot = v[0], v[1]
        return np.array([x_dot, self.system.epsilon*(1 - x**2)*x_dot - self.system.w_0*x])