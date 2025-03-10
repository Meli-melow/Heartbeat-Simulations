class VDPSystem:
    """
    An implementation of the Van der Pol model.
    """

    def __init__(self, epsilon: float, w_0: float):
        """
        Parameters
        ----------
        epsilon : float
            Value of the epsilon parameter.
        w_0 : float
            Value of the pulsation parameter.
        """

        self.epsilon = epsilon
        self.w_0 = w_0