class VDPSystem:
    """
    An implementation of the Van der Pol model.
    """
    def __init__(self, epsilon: float, w_0: float):
        """
        Fields : parameters of the initialized Van der Pol model.
        ------
        epsilon: float
            Value of the epsilon parameter.
        w_0: float
            Value of the pulsation parameter
        """

        self.epsilon = epsilon
        self.w_0 = w_0

    def to_string_system(self):
        """
        Print in the command line parameters or the initialized model.
        """

        print(f"System parameters : epsilon = {self.epsilon}, pulsation = {self.w_0}")