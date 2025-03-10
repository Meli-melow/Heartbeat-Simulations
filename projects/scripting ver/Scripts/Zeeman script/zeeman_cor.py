import math


def zeeman_no_peacemaker(v, eps: float, tension: float, x_s: float):
    """
    Compute the Zeeman model for a given two-dimensional vector describing the length of the length of myocardial cells
    and a given the global cardiac electrochemical activity.

     Parameters
     ----------
     v : array_like
          A two-dimensional vector. It is associated to the length of ...
          Second component is the derivative.
     eps : float
          The epsilon value of the zeeman system.
     tension : float
          Value of the tension parameter.
     x_s : float
          The average length of the myocardial cells during systole.

     Returns
     -------
     array_like[float]
          A vector computing the Zeeman model.
    """
    x, y = v
    return [-(y + x ** 3 - tension * x) / eps, x - x_s]


def zeeman_constant_peacemaker(v, eps: float, tension: float, x_s: float, x_d: float):
    """
    Compute the Zeeman model for a given two-dimensional vector containing the length of ... and a given heart beat tension.
    Note that the peacemaker (u function) is here simplified as a constant function (value = 1).

     Parameters
     ----------
     v : array_like
          A two-dimensional vector.
          The first component computes to the length of myocardial cells.
          The second component is the global electrochemical activity of the heart.
     eps : float
          The epsilon value of the zeeman system.
     tension : float
          Value of the tension parameter.
     x_s : float
          The average length of the myocardial cells during systole.
     x_d : float
          The average length of the myocardial cells during diastole.

     Returns
     -------
     array_like[float]
          A two-dimensional vector implementing the Zeeman system.
    """

    x, y = v
    return [-(y + x ** 3 - tension * x) / eps, x - x_s + x_d - x_s]


def zeeman_initial_point():
    """
    Initialize initial point for the Zeeman system by typing coordinates in the command line.

    Returns
    -------
    tuple_like[float]
        The float evaluations of the chosen initial values for both x and y coordinates.
    """

    print("Choose x initial value : ")
    x0 = input()
    print("Choose y initial value : ")
    y0 = input()
    return float(eval(x0)), float(eval(y0))


def x_extrema(tension: float):
    """
    Compute extrema of x the length of myocardial cells.

    Parameters
    ----------
    tension : float
        Heart beat tension of the Zeeman model.

    Returns
    -------
    float
        Extrema of x.
    """

    return math.pow(tension / 3.0, 0.5)