""" FULL REPAIRED :) """

epsilon, w_0 = 2, 1

def f(v):                                       # changement de variable matriciel
    x, y = v
    return y, epsilon*(1 - x**2)*y -w_0*x                                          # changement y = x_point