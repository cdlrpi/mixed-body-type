import numpy as np

def rk4(func, tspan, y, *args):
    """
    Inputs
       y: state vector [[x], [xdot]]
        func: function to be called, must return a state vector [xdot, xddot]
    """
    
    k1 = func(y, *args)
    
    y2 = np.array([y[0] + 1/2*dt*k1[0], y[1] + 1/2*dt*k1[1]]) 
    k2 = func(y2, *args)

    y3 = np.array([y[0] + 1/2*dt*k2[0], y[1] + 1/2*dt*k2[1]]) 
    k3 = func(y2, *args)

    y4 = np.array([y[0] + dt*k3[0], y[1] + dt*k3[1]]) 
    k4 = func(y2, *args)

    y = y + dt*(k1 + 2*k2 + 3*k3 + k4)/6 

    return(y)

# def rk4( f, x0, t ):
#     """Fourth-order Runge-Kutta method to solve x' = f(x,t) with x(t[0]) = x0.
# 
#     USAGE:
#         x = rk4(f, x0, t)
# 
#     INPUT:
#         f     - function of x and t equal to dx/dt.  x may be multivalued,
#                 in which case it should a list or a NumPy array.  In this
#                 case f must return a NumPy array with the same dimension
#                 as x.
#         x0    - the initial condition(s).  Specifies the value of x when
#                 t = t[0].  Can be either a scalar or a list or NumPy array
#                 if a system of equations is being solved.
#         t     - list or NumPy array of t values to compute solution at.
#                 t[0] is the the initial condition point, and the difference
#                 h=t[i+1]-t[i] determines the step size h.
# 
#     OUTPUT:
#         x     - NumPy array containing solution values corresponding to each
#                 entry in t array.  If a system is being solved, x will be
#                 an array of arrays.
#     """
# 
#     n = len( t )
#     x = numpy.array( [ x0 ] * n )
#     for i in xrange( n - 1 ):
#         h = t[i+1] - t[i]
#         k1 = h * f( x[i], t[i] )
#         k2 = h * f( x[i] + 0.5 * k1, t[i] + 0.5 * h )
#         k3 = h * f( x[i] + 0.5 * k2, t[i] + 0.5 * h )
#         k4 = h * f( x[i] + k3, t[i+1] )
#         x[i+1] = x[i] + ( k1 + 2.0 * ( k2 + k3 ) + k4 ) / 6.0
# 
#     return x
