# rotational motion of artificial satellite
# referece: T. R. Kane, P. W. Likins, and D. A. Levinson
#   "Spacecraft Dynamics" section 3.2

import numpy as np
from scipy.integrate import odeint

def SatRot(y,t,K1,K2,e=0):
    """
    y = initial condition for dependent variables
    y[0:4] = Euler parameters e1,e2,e3,e4
    y[4:7] = angular velocity / orbital ang. vel.
    y[7] = true anomaly / rad
    t = time sequence at which to evaluate y
        unit of time is (orbital period)/(2*pi)
    K1,K2 = (I2-I3)/I1, (I3-I1)/I2 (moments of inertia)
    e = eccentricity of orbit
    """
    a = 1/(1 - e**2)**1.5
    b = 3*a**2
    K3 = -(K1 + K2)/(1 + K1*K2)

    def difeq(y,t):
        e1,e2,e3,e4 = y[:4] # euler parameters
        u1,u2,u3 = y[4:7] # angular velocity
        phi = y[7] # true anomaliy
        e_cos = 1 + e*np.cos(phi)
        dphi = a*e_cos**2
        r3 = b*e_cos**3

        de1 = (e2*u3 - e3*u2 + e4*u1 + e2*dphi)/2
        de2 = (e3*u1 - e1*u3 + e4*u2 - e1*dphi)/2
        de3 = (e1*u2 - e2*u1 + e4*u3 - e4*dphi)/2
        de4 = (-e1*u1 - e2*u2 - e3*(u3 - dphi))/2

        C11 = 1 - 2*(e2**2 + e3**2)
        C12 = 2*(e1*e2 - e3*e4)
        C13 = 2*(e1*e3 + e2*e4)

        du1 = K1*(u2*u3 - r3*C12*C13)
        du2 = K2*(u3*u1 - r3*C11*C13)
        du3 = K3*(u1*u2 - r3*C11*C12)

        return de1,de2,de3,de4,du1,du2,du3,dphi

    return odeint(difeq, y, t)
