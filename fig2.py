# plot instability region
# reference: T. R. Kane, P. W. Likins and D. A. Levinson
#   "Spacecraft Dynamics" section 3.2

import numpy as np
import matplotlib.pyplot as plt
from numpy.polynomial import Polynomial

def KLL_fig325(e, phi, s1, s2, J1=1e-3, J2=2, N=200):
    """ plot figure 3.2.5 of Kane Likins Levinson
        (but ordinate and abscissa are interchanged)
    e = eccentricity
    phi = true anomaly
    s1,J1 = lower left corner of plot region
    s2,J2 = upper right corner of plot region
    J = I3/I1, s = spin freq. / orbital freq.
    """
    e2 = 1 - e**2
    e_cos = 1 + e*np.cos(phi)
    dphi = e_cos**2/e2**1.5
    r3 = 3*(e_cos/e2)**3
    dphi2 = dphi**2
    b = 4*dphi*r3

    def f1(s):
        t = r3/s
        a = np.sqrt(t**2 - 4*(dphi2 + t*(dphi - s)))
        J1 = (dphi - (t-a)/2)/s
        J2 = (dphi - (t+a)/2)/s
        return J1,J2

    def f2(J):
        s1 = dphi/J
        s2 = ((1-J)*r3/dphi + dphi)/J
        return s1,s2

    def f3(s):
        J1,J2 = [],[]
        for s in s:
            t = r3/s
            a = dphi2 + t*(dphi - s)
            p1 = b*(dphi - s)/s - 2*t*a
            p2 = t**2 + 2*a - 4*dphi*(dphi + t)
            p = Polynomial([a**2, p1, p2, -2*t, 1])
            z = p.roots()
            J = (dphi - np.real(z[np.isreal(z)]))/s
            J1.append(np.max(J))
            J2.append(np.min(J))

        return J1,J2

    x = np.linspace(s1, s2, N)
    x = x[x!=0]
    y1,y2 = f1(x)
    plt.fill_between(x, y1, y2, color='r', hatch='/', alpha=0.4)

    y = np.linspace(J1, J2, N)
    x1,x2 = f2(y)
    plt.fill_betweenx(y, x1, x2, color='b', hatch='-', alpha=0.4)

    y1,y2 = f3(x)
    plt.fill_between(x, y1, y2, color='g', hatch='|', alpha=0.4)


e = [0, 0.5, 0.5]
phi = [0, 0, np.pi]
label = [r'$(0,0)$',
         r'$(\frac{1}{2},0)$',
         r'$(\frac{1}{2},\pi)$']

s1,s2 = -4,6

plt.figure(figsize=(5, 10))

for i,(e,phi) in enumerate(zip(e,phi)):
    plt.subplot(3,1,i+1)
    KLL_fig325(e, phi, s1, s2)
    plt.ylabel(r'$J$')
    plt.yticks([0, 0.5, 1, 1.5, 2])
    plt.axis([s1, s2, 0, 2])
    plt.text(3.7, 1.85, r'$(e,\phi)=$' + label[i])
    if i<2: plt.tick_params(labelbottom=False)

plt.xlabel(r'$s$')

plt.tight_layout()
plt.savefig('fig2.png')
plt.show()
