# plot instability region
# reference: T. R. Kand, P. W. Likins and D. A. Levinson
#   "Spacecraft Dynamics" figure 3.5.3

import numpy as np
import matplotlib.pyplot as plt

plt.figure(figsize=(5,5))

K1 = np.linspace(-1,1,100)

# K3 = -1/3 (resonance with eccentricity)
K2 = (1-3*K1)/(3-K1)
plt.plot(K1, K2, ':')

# K3 > 0
plt.fill_between(K1,-1,-K1, color='g', alpha=0.4, hatch='/')

# 1 - K1*K2 + 3*K2 < 0
K2 = 1/(K1-3)
plt.fill_between(K1,-1,K2, color='r', alpha=0.4, hatch='|')

# K1*K2 > 0
plt.fill([0,-1,-1,0], [0,0,-1,-1], color='b', alpha=0.4, hatch='-')
plt.fill([0,1,1,0], [0,0,1,1], color='b', alpha=0.4, hatch='-')

K1 = np.linspace(0,1,100)

# (1 - K1*K2 + 3*K2)^2 + 16*K1*K2 < 0
a = (3-K1)**2
b = 3 + 7*K1
d = 4*np.sqrt(3*K1*(K1+1))
K21 = (-b + d)/a
K22 = (-b - d)/a
plt.fill_between(K1, K22, K21, color='y', alpha=0.4, hatch='o')

plt.xlabel(r'$K_1$')
plt.ylabel(r'$K_2$')
plt.xticks([-1,0,1])
plt.yticks([-1,0,1])
plt.axis([-1,1,-1,1])
plt.tight_layout()
plt.savefig('fig3.png')
plt.show()
