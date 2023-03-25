import numpy as np
import matplotlib.pyplot as plt
from SatRot import SatRot

psi = 0.1
K1 = [0.8, 0.2, 0.3,-0.7,-0.5]
K2 = [0.3, 0.5,-0.6,-0.2, 0.7]
s = [2, 0, 2, 0, 2]
e = [0, 0.1, 0.2, 0.3, 0.5]

t = np.linspace(0, 6, 256) # time / orbital period
y0 = np.zeros(8)
y0[1],y0[3] = np.sin(psi/2), np.cos(psi/2)

plt.figure(figsize=(5,10))

for i,(K1,K2,s,e) in enumerate(zip(K1,K2,s,e)):
    y0[6] = s
    y = SatRot(y0, 2*np.pi*t, K1, K2, e)
    e1,e2 = y[:,:2].T
    psi = np.arccos(1 - 2*(e1**2 + e2**2))

    plt.subplot(5,1,i+1)
    plt.plot(t, psi)
    plt.axis([t[0], t[-1], 0, 3.2])
    plt.ylabel(r'$\psi$ / rad')
    plt.text(t[0]+0.05, 2.85,
             r'$(K_1,K_2,s,e)=(%g,%g,%g,%g)$'%(K1,K2,s,e))
    if i<4: plt.tick_params(labelbottom=False)

plt.xlabel(r'$\Omega t/(2\pi)$')
plt.tight_layout()
plt.savefig('fig1.eps')
plt.show()
