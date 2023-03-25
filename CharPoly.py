import sympy as sp

# derive eq.(3.5.81) of Kane Likins Levinson
K1,K2 = sp.symbols('K1,K2')
half = sp.S(1)/2
A = sp.Matrix([
        [0, 1, half, 0],
        [-1, 0, 0, half],
        [0, 0, 0, K1],
        [0, -6*K2, K2, 0]])

print(A.charpoly())

# derive eq.(3.2.52) of Kane Likins Levinson
# W = \dot\phi
# Q = \dot\phi - (1-K1)*s
# M = K1*mu/r^3
W,Q,M = sp.symbols('Omega, Q, mu')
A = sp.Matrix([
        [0, W, half, 0],
        [-W, 0, 0, half],
        [0, 0, 0, Q],
        [0, 6*M, -Q, 0]])

print(A.charpoly())
