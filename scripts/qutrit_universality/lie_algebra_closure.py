
import sympy as sp

z = sp.symbols('z')
Phi = z**6 + z**3 + 1
def reduce_exact(expr):
    poly = sp.Poly(sp.expand(expr), z, domain=sp.QQ)
    modulus = sp.Poly(Phi, z, domain=sp.QQ)
    return sp.rem(poly, modulus).as_expr()

def reduce_matrix(M):
    return M.applyfunc(reduce_exact)

def zp(k):
    return z**(k % 9)

H = sp.Matrix([
    [zp(3*j*k) for k in range(3)]
    for j in range(3)
]) / sp.sqrt(3)

Hdag = sp.Matrix([
    [zp(-3*j*k) for k in range(3)]
    for j in range(3)
]) / sp.sqrt(3)

T = sp.diag(1, z**8, z)
Tdag = sp.diag(1, z, z**8)

# Group elements
R = reduce_matrix(H*T*Hdag*Tdag)
S = reduce_matrix(H*Tdag*Hdag*T)

# Exact complex conjugation: conjugate(zeta) = zeta^8
def conjugate_exact(expr):
    return reduce_exact(expr.subs(z, z**8))

def dagger_exact(M):
    return M.T.applyfunc(conjugate_exact)

def comm(A, B):
    return reduce_matrix(A*B - B*A)

A = reduce_matrix(R - dagger_exact(R))
B = reduce_matrix(S - dagger_exact(S))

L1 = A
L2 = B
L3 = comm(A, B)
L4 = comm(B, L3)
L5 = comm(L3, A)
L6 = comm(L3, L4)
L7 = comm(L4, A)
L8 = comm(L4, B)

L = [L1, L2, L3, L4, L5, L6, L7, L8]

# Use the first eight row-major entries.
entries = [
    (0,0), (0,1), (0,2),
    (1,0), (1,1), (1,2),
    (2,0), (2,1)
]

M = sp.Matrix([
    [L[j][r,c] for j in range(8)]
    for r,c in entries
])

# Exact determinant
det_M = reduce_exact(M.det())

print("Exact determinant:")
print(sp.factor(det_M))

print("\nIs the determinant identically zero?",
      det_M == 0)
