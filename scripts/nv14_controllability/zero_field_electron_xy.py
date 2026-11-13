import numpy as np
I=[[1,0,0],[0,1,0],[0,0,1]]
S_x=(1/np.sqrt(2))*np.array([[0,1,0],[1,0,1],[0,1,0]])
S_y=(1/np.sqrt(2))*np.array([[0,-1j,0],[1j,0,-1j],[0,1j,0]])
S_z=np.array([[1,0,0],[0,0,0],[0,0,-1]])
I_x=(1/np.sqrt(2))*np.array([[0,1,0],[1,0,1],[0,1,0]])
I_y=(1/np.sqrt(2))*np.array([[0,-1j,0],[1j,0,-1j],[0,1j,0]])
I_z=np.array([[1,0,0],[0,0,0],[0,0,-1]])
S_zsquared=np.kron(S_z@S_z,I)
BSz=np.kron(S_z,I)
BIz=np.kron(I,I_z)
QIz2=np.kron(I,I_z@I_z)
SAI1=np.kron(S_x,I_x)
SAI2=np.kron(S_y,I_y)
SAI3=np.kron(S_z,I_z)
# Numerical coefficients (GHz)
D = 2.87028
Q = -0.00494588
A_perp = -0.002635
A_parallel = -0.00216519
# Drift Hamiltonian
H0 = (D*S_zsquared+ Q*QIz2+ A_perp*(SAI1+SAI2)+ A_parallel*SAI3)
#control Hamiltonian
H_x=np.kron(S_x,I)
H_y=np.kron(S_y,I)
l=[H0,H_x,H_y]
new_l=[-1j*x for x in l]
def commutator(A, B):
    return A @ B - B @ A

def inner(A, B):
    return np.trace(A.conj().T @ B).real

# Gram-Schmidt independence check
def add_independent(B, V, M, depth, tol=1e-12):
    R = M.copy()
    original_norm = np.sqrt(inner(R, R))

    if original_norm < tol:
        return False

    for _ in range(2):
        for q in V:
            R -= inner(q, R)*q

    residual_norm = np.sqrt(inner(R, R))
    if residual_norm > tol*original_norm:
        B.append((depth, M.copy()))
        V.append(R/residual_norm)
        return True
    return False

B = []
V = []
for X in new_l:
    add_independent(B, V, X, 0)
depth_zero = [M for depth, M in B if depth == 0]
print("Depth 0:", len(depth_zero))
new_level = depth_zero.copy()
current_depth = 1
max_dimension = 81
# Lie algebra
while new_level:
    next_level = []
    for X in new_level:
        for A in depth_zero:
            C = commutator(X, A)
            if add_independent(B, V, C, current_depth):
                next_level.append(C)
    print(f"Depth {current_depth}: "f"{len(next_level)} new generators, "f"total dimension = {len(B)}")
    if len(B) == max_dimension:
        break
    if len(next_level) == 0:
        break
    new_level = next_level
    current_depth += 1

print("\nFinal Lie algebra dimension:", len(B))

if len(B) == 81:
    print("Full u(9) controllability!")
elif len(B) == 80:
    print("Full su(9) controllability!")
else:
    print("not controllable")
