
import numpy as np
from scipy.linalg import expm

T = 1e-6
N = 198
N_total = N + 2
dt = T / N_total

Omega_max = 2 * np.pi * 10e6

D = 2.87028
carrier_frequency = D
gamma_e = 28.03
B = 0.001

rng = np.random.default_rng()

def trace(A):
    return np.trace(A)

def dagger(A):
    return np.conjugate(A).T

def initial_pulse():
    raw = rng.normal(size=N)
    kernel = np.ones(9) / 9
    smooth = np.convolve(raw, kernel, mode="same")
    t = (np.arange(N) + 0.5) / N
    envelope = np.sin(np.pi * t)**2
    smooth *= envelope
    smooth /= np.max(np.abs(smooth))
    return 0.2 * Omega_max * smooth

def expand_controls(pulses):
    return [np.concatenate(([0.0], pulse, [0.0])) for pulse in pulses]

S_z = np.array([[1, 0, 0], [0, 0, 0], [0, 0, -1]], dtype=complex)

X_plus = (1 / np.sqrt(2)) * np.array([[0, 1, 0], [1, 0, 0], [0, 0, 0]], dtype=complex)
Y_plus = (1 / np.sqrt(2)) * np.array([[0, -1j, 0], [1j, 0, 0], [0, 0, 0]], dtype=complex)
X_minus = (1 / np.sqrt(2)) * np.array([[0, 0, 0], [0, 0, 1], [0, 1, 0]], dtype=complex)
Y_minus = (1 / np.sqrt(2)) * np.array([[0, 0, 0], [0, 0, -1j], [0, 1j, 0]], dtype=complex)

operators = [X_plus, Y_plus, X_minus, Y_minus]

detuning = D - carrier_frequency
H_0 = 2 * np.pi * (detuning * 1e9 * (S_z @ S_z) + gamma_e * B * 1e9 * S_z)

omega = np.exp(2j * np.pi / 3)
F_3 = (1 / np.sqrt(3)) * np.array([[1, 1, 1], [1, omega, omega**2], [1, omega**2, omega]], dtype=complex)

def project_controls(pulses):
    projected = [pulse.copy() for pulse in pulses]
    for i_x, i_y in ((0, 1), (2, 3)):
        amplitude = np.sqrt(projected[i_x]**2 + projected[i_y]**2)
        scale = np.maximum(1.0, amplitude / Omega_max)
        projected[i_x] /= scale
        projected[i_y] /= scale
    return projected

def total_hamiltonian(pulses, i):
    H = H_0.copy()
    for k in range(4):
        H += pulses[k][i] * operators[k]
    return H

def calculate_propagator(pulses):
    full_pulses = expand_controls(pulses)
    U = np.eye(3, dtype=complex)
    for i in range(N_total):
        H_total = total_hamiltonian(full_pulses, i)
        U_i = expm(-1j * H_total * dt)
        U = U_i @ U
    return U

def fidelity(U_target, U):
    overlap = trace(dagger(U_target) @ U)
    return abs(overlap)**2 / 9

def derivative_exponential(H_total, operator, dt):
    X = -1j * H_total * dt
    A = -1j * operator * dt
    zeros_3 = np.zeros((3, 3), dtype=complex)
    M = np.block([[X, A], [zeros_3, X]])
    exponential = expm(M)
    return exponential[:3, 3:]

def compute_gradient(U_target, pulses):
    full_pulses = expand_controls(pulses)
    V = []
    dV = [[] for _ in range(4)]

    for i in range(N_total):
        H_total = total_hamiltonian(full_pulses, i)
        V_i = expm(-1j * H_total * dt)
        V.append(V_i)
        for k in range(4):
            dV[k].append(derivative_exponential(H_total, operators[k], dt))

    U_before = [np.eye(3, dtype=complex)]
    for i in range(N_total):
        U_before.append(V[i] @ U_before[i])

    U_after = [None] * (N_total + 1)
    U_after[N_total] = np.eye(3, dtype=complex)
    for i in range(N_total - 1, -1, -1):
        U_after[i] = U_after[i + 1] @ V[i]

    U = U_before[N_total]
    overlap = trace(dagger(U_target) @ U)
    F = abs(overlap)**2 / 9
    target_dagger = dagger(U_target)
    gradients = [np.zeros(N) for _ in range(4)]

    for j in range(N):
        i = j + 1
        for k in range(4):
            dU = U_after[i + 1] @ dV[k][i] @ U_before[i]
            d_overlap = trace(target_dagger @ dU)
            gradients[k][j] = 2 * np.real(np.conjugate(overlap) * d_overlap) / 9

    return gradients, F

pulses = project_controls([initial_pulse() for _ in range(4)])
U_initial = calculate_propagator(pulses)
F = fidelity(F_3, U_initial)
print(f"Initial fidelity: {F:.10f}")

gradients, F = compute_gradient(F_3, pulses)
directions = [gradient.copy() for gradient in gradients]

step = 0.05
min_step = 1e-8
max_iterations = 300

for iteration in range(1, max_iterations + 1):
    direction_norm = np.sqrt(np.mean(sum(direction**2 for direction in directions)))
    if direction_norm < 1e-30:
        print("Convergence threshold reached.")
        break

    normalized_directions = [direction / direction_norm for direction in directions]

    trial_step = step
    success = False

    while trial_step >= min_step:
        trial_pulses = project_controls([
            pulses[k] + trial_step * Omega_max * normalized_directions[k]
            for k in range(4)
        ])
        U_trial = calculate_propagator(trial_pulses)
        F_trial = fidelity(F_3, U_trial)

        if F_trial > F:
            pulses = trial_pulses
            F = F_trial
            success = True
            step = min(trial_step * 1.2, 0.5)
            break

        trial_step *= 0.5

    if not success:
        directions = [gradient.copy() for gradient in gradients]
        step *= 0.2
        if step < min_step:
            print("Line search could not find an improving step.")
            break
        continue

    new_gradients, F = compute_gradient(F_3, pulses)
    denominator = sum(np.sum(gradient**2) for gradient in gradients)

    if denominator > 1e-30:
        numerator = sum(np.sum(new_gradients[k] * (new_gradients[k] - gradients[k])) for k in range(4))
        beta = max(0.0, numerator / denominator)
    else:
        beta = 0.0

    new_directions = [new_gradients[k] + beta * directions[k] for k in range(4)]

    if sum(np.sum(new_gradients[k] * new_directions[k]) for k in range(4)) <= 0:
        new_directions = [gradient.copy() for gradient in new_gradients]

    gradients = new_gradients
    directions = new_directions

    print(f"Iteration {iteration:04d} | Fidelity: {F:.10f} | Step: {trial_step:.5e}")

    if F >= 0.999:
        print("Target fidelity reached.")
        break

U_final = calculate_propagator(pulses)
F_final = fidelity(F_3, U_final)
full_pulses = expand_controls(pulses)
unitarity_error = np.linalg.norm(dagger(U_final) @ U_final - np.eye(3))

print("\nFinal results")
print(f"Final fidelity: {F_final:.10f}")
print(f"Unitarity error: {unitarity_error:.3e}")



