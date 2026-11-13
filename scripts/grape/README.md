
# GRAPE-Based Quantum Gate Optimization

## Overview

This directory contains two implementations of Gradient Ascent Pulse Engineering (GRAPE) for quantum gate optimization in the three-level electronic spin system of a nitrogen-vacancy (NV) centre in diamond.

The objective is to optimize time-dependent control pulses that implement two target qutrit gates:

- The qutrit Fourier transform, $H_3$.
- The diagonal qutrit gate, $T_3$.

Both implementations use four piecewise-constant control channels corresponding to the two quadratures of the two electronic-spin transitions.

## Physical Model

The controlled system is described in a rotating frame, with the rotating-wave approximation (RWA).

The drift Hamiltonian is

$$
H_0 = 2\pi\left[
(D-\omega_{\mathrm{rf}})S_z^2
+\gamma_e B S_z
\right],
$$

where:

- $D$ is the zero-field splitting.
- $\omega_{\mathrm{rf}}$ is the carrier frequency.
- $\gamma_e$ is the electronic gyromagnetic ratio.
- $B$ is the applied magnetic field.
- $S_z$ is the spin-1 $z$-component operator.

The control Hamiltonian is

$$
\begin{aligned}
H_c(t) ={}&
u_{x,+}(t)X_+ + u_{y,+}(t)Y_+\\
&+u_{x,-}(t)X_- + u_{y,-}(t)Y_-.
\end{aligned}
$$

The four control channels are $X_+$, $Y_+$, $X_-$, and $Y_-$, corresponding to the two quadratures of the two allowed transitions.

## Optimization Method

The total gate duration is divided into short time intervals, with constant control amplitudes within each interval.

The total propagator is calculated as

$$
U(T)=U_N U_{N-1}\cdots U_1,
$$

where each time-slice propagator is

$$
U_k=\exp(-iH_k\Delta t).
$$

The objective function is the gate fidelity

$$
F=\frac{1}{9}
\left|\operatorname{Tr}
\left(U_{\mathrm{target}}^\dagger U(T)\right)
\right|^2.
$$

The optimization includes the following components:

- Random initialization of the control pulses.
- Smoothing of the initial pulses using a moving-average filter.
- A $\sin^2(\pi t)$ temporal envelope.
- Amplitude constraints on both transition-control pairs.
- Analytical gradient calculations using a block-matrix exponential.
- A Polak–Ribière nonlinear conjugate-gradient method.
- Adaptive step-size selection with a backtracking line search.

For each transition, the amplitude constraint is

$$
\sqrt{u_{x,\pm}^2+u_{y,\pm}^2}
\leq \Omega_{\max}.
$$

The optimization terminates when the target fidelity reaches the specified threshold or when the line search can no longer find an improving step.

## Scripts

### 1. Fourier Gate Optimization

`fourier_gate.py`

Optimizes the four control channels to synthesize the qutrit Fourier transform $H_3$.

Run from the repository root:

```bash
python scripts/grape/fourier_gate.py
```

### 2. Diagonal Gate Optimization

`t_gate.py`

Uses the same GRAPE framework to optimize the control pulses for the diagonal qutrit gate $T_3$.

Run from the repository root:

```bash
python scripts/grape/t_gate.py
```

## Implementation Notes

The two scripts use independent random initial pulses. No fixed random seed is specified in the current implementations, so separate executions may produce different optimization trajectories.

The final output includes the achieved gate fidelity and a numerical check of the unitarity of the resulting propagator.

## Dependencies

- NumPy
- SciPy

Install the required packages from the repository root:

```bash
pip install -r requirements.txt
```

## Future Development

An extension of this implementation to open quantum systems is planned. The intended next step is to incorporate environmental effects, including relaxation and dephasing, through the Lindblad master equation and investigate their influence on quantum gate optimization.
