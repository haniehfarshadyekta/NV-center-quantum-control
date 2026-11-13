# NV-Center Quantum Control and Controllability

**Ongoing BSc research project**  
**Researcher:** Hanieh Farshad Yekta  
**Supervisor:** Dr. Maryam Mirkamali  
**Programme:** BSc in Physics, Sharif University of Technology

This repository contains Python implementations developed as part of an ongoing research project on qutrit-gate universality, numerical quantum gate synthesis, and Lie-algebraic controllability of NV-centre spin systems. The current code is organized into three connected research areas.

> **Project status:** These are research scripts under active development, not a released or fully validated software package. The open-quantum-system extension using Lindblad dynamics is a planned next step and is not implemented in the current repository.

## Research overview

### 1. Qutrit-gate universality

The qutrit-universality calculation studies a pair of three-level gates, the qutrit Fourier transform $H_3$ and a diagonal gate $T_3$. From these gates, the code constructs auxiliary group elements and forms matrices from their anti-Hermitian parts and nested commutators. SymPy is used to carry out the algebra in exact arithmetic, reducing expressions modulo the cyclotomic relation for a ninth root of unity. The script then evaluates an exact determinant to test the linear independence of eight candidate traceless matrices. This calculation supports the Lie-algebraic part of the universality analysis.

### 2. GRAPE-based quantum gate synthesis

The two GRAPE scripts optimize control pulses for two target qutrit gates: the Fourier gate $H_3$ and the diagonal gate $T_3$. The controlled system is the three-level electronic spin of an NV centre.Each transition is driven through two quadratures, giving four piecewise-constant control channels: $X_{+}$, $Y_{+}$, $X_{-}$, and $Y_{-}$.
The implementation includes:

- propagation through a sequence of piecewise-constant Hamiltonians;
- a rotating-frame Hamiltonian and the rotating-wave approximation;
- amplitude projection for the two transition-control pairs;
- gate fidelity based on the trace overlap with the target unitary;
- analytical derivatives of matrix exponentials using a block-matrix exponential;
- a Polak–Ribière nonlinear conjugate-gradient update and an adaptive line search.

The initial control pulses are generated randomly and smoothed before optimization. No fixed random seed is set in the current scripts, so separate runs may follow different optimization trajectories.

### 3. Lie-algebraic controllability of the coupled NV–¹⁴N system

The controllability scripts model the NV electronic spin and the nitrogen-14 nuclear spin as two spin-1 subsystems. Their tensor-product Hilbert space has dimension nine. The drift Hamiltonian includes the electronic zero-field splitting, nuclear quadrupole term, and transverse and longitudinal hyperfine couplings. The finite-field scripts additionally include longitudinal Zeeman terms.

For each control configuration, the code generates nested commutators of the drift and control generators. It uses the Hilbert–Schmidt inner product and a re-orthogonalized Gram–Schmidt procedure to identify linearly independent matrices and estimate the dimension of the generated dynamical Lie algebra.

The four scripts compare:

- zero magnetic field with electronic $S_x$ and $S_y$ controls;
- zero magnetic field with electronic $S_x$ and nuclear $I_x$ controls;
- a finite longitudinal magnetic field with electronic $S_x$ control, including the nuclear Zeeman term;
- the same finite-field control configuration with the nuclear Zeeman term neglected.

## Repository structure

```text
NV-center-quantum-control/
├── README.md
├── ORIGINAL_FILES.md
├── requirements.txt
├── .gitignore
└── scripts/
    ├── qutrit_universality/
    │   ├── README.md
    │   └── lie_algebra_closure.py
    ├── grape/
    │   ├── README.md
    │   ├── fourier_gate.py
    │   └── t_gate.py
    └── nv14_controllability/
        ├── README.md
        ├── field_electron_x.py
        ├── field_electron_x_no_nuclear_zeeman.py
        ├── zero_field_electron_x_nuclear_x.py
        └── zero_field_electron_xy.py
```

## Script guide

| Script | Research question / calculation |
|---|---|
| `scripts/qutrit_universality/lie_algebra_closure.py` | Exact symbolic construction of qutrit group elements and a determinant-based linear-independence check for nested-commutator matrices. |
| `scripts/grape/fourier_gate.py` | GRAPE optimization with the qutrit Fourier transform as the target gate. |
| `scripts/grape/t_gate.py` | GRAPE optimization with the diagonal qutrit $T_3$ gate as the target. |
| `scripts/nv14_controllability/zero_field_electron_xy.py` | Coupled  NV–¹⁴N system at zero field with electronic $S_x,S_y$ controls. |
| `scripts/nv14_controllability/zero_field_electron_x_nuclear_x.py` | Coupled system at zero field with electronic $S_x$ and nuclear $I_x$ controls. |
| `scripts/nv14_controllability/field_electron_x.py` | Finite longitudinal field, electronic $S_x$ control, and both electronic and nuclear Zeeman terms. |
| `scripts/nv14_controllability/field_electron_x_no_nuclear_zeeman.py` | Same finite-field configuration with the nuclear Zeeman contribution set to zero. |

## Installation

The scripts use NumPy, SciPy, and SymPy. From the repository root, create and activate a virtual environment, then install the dependencies.

**Windows (PowerShell):**

```powershell
python -m venv .venv
.venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
```

**macOS / Linux:**

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

## Running the scripts

Run a script from the repository root. For example:

```bash
python scripts/qutrit_universality/lie_algebra_closure.py
python scripts/grape/fourier_gate.py
python scripts/grape/t_gate.py
python scripts/nv14_controllability/zero_field_electron_xy.py
python scripts/nv14_controllability/zero_field_electron_x_nuclear_x.py
python scripts/nv14_controllability/field_electron_x.py
python scripts/nv14_controllability/field_electron_x_no_nuclear_zeeman.py
```

The GRAPE scripts perform iterative numerical optimizations and may take longer to run than the Lie-algebra scripts.

## Numerical and reproducibility notes

- The GRAPE scripts initialize their pulse arrays using `numpy.random.default_rng()` without a fixed seed. Results and optimization trajectories may therefore differ between runs. A seed should be set explicitly when reproducible pulse initializations are required.
- The controllability scripts use a numerical tolerance (`1e-12`) when testing linear independence. The computed dimension can depend on numerical conditioning and tolerance, so it should be checked for stability when parameters or precision are changed.
- The parameter values and units are defined in the individual scripts. Check them before adapting the model to a different sample or experimental configuration.
- The scripts currently run as standalone files; they have not yet been refactored into a shared Python package or accompanied by a formal automated test suite.
- Before treating numerical outputs as final research results, validate the Hamiltonian conventions and units, check matrix properties, and compare analytical GRAPE gradients with finite-difference estimates.

## Planned development

A planned extension is to incorporate environmental effects into the GRAPE optimization by evolving density matrices with a Lindblad master equation. This would allow relaxation and dephasing to be included directly in the pulse optimization. **No Lindblad-based OQS implementation is included yet.**

