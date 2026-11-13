
## Script

`lie_algebra_closure.py`

The script defines the qutrit Fourier gate $H_3$ and the diagonal gate $T_3$, then constructs the auxiliary matrices

$$
R = H_3 T_3 H_3^\dagger T_3^\dagger,
$$

$$
S = H_3 T_3^\dagger H_3^\dagger T_3.
$$

It forms the anti-Hermitian parts of these matrices and a collection of nested commutators. SymPy performs the calculation exactly, using the cyclotomic relation for the ninth root of unity to reduce powers and simplify expressions.

The resulting eight candidate matrices are assembled into an $8 \times 8$ matrix, and its determinant is evaluated exactly. The determinant test is used to assess whether the selected matrices are linearly independent.

## Run

From the repository root:

```bash
python scripts/qutrit_universality/lie_algebra_closure.py
```

The script prints the exact determinant and whether it is identically zero. This calculation is part of the broader qutrit universality analysis.
