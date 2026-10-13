import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigh
import scipy.constants as const

# Physical Constants
h = const.h
hbar = const.hbar
m_amu = const.u
L = 1.0

# Exact analytical ground state energy for comparison
E_analytical_ground = (1**2 * h**2) / (8 * m_amu * L**2)

# Test a range of grid resolutions (N)
N_values = [50, 100, 200, 500, 1000, 1500, 2000]
percent_errors = []

print("Running grid convergence study...")
for N in N_values:
    dx = L / (N + 1)
    x = np.linspace(dx, L - dx, N)
    
    # Hamiltonian matrix construction
    k_factor = (hbar**2) / (2 * m_amu * dx**2)
    main_diag = 2 * np.ones(N) * k_factor
    off_diag = -1 * np.ones(N - 1) * k_factor
    
    T = np.diag(main_diag) + np.diag(off_diag, k=1) + np.diag(off_diag, k=-1)
    V = np.zeros((N, N))
    H = T + V
    
    # Solve eigenvalues
    energies, _ = eigh(H)
    E_numerical = energies[0]
    
    # Compute error percentage
    err = abs(E_numerical - E_analytical_ground) / E_analytical_ground * 100
    percent_errors.append(err)
    print(f"N = {N:4d} | Error = {err:.6f}%")

# Plot Convergence Curve
plt.figure(figsize=(8, 5))
plt.plot(N_values, percent_errors, 'bo-', linewidth=2, markersize=6)
plt.yscale('log')
plt.title("Finite Difference Method: Grid Convergence Analysis", fontsize=12)
plt.xlabel("Number of Internal Grid Points ($N$)", fontsize=10)
plt.ylabel("Ground State Energy Error (%) [Log Scale]", fontsize=10)
plt.grid(True, which="both", ls="--")

# Save high-res figure for portfolio
plt.savefig("grid_convergence_error.png", dpi=300, bbox_inches="tight")
plt.show()