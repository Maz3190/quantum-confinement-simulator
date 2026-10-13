import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigh
import scipy.constants as const

# Physical Constants (SI Units)
h = const.h
hbar = const.hbar
m_amu = const.u  # 1 atomic mass unit in kg
eV = const.e      # 1 eV in Joules

# 1. Box and Grid Parameters
N = 500       # Number of internal grid points
L = 1.0       # Box length in meters
dx = L / (N + 1)
x = np.linspace(dx, L - dx, N)

# 2. Construct Hamiltonian Matrix (T + V)
k_factor = (hbar**2) / (2 * m_amu * dx**2)
main_diag = 2 * np.ones(N) * k_factor
off_diag = -1 * np.ones(N - 1) * k_factor

T = np.diag(main_diag) + np.diag(off_diag, k=1) + np.diag(off_diag, k=-1)
V = np.zeros((N, N))  # Infinite potential well interior is zero
H = T + V

# 3. Solve the Matrix Eigenvalue Problem (Find all eigenvalues/eigenvectors)
energies, wavefunctions = eigh(H)

# 4. Compare Ground State with Analytical Theory
n_quantum = 1
E_analytical = (n_quantum**2 * h**2) / (8 * m_amu * L**2)
E_numerical = energies[0]
percent_error = abs(E_numerical - E_analytical) / E_analytical * 100

print(f"Numerical Ground State Energy: {E_numerical:.5e} J")
print(f"Analytical Ground State Energy: {E_analytical:.5e} J")
print(f"Precision Error: {percent_error:.6f}%")

# Print first 3 energy levels in eV
energies_eV = energies[:3] / eV
print(f"First 3 Energy Levels (eV): {energies_eV}")

# 5. Plot the First 3 Wavefunctions
plt.figure(figsize=(9, 6))

for n in range(3):
    # Normalize wavefunctions for clean plotting presentation
    psi = wavefunctions[:, n]
    plt.plot(x, psi, label=f"n = {n+1} ($E$ = {energies_eV[n]:.3e} eV)")

plt.title("Particle in a 1D Infinite Potential Well (First 3 States)", fontsize=12)
plt.xlabel("Position ($x$) [m]", fontsize=10)
plt.ylabel("Wavefunction Amplitude ($\psi$)", fontsize=10)
plt.legend()
plt.grid(True)
plt.show()