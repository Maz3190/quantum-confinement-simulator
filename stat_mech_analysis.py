import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as const

# 1. Physical Constants & System Parameters
h = const.h
kB = const.k
m = const.u  # 1 atomic mass unit in kg
L = 1.0      # 1 meter box length

# Calculate ground state energy E1 = h^2 / (8 m L^2)
E1 = (h**2) / (8 * m * L**2)

# 2. Professor's Suggested Temperature Range in femtoKelvin (fK)
fK_values = np.array([0.1, 1.0, 10.0, 100.0, 1000.0, 10000.0, 100000.0, 1000000.0])
T_values = fK_values * 1e-15  # Convert fK to Kelvin

q_exact_list = []
q_approx_list = []
percent_errors = []

print(f"{'T (fK)':<12} | {'n_char':<10} | {'N_max':<10} | {'q_exact':<15} | {'q_approx':<15} | {'Error (%)':<10}")
print("-" * 80)

for T, fK in zip(T_values, fK_values):
    thermal_energy = kB * T
    
    # Characteristic quantum number where n^2 * E1 = kT
    n_char = np.sqrt(thermal_energy / E1)
    
    # Dynamic Truncation: 10x n_char ensures full Gaussian convergence (minimum N=10)
    N_max = max(10, int(np.ceil(10 * n_char)))
    
    # Exact discrete sum: q_exact = sum_{n=1}^{N_max} exp(-n^2 E1 / kT)
    n_array = np.arange(1, N_max + 1)
    terms = np.exp(-(n_array**2 * E1) / thermal_energy)
    q_exact = np.sum(terms)
    
    # Continuous classical approximation: q_approx = L / Lambda
    Lambda = np.sqrt((h**2) / (2 * np.pi * m * thermal_energy))
    q_approx = L / Lambda
    
    # Calculate percentage error
    error = abs(q_exact - q_approx) / q_approx * 100
    
    q_exact_list.append(q_exact)
    q_approx_list.append(q_approx)
    percent_errors.append(error)
    
    print(f"{fK:<12.1f} | {n_char:<10.2f} | {N_max:<10} | {q_exact:<15.5e} | {q_approx:<15.5e} | {error:<10.2f}%")

# 3. Generate and Save the Updated Breakdown Plot
plt.figure(figsize=(10, 6))
plt.loglog(fK_values, q_exact_list, 'o-', label='Exact Discrete Sum ($q_{exact}$)', linewidth=2, color='blue')
plt.loglog(fK_values, q_approx_list, 's--', label='Classical Integral ($L/\Lambda$)', linewidth=2, color='orange')
plt.xlabel('Temperature $T$ (fK, log scale)', fontsize=12)
plt.ylabel('Translational Partition Function $q_t$ (log scale)', fontsize=12)
plt.title('Quantum vs. Classical Partition Function Breakdown Analysis (Dynamic Truncation)', fontsize=14)
plt.legend()
plt.grid(True, which="both", ls="--")
plt.savefig('stat_mech_breakdown.png', dpi=300)
plt.show()