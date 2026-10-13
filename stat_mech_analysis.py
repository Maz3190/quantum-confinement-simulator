import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as const

print("Starting statistical mechanics analysis...")

# Physical Constants
h = const.h
k_B = const.k
m_amu = const.u  # 1 amu in kg
L = 1.0          # 1 meter box

def thermal_de_broglie(T, m):
    return h / np.sqrt(2 * np.pi * m * k_B * T)

def q_discrete(T, m, L, max_n=1500):
    n = np.arange(1, max_n + 1)
    E_n = (n**2 * h**2) / (8 * m * L**2)
    terms = np.exp(-E_n / (k_B * T))
    return np.sum(terms)

def q_continuous(T, m, L):
    Lambda = thermal_de_broglie(T, m)
    return L / Lambda

try:
    # Temperature range spanning ultracold to ambient (1 fK to 1 K)
    temperatures = np.logspace(-15, 0, 200)

    q_exact_list = []
    q_approx_list = []
    lambda_over_L_list = []
    percent_errors = []

    print("Computing partition functions across temperature range...")
    for T in temperatures:
        Lambda = thermal_de_broglie(T, m_amu)
        ratio = Lambda / L
        lambda_over_L_list.append(ratio)
        
        q_ex = q_discrete(T, m_amu, L, max_n=1500)
        q_ap = q_continuous(T, m_amu, L)
        
        q_exact_list.append(q_ex)
        q_approx_list.append(q_ap)
        
        err = abs(q_ap - q_ex) / q_ex * 100
        percent_errors.append(err)

    temperatures = np.array(temperatures)
    lambda_over_L_list = np.array(lambda_over_L_list)
    percent_errors = np.array(percent_errors)

    cutoff_mask = lambda_over_L_list >= 0.1
    T_quantum_boundary = temperatures[cutoff_mask][-1] if np.any(cutoff_mask) else 0.0
    print(f"[SUCCESS] Quantum regime threshold detected at T <= {T_quantum_boundary:.3e} K")

    # Plotting
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8), sharex=True)

    ax1.plot(temperatures, q_exact_list, 'b-', label='Exact Discrete Sum ($q_{exact}$)')
    ax1.plot(temperatures, q_approx_list, 'r--', label='Continuous Integral ($q_{approx} = L/\Lambda$)')
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_ylabel('Partition Function ($q_T$)')
    ax1.set_title('Partition Function: Discrete Quantum Sum vs Classical Integral')
    ax1.legend()
    ax1.grid(True, which="both", ls="--")

    ax2.plot(temperatures, percent_errors, 'g-', label='Percentage Error (%)')
    ax2.axvline(x=T_quantum_boundary, color='orange', linestyle=':', label='Threshold $\Lambda/L = 0.1$')
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.set_xlabel('Temperature ($T$) [K]')
    ax2.set_ylabel('Percentage Error (%)')
    ax2.set_title('Classical Approximation Error vs Temperature')
    ax2.legend()
    ax2.grid(True, which="both", ls="--")

    plt.tight_layout()
    plt.savefig("stat_mech_breakdown.png", dpi=300, bbox_inches="tight")
    print("[SUCCESS] Figure saved as stat_mech_breakdown.png!")
    
    plt.show()

except Exception as e:
    print(f"[ERROR OCCURRED]: {e}")