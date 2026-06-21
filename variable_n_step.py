import math

EPS = 1e-9
N_TOTAL = int(1e9)
c = 299792.458
beta_fixed = 0.05

# =============================================
# Helpers
# =============================================
def is_scalar(x):
    return isinstance(x, (int, float))

def to_list(x):
    if is_scalar(x):
        return [float(x)]
    return [float(xi) for xi in x]

def cumulative_trapezoid_pure(y, x):
    n = len(y)
    if n == 0:
        return [0.0]
    cum = [0.0]
    for i in range(1, n):
        dx = x[i] - x[i-1]
        avg_y = (y[i] + y[i-1]) / 2.0
        cum.append(cum[-1] + avg_y * dx)
    return cum

# =============================================
# Isotropic Sector
# =============================================
def H_base(z):
    z = float(z)
    H0_base = 70.0
    Om = 0.3
    Ol = 0.7
    return H0_base * math.sqrt(Om * (1 + z)**3 + Ol)

def H_eff(z):
    z = float(z)
    return H_base(z) / (1 + beta_fixed * z / (1 + z))

# =============================================
# χ-Mitigation
# =============================================
_chi_cache = {}

def chi_mitigation(z_input):
    z_list = to_list(z_input)
    if not z_list:
        return 0.0
    
    sorted_z = sorted(z_list)
    key = tuple(sorted_z)
    
    if key in _chi_cache:
        chi_sorted = _chi_cache[key]
    else:
        integrand = [c / H_eff(zz) for zz in sorted_z]
        chi_sorted = cumulative_trapezoid_pure(integrand, sorted_z)
        _chi_cache[key] = chi_sorted
    
    if is_scalar(z_input):
        z_val = float(z_input)
        for i in range(len(sorted_z)-1):
            if sorted_z[i] <= z_val <= sorted_z[i+1]:
                t = (z_val - sorted_z[i]) / (sorted_z[i+1] - sorted_z[i])
                return chi_sorted[i] + t * (chi_sorted[i+1] - chi_sorted[i])
        if z_val <= sorted_z[0]:
            return chi_sorted[0]
        return chi_sorted[-1]
    
    return [chi_mitigation(z) for z in z_list]

# =============================================
# Torque Functions
# =============================================
def Omega_zt_adaptable(z, n_step=None):
    z = float(z)
    tau_v = 0.8
    omega_max = 1.13
    damping = 1 - math.exp(-1.0 / (1 + z) / tau_v)
    
    if n_step is not None:
        n_step = float(n_step)
        y_n = min(n_step * EPS, 1.0)
        freq_mod = 1.0 + 5.0 * math.exp(-z / 2.0)
        control_mod = 0.5 * math.sin(2 * math.pi * (y_n * 10))
        phase = 2 * math.pi * y_n * freq_mod + control_mod * (math.pi / 4) + 0.1047
        damping = damping * (0.5 + 0.5 * math.sin(phase))
    else:
        damping = damping * 0.5
        
    return omega_max * damping

def delta_torque_adaptable(z, n_step=None):
    z = float(z)
    omega = Omega_zt_adaptable(z, n_step)
    chi = chi_mitigation(z)
    
    A_base = 9.8
    chi_scale = 4500.0
    lambda_decay = 5.8
    A_mitigated = A_base * math.exp(-chi / chi_scale)
    
    return A_mitigated * omega * math.exp(-z / lambda_decay)

def H_total_eff(z, n_step=None):
    return H_eff(z) + delta_torque_adaptable(z, n_step)

# =============================================
# Generate Table
# =============================================
if __name__ == "__main__":
    print("Code executed successfully.")
    print("Here are the results for H_total_eff(z, n_step), delta_torque, and Omega_zt at several redshifts:")
    print("z      n_step     Omega_zt   delta_torque H_total_eff ")
    print("-----------------------------------------------------------------")
    
    test_cases = [
        (0.0, 0), (0.0, 10), (0.0, 100), (0.0, 1000), (0.0, 10000), 
        (0.0, 100000), (0.0, None),
        (0.5, 0), (0.5, 10), (0.5, 100), (0.5, 1000), (0.5, 10000), 
        (0.5, 100000), (0.5, None),
        (1.0, 0), (1.0, 10), (1.0, 100), (1.0, 1000), (1.0, 10000), 
        (1.0, 100000), (1.0, None),
        (2.0, 0), (2.0, 10), (2.0, 100), (2.0, 1000), (2.0, 10000), 
        (2.0, 100000), (2.0, None),
        (5.0, 0), (5.0, 10), (5.0, 100), (5.0, 1000), (5.0, 10000), 
        (5.0, 100000), (5.0, None)
    ]
    
    for z, n_step in test_cases:
        omega = Omega_zt_adaptable(z, n_step)
        torque = delta_torque_adaptable(z, n_step)
        htot = H_total_eff(z, n_step)
        n_str = str(n_step) if n_step is not None else "None"
        print(f"{z:<6.1f} {n_str:<10} {omega:<10.5f} {torque:<12.5f} {htot:<10.4f}")
