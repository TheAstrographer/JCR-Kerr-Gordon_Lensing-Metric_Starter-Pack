import math

EPS = 1e-9
N_TOTAL = int(1e9)
c = 299792.458
beta_fixed = 0.05

# =============================================
# Python
# =============================================
def is_scalar(x):
    return isinstance(x, (int, float))

def to_list(x):
    if is_scalar(x):
        return [float(x)]
    return [float(xi) for xi in x]

def cumulative_trapezoid_pure(y, x):
    """Manual cumulative trapezoidal integration"""
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
# Fixed Isotropic Conformal Sector
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

def mu_lens(z):
    z = float(z)
    return 1 + beta_fixed * z / (1 + z)

# =============================================
# χ-Mitigation (Comoving Distance)
# =============================================
_chi_cache = {}  # Simple cache for repeated calls

def chi_mitigation(z_input):
    """Comoving distance χ(z) - pure Python"""
    z_list = to_list(z_input)
    if not z_list:
        return 0.0
    
    # Sort for integration
    sorted_z = sorted(z_list)
    key = tuple(sorted_z)
    
    if key in _chi_cache:
        chi_sorted = _chi_cache[key]
    else:
        integrand = [c / H_eff(zz) for zz in sorted_z]
        chi_sorted = cumulative_trapezoid_pure(integrand, sorted_z)
        _chi_cache[key] = chi_sorted
    
    # Interpolate back
    if is_scalar(z_input):
        z_val = float(z_input)
        # Simple linear interpolation
        for i in range(len(sorted_z)-1):
            if sorted_z[i] <= z_val <= sorted_z[i+1]:
                t = (z_val - sorted_z[i]) / (sorted_z[i+1] - sorted_z[i])
                return chi_sorted[i] + t * (chi_sorted[i+1] - chi_sorted[i])
        # Edge cases
        if z_val <= sorted_z[0]:
            return chi_sorted[0]
        return chi_sorted[-1]
    
    # Return list for array input
    return [chi_mitigation(z) for z in z_list]  # recursive for simplicity

# =============================================
# χ-Mitigated Anisotropic Torque
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
    A_mitigated = A_base * math.exp(-chi / chi_scale)
    
    lambda_decay = 5.8
    return A_mitigated * omega * math.exp(-z / lambda_decay)

def H_total_eff(z, n_step=None):
    """Full effective Hubble parameter"""
    return H_eff(z) + delta_torque_adaptable(z, n_step)
