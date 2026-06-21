import math

EPS = 1e-9
N_TOTAL = int(1e9)
c = 299792.458
beta_fixed = 0.05
A_prefactor = 0.7265

# =============================================
# Pure Python
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
        dx = x[i] - x[i - 1]
        avg_y = (y[i] + y[i - 1]) / 2.0
        cum.append(cum[-1] + avg_y * dx)
    return cum


# =============================================
# Isotropic Conformal Sector
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
    z_list = to_list(z_input)
    if not z_list:
        return 0.0

    # Always integrate from 0
    zs = sorted(set([0.0] + z_list))
    key = tuple(zs)

    if key in _chi_cache:
        chi_sorted = _chi_cache[key]
    else:
        integrand = [c / H_eff(zz) for zz in zs]
        chi_sorted = cumulative_trapezoid_pure(integrand, zs)
        _chi_cache[key] = chi_sorted

    # Interpolation (same as before)
    if is_scalar(z_input):
        z_val = float(z_input)
        for i in range(len(zs)-1):
            if zs[i] <= z_val <= zs[i+1]:
                t = (z_val - zs[i]) / (zs[i+1] - zs[i])
                return chi_sorted[i] + t * (chi_sorted[i+1] - chi_sorted[i])
        return chi_sorted[0] if z_val <= zs[0] else chi_sorted[-1]

    return [chi_mitigation(z) for z in z_list]

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

    A_mitigated = A_base * A_prefactor * math.exp(-chi / chi_scale)

    lambda_decay = 5.8
    return A_mitigated * omega * math.exp(-z / lambda_decay)


def H_total_eff(z, n_step=None):
    """Full effective Hubble parameter"""
    return H_eff(z) + delta_torque_adaptable(z, n_step)


# =============================================
# Main Execution / Demo
# =============================================
if __name__ == "__main__":
    print("=== Joshua Christopher Ryan's Kerr-Gordon Metric (Pure Python) ===")
    print(f"Fixed β = {beta_fixed} | χ-mitigation Active\n")

    h_iso_today = H_eff(0)
    torque_today = delta_torque_adaptable(0, N_TOTAL - 1)
    h_total_today = H_total_eff(0, N_TOTAL - 1)

    print(f"Today (z=0):")
    print(f"   Isotropic H     = {h_iso_today:.2f} km/s/Mpc")
    print(f"   Mitigated Torque = {torque_today:.2f} km/s/Mpc")
    print(f"   Total Effective H0 = {h_total_today:.2f} km/s/Mpc\n")

    # Calibration Table
    print("=== Combined Table (χ-Mitigated Torque) ===")
    regimes = {
        -10.0: (22025, 0),
        -8.0: (2979, 1000),
        -6.0: (402, 20000),
        -4.0: (53.6, 300000),
        -2.0: (6.39, 5e7),
        -1.0: (1.72, 6e8),
        -0.2: (0.221, 9.5e8),
        0.0: (0.0, N_TOTAL - 1)
    }

    print("ln(a)   z        H_iso    Torque   H_total    Regime")
    print("-" * 78)

    regime_names = {
        22025: "Early",
        2979: "UV",
        402: "Winding",
        53.6: "Ramp-up",
        6.39: "Mid-z",
        1.72: "Pantheon+",
        0.221: "SH0ES",
        0.0: "Today"
    }

    for lna, (z_val, n_sample) in regimes.items():
        h_iso = H_eff(z_val)
        torque = delta_torque_adaptable(z_val, n_sample)
        h_tot = H_total_eff(z_val, n_sample)
        regime_name = regime_names.get(z_val, "")
        print(f"{lna:6.1f} {z_val:8.2f} {h_iso:8.2f} {torque:8.2f} {h_tot:8.2f}   {regime_name}")
