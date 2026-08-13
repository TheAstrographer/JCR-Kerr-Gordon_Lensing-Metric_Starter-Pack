# ============================================================
# Cosmological Epoch Numerical Value Tables
# ============================================================

import math

# ------------------------------------------------------------
# Cosmological parameters (Planck 2018 approximate)
# ------------------------------------------------------------
H0 = 67.66          # km/s/Mpc
Om = 0.3111         # matter density
Ol = 1.0 - Om       # dark energy density (flat)

# Convert H0 to 1/Gyr
# 1 Mpc / (km/s) ≈ 977.8 Gyr
H0_Gyr = H0 / 977.8

# ------------------------------------------------------------
# Lookback time and age
# ------------------------------------------------------------
def E(z):
    """Hubble function E(z) = H(z)/H0"""
    return math.sqrt(Om * (1 + z)**3 + Ol)

def integrate(f, a, b, n=1000):
    """Simple trapezoidal integration"""
    if a == b:
        return 0.0
    h = (b - a) / n
    s = 0.5 * (f(a) + f(b))
    for i in range(1, n):
        s += f(a + i * h)
    return s * h

def lookback_time(z):
    """Lookback time in Gyr from z=0 to redshift z"""
    def integrand(zp):
        return 1.0 / ((1 + zp) * E(zp))
    return integrate(integrand, 0.0, z) / H0_Gyr

def age_of_universe(z):
    """Age of the Universe at redshift z (Gyr)"""
    z_max = 5000.0
    def integrand(zp):
        return 1.0 / ((1 + zp) * E(zp))
    return integrate(integrand, z, z_max) / H0_Gyr

# ------------------------------------------------------------
# Epoch definitions
# ------------------------------------------------------------
epochs = [
    ("Local Universe (SH0ES regime)",           0.00,   0.15),
    ("Matter-Dark Energy Equality",             0.30,   0.40),
    ("Reionization / First Galaxies",           6.0,   15.0),
    ("Dark Ages (approximate)",                20.0, 1100.0),
    ("Recombination / CMB last scattering",  1080.0, 1100.0),
]

# ------------------------------------------------------------
# Print Table 1: Redshift Ranges
# ------------------------------------------------------------
print("=" * 70)
print("1. Redshift Ranges (as displayed in the chart)")
print("=" * 70)
print(f"{'Epoch':<45} {'z_min':>10} {'z_max':>10}")
print("-" * 70)
for name, zmin, zmax in epochs:
    print(f"{name:<45} {zmin:10.2f} {zmax:10.2f}")

# ------------------------------------------------------------
# Print Table 2: Cosmic Times
# ------------------------------------------------------------
print()
print("=" * 95)
print("2. Corresponding Cosmic Times (Planck 2018 inspired parameters)")
print("=" * 95)
print(f"{'Epoch':<40} {'t_lb(zmin)':>12} {'t_lb(zmax)':>12} {'Age(zmin)':>12}")
print(f"{'':40} {'(Gyr)':>12} {'(Gyr)':>12} {'(Gyr)':>12}")
print("-" * 95)

for name, zmin, zmax in epochs:
    lb_min = lookback_time(zmin)
    lb_max = lookback_time(zmax)
    age_min = age_of_universe(zmin)
    print(f"{name:<40} {lb_min:12.2f} {lb_max:12.2f} {age_min:12.2f}")

print()
print("Notes:")
print("  - H0 = 67.66 km/s/Mpc, Om = 0.3111, flat LambdaCDM")
print("  - Lookback time = time from today back to that redshift")
print("  - Age = age of the Universe at that redshift")
print("  - Pure Python trapezoidal integration (no external libraries)")
print("=" * 95)
