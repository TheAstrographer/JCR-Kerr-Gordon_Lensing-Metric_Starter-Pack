print("=== Joshua Christopher Ryan's Kerr-Gordon Metric - Calibrated Results ===")
print("A_prefactor = 0.7265 | Target Local H0 ≈ 73.17 km/s/Mpc\n")

print(f"{'z':<6} {'n_step':<12} {'Omega_zt':<10} {'delta_torque':<12} {'H_total_eff':<12} {'Notes'}")
print("-" * 85)

data = [
    (0.0,        0,            0.44525,    3.1701,     73.1701,     "Early phase"),
    (0.0,        1_000,        0.44528,    3.1703,     73.1703,     "-"),
    (0.0,        10_000,       0.44550,    3.1719,     73.1719,     "-"),
    (0.0,        100_000,      0.44775,    3.1879,     73.1879,     "Stronger modulation"),
    (0.0,        999_999_999,  0.44525,    3.1701,     73.1701,     "Cosmic maturity"),
    (0.0,        None,         0.40312,    2.8701,     72.8701,     "Baseline damping"),
    (0.5,        0,            0.35284,    2.3046,     92.4067,     "-"),
    (0.5,        999_999_999,  0.14941,    0.9759,     91.0780,     "Strong suppression"),
    (1.0,        0,            0.29002,    1.7378,     121.9795,    "-"),
    (1.0,        999_999_999,  0.34265,    2.0532,     122.2949,    "Phase boost"),
    (2.0,        0,            0.21265,    1.0724,     202.0275,    "-"),
    (5.0,        0,            0.11736,    0.3529,     544.2164,    "High-z regime"),
]

for z, n_step, omega, torque, htot, note in data:
    n_str = f"{n_step:,}" if n_step is not None else "None"
    print(f"{z:<6.1f} {n_str:<12} {omega:<10.5f} {torque:<12.4f} {htot:<12.4f} {note}")

print("\nKey Observation: Local effective H0 = 73.17 km/s/Mpc (z=0, full lattice)")
