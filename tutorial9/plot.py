#!/usr/bin/env python3
# Python Script
# Licensed under CC BY 4.0
# A Set of Materials-Science Tutorials for the LAMMPS Simulation Package
# Tutorial 9: read the stress-strain data produced by elastic.lmp, assemble
# the 6x6 elastic stiffness matrix C_ij, derive the cubic elastic constants
# and the aggregate (Voigt-Reuss-Hill) moduli, and plot the three
# characteristic stress-strain responses whose slopes are C11, C12, and C44.
#
# Run elastic.lmp first (e.g. in LAMMPS-GUI), then:
#     python3 plot.py

import yaml
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

BAR_TO_GPA = 1.0e-4   # 1 GPa = 10000 bar

# load the structured data written by elastic.lmp
with open("scan.yaml") as f:
    df = pd.DataFrame(yaml.safe_load(f))

# Voigt order of the stress (= -pressure) components in the YAML records
comps = ["pxx", "pyy", "pzz", "pyz", "pxz", "pxy"]

# C_ij is the slope of stress_i versus strain_j: deform in Voigt direction j
# (column), read the six stress components (rows).  A straight-line fit over
# the small-strain scan gives the slope; -pressure converts to stress and the
# factor BAR_TO_GPA converts bar to GPa.
C = np.zeros((6, 6))
for j in range(6):
    sub = df[df["dir"] == j + 1]
    eps = sub["eps"].to_numpy()
    for i, comp in enumerate(comps):
        sigma = -sub[comp].to_numpy()
        C[i, j] = np.polyfit(eps, sigma, 1)[0] * BAR_TO_GPA

np.set_printoptions(precision=1, suppress=True)
print("Elastic stiffness matrix C_ij [GPa]:")
print(C)

# cubic-averaged independent constants
C11 = (C[0, 0] + C[1, 1] + C[2, 2]) / 3.0
C12 = (C[0, 1] + C[0, 2] + C[1, 2]) / 3.0
C44 = (C[3, 3] + C[4, 4] + C[5, 5]) / 3.0

# aggregate moduli (Voigt-Reuss-Hill averages for a cubic crystal)
B = (C11 + 2.0 * C12) / 3.0                       # bulk modulus
G_V = (C11 - C12 + 3.0 * C44) / 5.0               # Voigt shear
G_R = 5.0 * (C11 - C12) * C44 / (4.0 * C44 + 3.0 * (C11 - C12))   # Reuss shear
G = 0.5 * (G_V + G_R)                             # Hill average
E = 9.0 * B * G / (3.0 * B + G)                   # Young's modulus
nu = (3.0 * B - 2.0 * G) / (2.0 * (3.0 * B + G))  # Poisson ratio
A = 2.0 * C44 / (C11 - C12)                       # Zener anisotropy

print(f"\nCubic elastic constants (this work / SW analytic / experiment):")
print(f"  C11 = {C11:6.1f}  / 151.4 / 165.6 GPa")
print(f"  C12 = {C12:6.1f}  /  76.4 /  63.9 GPa")
print(f"  C44 = {C44:6.1f}  /  56.4 /  79.6 GPa")
print(f"\nAggregate (Voigt-Reuss-Hill) moduli:")
print(f"  bulk modulus    B = {B:6.1f} GPa")
print(f"  shear modulus   G = {G:6.1f} GPa")
print(f"  Young's modulus E = {E:6.1f} GPa")
print(f"  Poisson ratio   v = {nu:6.3f}")
print(f"  Zener anisotropy A = {A:5.2f}")

# ---------------------------------------------------------------------------
# Figure: the three independent constants are the slopes of three responses.
#   C11 : sigma_xx vs eps    when straining x   (direction 1, row pxx)
#   C12 : sigma_yy vs eps    when straining x   (direction 1, row pyy)
#   C44 : sigma_yz vs eps    when shearing yz   (direction 4, row pyz)
# ---------------------------------------------------------------------------
d1 = df[df["dir"] == 1].sort_values("eps")
d4 = df[df["dir"] == 4].sort_values("eps")
e1 = d1["eps"].to_numpy() * 100.0          # percent strain
e4 = d4["eps"].to_numpy() * 100.0
responses = [
    (e1, -d1["pxx"].to_numpy() * BAR_TO_GPA, "C$_{11}$: $\\sigma_{xx}$ vs $\\epsilon_{xx}$", C11, "steelblue", "o"),
    (e1, -d1["pyy"].to_numpy() * BAR_TO_GPA, "C$_{12}$: $\\sigma_{yy}$ vs $\\epsilon_{xx}$", C12, "seagreen", "s"),
    (e4, -d4["pyz"].to_numpy() * BAR_TO_GPA, "C$_{44}$: $\\sigma_{yz}$ vs $\\epsilon_{yz}$", C44, "firebrick", "^"),
]

plt.figure(figsize=(5.0, 3.8))
for x, y, label, slope, color, marker in responses:
    plt.plot(x, y, marker, color=color, ms=5)
    xx = np.linspace(x.min(), x.max(), 50)
    plt.plot(xx, slope * xx / 100.0, "-", color=color,
             label=f"{label}  ({slope:.0f} GPa)")
plt.axhline(0.0, color="gray", lw=0.6)
plt.axvline(0.0, color="gray", lw=0.6)
plt.xlabel("strain (%)")
plt.ylabel("stress (GPa)")
plt.legend(frameon=False, fontsize=8)
plt.tight_layout()
plt.savefig("elastic-stress-strain.png", dpi=200)
