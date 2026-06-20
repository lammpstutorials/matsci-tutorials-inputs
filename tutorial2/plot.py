#!/usr/bin/env python3
# Python Script
# Licensed under CC BY 4.0
# A Set of Materials-Science Tutorials for the LAMMPS Simulation Package
# Tutorial 2: read the energy-volume data produced by scan.lmp, fit a
# third-order Birch-Murnaghan equation of state, and plot the result.
#
# scan.lmp writes ev.yaml, a structured file that loads directly with pandas
# (no custom parser needed). Run scan.lmp first (e.g. in LAMMPS-GUI), then:
#     python3 plot.py

import yaml
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# load the structured data written by scan.lmp
with open("ev.yaml") as f:
    df = pd.DataFrame(yaml.safe_load(f))
V = df["volume_per_atom"].to_numpy()
E = df["cohesive_energy"].to_numpy()


# third-order Birch-Murnaghan equation of state
def birch_murnaghan(V, E0, V0, B0, B0p):
    t = (V0 / V) ** (2.0 / 3.0) - 1.0
    return E0 + 9.0 * V0 * B0 / 16.0 * (t**3 * B0p + t**2 * (6.0 - 4.0 * (V0 / V) ** (2.0 / 3.0)))


p0 = [E.min(), V[np.argmin(E)], 0.5, 4.0]   # initial guess (B0 in eV/Ang^3)
popt, _ = curve_fit(birch_murnaghan, V, E, p0=p0)
E0, V0, B0, B0p = popt

a0_eq = (4.0 * V0) ** (1.0 / 3.0)            # fcc: 4 atoms per conventional cell
B0_GPa = B0 * 160.21766208                   # eV/Ang^3 -> GPa
print(f"Equilibrium lattice constant a0 = {a0_eq:.4f} Ang")
print(f"Minimum cohesive energy      E0 = {E0:.4f} eV/atom")
print(f"Bulk modulus                 B0 = {B0_GPa:.1f} GPa")

# plot the energy-volume curve with the fit
vfit = np.linspace(V.min(), V.max(), 300)
plt.figure(figsize=(5.0, 3.8))
plt.plot(V, E, "o", color="steelblue", label="LAMMPS")
plt.plot(vfit, birch_murnaghan(vfit, *popt), "-", color="firebrick",
         label="Birch--Murnaghan fit")
plt.axvline(V0, ls="--", color="gray", lw=1)
plt.xlabel(r"volume per atom ($\mathrm{\AA}^3$)")
plt.ylabel("cohesive energy (eV/atom)")
plt.legend(frameon=False)
plt.tight_layout()
plt.savefig("ev-curve.png", dpi=200)
