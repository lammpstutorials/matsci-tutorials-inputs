#!/usr/bin/env python3
# Python Script
# Licensed under CC BY 4.0
# A Set of Materials-Science Tutorials for the LAMMPS Simulation Package
# Tutorial 10: read the gamma-line data produced by gsfe.lmp, convert the
# energies to a generalized stacking fault energy per unit area, identify the
# intrinsic (ISF) and unstable (USF) stacking fault energies, and plot the
# gamma-line.
#
# Run gsfe.lmp first (e.g. in LAMMPS-GUI), then:
#     python3 plot.py

import yaml
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

EV_PER_A2_TO_MJM2 = 16021.766   # 1 eV/Ang^2 = 16021.766 mJ/m^2

with open("gsfe.yaml") as f:
    df = pd.DataFrame(yaml.safe_load(f))

shift = df["shift"].to_numpy()
area = df["area"].to_numpy()
pe = df["pe"].to_numpy()

# gamma is the excess energy (relative to the unfaulted slab, shift = 0) per
# unit fault area; there is a single fault plane, so no factor of two
gamma = (pe - pe[0]) / area * EV_PER_A2_TO_MJM2

# the intrinsic stacking fault is the local minimum near one third of the
# period (one Shockley partial); the unstable SF is the maximum before it
period = shift[-1]                      # last point closes one full period
b_partial = period / 3.0                # a0/sqrt(6)
isf_region = (shift > 0.2 * period) & (shift < 0.55 * period)
usf_region = (shift > 0.0) & (shift < b_partial)
i_isf = np.where(isf_region)[0][np.argmin(gamma[isf_region])]
i_usf = np.where(usf_region)[0][np.argmax(gamma[usf_region])]


def refine(i):
    """parabolic interpolation through three points around index i"""
    x0, x1, x2 = shift[i - 1], shift[i], shift[i + 1]
    y0, y1, y2 = gamma[i - 1], gamma[i], gamma[i + 1]
    a = ((y2 - y1) / (x2 - x1) - (y1 - y0) / (x1 - x0)) / (x2 - x0)
    b = (y1 - y0) / (x1 - x0) - a * (x0 + x1)
    xv = -b / (2 * a)
    yv = y1 + a * (xv - x1) ** 2 + (2 * a * x1 + b) * (xv - x1)
    return xv, yv


x_isf, g_isf = refine(i_isf)
x_usf, g_usf = refine(i_usf)

print(f"period (a0*sqrt(6)/2)      = {period:.3f} Ang")
print(f"Shockley partial b_p       = {b_partial:.3f} Ang")
print(f"intrinsic SF energy  (ISF) = {g_isf:6.1f} mJ/m^2  at shift {x_isf:.3f} Ang")
print(f"unstable  SF energy  (USF) = {g_usf:6.1f} mJ/m^2  at shift {x_usf:.3f} Ang")

# plot the gamma-line
plt.figure(figsize=(5.0, 3.8))
plt.plot(shift, gamma, "o-", color="steelblue", ms=4, label=r"$\gamma$-line")
plt.plot(x_isf, g_isf, "v", color="seagreen", ms=9,
         label=f"ISF = {g_isf:.0f} mJ/m$^2$")
plt.plot(x_usf, g_usf, "^", color="firebrick", ms=9,
         label=f"USF = {g_usf:.0f} mJ/m$^2$")
plt.axvline(b_partial, ls="--", color="gray", lw=1)
plt.text(b_partial, plt.ylim()[1] * 0.05, r"  $b_p=a_0/\sqrt{6}$",
         color="gray", fontsize=8, ha="left")
plt.xlabel(r"displacement along $[11\bar 2]$ ($\mathrm{\AA}$)")
plt.ylabel(r"stacking fault energy $\gamma$ (mJ/m$^2$)")
plt.legend(frameon=True, facecolor="white", edgecolor="black", framealpha=1.0,
           fontsize=9)
plt.tight_layout()
plt.savefig("gsfe-curve.png", dpi=200)
