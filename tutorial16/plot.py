#!/usr/bin/env python3
# Licensed under CC BY 4.0
# A Set of Materials-Science Tutorials for the LAMMPS Simulation Package
# Tutorial 16: self-diffusion in liquid aluminum
#
# Reads the mean-squared-displacement files msd_<T>.dat written by
# diffusion.lmp, fits the diffusive (linear) part of each curve to obtain the
# self-diffusion coefficient D(T) via the Einstein relation, and fits the
# temperature dependence to an Arrhenius law to extract the activation energy.

import glob
import re
import numpy as np
import matplotlib.pyplot as plt

DT = 0.002          # ps per step, must match the timestep in diffusion.lmp
FIT_TMIN = 20.0     # ps; skip the early ballistic/sub-diffusive part
KB = 8.617333e-5    # Boltzmann constant in eV/K
# 1 Ang^2/ps = 1e-8 m^2/s = 1e-4 cm^2/s
A2PS_TO_M2S = 1.0e-8

# --- collect (T, D) from every msd_<T>.dat file ---
temps, diff = [], []
fits = {}
for fname in sorted(glob.glob("msd_*.dat")):
    T = float(re.search(r"msd_(\d+)", fname).group(1))
    step, msd = np.loadtxt(fname, comments="#", unpack=True)
    t = step * DT
    sel = t >= FIT_TMIN
    slope, intercept = np.polyfit(t[sel], msd[sel], 1)
    D = slope / 6.0                      # D = slope / (2*dim), dim = 3
    temps.append(T)
    diff.append(D)
    fits[T] = (t, msd, slope, intercept, D)
    print(f"T = {T:6.0f} K   D = {D:7.4f} Ang^2/ps "
          f"= {D*A2PS_TO_M2S*1e9:5.2f} x 1e-9 m^2/s")

order = np.argsort(temps)
temps = np.array(temps)[order]
diff = np.array(diff)[order]

# --- Figure 1: MSD(t) with the fitted diffusive slope for each temperature ---
fig, ax = plt.subplots(figsize=(5.0, 3.8))
colors = plt.cm.viridis(np.linspace(0.1, 0.85, len(temps)))
for T, c in zip(temps, colors):
    t, msd, slope, intercept, D = fits[T]
    ax.plot(t, msd, color=c, lw=1.2, label=f"{T:.0f} K")
    tf = t[t >= FIT_TMIN]
    ax.plot(tf, slope * tf + intercept, color="k", ls="--", lw=0.8)
ax.set_xlabel("time (ps)")
ax.set_ylabel(r"MSD ($\mathrm{\AA}^2$)")
ax.set_title("Mean-squared displacement, liquid Al")
ax.legend(frameon=False, fontsize=8, loc="upper left")
fig.tight_layout()
fig.savefig("diffusion-msd.png", dpi=200)

# --- Figure 2: Arrhenius plot, ln D vs 1/T, slope = -Ea/kB ---
invT = 1.0 / temps
lnD = np.log(diff)
a_slope, a_int = np.polyfit(invT, lnD, 1)
Ea = -a_slope * KB                       # eV
D0 = np.exp(a_int)
print(f"\nArrhenius fit:  D0 = {D0:.3f} Ang^2/ps, "
      f"Ea = {Ea:.3f} eV ({Ea*96.485:.1f} kJ/mol)")

fig, ax = plt.subplots(figsize=(5.0, 3.8))
ax.plot(1000.0 / temps, diff * A2PS_TO_M2S * 1e9, "o", color="firebrick",
        ms=7, label="LAMMPS")
xfit = np.linspace(invT.min(), invT.max(), 50)
ax.plot(1000.0 * xfit, np.exp(a_slope * xfit + a_int) * A2PS_TO_M2S * 1e9,
        "k--", lw=1.0, label=f"$E_a$ = {Ea:.2f} eV")
ax.set_yscale("log")
ax.set_ylim((diff.min() * A2PS_TO_M2S * 1e9) * 0.9,
            (diff.max() * A2PS_TO_M2S * 1e9) * 1.18)
ax.set_xlabel(r"$1000/T$ (K$^{-1}$)")
ax.set_ylabel(r"$D$ ($10^{-9}\,\mathrm{m^2/s}$)")
ax.set_title("Arrhenius plot of liquid-Al self-diffusion")
ax.legend(frameon=False, fontsize=9)
fig.tight_layout()
fig.savefig("diffusion-arrhenius.png", dpi=200)
print("wrote diffusion-msd.png and diffusion-arrhenius.png")
