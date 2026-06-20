#!/usr/bin/env python3
# Licensed under CC BY 4.0
# A Set of Materials-Science Tutorials for the LAMMPS Simulation Package
# Tutorial 9: phonon dispersion from MD -- vibrational density of states
#
# Reads the velocity autocorrelation function vacf.dat written by
# vacfdos.lmp, Fourier-transforms it into the vibrational density of states
# g(nu), and integrates g(nu) with the harmonic-oscillator weight to obtain
# the quantum heat capacity Cv(T).  Classical MD itself always gives the
# Dulong-Petit value Cv = 3R (equipartition); the quantum curve below comes
# from the measured spectrum, not from the classical dynamics.

import numpy as np
import matplotlib.pyplot as plt

DT = 0.001          # ps per step
H = 4.135668e-3     # Planck constant (eV/THz = eV ps)
KB = 8.617333e-5    # Boltzmann constant (eV/K)
R = 8.314463        # gas constant (J/mol/K)

step, vacf = np.loadtxt("vacf.dat", comments="#", unpack=True)
t = (step - step[0]) * DT
c = vacf / vacf[0]                       # normalized VACF

# density of states: cosine transform of the VACF with a Hann window to
# suppress truncation ripples;  nu in THz
window = 0.5 * (1.0 + np.cos(np.pi * t / t[-1]))
nu = np.linspace(0.0, 12.0, 600)
g = np.array([np.trapezoid(c * window * np.cos(2 * np.pi * f * t), t)
              for f in nu])
g = np.clip(g, 0.0, None)
g /= np.trapezoid(g, nu)                 # normalize: integral g(nu) dnu = 1

# quantum harmonic heat capacity from the spectrum:
#   Cv(T) = 3R * Int g(nu) x^2 e^x / (e^x - 1)^2 dnu,  x = h nu / kB T
Ts = np.linspace(5.0, 500.0, 100)
cv = []
for T in Ts:
    x = H * nu[1:] / (KB * T)
    w = x**2 * np.exp(x) / np.expm1(x) ** 2
    cv.append(3.0 * R * np.trapezoid(g[1:] * w, nu[1:]))
cv = np.array(cv)
i300 = np.argmin(np.abs(Ts - 300.0))
print(f"Cv(300 K) = {cv[i300]:.1f} J/(mol K)   "
      f"(Dulong-Petit 3R = {3*R:.1f}, experiment Cu ~ 23.8)")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8.0, 3.4))
ax1.plot(nu, g, color="firebrick", lw=1.4)
ax1.set_xlabel(r"frequency $\nu$ (THz)")
ax1.set_ylabel(r"$g(\nu)$ (1/THz)")
ax1.set_title("Vibrational density of states")
ax1.set_xlim(0, 10)
ax2.plot(Ts, cv, color="steelblue", lw=1.4, label="quantum, from $g(\\nu)$")
ax2.axhline(3 * R, color="gray", ls="--", lw=1.0,
            label="classical MD (Dulong-Petit)")
ax2.set_xlabel("temperature (K)")
ax2.set_ylabel(r"$C_V$ (J mol$^{-1}$ K$^{-1}$)")
ax2.set_title("Harmonic heat capacity")
ax2.legend(frameon=False, fontsize=9)
ax2.set_ylim(0, 27)
fig.tight_layout()
fig.savefig("dos-cv.png", dpi=200)
print("wrote dos-cv.png")
