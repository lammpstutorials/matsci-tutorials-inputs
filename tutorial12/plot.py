#!/usr/bin/env python3
# Python Script
# Licensed under CC BY 4.0
# A Set of Materials-Science Tutorials for the LAMMPS Simulation Package
# Tutorial 12: plot the heating/cooling hysteresis of aluminum from the
# volume-per-atom data written by hysteresis.lmp (part 1).  The volume jumps
# at melting on heating and would jump back at freezing on cooling, but the
# defect-free crystal superheats and the melt supercools, so the two branches
# do not coincide -- the loop brackets the true melting point.
#
# Run hysteresis.lmp first, then:  python3 plot.py

import numpy as np
import matplotlib.pyplot as plt

heat = np.loadtxt("heat.dat")          # columns: step  temperature  vol/atom
cool = np.loadtxt("cool.dat")
Th, Vh = heat[:, 1], heat[:, 2]
Tc, Vc = cool[:, 1], cool[:, 2]

# the melting jump on heating is the single largest step in volume (the data
# are in time order, i.e. increasing temperature, so consecutive differences
# isolate the discontinuity)
i = np.argmax(np.diff(Vh))
T_melt = 0.5 * (Th[i] + Th[i + 1])
print(f"heating: volume jumps (melts) near {T_melt:.0f} K")

# thermal expansion, for free: the gentle slope of the solid branch IS the
# thermal expansion of the crystal.  Convert volume/atom to the lattice
# constant a = (4 V)^(1/3) and fit a straight line between 300 K and 900 K;
# the slope over a gives the linear expansion coefficient alpha.
sol = (Th > 300.0) & (Th < 900.0)
aT = (4.0 * Vh[sol]) ** (1.0 / 3.0)
slope, icpt = np.polyfit(Th[sol], aT, 1)
a300 = slope * 300.0 + icpt
alpha = slope / a300
print(f"solid branch: a(300 K) = {a300:.4f} A, "
      f"alpha = {alpha*1e6:.1f}e-6 /K (experiment: 23.1e-6 /K)")
print("cooling: the melt supercools strongly; in a small, defect-free cell it may")
print("         not recrystallize at all.  Both transitions are nucleation-limited")
print("         and only bracket Tm -- part 2 (two-phase) gives the true value.")

plt.figure(figsize=(5.0, 3.8))
plt.plot(Th, Vh, ".", ms=3, color="firebrick", label="heating")
plt.plot(Tc, Vc, ".", ms=3, color="steelblue", label="cooling")
Tf = np.linspace(300.0, 900.0, 20)
plt.plot(Tf, (slope * Tf + icpt) ** 3 / 4.0, "k--", lw=1.0,
         label="thermal expansion fit")
plt.axvline(T_melt, ls="--", color="gray", lw=1)
plt.text(T_melt, Vh.min() + 0.05 * (Vh.max() - Vh.min()),
         f"  melt ~{T_melt:.0f} K", color="gray", fontsize=8)
plt.xlabel("temperature (K)")
plt.ylabel(r"volume per atom ($\mathrm{\AA}^3$)")
plt.legend(frameon=True, facecolor="white", edgecolor="black", framealpha=1.0)
plt.tight_layout()
plt.savefig("hysteresis.png", dpi=200)
