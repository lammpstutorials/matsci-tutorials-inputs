#!/usr/bin/env python3
# Licensed under CC BY 4.0
# A Set of Materials-Science Tutorials for the LAMMPS Simulation Package
# Tutorial 15, part 2: crystalline fraction vs temperature -> melting point.
#
# velocity.lmp writes fcc_vs_T.dat with three columns: temperature, number of
# fcc-coordinated atoms in the inherent structure, and the total atom count.
# Below T_m the solid grows (fraction -> 1); above T_m it melts (-> 0).  The
# melting point is where the curve crosses one half.

import numpy as np
import matplotlib.pyplot as plt

T, nfcc, ntot = np.loadtxt("fcc_vs_T.dat", unpack=True)
frac = nfcc / ntot
order = np.argsort(T)
T, frac = T[order], frac[order]

# linear interpolation to fraction = 1/2 (the curve decreases with T)
Tm = np.interp(0.5, frac[::-1], T[::-1])

fig, ax = plt.subplots(figsize=(5.0, 3.6))
ax.plot(T, frac, "o-", color="steelblue", lw=2, ms=8)
ax.axhline(0.5, ls=":", color="0.5")
ax.axvline(Tm, ls="--", color="firebrick")
ax.text(Tm + 3, 0.55, rf"$T_m \approx$ {Tm:.0f} K", color="firebrick")
ax.set_xlabel("temperature (K)")
ax.set_ylabel("crystalline fraction")
ax.set_ylim(-0.05, 1.05)
fig.tight_layout()
fig.savefig("coexist.png", dpi=150)
print(f"melting temperature Tm = {Tm:.0f} K")
