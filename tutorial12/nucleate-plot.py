#!/usr/bin/env python3
# Python Script
# Licensed under CC BY 4.0
# A Set of Materials-Science Tutorials for the LAMMPS Simulation Package
# Tutorial 12 (advanced): plot the seeded heating/cooling run from
# nucleate.lmp.  The seeded melt refreezes on cooling (volume drops back to
# the crystal branch and the fcc-atom count jumps up), unlike the unseeded
# run -- but only at deep undercooling, well below Tm.
#
# Run hysteresis.lmp and nucleate.lmp first, then:  python3 nucleate-plot.py

import os
import numpy as np
import matplotlib.pyplot as plt

# columns: step  temperature  vol/atom  N(fcc)  dev_avg  dev_max
h = np.loadtxt("seeded-heat.dat")[1:]          # drop the first averaging block
c = np.loadtxt("seeded-cool.dat")[1:]          # (it straddles the fix switch)

# The heating branch and the unseeded cooling come from the part-1 run
# (hysteresis.lmp -> heat.dat, cool.dat); the seeded cooling comes from this
# run.  Heating is essentially the same with or without the small seed.
hu = np.loadtxt("heat.dat")            # unseeded heating: step temperature vol/atom
u  = np.loadtxt("cool.dat")            # unseeded cooling

# vertical guide lines: melting on heating (largest volume rise, unseeded run)
# and the seeded crystallization on cooling (largest drop)
im = int(np.argmax(np.diff(hu[:, 2])))
T_melt = 0.5 * (hu[im, 1] + hu[im + 1, 1])
ifz = int(np.argmin(np.diff(c[:, 2])))
T_freeze = 0.5 * (c[ifz, 1] + c[ifz + 1, 1])

def guides(a):
    a.axvline(T_melt,   ls="--", color="0.4", lw=1)
    a.axvline(T_freeze, ls="--", color="0.4", lw=1)

fig, ax = plt.subplots(1, 2, figsize=(8.4, 3.6))

# (a) volume per atom: heating, seeded cooling, and the unseeded cooling
ax[0].plot(hu[:, 1], hu[:, 2], ".", ms=3, color="firebrick", label="heating")
ax[0].plot(c[:, 1], c[:, 2], ".", ms=3, color="seagreen", label="cooling, seeded")
ax[0].plot(u[:, 1], u[:, 2], ".", ms=2, color="steelblue", label="cooling, no seed")
guides(ax[0])
ax[0].set_xlabel("temperature (K)")
ax[0].set_ylabel(r"volume per atom ($\mathrm{\AA}^3$)")
ax[0].legend(frameon=True, facecolor="white", edgecolor="black", framealpha=1.0,
             fontsize=8)
ylo, yhi = ax[0].get_ylim()
ax[0].text(T_freeze, ylo + 0.04 * (yhi - ylo), f" freeze\n ~{T_freeze:.0f} K",
           color="0.4", fontsize=8)
ax[0].text(T_melt, ylo + 0.04 * (yhi - ylo), f" melt\n ~{T_melt:.0f} K",
           color="0.4", fontsize=8)

# (b) number of fcc-coordinated atoms: the seed (flat line) then regrowth
ax[1].plot(h[:, 1], h[:, 3], ".", ms=3, color="firebrick", label="heating")
ax[1].plot(c[:, 1], c[:, 3], ".", ms=3, color="seagreen", label="cooling, seeded")
guides(ax[1])
ax[1].set_xlabel("temperature (K)")
ax[1].set_ylabel("number of fcc-coordinated atoms")
ax[1].legend(frameon=True, facecolor="white", edgecolor="black", framealpha=1.0,
             fontsize=8)

plt.tight_layout()
plt.savefig("seeded.png", dpi=200)

# report where the melt refreezes on cooling (largest jump in N_fcc as T falls)
cb = c[np.argsort(-c[:, 1])]
i = int(np.argmax(np.diff(cb[:, 3])))
print(f"seeded cooling: fcc count jumps near {0.5*(cb[i,1]+cb[i+1,1]):.0f} K")
print(f"seed stays within {h[:,4].max():.2f} A (avg) of its target over the run")
print("the seed cures the supercooling, but freezes far below Tm (Gibbs-Thomson)")
