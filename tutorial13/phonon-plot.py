#!/usr/bin/env python3
# Python Script -- Licensed under CC BY 4.0
# A Set of Materials-Science Tutorials for the LAMMPS Simulation Package
# Tutorial 13: plot the phonon dispersion of fcc copper.
#
# Reads disp.dat (computed by the phana post-processor from the fix-phonon
# output) and overlays experimental neutron-scattering points (Cu_exp.dat).
# Run phonon.lmp, then phana < dispersion.in, then:  python3 phonon-plot.py

import numpy as np
import matplotlib.pyplot as plt

# computed dispersion: columns qx qy qz qr f1 f2 f3 (THz); qr is the
# cumulative distance along the path, the three f are the phonon branches
d = np.loadtxt("disp.dat")
qr = d[:, 3]
branches = d[:, 4:7]

# experimental points: qr  frequency (THz)
exp = np.loadtxt("Cu_exp.dat")

# high-symmetry points along the Gamma-X-W-X'-Gamma-L path
ticks = [0.0, 0.70711, 1.06066, 1.41421, 2.63896, 3.50498]
labels = [r"$\Gamma$", "X", "W", "X'", r"$\Gamma$", "L"]

fig, ax = plt.subplots(figsize=(6.2, 4.2))
for b in range(3):
    ax.plot(qr, branches[:, b], "-", color="firebrick", lw=1.6,
            label="MD (fix phonon)" if b == 0 else None)
ax.plot(exp[:, 0], exp[:, 1], "o", ms=4, mfc="white", mec="steelblue",
        mew=1.2, label="neutron exp. (298 K)")

for t in ticks:
    ax.axvline(t, color="0.6", lw=0.7)
ax.set_xticks(ticks)
ax.set_xticklabels(labels)
ax.set_xlim(0, ticks[-1])
ax.set_ylim(0, 8.0)
ax.set_ylabel("frequency (THz)")
ax.set_xlabel("wave vector")
ax.legend(frameon=True, facecolor="white", edgecolor="black",
          framealpha=1.0, loc="upper right", fontsize=9)
plt.tight_layout()
plt.savefig("phonon.png", dpi=200)
print("wrote phonon.png")
