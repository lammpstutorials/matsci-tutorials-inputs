#!/usr/bin/env python3
# Licensed under CC BY 4.0
# A Set of Materials-Science Tutorials for the LAMMPS Simulation Package
# Appendix (multi-partition): plot the NEB minimum-energy path for vacancy
# migration in fcc aluminum.
#
# mep.dat holds two columns (normalized reaction coordinate, energy relative to
# the initial state in eV), extracted from the final line of the NEB output:
#   tail -1 neb.out | awk '{n=(NF-9)/2; e0=$11; for(i=0;i<n;i++) \
#       printf "%.5f %.5f\n",$(10+2*i),$(11+2*i)-e0}'  > mep.dat

import numpy as np
import matplotlib.pyplot as plt

rc, de = np.loadtxt("mep.dat", comments="#", unpack=True)
ibar = np.argmax(de)
barrier = de[ibar]
print(f"migration barrier E_m = {barrier:.3f} eV at reaction coordinate {rc[ibar]:.2f}")

fig, ax = plt.subplots(figsize=(5.0, 3.6))
ax.plot(rc, de, "o-", color="firebrick", ms=6, lw=1.3)
ax.axhline(barrier, ls=":", color="gray", lw=0.8)
ax.annotate(f"$E_m$ = {barrier:.2f} eV",
            xy=(rc[ibar], barrier), xytext=(rc[ibar], barrier * 0.62),
            ha="center", fontsize=10,
            arrowprops=dict(arrowstyle="->", color="gray"))
ax.set_xlabel("reaction coordinate")
ax.set_ylabel(r"energy relative to initial (eV)")
ax.set_title("NEB minimum-energy path: vacancy hop in fcc Al")
ax.set_xlim(0, 1)
ax.set_ylim(bottom=-0.03)
fig.tight_layout()
fig.savefig("neb-barrier.png", dpi=200)
print("wrote neb-barrier.png")
