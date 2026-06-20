#!/usr/bin/env python3
# Python Script
# Licensed under CC BY 4.0
# A Set of Materials-Science Tutorials for the LAMMPS Simulation Package
# Tutorial 10: diagonalize the Gamma-point dynamical matrix written by
# dynmat.lmp and report the phonon frequencies.
#
# dynmat.dat holds the mass-weighted dynamical matrix D (eskm units, 1/fs^2),
# printed three numbers (the three beta components) per line, in the order
# j fastest, then alpha, then i.  That is exactly the row-major layout of a
# (N,3,N,3) array, so it reshapes straight into the 3N x 3N matrix.
#
# Run dynmat.lmp first, then:  python3 phonon-gamma.py

import numpy as np
import matplotlib.pyplot as plt

d = np.loadtxt("dynmat.dat")
N = int(round((d.shape[0] / 3) ** 0.5))      # 3*N^2 lines  ->  N atoms
D = d.reshape(N, 3, N, 3).reshape(3 * N, 3 * N)
D = 0.5 * (D + D.T)                          # symmetrize (finite-diff noise)

eigval = np.linalg.eigvalsh(D)
# frequency in THz: nu = sqrt(eigenvalue) / (2 pi); keep the sign of tiny
# negative acoustic eigenvalues so they show up as ~0 rather than NaN
nu = np.sign(eigval) * np.sqrt(np.abs(eigval)) / (2.0 * np.pi)
nu.sort()

exp_optical = 15.53          # experimental Si Gamma optical mode (THz)

print("%d atoms in the cell -> %d phonon modes at Gamma" % (N, 3 * N))
print("frequencies (THz):")
for i, f in enumerate(nu):
    kind = "acoustic" if abs(f) < 1.0 else "optical"
    print("  mode %d:  %7.3f  (%s)" % (i + 1, f, kind))
print("optical mode: %.2f THz   (experiment: %.2f THz)" % (nu.max(), exp_optical))

# level diagram: every mode as a horizontal line, acoustic near 0, optical high
fig, ax = plt.subplots(figsize=(3.6, 4.2))
for f in nu:
    ax.hlines(f, 0.15, 0.85, color="goldenrod", lw=2.5)
ax.axhline(exp_optical, ls="--", color="0.4", lw=1.2)
ax.text(0.5, exp_optical - 0.9, "experiment", ha="center", color="0.4", fontsize=9)
ax.text(0.5, nu.max() + 0.4, "3 optical", ha="center", fontsize=9)
ax.text(0.5, 0.6, "3 acoustic", ha="center", fontsize=9)
ax.set_xticks([])
ax.set_xlim(0, 1)
ax.set_ylim(-1.0, 20.0)
ax.set_ylabel("frequency (THz)")
ax.set_title("Si phonons at $\\Gamma$ (SW)")
plt.tight_layout()
plt.savefig("gamma-levels.png", dpi=200)
print("wrote gamma-levels.png")
