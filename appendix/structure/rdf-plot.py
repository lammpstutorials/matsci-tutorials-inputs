#!/usr/bin/env python3
# Licensed under CC BY 4.0
# A Set of Materials-Science Tutorials for the LAMMPS Simulation Package
# Appendix: structure-analysis toolkit -- radial distribution functions
#
# Reads the time-averaged g(r) files rdf_<material>_<T>.dat written by
# rdf-al.lmp, rdf-fe.lmp, rdf-mg.lmp, and rdf-si.lmp and draws one panel
# per lattice with the three temperatures overlaid.  Each crystal has its
# own fingerprint of coordination shells; heating broadens the shells, and
# in the liquid only short-range order survives: a broad first peak, a
# weak second, and g -> 1 beyond.

import glob
import re
import numpy as np
import matplotlib.pyplot as plt

PANELS = [("al", "fcc aluminum"), ("fe", "bcc iron"),
          ("mg", "hcp magnesium"), ("si", "diamond silicon")]
COLORS = ["steelblue", "darkgoldenrod", "firebrick"]


def read_rdf(fname):
    """Return (r, g) from a fix ave/time mode-vector file (last block)."""
    rows = []
    with open(fname) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) == 2:        # block header: timestep, nrows
                rows = []
            else:
                rows.append((float(parts[1]), float(parts[2])))
    return np.array(rows).T


fig, axes = plt.subplots(2, 2, figsize=(8.0, 6.0), sharex=True)
for ax, (tag, title) in zip(axes.flat, PANELS):
    files = sorted(glob.glob(f"rdf_{tag}_*.dat"),
                   key=lambda f: float(re.search(r"_(\d+)\.dat", f).group(1)))
    for fname, c in zip(files, COLORS):
        T = re.search(r"_(\d+)\.dat", fname).group(1)
        r, g = read_rdf(fname)
        ax.plot(r, g, color=c, lw=1.2, label=f"{T} K")
    ax.set_title(title, fontsize=10)
    ax.axhline(1.0, color="gray", lw=0.5, ls=":")
    ax.legend(frameon=False, fontsize=8)
    ax.set_xlim(1.5, 8.5)
for ax in axes[1]:
    ax.set_xlabel(r"$r$ ($\mathrm{\AA}$)")
for ax in axes[:, 0]:
    ax.set_ylabel(r"$g(r)$")
fig.tight_layout()
fig.savefig("rdf-lattices.png", dpi=200)
print("wrote rdf-lattices.png")
