#!/usr/bin/env python3
# Licensed under CC BY 4.0
# A Set of Materials-Science Tutorials for the LAMMPS Simulation Package
# Tutorial 17: chemical short-range order in a high-entropy alloy
#
# Side-by-side comparison of the first-shell Warren-Cowley matrices produced
# by the three potentials at 1300 K (sro-2nnmeam.lmp, sro-sharifi.lmp,
# sro-hot.lmp).  Each panel shows the full 4x4 alpha matrix averaged over the
# tail of the corresponding trajectory: red = attraction (alpha < 0, more
# pairs than random), blue = avoidance (alpha > 0), white = ideal solution.

import numpy as np
import matplotlib.pyplot as plt

RCUT = 3.0
NTYPES = 4
ELEMENTS = ["Co", "Cr", "Fe", "Ni"]
TAIL = 10
RUNS = [("anneal-2nnmeam.dump", "2NN-MEAM (Wang 2023)"),
        ("anneal-sharifi.dump", "MEAM (Sharifi-Wick 2025)"),
        ("anneal-hot.dump", "EAM (Farkas-Caro 2018)")]


def read_dump(fname):
    frames = []
    with open(fname) as f:
        lines = f.read().splitlines()
    i = 0
    while i < len(lines):
        if not lines[i].startswith("ITEM: TIMESTEP"):
            i += 1
            continue
        natoms = int(lines[i + 3])
        box = np.array([[float(v) for v in lines[i + 5 + k].split()[:2]]
                        for k in range(3)])
        cols = lines[i + 8].split()[2:]
        data = np.array([[float(v) for v in lines[i + 9 + k].split()]
                         for k in range(natoms)])
        col = {c: data[:, j] for j, c in enumerate(cols)}
        xyz = np.stack([col["x"], col["y"], col["z"]], axis=1)
        frames.append((col["type"].astype(int), xyz, box[:, 1] - box[:, 0]))
        i += 9 + natoms
    return frames


def warren_cowley(types, xyz, blen):
    conc = np.array([np.mean(types == t) for t in range(1, NTYPES + 1)])
    pairs = np.zeros((NTYPES, NTYPES))
    total = np.zeros(NTYPES)
    for i in range(len(types)):
        d = xyz - xyz[i]
        d -= blen * np.round(d / blen)
        r = np.sqrt((d * d).sum(axis=1))
        r[i] = 1e9
        for j in np.where(r < RCUT)[0]:
            pairs[types[i] - 1, types[j] - 1] += 1
            total[types[i] - 1] += 1
    alpha = 1.0 - (pairs / total[:, None]) / conc[None, :]
    return 0.5 * (alpha + alpha.T)


fig, axes = plt.subplots(1, 3, figsize=(10.0, 3.6))
for ax, (fname, title) in zip(axes, RUNS):
    frames = read_dump(fname)
    A = np.array([warren_cowley(t, x, b) for t, x, b in frames[-TAIL:]]).mean(axis=0)
    im = ax.imshow(A, cmap="RdBu", vmin=-0.42, vmax=0.42)
    ax.set_xticks(range(NTYPES), ELEMENTS)
    ax.set_yticks(range(NTYPES), ELEMENTS)
    ax.set_title(title, fontsize=10)
    for i in range(NTYPES):
        for j in range(NTYPES):
            ax.text(j, i, f"{A[i, j]:+.2f}", ha="center", va="center",
                    fontsize=8.5,
                    color="white" if abs(A[i, j]) > 0.25 else "black")
cb = fig.colorbar(im, ax=axes, shrink=0.85, pad=0.02)
cb.set_label(r"Warren-Cowley $\alpha$")
fig.savefig("compare-alpha.png", dpi=200, bbox_inches="tight")
print("wrote compare-alpha.png")
