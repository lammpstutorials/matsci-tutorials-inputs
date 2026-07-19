#!/usr/bin/env python3
# Licensed under CC BY 4.0
# A Set of Materials-Science Tutorials for the LAMMPS Simulation Package
# Tutorial 17: chemical short-range order in a high-entropy alloy
#
# Reads the trajectory anneal.dump written by sro.lmp and computes the
# first-shell Warren-Cowley short-range-order parameters
#     alpha_AB = 1 - P_AB / c_B ,
# where P_AB is the probability that a first-neighbor of an A atom is of type
# B and c_B is the concentration of B.  alpha = 0 for a random alloy,
# alpha > 0 means A avoids B in its first shell, alpha < 0 means A attracts B.
# Plots the evolution of selected pairs during the hybrid MD/MC run and prints
# the full alpha matrix averaged over the converged tail of the trajectory.

import numpy as np
import matplotlib.pyplot as plt

DUMPFILE = "anneal.dump"
EFILE = "energy.dat"
RCUT = 3.0          # A; first-neighbor shell (fcc a0=3.56 -> d_NN = 2.52 A)
NTYPES = 4
ELEMENTS = ["Co", "Cr", "Fe", "Ni"]
TAIL = 20           # number of final frames averaged for the alpha matrix


def read_dump(fname):
    """Return a list of (step, types, xyz, box_lengths) tuples."""
    frames = []
    with open(fname) as f:
        lines = f.read().splitlines()
    i = 0
    while i < len(lines):
        if not lines[i].startswith("ITEM: TIMESTEP"):
            i += 1
            continue
        step = int(lines[i + 1])
        natoms = int(lines[i + 3])
        box = np.array([[float(v) for v in lines[i + 5 + k].split()[:2]]
                        for k in range(3)])
        cols = lines[i + 8].split()[2:]
        data = np.array([[float(v) for v in lines[i + 9 + k].split()]
                         for k in range(natoms)])
        col = {c: data[:, j] for j, c in enumerate(cols)}
        xyz = np.stack([col["x"], col["y"], col["z"]], axis=1)
        frames.append((step, col["type"].astype(int), xyz, box[:, 1] - box[:, 0]))
        i += 9 + natoms
    return frames


def warren_cowley(types, xyz, blen):
    """First-shell Warren-Cowley matrix, symmetrized."""
    conc = np.array([np.mean(types == t) for t in range(1, NTYPES + 1)])
    pairs = np.zeros((NTYPES, NTYPES))
    total = np.zeros(NTYPES)
    for i in range(len(types)):
        d = xyz - xyz[i]
        d -= blen * np.round(d / blen)        # periodic boundaries
        r = np.sqrt((d * d).sum(axis=1))
        r[i] = 1e9
        for j in np.where(r < RCUT)[0]:
            pairs[types[i] - 1, types[j] - 1] += 1
            total[types[i] - 1] += 1
    alpha = 1.0 - (pairs / total[:, None]) / conc[None, :]
    return 0.5 * (alpha + alpha.T)


frames = read_dump(DUMPFILE)
steps = np.array([f[0] for f in frames])
alphas = np.array([warren_cowley(t, x, b) for _, t, x, b in frames])

# --- full alpha matrix, averaged over the converged tail ---
amean = alphas[-TAIL:].mean(axis=0)
print(f"Warren-Cowley alpha (first shell, average of last {TAIL} frames):")
print("      " + "".join(f"{e:>8}" for e in ELEMENTS))
for i, e in enumerate(ELEMENTS):
    print(f"  {e:>4}" + "".join(f"{amean[i, j]:+8.3f}" for j in range(NTYPES)))

# --- Figure 1: evolution of selected pairs during the MD/MC run ---
# faint markers: per-snapshot values; solid lines: running average
fig, ax = plt.subplots(figsize=(5.0, 3.8))
sel = [((1, 1), "Cr-Cr", "firebrick"), ((1, 2), "Cr-Fe", "steelblue"),
       ((0, 3), "Co-Ni", "darkgoldenrod")]
NRUN = 7                                  # running-average window (frames)
ker = np.ones(NRUN) / NRUN
for (i, j), lbl, c in sel:
    y = alphas[:, i, j]
    ax.plot(steps, y, "o", ms=3.5, color=c, alpha=0.35)
    avg = np.convolve(y, ker, mode="valid")
    ax.plot(steps[NRUN - 1:], avg, "-", lw=1.8, color=c, label=lbl)
ax.axhline(0.0, color="gray", lw=0.8, ls=":")
ax.annotate("random", (steps[2], 0.004), color="gray", fontsize=8)
ax.set_xlabel("time step")
ax.set_ylabel(r"Warren-Cowley $\alpha$ (first shell)")
ax.set_title("Development of chemical short-range order")
ax.legend(frameon=False, fontsize=9)
fig.tight_layout()
fig.savefig("sro-alpha.png", dpi=200)

# --- Figure 2: potential energy during the run ---
step_e, pe = np.loadtxt(EFILE, comments="#", unpack=True)
fig, ax = plt.subplots(figsize=(5.0, 3.8))
ax.plot(step_e, 1000.0 * (pe - pe[0]), color="seagreen", lw=1.0)
ax.set_xlabel("time step")
ax.set_ylabel(r"$\Delta E$ (meV/atom)")
ax.set_title("Potential energy during the MD/MC run")
fig.tight_layout()
fig.savefig("sro-energy.png", dpi=200)
print("wrote sro-alpha.png and sro-energy.png")
