#!/usr/bin/env python3
# Python Script -- Licensed under CC BY 4.0
# A Set of Materials-Science Tutorials for the LAMMPS Simulation Package
# Appendix (HPC offload): thermal conductivity of germanium with a classical
# (Tersoff) and a machine-learned (SNAP) potential, from the reverse-NEMD output
# of kappa_ge.lmp.  For each potential it fits the steady-state temperature
# gradient at every cell length, forms the apparent kappa(L), and extrapolates
# 1/kappa vs 1/L to the bulk (infinite-length) value.
#
# The extrapolation is done twice per potential: once with all cells and equal
# weights (solid line) and once with the smallest cell excluded (dashed line).
# The spread between the two bulk values is the systematic uncertainty of the
# extrapolation itself: the smallest cell has the largest leverage on the fit,
# the largest finite-size deviation from the asymptotic 1/kappa-1/L law, and
# (since errors in kappa propagate as d(1/kappa) = dk/kappa^2) the least
# reliable 1/kappa.  A reliability-weighted fit is also reported on stdout.
#
# Run kappa_ge.lmp for both potentials first, then:  python3 plot_ge.py

import yaml
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

CONV = 1602.1766          # eV/(Ang*ps*K) -> W/(m*K)
KAPPA_EXP = 60.0          # experimental Ge at 300 K (W/m/K)
TAGS = {"tersoff": ("Tersoff (classical)", "steelblue"),
        "snap":    ("SNAP (machine-learned)", "firebrick")}


def read_profile(fname):
    blocks, cur = [], None
    for line in open(fname):
        if line.startswith("#"):
            continue
        p = line.split()
        if len(p) == 3:
            if cur:
                blocks.append(cur)
            cur = []
        elif len(p) == 4 and cur is not None:
            cur.append((float(p[1]), float(p[3])))
    if cur:
        blocks.append(cur)
    coord = np.array([c for c, _ in blocks[0]])
    # discard the first few blocks as transient, average the rest
    keep = blocks[len(blocks) // 3:]
    T = np.array([[v for _, v in b] for b in keep]).mean(axis=0)
    return coord, T


def gradient(zfrac, T, L):
    zA = zfrac * L
    rise = (zfrac > 0.10) & (zfrac < 0.45)
    fall = (zfrac > 0.60) & (zfrac < 0.95)
    s1 = np.polyfit(zA[rise], T[rise], 1)[0]
    s2 = np.polyfit(zA[fall], T[fall], 1)[0]
    return 0.5 * (abs(s1) + abs(s2))


def analyze(tag):
    recs = sorted(yaml.safe_load(open(f"kappa_{tag}.yaml")),
                  key=lambda r: r["ncells"])
    L_nm, kappa = [], []
    for r in recs:
        z, T = read_profile(f"profile_{tag}_{r['ncells']}.dat")
        dTdz = gradient(z, T, r["length_A"])
        flux = r["edelta_eV"] / (2.0 * r["area_A2"] * r["time_ps"])
        kappa.append(flux / dTdz * CONV)
        L_nm.append(r["length_A"] / 10.0)
    return np.array(L_nm), np.array(kappa)


fig, ax = plt.subplots(figsize=(6.0, 4.2))
handles = []
for tag, (label, color) in TAGS.items():
    try:
        L, k = analyze(tag)
    except FileNotFoundError:
        print(f"(no data for {tag}; run kappa_ge.lmp -v tag {tag} first)")
        continue
    # 1/kappa = 1/kappa_bulk * (1 + l_mfp/L)  ->  linear in 1/L
    # (L is sorted ascending, so index 0 is the smallest cell)
    invL, invk = 1.0 / L, 1.0 / k
    s_all, i_all = np.polyfit(invL, invk, 1)
    s_big, i_big = np.polyfit(invL[1:], invk[1:], 1)
    # weights 1/sigma(1/kappa) ~ kappa^2 for a size-independent error in kappa
    s_wgt, i_wgt = np.polyfit(invL, invk, 1, w=k**2)
    print(f"{label}: kappa(L) = {k.min():.0f}-{k.max():.0f} W/m/K; bulk = "
          f"{1.0/i_all:.0f} (all cells) / {1.0/i_big:.0f} (smallest excluded) / "
          f"{1.0/i_wgt:.0f} (weighted) W/m/K")
    ax.plot(invL, invk, "o", color=color, ms=6)
    handles.append(Line2D([], [], marker="o", ls="", color=color, ms=6,
                          label=label))
    xx = np.linspace(0, invL.max() * 1.05, 50)
    ax.plot(xx, s_all * xx + i_all, "-", color=color, lw=1.5)
    ax.plot(xx, s_big * xx + i_big, "--", color=color, lw=1.2)
    ax.plot(0, i_all, "s", color=color, ms=8, mfc="white", mew=1.6,
            clip_on=False, zorder=5)
    ax.plot(0, i_big, "D", color=color, ms=7, mfc="white", mew=1.6,
            clip_on=False, zorder=5)

handles += [Line2D([], [], ls="-", color="black", lw=1.5,
                   label="fit to all cells"),
            Line2D([], [], ls="--", color="black", lw=1.2,
                   label="smallest cell excluded")]
ax.axhline(1.0 / KAPPA_EXP, color="0.5", lw=0.8, ls=":")
ax.text(ax.get_xlim()[1], 1.0 / KAPPA_EXP, f" exp. {KAPPA_EXP:.0f} W/m/K",
        va="bottom", ha="right", fontsize=9, color="0.4")
ax.set_xlabel(r"inverse cell length $1/L$ (nm$^{-1}$)")
ax.set_ylabel(r"inverse conductivity $1/\kappa$ (m$\cdot$K/W)")
ax.set_xlim(left=0)
ax.legend(handles=handles, frameon=True, edgecolor="black", framealpha=1.0,
          fontsize=9)
plt.tight_layout()
plt.savefig("ge-kappa.png", dpi=200)
print("wrote ge-kappa.png  (open squares/diamonds = extrapolated bulk values)")
