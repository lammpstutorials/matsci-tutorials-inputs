#!/usr/bin/env python3
# Python Script
# Licensed under CC BY 4.0
# A Set of Materials-Science Tutorials for the LAMMPS Simulation Package
# Tutorial 8: from the reverse-NEMD output of kappa.lmp, fit the steady-state
# temperature gradient for each cell length and compute the (size-dependent)
# thermal conductivity kappa(L).  Plots the temperature profile and the
# strong finite-size dependence of kappa.
#
# Run kappa.lmp first (e.g. in LAMMPS-GUI), then:  python3 plot.py

import yaml
import numpy as np
import matplotlib.pyplot as plt

CONV = 1602.1766          # eV/(Ang*ps*K) -> W/(m*K)
KAPPA_EXP = 148.0         # experimental Si at 300 K (W/m/K)


def read_profile(fname):
    """average all time blocks of a fix ave/chunk file; return z-fraction, T[K]"""
    blocks, cur = [], None
    for line in open(fname):
        if line.startswith("#"):
            continue
        p = line.split()
        if len(p) == 3:                       # block header: timestep nchunks total
            if cur:
                blocks.append(cur)
            cur = []
        elif len(p) == 4 and cur is not None:
            cur.append((float(p[1]), float(p[3])))
    if cur:
        blocks.append(cur)
    coord = np.array([c for c, _ in blocks[0]])
    T = np.array([[v for _, v in b] for b in blocks]).mean(axis=0)
    return coord, T


def gradient(zfrac, T, L):
    """slope magnitude of the two linear branches away from the swap slabs"""
    zA = zfrac * L
    rise = (zfrac > 0.10) & (zfrac < 0.45)
    fall = (zfrac > 0.60) & (zfrac < 0.95)
    s1 = np.polyfit(zA[rise], T[rise], 1)[0]
    s2 = np.polyfit(zA[fall], T[fall], 1)[0]
    return 0.5 * (abs(s1) + abs(s2))


recs = sorted(yaml.safe_load(open("kappa.yaml")), key=lambda r: r["ncells"])
L_nm, kappa = [], []
for r in recs:
    z, T = read_profile(f"profile_{r['ncells']}.dat")
    dTdz = gradient(z, T, r["length_A"])                          # K/Ang
    flux = r["edelta_eV"] / (2.0 * r["area_A2"] * r["time_ps"])   # eV/Ang^2/ps
    k = flux / dTdz * CONV                                        # W/m/K
    L_nm.append(r["length_A"] / 10.0)
    kappa.append(k)
    print(f"L = {r['length_A']/10:6.1f} nm   dT/dz = {dTdz*10:5.2f} K/nm   "
          f"kappa = {k:6.1f} W/m/K")

L_nm = np.array(L_nm)
kappa = np.array(kappa)
print(f"\nkappa grows by {kappa[-1]/kappa[0]:.1f}x over this length range and is still "
      f"rising:\nthe cells are far shorter than the silicon phonon mean free path "
      f"(~300 nm),\nso transport is quasi-ballistic and kappa is severely "
      f"size-limited.\nexperiment (Si, 300 K): ~{KAPPA_EXP:.0f} W/m/K.")

# --- figure 1: temperature profile for the shortest cell -------------------
r0 = recs[0]
z, T = read_profile(f"profile_{r0['ncells']}.dat")
zA = z * r0["length_A"] / 10.0      # nm
plt.figure(figsize=(5.0, 3.8))
plt.plot(zA, T, "o", color="steelblue", ms=4)
for sel in [(z > 0.10) & (z < 0.45), (z > 0.60) & (z < 0.95)]:
    fit = np.polyfit(zA[sel], T[sel], 1)
    plt.plot(zA[sel], np.polyval(fit, zA[sel]), "-", color="firebrick", lw=2)
plt.xlabel("position along heat-flow direction (nm)")
plt.ylabel("temperature (K)")
plt.title(f"L = {r0['length_A']/10:.1f} nm  (hot slab in the middle)", fontsize=10)
plt.tight_layout()
plt.savefig("kappa-profile.png", dpi=200)

# --- figure 2: finite-size dependence of kappa -----------------------------
plt.figure(figsize=(5.0, 3.8))
plt.plot(L_nm, kappa, "o-", color="steelblue", ms=6)
plt.xlabel("cell length $L$ (nm)")
plt.ylabel(r"apparent thermal conductivity $\kappa$ (W/m/K)")
plt.tight_layout()
plt.savefig("kappa-size.png", dpi=200)
