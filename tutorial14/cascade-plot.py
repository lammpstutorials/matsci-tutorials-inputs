#!/usr/bin/env python3
# Python Script -- Licensed under CC BY 4.0
# A Set of Materials-Science Tutorials for the LAMMPS Simulation Package
# Tutorial 14: plot the number of Frenkel pairs during the cascade.
#
# Reads cascade.dat (written by compute frenkel through fix ave/time in
# cascade.lmp).  Columns are:  step  time(ps)  N_vacancy  N_interstitial.
# Run cascade.lmp, then:  python3 cascade-plot.py

import numpy as np
import matplotlib.pyplot as plt

d = np.loadtxt("cascade.dat")
t = d[:, 1]
nvac = d[:, 2]
nint = d[:, 3]

ipeak = int(np.argmax(nvac))

fig, ax = plt.subplots(figsize=(6.2, 4.2))
ax.plot(t, nvac, "-", color="firebrick", lw=1.6, label="vacancies")
ax.plot(t, nint, "--", color="steelblue", lw=1.4, label="interstitials")

# the cascade fully recombines within ~2.5 ps and the count is then flat to the
# end of the run (15 ps), so zoom in on the active window
tmax = 6.0

# mark the ballistic peak and the surviving plateau
ax.annotate(f"peak {int(nvac[ipeak])}", xy=(t[ipeak], nvac[ipeak]),
            xytext=(t[ipeak] + 0.10 * tmax, nvac[ipeak]), va="center",
            fontsize=9, arrowprops=dict(arrowstyle="->", color="0.4"))
ax.axhline(nvac[-1], color="0.6", lw=0.7, ls=":")
ax.text(tmax, nvac[-1] + 0.03 * nvac[ipeak], f"{int(nvac[-1])} surviving",
        ha="right", va="bottom", fontsize=9)

ax.set_xlim(0, tmax)
ax.set_ylim(0, nvac[ipeak] * 1.12)
ax.set_xlabel("time (ps)")
ax.set_ylabel("number of Frenkel defects")
ax.legend(frameon=True, facecolor="white", edgecolor="black",
          framealpha=1.0, loc="upper right", fontsize=9)
plt.tight_layout()
plt.savefig("cascade.png", dpi=200)
print("wrote cascade.png")
