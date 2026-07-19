#!/usr/bin/env python3
# Python Script
# Licensed under CC BY 4.0
# A Set of Materials-Science Tutorials for the LAMMPS Simulation Package
# Tutorial 10: overlay the gamma-lines of two fcc metals to show how the
# stacking fault energy depends on the material (and the potential).
#
# Produce the two data sets by running gsfe.lmp twice and renaming the output:
#   lmp -in gsfe.lmp                                  # aluminum (default)
#   mv gsfe.yaml gsfe-al.yaml
#   lmp -in gsfe.lmp -var pot Cu_mishin1.eam.alloy -var elem Cu -var lat 3.61493
#   mv gsfe.yaml gsfe-cu.yaml
# then:  python3 compare.py
# (The shipped solution/ already contains gsfe.yaml for Al and gsfe-cu.yaml.)

import os
import yaml
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

EV_PER_A2_TO_MJM2 = 16021.766

# (file, label, color); fall back to the single-run name for aluminum
al_file = "gsfe-al.yaml" if os.path.exists("gsfe-al.yaml") else "gsfe.yaml"
datasets = [
    (al_file, "Al (Zhou)", "steelblue"),
    ("gsfe-cu.yaml", "Cu (Mishin)", "firebrick"),
]

plt.figure(figsize=(5.2, 3.8))
for fname, label, color in datasets:
    df = pd.DataFrame(yaml.safe_load(open(fname)))
    shift = df["shift"].to_numpy()
    gamma = (df["pe"].to_numpy() - df["pe"][0]) / df["area"].to_numpy() * EV_PER_A2_TO_MJM2
    period = shift[-1]
    # normalize the displacement by the period so the two metals share an x-axis
    x = shift / period
    isf = gamma[(x > 0.2) & (x < 0.55)].min()
    usf = gamma[(x > 0.0) & (x < 0.34)].max()
    plt.plot(x, gamma, "o-", color=color, ms=4,
             label=f"{label}:\nISF {isf:.0f}, USF {usf:.0f} mJ/m$^2$")
    print(f"{label:12s}  period={period:.3f} Ang  ISF={isf:5.1f}  USF={usf:5.1f} mJ/m^2")

plt.axvline(1.0 / 3.0, ls="--", color="gray", lw=1)
plt.text(1.0 / 3.0, plt.ylim()[1] * 0.04, r"  $b_p$", color="gray", fontsize=9)
plt.xlabel(r"displacement along $[11\bar 2]$ (in units of the period $a_0\sqrt{6}/2$)")
plt.ylabel(r"stacking fault energy $\gamma$ (mJ/m$^2$)")
plt.legend(frameon=True, facecolor="white", edgecolor="black", framealpha=1.0,
           fontsize=8.5)
plt.tight_layout()
plt.savefig("gsfe-compare.png", dpi=200)
