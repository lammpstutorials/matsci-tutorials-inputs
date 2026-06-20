#!/usr/bin/env python3
# Python Script
# Licensed under CC BY 4.0
# A Set of Materials-Science Tutorials for the LAMMPS Simulation Package
# Tutorial 6: plot the indentation load-displacement curves written by
# indent.lmp.  Each CSV file loads directly with pandas (no parser needed).
#
# Run indent.lmp once with shape=sphere and once with shape=plane (e.g. in
# LAMMPS-GUI), then:   python3 plot.py

import pandas as pd
import matplotlib.pyplot as plt

# metal-units force (eV/Angstrom) -> nN
EV_PER_ANG_TO_NN = 1.602176634

curves = [("load-sphere.csv", "sphere", "firebrick"),
          ("load-plane.csv", "flat punch", "steelblue")]

plt.figure(figsize=(5.0, 3.8))
for fname, label, color in curves:
    try:
        d = pd.read_csv(fname)
    except FileNotFoundError:
        continue
    plt.plot(d["depth_Ang"], d["load_eV_per_Ang"] * EV_PER_ANG_TO_NN,
             "-", color=color, lw=1.4, label=label)

plt.xlabel(r"indenter depth ($\mathrm{\AA}$)")
plt.ylabel("load (nN)")
plt.legend(frameon=False)
plt.tight_layout()
plt.savefig("load-displacement.png", dpi=200)
print("wrote load-displacement.png")
