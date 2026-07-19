#!/usr/bin/env python3
# Python Script
# Licensed under CC BY 4.0
# A Set of Materials-Science Tutorials for the LAMMPS Simulation Package
# Tutorial 5: plot the tension and compression stress-strain curves written
# by deform.lmp. Each CSV file loads directly with pandas (no parser needed).
#
# Run deform.lmp once with mode=tension and once with mode=compress
# (e.g. in LAMMPS-GUI), then:   python3 plot.py

import pandas as pd
import matplotlib.pyplot as plt

tension = pd.read_csv("stress-tension.csv")
compression = pd.read_csv("stress-compress.csv")

plt.figure(figsize=(5.0, 3.8))
plt.plot(tension["strain"], tension["sxx_GPa"], "-", color="firebrick", label="tension")
plt.plot(compression["strain"], compression["sxx_GPa"], "-", color="steelblue",
         label="compression")
plt.axhline(0.0, color="gray", lw=0.6)
plt.axvline(0.0, color="gray", lw=0.6)
plt.xlabel("engineering strain")
plt.ylabel(r"stress $\sigma_{xx}$ (GPa)")
plt.legend(frameon=False)
plt.tight_layout()
plt.savefig("stress-strain.png", dpi=200)

for name, df in [("tension", tension), ("compression", compression)]:
    idx = df["sxx_GPa"].abs().idxmax()
    print(f"{name}: peak |sigma_xx| = {df['sxx_GPa'][idx]:.2f} GPa "
          f"at strain = {df['strain'][idx]:.3f}")
