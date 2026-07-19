#!/usr/bin/env python3
# Python Script
# Licensed under CC BY 4.0
# A Set of Materials-Science Tutorials for the LAMMPS Simulation Package
# Tutorial 7: plot the grain-boundary fracture stress-strain curve written by
# fracture.lmp.  The CSV file loads directly with pandas (no parser needed).
#
# Run fracture.lmp first (e.g. in LAMMPS-GUI), then:   python3 plot.py

import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_csv("stress-fracture.csv")
peak = data.loc[data["stress_GPa"].idxmax()]

plt.figure(figsize=(5.0, 3.8))
plt.plot(data["strain"], data["stress_GPa"], "-o", color="firebrick",
         markersize=3, lw=1.2)
plt.plot(peak["strain"], peak["stress_GPa"], "o", color="black", markersize=6,
         label=f"peak {peak['stress_GPa']:.1f} GPa at "
               f"strain {peak['strain']:.3f}")
plt.axhline(0.0, color="gray", lw=0.6)
plt.xlabel("engineering strain")
plt.ylabel(r"stress $\sigma_{yy}$ (GPa)")
plt.legend(frameon=False)
plt.tight_layout()
plt.savefig("fracture-stress-strain.png", dpi=200)

print(f"peak stress {peak['stress_GPa']:.2f} GPa at strain {peak['strain']:.4f}")
print("the stress drops abruptly to zero when the grain boundary cleaves")
