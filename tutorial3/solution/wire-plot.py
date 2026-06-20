#!/usr/bin/env python3
# Licensed under CC BY 4.0
# A Set of Materials-Science Tutorials for the LAMMPS Simulation Package
# Tutorial 3: compare the tensile stress-strain curves of the bulk crystal
# (deform.lmp, mode=tension), the pure nanowire with the same Zhou
# potential (wire.lmp), and the pure and alloyed nanowires with the
# Apostol-Mishin Al-Cu potential (wire-alloy.lmp).
# Run the inputs first, then:  python3 wire-plot.py

import pandas as pd
import matplotlib.pyplot as plt

CURVES = [("stress-tension.csv", "bulk, Zhou EAM", "firebrick"),
          ("stress-wire.csv", "wire, Zhou EAM", "darkgoldenrod"),
          ("stress-wire-0.csv", "wire, Al-Cu ADP", "steelblue"),
          ("stress-wire-10.csv", "wire +10% Cu, Al-Cu ADP", "seagreen")]

plt.figure(figsize=(5.4, 3.8))
for fname, label, color in CURVES:
    df = pd.read_csv(fname)
    plt.plot(df["strain"], df["sxx_GPa"], "-", lw=1.1, color=color,
             label=label)
    i = df["sxx_GPa"].idxmax()
    print(f"{label:26s} yield {df.sxx_GPa[i]:.2f} GPa at strain {df.strain[i]:.3f}")
plt.axhline(0.0, color="gray", lw=0.6)
plt.xlim(0, 0.7)
plt.xlabel("engineering strain")
plt.ylabel(r"stress $\sigma_{xx}$ (GPa)")
plt.legend(frameon=False, fontsize=8)
plt.tight_layout()
plt.savefig("wire-stress-strain.png", dpi=200)
print("wrote wire-stress-strain.png")
