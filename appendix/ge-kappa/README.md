# Germanium thermal conductivity: classical vs. machine-learned (HPC offload)

These inputs compute the lattice thermal conductivity of crystalline germanium
by reverse non-equilibrium molecular dynamics (the same Mueller-Plathe method as
Tutorial 8), for a **classical** Tersoff potential and a **machine-learned**
SNAP potential, and compare their finite-size behaviour.  They are deliberately
**not** sized for an interactive run: the cells are large, the runs are long, and
the SNAP evaluations are ~100x more expensive than Tersoff, so this is meant to
be submitted to an HPC cluster with MPI.

## Files

| file | purpose |
|---|---|
| `kappa_ge.lmp`        | reverse-NEMD conductivity scan (parameterized) |
| `pot_ge_tersoff.mod`  | `pair_style tersoff` + `Ge.tersoff` |
| `pot_ge_snap.mod`     | `pair_style snap` + `Ge_Zuo_JPCA2020` |
| `plot_ge.py`          | fit gradients, kappa(L), and 1/kappa-vs-1/L extrapolation |
| `run.sbatch`          | example SLURM batch script |
| `Ge.tersoff`, `Ge_Zuo_JPCA2020.*` | potential files (also in LAMMPS `potentials/`) |

## Running

The script loops over a list of cell lengths (`nz`, in conventional cells along
the heat-flow z axis) for one potential at a time.  Run the cheap classical
scan over a wide length range, and the expensive SNAP scan over a shorter range:

    mpirun -np 256 lmp -in kappa_ge.lmp -v tag tersoff -v pot pot_ge_tersoff.mod \
                       -v a0 5.658 -v ncross 6 -v nz 20 40 80 120 160 200

    mpirun -np 512 lmp -in kappa_ge.lmp -v tag snap    -v pot pot_ge_snap.mod \
                       -v a0 5.760 -v ncross 6 -v nz 20 40 60 80

`a0` is only a starting guess (an NPT step relaxes the box to each potential's
own 300 K lattice constant).  `ncross` is the cross-section in unit cells
(6 => a 6x6 bar).  Decompose the domain along the long axis for good scaling,
e.g. `-pk omp 1` or a `processors 1 1 P` grid.

Each potential writes `kappa_<tag>.yaml` and `profile_<tag>_<nz>.dat`.  When both
have finished, analyze them together:

    python3 plot_ge.py

This prints the apparent kappa(L) and the extrapolated bulk values for each
potential and writes `ge-kappa.png` (1/kappa vs 1/L).  Each potential is fitted
twice — all cells with equal weights (solid line, open square intercept) and
with the smallest cell excluded (dashed line, open diamond) — because the
extrapolation, not the simulations, dominates the uncertainty (see below); a
reliability-weighted fit is also printed.

## Results (shipped with this folder)

The output of the runs reported in the appendix is included, so the figure can
be reproduced with `plot_ge.py` alone: `kappa_{tersoff,snap}.yaml`,
`profile_<tag>_<nz>.dat` (steady-state temperature profiles), and the LAMMPS
logs `kappa_ge.log` (Tersoff, HPC), `kappa_t20.log` (the smallest Tersoff
cell, run on a workstation), and `kappa_ml.log` (SNAP).

Germanium's phonon mean free path is much shorter than silicon's, but still
longer than affordable cells, so kappa(L) rises with L for both potentials and
must be extrapolated — and the extrapolation, not the simulation statistics,
dominates the uncertainty.  Tersoff: kappa(L) = 4.5-21 W/m/K over L = 11-111 nm;
the bulk value is 35 W/m/K when all cells are fitted at equal weight, 58 W/m/K
with the smallest cell excluded, and in between for reliability-weighted fits —
all with R^2 >= 0.97.  That is the magnitude of the experimental ~60 W/m/K, but
no precision determination.  SNAP: kappa(L) = 1.6-3.8 W/m/K over L = 12-46 nm,
extrapolating to 6-13 W/m/K under the same fitting choices — a factor of 5-10
*below* experiment in every variant: nothing in the SNAP training data
constrains anharmonic phonon-phonon scattering, so the model scatters phonons
far too strongly, even though its static properties (elastic constants etc.)
are good.  See the appendix for the discussion.

## Resource notes

Measured on 64 MPI ranks (one node), except where noted:

- Tersoff: the five-length HPC scan (up to 57600 atoms) took ~35 core-hours
  (about half an hour wall time); the smallest cell (`nz 20`, 5760 atoms)
  adds ~1 core-hour and also runs on a desktop machine.
- SNAP: the four-length scan (up to 23040 atoms) took ~1800 core-hours
  (about 28 hours wall time).  If memory or time is tight, reduce `ncross`
  to 5 or the production `run` length in `kappa_ge.lmp` (currently 400000
  steps).
