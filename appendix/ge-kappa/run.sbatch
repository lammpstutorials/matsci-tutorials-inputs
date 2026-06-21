#!/bin/bash
#SBATCH --job-name=ge-kappa
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=64
#SBATCH --time=24:00:00
# Example SLURM script for the germanium RNEMD conductivity scan.
# Adjust the partition, account, module names, and the lmp path for your site.

module load lammps        # or load your own MPI + LAMMPS build (needs ML-SNAP)

# classical Tersoff: cheap, wide length range
srun lmp -in kappa_ge.lmp -v tag tersoff -v pot pot_ge_tersoff.mod \
         -v a0 5.658 -v ncross 6 -v nz 40 80 120 160 200

# machine-learning SNAP: ~100x heavier, shorter range
srun lmp -in kappa_ge.lmp -v tag snap -v pot pot_ge_snap.mod \
         -v a0 5.760 -v ncross 6 -v nz 20 40 60 80

# analyze on a login node afterwards:  python3 plot_ge.py
