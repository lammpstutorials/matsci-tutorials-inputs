# Appendix: machine-learned vs. classical potentials

Inputs that reproduce Table 2 of the article (lattice constant, cohesive and
vacancy-formation energy, elastic constants, and cost for silicon and tantalum
with several interatomic potentials).  All files are licensed CC BY 4.0 except
the potential files (`*.sw`, `*.tersoff`, `*.meam.spline`, `*.snap*`), which are
redistributed from the LAMMPS `potentials/` directory under their original
licenses.  The SNAP models are from Zuo *et al.*, J. Phys. Chem. A **124**, 731
(2020).

The pair style for each potential is a one-line include file:

| file | potential |
|---|---|
| `pot_si_sw.mod`      | silicon, Stillinger–Weber (`pair_style sw`) |
| `pot_si_tersoff.mod` | silicon, Tersoff 1988 (`pair_style tersoff`) |
| `pot_si_meam.mod`    | silicon, MEAM spline (`pair_style meam/spline`) |
| `pot_si_snap.mod`    | silicon, SNAP (`pair_style snap`) |
| `pot_ta_snap.mod`    | tantalum, SNAP + ZBL (`pair_style hybrid/overlay`) |

## Elastic constants, lattice constant, moduli

`in.elastic` is the standard LAMMPS `examples/ELASTIC` driver (it `include`s
`init.mod`, `potential.mod`, and `displace.mod`), extended to print the relaxed
lattice constant `RESULT_A0` and energy.  It is parameterized by command-line
variables.  For each silicon potential (diamond, 3x3x3 cell):

    lmp -in in.elastic -var lat diamond -var a 5.43 -var mass 28.0855 \
        -var nx 3 -var ny 3 -var nz 3 -var pot pot_si_sw.mod

and for tantalum (bcc, 4x4x4 cell):

    lmp -in in.elastic -var lat bcc -var a 3.31 -var mass 180.95 \
        -var nx 4 -var ny 4 -var nz 4 -var pot pot_ta_snap.mod

The Voigt–Reuss–Hill bulk modulus *B*, shear modulus *G*, Young's modulus, and
Poisson ratio *nu* in the table are computed from C11, C12, C44 with the usual
cubic averaging formulas.

## Cohesive and vacancy-formation energy

`vac.lmp` relaxes a perfect cell and a one-vacancy cell and prints both the
cohesive energy per atom and the vacancy-formation energy, e.g.

    lmp -in vac.lmp -var lat diamond -var a 5.431 -var mass 28.0855 \
        -var n 4 -var pot pot_si_sw.mod

The cohesive energy is meaningful only for the classical potentials; an ML
model's energy zero is its DFT reference, not the isolated atom.

## Cost

`timing.lmp` runs a short serial molecular-dynamics trajectory; divide the
reported loop time by (atoms x steps) to get the per-atom-step cost:

    OMP_NUM_THREADS=1 lmp -in timing.lmp -var lat diamond -var a 5.43 \
        -var mass 28.0855 -var n 4 -var pot pot_si_snap.mod

## Note on the potential files

The pair-style include files reference the potential files by name, so run from
this directory (or set the `LAMMPS_POTENTIALS` environment variable to the LAMMPS
`potentials/` directory, where the same files are installed).
