# Appendix: pairwise vs. many-body potentials for metals

Inputs that reproduce the three tables of the "Pairwise versus many-body
potentials" appendix: crystal-structure stability, the elastic constants
(Cauchy relation), and the vacancy-formation energy, for aluminum and iron with
five parameterizations — the Lennard-Jones form in two flavors (a *tuned* LJ and
the OpenKIM *universal* LJ), Morse, EAM, and MEAM.

All files are licensed CC BY 4.0 except the potential files, which are
redistributed from the LAMMPS `potentials/` directory under their original
licenses: `Al_zhou.eam.alloy` (Zhou *et al.*), `Fe_mm.eam.fs` (Mendelev *et
al.*), and the Jelinek *et al.* MEAM (`library.meam` + `AlSiMgCuFe.meam`).

All five parameterizations run with the built-in `lj/cut`, `morse`, `eam`, and
`meam` pair styles, so **no optional packages are required**. The pair style for
each case is a one-line include file:

| file | potential |
|---|---|
| `pot_al_lj.mod`, `pot_fe_lj.mod`       | LJ tuned to the lattice constant (`lj/cut`) |
| `pot_al_ulj.mod`, `pot_fe_ulj.mod`     | OpenKIM "universal" LJ (`lj/cut` + `shift yes`) |
| `pot_al_morse.mod`, `pot_fe_morse.mod` | Morse, Girifalco-Weaver (`morse`) |
| `pot_al_eam.mod`, `pot_fe_eam.mod`     | EAM (`eam/alloy`, `eam/fs`) |
| `pot_al_meam.mod`, `pot_fe_meam.mod`   | MEAM, Jelinek (`meam`) |

The *tuned* LJ has sigma/epsilon chosen so the relaxed fcc crystal reproduces the
experimental lattice constant and cohesive energy. The Morse parameters are those
of Girifalco & Weaver, Phys. Rev. **114**, 687 (1959).

## The "universal" LJ

`pot_*_ulj.mod` reproduce the Elliott-Akerson "universal" shifted Lennard-Jones
model distributed through OpenKIM (`LennardJones612_UniversalShifted__MO_959249795837_003`),
whose parameters are fixed by the isolated dimer (minimum at the sum of covalent
radii, well depth at the bond dissociation energy). Because that model is just a
shifted 12-6 LJ, it is reproduced exactly by `lj/cut` + `pair_modify shift yes`,
so the shipped inputs need no KIM package. The effective sigma/epsilon/cutoff were
read off a dimer energy scan (`dimer.lmp`, which *does* require the KIM package);
the resulting `lj/cut` rows match `pair_style kim` to better than 0.03 % on the
fcc cohesive energy.

## Crystal-structure stability

`structure.lmp` relaxes one crystal structure (`fcc`, `bcc`, or `hcp`) to its
own equilibrium and prints the cohesive energy per atom. The lowest energy over
the three structures is the predicted ground state, e.g. aluminum in bcc:

    lmp -in structure.lmp -var lat bcc -var a 3.20 -var mass 26.9815 \
        -var n 8 -var pot pot_al_eam.mod

## Elastic constants (Cauchy relation)

`in.elastic` is the standard LAMMPS `examples/ELASTIC` driver (it `include`s
`init.mod`, `potential.mod`, and `displace.mod`). For fcc aluminum:

    lmp -in in.elastic -var lat fcc -var a 4.05 -var mass 26.9815 \
        -var nx 5 -var ny 5 -var nz 5 -var pot pot_al_morse.mod

Both LJ parameterizations and Morse give the equal pair `C12all = C44all` (the
Cauchy relation); EAM and MEAM give C12 ~ 2*C44.

## Vacancy formation

`vac.lmp` relaxes a perfect cell and a one-vacancy cell at fixed box and prints
the cohesive energy per atom and the vacancy-formation energy. Pass each
potential's own equilibrium lattice constant (from `structure.lmp`):

    lmp -in vac.lmp -var lat fcc -var a 4.0817 -var mass 26.9815 \
        -var n 5 -var pot pot_al_eam.mod

## Reproduce everything

`run.sh` runs all cases and prints the numbers that go into the three tables.
Every calculation is a static minimization that finishes in seconds:

    sh run.sh /path/to/lmp
