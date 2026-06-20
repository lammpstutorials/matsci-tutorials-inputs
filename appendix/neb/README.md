# Vacancy migration barrier in fcc aluminum by NEB (multi-partition)

These inputs compute the migration energy barrier for a vacancy hop in fcc
aluminum with the nudged elastic band (NEB) method, using the same Mendelev
Finnis-Sinclair Al potential (`Al_mm.eam.fs`) as Tutorial 13.  NEB is a
**multi-partition** calculation: each image of the band is a separate replica,
so it must be run from the command line with MPI and cannot be run inside
LAMMPS-GUI (which hosts a single LAMMPS instance).

## Files

| file | purpose |
|---|---|
| `in.neb`        | the NEB run (build cell + vacancy, relax, then `neb`) |
| `final.neb`     | final-state file: the migrating atom's target position |
| `Al_mm.eam.fs`  | Mendelev Al EAM potential (also in LAMMPS `potentials/`) |
| `plot.py`       | plot the minimum-energy path from `mep.dat` |
| `mep.dat`       | extracted reaction coordinate vs energy (see below) |
| `neb-barrier.png`, `neb-path.png` | reference figures |

## Running

NEB runs one replica per band image with the `-partition` command-line flag.
For a 12-image band, 12 MPI ranks (one per replica):

    mpirun -np 12 lmp -partition 12x1 -in in.neb

The forward and reverse barriers are printed on the final line of the NEB
output as the `EBF` and `EBR` columns; that line also lists, for each replica,
its reaction coordinate (`RD`) and potential energy (`PE`).  Extract the
minimum-energy path and plot it with:

    tail -1 log.lammps | awk '{n=(NF-9)/2; e0=$11; \
        for(i=0;i<n;i++) printf "%.5f %.5f\n",$(10+2*i),$(11+2*i)-e0}' > mep.dat
    python3 plot.py

## Visualizing the band

`in.neb` also illustrates two multi-partition features (uncomment the block
near the `neb` command; it needs the GRAPHICS package).  `fix graphics/replica`
gathers the migrating atom's position from *every* replica, and the `partition
yes 1` prefix runs a `dump image` on a single partition, so one picture
(`neb-path.png`) shows the whole migration path at once.

## What to expect

The path is symmetric (the two end states are equivalent sites) and peaks at
the saddle.  With the Mendelev potential the migration barrier is
**E_m = 0.63 eV**, close to the experimental value for aluminum (~0.61 eV).
This barrier is the activation energy for vacancy diffusion in the *solid*,
the crystalline counterpart of the liquid-diffusion activation energy measured
in Tutorial 13.
