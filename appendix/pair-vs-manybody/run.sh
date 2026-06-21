#!/bin/sh
# Reproduce the three tables of the pairwise-vs-many-body appendix.
# Usage: sh run.sh [lmp-executable]   (defaults to "lmp" in PATH)
# Every run is a static minimization and finishes in seconds.
LMP="${1:-lmp}"

echo "===== Table 1: crystal-structure stability (Ecoh per atom) ====="
for el in "al 26.9815 fcc:4.05 bcc:3.20 hcp:2.86" "fe 55.847 fcc:3.60 bcc:2.86 hcp:2.55"; do
  set -- $el; e=$1; m=$2; shift 2
  for pot in lj ulj morse eam meam; do
    for ls in "$@"; do
      lat=${ls%:*}; a=${ls#*:}
      "$LMP" -in structure.lmp -var lat $lat -var a $a -var mass $m \
             -var n 8 -var pot pot_${e}_${pot}.mod -log none 2>/dev/null | grep '^RESULT'
    done
  done
done

echo "===== Table 2: elastic constants of fcc aluminum (Cauchy: C12 = C44) ====="
# seed each relaxation near that potential's own a0 (the universal LJ sits at 3.33)
for pa in "lj:4.05" "ulj:3.33" "morse:4.05" "eam:4.05" "meam:4.05"; do
  pot=${pa%:*}; a=${pa#*:}
  echo "--- $pot ---"
  "$LMP" -in in.elastic -var lat fcc -var a $a -var mass 26.9815 \
         -var nx 6 -var ny 6 -var nz 6 -var pot pot_al_${pot}.mod -log none 2>/dev/null \
    | grep -E 'C11all|C12all|C44all'
done

echo "===== Table 3: vacancy formation in fcc aluminum (each at its own a0) ====="
# equilibrium a0 from Table 1 (fcc aluminum)
for pa in "lj:4.0602" "ulj:3.3282" "morse:4.0398" "eam:4.0817" "meam:4.0500"; do
  pot=${pa%:*}; a=${pa#*:}
  "$LMP" -in vac.lmp -var lat fcc -var a $a -var mass 26.9815 \
         -var n 5 -var pot pot_al_${pot}.mod -log none 2>/dev/null | grep -E '^EVAC' \
    | sed "s/^/[$pot] /"
done
