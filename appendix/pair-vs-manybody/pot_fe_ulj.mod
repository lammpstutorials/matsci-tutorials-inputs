# Elliott-Akerson "universal" LJ for Fe, reproduced with lj/cut + shift yes.
# Extracted from OpenKIM model LennardJones612_UniversalShifted__MO_959249795837_003
# via a dimer scan (dimer.lmp): r_min = 2*r_cov(Fe) = 2.64 A; sigma = 2.352 A;
# eps = 1.184 eV; cutoff = 4*sigma.  Reproduces pair_style kim to < 0.03%.
pair_style lj/cut 9.408
pair_coeff 1 1 1.184 2.352
pair_modify shift yes
