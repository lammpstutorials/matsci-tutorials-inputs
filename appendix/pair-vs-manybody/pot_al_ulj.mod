# Elliott-Akerson "universal" LJ for Al, reproduced with lj/cut + shift yes.
# Extracted from OpenKIM model LennardJones612_UniversalShifted__MO_959249795837_003
# via a dimer scan (dimer.lmp): r_min = 2*r_cov(Al) = 2.42 A; sigma = r_min/2^(1/6)
# = 2.156 A; eps = dimer well depth = 2.700 eV; cutoff = 4*sigma.  Reproduces
# pair_style kim to < 0.03% on the fcc cohesive energy.
pair_style lj/cut 8.624
pair_coeff 1 1 2.700 2.156
pair_modify shift yes
