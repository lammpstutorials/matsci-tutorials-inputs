# parameterized init for the ML-appendix elastic-constants sweep
variable up equal 1.0e-6
variable atomjiggle equal 1.0e-5
units metal
variable cfac equal 1.0e-4
variable cunits string GPa
variable etol equal 0.0
variable ftol equal 1.0e-10
variable maxiter equal 1000
variable maxeval equal 10000
variable dmax equal 1.0e-2
boundary p p p
lattice ${lat} ${a}
region box prism 0 ${nx} 0 ${ny} 0 ${nz} 0 0 0
create_box 1 box
create_atoms 1 box
mass 1 ${mass}
