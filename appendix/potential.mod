include ${pot}
neighbor 1.0 bin
neigh_modify once no every 1 delay 0 check yes
min_style cg
min_modify dmax ${dmax} line quadratic
thermo 50
thermo_style custom step temp pe press pxx pyy pzz pxy pxz pyz lx ly lz vol
thermo_modify norm no
