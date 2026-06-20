# Structure-analysis toolkit (appendix)

`hotsolid.lmp` compares common neighbor analysis (`compute cna/atom`) and
polyhedral template matching (`compute ptm/atom`) on a perfect fcc aluminum
crystal vibrating at 800 K (Mendelev potential, 2048 atoms, ~1 minute on a
few cores).  It prints the fcc-classified fraction for both methods
(CNA ~0.70, PTM ~0.996) and renders one image per method
(`cna-hot.png`, `ptm-hot.png`; blue = fcc, red = unclassified), which are
combined side by side into the figure of the appendix.

The `rdf-al.lmp` / `rdf-fe.lmp` / `rdf-mg.lmp` / `rdf-si.lmp` inputs compute
the radial distribution function of fcc Al, bcc Fe, hcp Mg (Mendelev EAM
potentials), and diamond Si (Stillinger-Weber) at three temperatures each;
the highest temperature lies above the superheating limit, so the crystal
melts during the run and the last curve is the liquid g(r).  `rdf-plot.py`
combines the twelve curves into the four-panel figure of the appendix
(each run a few minutes on a couple of cores).
