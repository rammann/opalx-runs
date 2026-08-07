# rectangle

COLLIMATOR with `APERTURE="RECTANGLE(0.02,0.04)"` (half-apertures a=0.01 m, b=0.02 m), body s = 0.5 .. 0.6 m. 0.1 GeV electrons, no space charge.

49 particles on a fixed transverse grid with momenta exactly (0, 0, BG0), so (x, y) never changes. Expected: **20 particles deleted** inside the collimator, **29 survive** -- exactly the grid points inside the aperture (>= 20% boundary clearance by construction).

`run_tests.py` asserts the survivor (x, y) set, the final particle count in the .stat, and that the count drop happens inside the collimator body.
