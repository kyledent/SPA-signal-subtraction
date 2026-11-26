# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Diamond Light Source
# render_masks.py -- micelle_masks.py's protein (blue) and micelle (orange) masks at 0.5, with the
# deposited map at the author contour inside (grey). Side and top views.
# Writes renders/masks_side.png, renders/masks_top.png.
#   python3 -I viz/micelle_masks.py      (first)
#   chimerax --nogui --offscreen --exit --script viz/render_masks.py
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(sys.argv[0])))
from chimerax.core.commands import run
import emdviz as E

run(session, f"open {E.DATA}/7w4a.cif")                                   # #1 (frame only)
cen, up, perp = E.membrane_frame(session)
run(session, "hide #1 models; set bgColor white; lighting soft; graphics silhouettes true width 1.2")
run(session, f"open {E.DATA}/emd_32299.map")                              # #2
run(session, f"volume #2 level {E.AUTHOR_LEVEL} style surface step 1 color #6b6b6b; surface dust #2 size 6")
run(session, f"open {E.VIZ}/masks/protein_mask.mrc")                      # #3
run(session, "volume #3 level 0.5 style surface step 1 color #4f7fc9 transparency 0.7")
run(session, f"open {E.VIZ}/masks/micelle_mask.mrc")                      # #4
run(session, "volume #4 level 0.5 style surface step 1 color #e39b36 transparency 0.35")
E.set_camera(session, cen, perp, up);  E.save(session, "masks_side.png", 1200)
E.set_camera(session, cen, -up, perp); E.save(session, "masks_top.png", 1200)
