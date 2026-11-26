# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Diamond Light Source
# render_overview.py -- EMD-32299 primary map at the author contour, coloured by PDB 7W4A chain
# (zone 2.2 A). Writes renders/side_0.png, side_90.png, top_cytoplasm.png.
#   chimerax --nogui --offscreen --exit --script viz/render_overview.py
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(sys.argv[0])))
from chimerax.core.commands import run
import emdviz as E

run(session, f"open {E.DATA}/emd_32299.map")
run(session, f"open {E.DATA}/7w4a.cif")
cen, up, perp = E.membrane_frame(session)
run(session, "set bgColor white; lighting soft; graphics silhouettes true width 1.5")
run(session, f"volume #1 level {E.AUTHOR_LEVEL} style surface step 1 color #d9d9d9")
run(session, E.COLOURS)
run(session, "color zone #1 near #2 distance 2.2")
run(session, "surface dust #1 size 6")
run(session, "hide #2 models")
E.set_camera(session, cen, perp, up);                 E.save(session, "side_0.png")
E.set_camera(session, cen, np.cross(up, perp), up);   E.save(session, "side_90.png")
E.set_camera(session, cen, -up, perp);                E.save(session, "top_cytoplasm.png")
