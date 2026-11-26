# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Diamond Light Source
# render_ligand.py -- revaprazan (8CK, /A:1101) density as mesh (step 1, zone 2.2 A).
# Writes renders/revaprazan_closeup.png (pocket residues within 4.5 A, labelled, slab +-9 A),
# revaprazan_density.png (author contour) and revaprazan_density_0p020.png (contour 0.020).
#   chimerax --nogui --offscreen --exit --script viz/render_ligand.py
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(sys.argv[0])))
from chimerax.core.commands import run
import emdviz as E

for c in [
    f"open {E.DATA}/emd_32299.map", f"open {E.DATA}/7w4a.cif",
    "set bgColor white; lighting simple; graphics silhouettes false",
    "hide #2 atoms; cartoon hide #2; show #2:8CK :<4.5 & ~H; style #2 stick",
    "color #2 & C #b0bec5; color #2:8CK & C #d81b60; color #2 byhetero",
    f"volume #1 level {E.AUTHOR_LEVEL} style mesh color #4f7fc9 squareMesh false step 1",
    "volume zone #1 nearAtoms #2:8CK range 2.2",
    "view #2:8CK :<4.5; view cofr false; clip near -9 far 9 position #2:8CK",
    "label #2:8CK :<4.5 & protein & ~:8CK residues text {0.name}{0.number} height 0.5 color #37474f",
]:
    run(session, c)
E.save(session, "revaprazan_closeup.png")
run(session, "label delete; hide #2:8CK :<4.5 & ~:8CK atoms; view #2:8CK; clip off")
E.save(session, "revaprazan_density.png", 1000)
run(session, "volume #1 level 0.020 step 1; volume zone #1 nearAtoms #2:8CK range 2.2")
E.save(session, "revaprazan_density_0p020.png", 1000)
