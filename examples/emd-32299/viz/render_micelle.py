# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Diamond Light Source
# render_micelle.py -- the micelle test as pictures: rows = map (deposited map low-passed 15 A as the
# control, then the unmasked half-map average unfiltered and low-passed 6/10/15/20 A), columns = the
# four render levels micelle_maps.py wrote to levels.json (solvent median + k * solvent sd for the
# half-map average; fractions of the protein mean for the deposited map, whose solvent is exactly 0).
# Side view, membrane normal vertical, cytoplasm up. Surface coloured by the 7W4A chain within the
# zone distance (alpha blue, beta orange); density the model does not explain stays grey.
# Writes renders/micelle/<map>_<level>.png and renders/micelle_hero.png.
#   python3 -I viz/micelle_maps.py   (first: writes maps/ and levels.json)
#   chimerax --nogui --offscreen --exit --script viz/render_micelle.py
#   python3 -I viz/montage_micelle.py
import json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(sys.argv[0])))
from chimerax.core.commands import run
import emdviz as E

L = json.load(open(os.path.join(E.VIZ, "levels.json")))["levels"]
ROWS = ["deposited_lp15", "unfiltered", "lp6", "lp10", "lp15", "lp20"]
ZONE = {"unfiltered": 2.5, "lp6": 4, "lp10": 6, "lp15": 8, "lp20": 10, "deposited_lp15": 8}
out = os.path.join(E.RENDERS, "micelle")
os.makedirs(out, exist_ok=True)

run(session, f"open {E.DATA}/7w4a.cif")                 # #1
cen, up, perp = E.membrane_frame(session)
run(session, "set bgColor white; lighting soft; graphics silhouettes true width 1")
run(session, E.COLOURS.replace("#2", "#1") + "; hide #1 models")
for r in ROWS:
    run(session, f"open {E.VIZ}/maps/{'avg_' if not r.startswith('deposited') else ''}{r}.mrc")   # #2
    for key, lev in L[r].items():
        run(session, f"volume #2 level {lev} style surface step 1 color #c8c8c8")
        run(session, f"color zone #2 near #1 distance {ZONE[r]}")
        E.set_camera(session, cen, perp, up)
        run(session, f"save {out}/{r}_{key}.png width 520 height 520 supersample 3")
    run(session, "close #2")

# hero: deposited map at the author contour (coloured by chain) inside the 15 A half-map average
# at k=10, transparent grey -- the micelle as the context for the high-resolution map
run(session, f"open {E.DATA}/emd_32299.map")              # #2
run(session, f"volume #2 level {E.AUTHOR_LEVEL} style surface step 1; color zone #2 near #1 distance 2.2; surface dust #2 size 6")
run(session, f"open {E.VIZ}/maps/avg_lp15.mrc")           # #3
run(session, f"volume #3 level {L['lp15']['10']} style surface step 1 color #9e9e9e transparency 0.65")
run(session, "graphics silhouettes true width 1.2")
E.set_camera(session, cen, perp, up)
run(session, f"save {E.RENDERS}/micelle_hero.png width 1400 height 1400 supersample 3")
E.set_camera(session, cen, -up, perp)
run(session, f"save {E.RENDERS}/micelle_hero_top.png width 1400 height 1400 supersample 3")
