# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Diamond Light Source
# emdviz.py -- shared helpers for the EMD-32299 ChimeraX scripts (imported inside ChimeraX).
#
# Layout assumed (the script folder sits beside the download):
#   <data>/emd_32299.map, emd_32299_half_map_{1,2}.map, emd_32299_msk_1.map, 7w4a.cif
#   <data>/viz/            these scripts
#   <data>/viz/renders/    PNG output
#   <data>/viz/maps/       derived maps (micelle_maps.py)
import os
import numpy as np
from chimerax.core.commands import run
from chimerax.geometry import Place

VIZ = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.dirname(VIZ)
RENDERS = os.path.join(VIZ, "renders")
MAPS = os.path.join(VIZ, "maps")
AUTHOR_LEVEL = 0.0298          # EMDB primary contour (author)

# colours: alpha subunit, beta subunit, revaprazan, glycans, lipid, ions
COLOURS = ("color #2/A #4f7fc9; color #2/B #e39b36; color #2:8CK #d81b60; "
           "color #2:NAG #43a047; color #2:PCW #9e9e9e; color #2:BEF,MG,NA #fdd835")


def model(session):
    from chimerax.atomic import AtomicStructure
    return [m for m in session.models if isinstance(m, AtomicStructure)][0]


def membrane_frame(session):
    """(centre, up, perp): longest inertia axis of 7W4A ~ membrane normal; 'up' points to the
    cytoplasmic headpiece (away from the beta ectodomain)."""
    m = model(session)
    xyz = m.atoms.scene_coords
    cen = xyz.mean(0)
    w, v = np.linalg.eigh(np.cov((xyz - cen).T))
    up = v[:, -1]
    beta = m.atoms[m.atoms.chain_ids == "B"].scene_coords.mean(0)
    if np.dot(beta - cen, up) > 0:
        up = -up
    return cen, up, v[:, 1]


def set_camera(session, cen, d, upv):
    d = d / np.linalg.norm(d)
    z = -d
    y = upv - np.dot(upv, z) * z
    y /= np.linalg.norm(y)
    x = np.cross(y, z)
    session.main_view.camera.position = Place(axes=np.array([x, y, z]), origin=cen - 400 * d)
    run(session, "view; view cofr false")


def save(session, name, w=1400, h=None):
    os.makedirs(RENDERS, exist_ok=True)
    run(session, f"save {os.path.join(RENDERS, name)} width {w} height {h or w} supersample 3")
