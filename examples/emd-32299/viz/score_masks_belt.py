#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Diamond Light Source
"""score_masks_belt.py -- where does the measured micelle belt land in micelle_masks.py's masks?
Belt = voxels 10-25 A from 7W4A atoms at heights -42..+2 A along the membrane normal (the FWHM band
from micelle_maps.py, levels.json; + = cytoplasm). Reports the fraction of that belt in the micelle
mask, in the protein mask (bad: it would be KEPT by subtraction), and in neither.
    python3 -I viz/score_masks_belt.py [masks_dir]
"""
import json, os, sys
import numpy as np, mrcfile
from scipy import ndimage
here = os.path.dirname(os.path.abspath(__file__)); data = os.path.dirname(here)
sys.path.insert(0, here)
from micelle_maps import read_atoms
md = sys.argv[1] if len(sys.argv) > 1 else f"{here}/masks"
P = mrcfile.read(f"{md}/protein_mask.mrc"); M = mrcfile.read(f"{md}/micelle_mask.mrc")
apix = 0.752; n = P.shape[0]
L = json.load(open(f"{here}/levels.json")); cen = np.array(L["centre"]); up = np.array(L["up"])
lo, hi = L["summary"]["lp15"]["band_A"]
xyz, res, ch = read_atoms(f"{data}/7w4a.cif")
occ = np.zeros(P.shape, bool); ijk = np.clip(np.rint(xyz / apix).astype(int), 0, n - 1); occ[ijk[:, 2], ijk[:, 1], ijk[:, 0]] = True
dist = ndimage.distance_transform_edt(~occ, sampling=apix)
g = np.arange(n) * apix; Z, Y, X = np.meshgrid(g, g, g, indexing="ij", sparse=True)
h = (X - cen[0]) * up[0] + (Y - cen[1]) * up[1] + (Z - cen[2]) * up[2]
belt = (dist > 10) & (dist <= 25) & (h >= lo) & (h < hi)
out = dict(belt_band_A=[lo, hi], belt_voxels=int(belt.sum()),
           belt_in_micelle_mask=float(np.mean(M[belt] >= 0.5)), belt_in_protein_mask=float(np.mean(P[belt] >= 0.5)),
           belt_in_neither=float(np.mean((M[belt] < 0.5) & (P[belt] < 0.5))),
           micelle_mask_in_belt_band=float(np.mean(((h >= lo - 5) & (h < hi + 5))[M >= 0.5])))
json.dump(out, open(f"{md}/belt_score.json", "w"), indent=1); print(json.dumps(out, indent=1))
