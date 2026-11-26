#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Diamond Light Source
"""sweep_masks.py -- choose micelle_masks.py's envelope level, local-CC threshold and dilation by
scoring against PDB 7W4A on EMD-32299 (the model is used to score only). Writes masks_sweep.tsv.
Scores: SENSITIVITY (protein atoms inside P) AND SPECIFICITY (P volume / model volume, fraction of P
> 6 A from atoms). A first version scored sensitivity only and so rewarded masks 3-4x too big.
    python3 -I viz/sweep_masks.py
"""
import itertools, json, os, subprocess, sys
here = os.path.dirname(os.path.abspath(__file__))
rows = []
for k_env, t, dil in itertools.product([10, 20, 40], [0.5, 0.6, 0.7], [2, 4]):
    out = f"{here}/masks_sweep/k{k_env}_t{t}_d{dil}"
    subprocess.run([sys.executable, "-I", f"{here}/micelle_masks.py", "--k_env", str(k_env), "--lcc_t", str(t),
                    "--dilate", str(dil), "--out", out], check=True, capture_output=True)
    r = json.load(open(f"{out}/masks.json")); s = r["model_score"]; v = r["volumes_A3"]
    rows.append((k_env, t, dil, v["envelope"], v["protein"], v["micelle"], r["band_A"][0], r["band_A"][1],
                 s["protein_atoms_in_P"], s["protein_atoms_touching_M"], s["P_volume_over_model_volume"],
                 s["P_frac_further_than_6A"], s["M_frac_within_3A_of_protein"], s["lipid_PCW_in_M"]))
    print("k_env %2d lcc_t %.1f dil %d | env %7.0f P %7.0f M %7.0f | band %+.0f..%+.0f | atoms in P %.3f, touching M %.3f | P/model %.2f, P>6A %.3f | M<=3A %.3f | PCW in M %.2f" % rows[-1], flush=True)
with open(f"{here}/masks_sweep.tsv", "w") as f:
    f.write("k_env\tlcc_t\tdilate_A\tenvelope_A3\tprotein_A3\tmicelle_A3\tband_lo_A\tband_hi_A\tatoms_in_P\tatoms_touching_M\tP_over_model\tP_far_6A\tM_within_3A\tPCW_in_M\n")
    for r in rows:
        f.write("\t".join(str(round(x, 4)) if isinstance(x, float) else str(x) for x in r) + "\n")
