#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Diamond Light Source
"""montage_micelle.py -- assemble renders/micelle/*.png (render_micelle.py) into one labelled grid,
renders/micelle_grid.png: rows = map, columns = render level (from levels.json).
    python3 -I viz/montage_micelle.py
"""
import json, os
from PIL import Image, ImageDraw, ImageFont

VIZ = os.path.dirname(os.path.abspath(__file__))
L = json.load(open(f"{VIZ}/levels.json"))
ROWS = [("deposited_lp15", "deposited map\n(masked, sharpened)\nlow-pass 15 Å"),
        ("unfiltered", "half-map average\nunfiltered"),
        ("lp6", "half-map average\nlow-pass 6 Å"), ("lp10", "half-map average\nlow-pass 10 Å"),
        ("lp15", "half-map average\nlow-pass 15 Å"), ("lp20", "half-map average\nlow-pass 20 Å")]
CELL, LW, TH = 360, 210, 70
font = lambda s: ImageFont.truetype("DejaVuSans.ttf", s)
img = Image.new("RGB", (LW + 4 * CELL, TH + len(ROWS) * CELL), "#fcfcfb")
d = ImageDraw.Draw(img)
d.text((LW, 8), "Threshold → (half-map rows: solvent median + k × solvent sd;  deposited row: fraction of its protein mean)",
       fill="#52514e", font=font(15))
for j, k in enumerate(["5", "10", "20", "40"]):
    d.text((LW + j * CELL + CELL // 2, 44), f"k = {k}", fill="#0b0b0b", font=font(18), anchor="mm")
for i, (r, label) in enumerate(ROWS):
    y = TH + i * CELL
    d.multiline_text((10, y + CELL // 2), label, fill="#0b0b0b", font=font(16), anchor="lm", spacing=4)
    for j, key in enumerate(L["levels"][r]):
        im = Image.open(f"{VIZ}/renders/micelle/{r}_{key}.png").convert("RGB").resize((CELL, CELL), Image.LANCZOS)
        img.paste(im, (LW + j * CELL, y))
        lev = L["levels"][r][key]
        sub = key.replace("xprot", " × protein") if r.startswith("deposited") else f"level {lev:.4f}"
        d.text((LW + j * CELL + CELL // 2, y + CELL - 16), sub + (f"  ({lev:.4f})" if r.startswith("deposited") else ""),
               fill="#52514e", font=font(13), anchor="mm")
img.save(f"{VIZ}/renders/micelle_grid.png")
print("wrote", f"{VIZ}/renders/micelle_grid.png", img.size)
