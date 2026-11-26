# Worked example: EMD-32299

EMD-32299 is the gastric H+,K+-ATPase (pig) with revaprazan, at 2.76 Å; its fitted model is PDB
7W4A (Tanaka *et al.* 2022, [doi:10.1021/acs.jmedchem.2c00338](https://doi.org/10.1021/acs.jmedchem.2c00338)).
The protein sits in a detergent belt. This example asks two things:
- can the belt be made visible?
- can it be masked cleanly enough for RELION's micelle options ([`docs`](../../docs/relion_micelle_options.md))?

## Reproduce

```bash
bash fetch.sh                                    # maps, half maps, mask, model (~320 MB)
python3 -I viz/micelle_maps.py                   # belt profile and threshold test   (~20 s)
python3 -I viz/micelle_masks.py                  # protein / micelle / envelope masks (~25 s)
python3 -I viz/score_masks_belt.py               # where the measured belt lands
python3 -I viz/sweep_masks.py                    # the 18-setting mask sweep          (~13 min)
chimerax --nogui --offscreen --exit --script viz/render_overview.py
chimerax --nogui --offscreen --exit --script viz/render_ligand.py
chimerax --nogui --offscreen --exit --script viz/render_micelle.py   # needs micelle_maps.py first
python3 -I viz/montage_micelle.py
chimerax --nogui --offscreen --exit --script viz/render_masks.py     # needs micelle_masks.py first
```

All commands run from this folder. Every path is taken from the script's own location.

**Requirements:**
- Python 3 with numpy, scipy, mrcfile, matplotlib and Pillow;
- UCSF ChimeraX (tested with 1.12; offscreen rendering, no display needed).

**Outputs:**
- figures go to `viz/renders/`;
- derived maps go to `viz/maps/` and masks to `viz/masks/`;
- the committed copies are in `figures/` and `results/`.

## Figures

| Figure | Script | What it shows |
|---|---|---|
| `side_0.png`, `side_90.png`, `top_cytoplasm.png` | `render_overview.py` | The deposited map at the author contour (0.0298), coloured by 7W4A chain within 2.2 Å: α blue, β orange, revaprazan magenta, glycans green. The membrane normal is vertical, cytoplasm up. |
| `revaprazan_closeup.png`, `revaprazan_density.png` | `render_ligand.py` | Revaprazan's density as mesh, with the pocket residues within 4.5 Å. |
| `micelle_grid.png` | `render_micelle.py`, `montage_micelle.py` | Rows are maps and columns are thresholds; unexplained density stays grey. |
| `micelle_profiles.png` | `micelle_maps.py` | Density 10–25 Å outside the model by height, and whether the belt stays closed as the threshold rises. |
| `micelle_hero.png`, `micelle_hero_top.png` | `render_micelle.py` | The deposited map inside the 15 Å half-map average, shown transparent. |
| `masks_side_top.png` | `render_masks.py` | The protein mask (blue) and the micelle mask (orange) at 0.5. |

## Findings

1. **The deposited map hides the belt.** It is masked (92% of voxels are exactly 0) and sharpened,
   and at the author contour the belt does not show. The deposited mask does contain 98% of the
   belt: low-passing the deposited map to 15 Å reveals it, clipped at the mask edge.
2. **Low-pass, not threshold.** The average of the unmasked half maps, low-passed to 6 Å or more,
   shows a closed belt, about 44 Å thick (FWHM, −42 to +2 Å along the membrane normal), around the
   bound lipid and the drug. Unfiltered, the belt is speckle at any threshold.

   | Half-map average | Belt peak above solvent | Belt still closed at |
   |---|---|---|
   | unfiltered | 4.9 sd | only speckle |
   | low-pass 6 Å | 24 sd | up to 12 sd |
   | low-pass 10 Å | 40 sd | up to 20 sd |
   | low-pass 15–20 Å | 49 sd | up to 20–30 sd |

3. **Order separates protein from micelle; density does not.** The local half-map correlation,
   band-passed 8–4 Å, has a median of 0.95 on protein, 0.26 in the belt, 0.07 around the soluble
   domains and 0.01 in solvent. The masks built on it are scored in
   [`docs`](../../docs/relion_micelle_options.md#4-the-masks) and in `results/`.

## Notes on the measurement

`micelle_maps.py` measures a continuous profile: the mean density in a shell 10–25 Å outside the
model, by height, relative to bulk solvent.

**Counting thresholded density does not work.** Several earlier designs counted thresholded
"unexplained" density, and every one saturated: each reported a full belt at every height and
setting.
- **The atom-distance criterion fails after a low-pass.** "Further than 3 Å from atoms" stops
  meaning anything once a 10–20 Å low-pass spreads the protein's own density 5–10 Å beyond its
  atoms.
- **Residuals inside the protein fill every sector.** Any residual inside the protein fills every
  azimuthal sector of every slice that cuts through it.

**Scoring masks on coverage alone rewards masks that are too big.** Scoring a mask only by how
many protein atoms it covers rewarded oversized masks. The mask sweep therefore also scores how
large the mask is relative to the model, and how much of it lies more than 6 Å from any atom.
