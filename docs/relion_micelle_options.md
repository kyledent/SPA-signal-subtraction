# Handling a detergent micelle inside a RELION refinement

How can the signal of a detergent micelle, nanodisc or amphipol belt be removed or tamed inside a
RELION single-particle refinement? This note lists what RELION already provides, what the
literature shows, and a test plan.

**Sources:**
- **RELION.** Read at upstream branch `ver5.1`, commit `afdd8167`, and each option checked as
  present in tag `5.0.0`; line numbers are from `ver5.1`.
- **Literature.** Papers read in full.
- **Measurements.** From the worked example in [`examples/emd-32299`](../examples/emd-32299).

## Summary

- **RELION 5 already has five mechanisms** that act on the micelle; four of them need no code
  change. The fifth, multi-body refinement with a fixed micelle body, reduces to a one-off
  subtraction.
- **No published study shows that subtracting a micelle from the particle images improves
  resolution.** The clearest gains for membrane proteins in detergent come from regularisation
  (Blush in RELION; non-uniform refinement in cryoSPARC).
- **Measure on your own data, cheapest options first.** Start with the in-refinement reference
  options (`--lowpass_mask`, `--solvent_mask2`, `--blush`). Then try particle subtraction followed
  by a local re-refinement.
- **Every option needs a protein mask and a micelle mask** in the refinement's own frame, and the
  two must not overlap. [`micelle_masks.py`](../examples/emd-32299/viz/micelle_masks.py) builds
  them from the two half maps alone (§4).

## 1. What RELION 5 offers

| # | Mechanism | Acts on | What it does | Notes |
|---|---|---|---|---|
| A | **Particle subtraction**: `relion_particle_subtract` (GUI: "Particle subtraction") | Images, once, after a refinement | Projects `reference × (1 − keep_mask)` with each particle's orientation and CTF, and subtracts it. `--mask` is what is **kept** (`particle_subtractor.cpp` l.29, 224) | MPI ranks are split by half-set (l.65–68, 475). Needs `rlnImageOriginalName`. `--revert` restores the original images, which CTF refinement and Bayesian polishing need. Options `--recenter_on_mask` and `--new_box`. Re-refine with **local** searches afterwards. |
| B1 | Tight solvent mask: `--solvent_mask protein --flatten_solvent` | Reference, every iteration | Zeroes the micelle in the reference | The usual practice. The images still contain the micelle, so it acts as unmodelled signal. |
| B2 | `--solvent_mask2 M` (only via "Additional arguments") | Reference, every iteration | Sets region M to **its own mean density** (`softMaskOutsideMap(..., invert=true)`, `mask.cpp`) | **Ordering trap:** `solventFlatten()` applies the main solvent mask **first**, so if the main mask excludes the micelle, M's "mean" is 0. Use `--solvent_mask envelope --solvent_mask2 micelle`. |
| B3 | **`--lowpass_mask M --lowpass 15`** (not in the GUI; added in 2018, commit `478df5df`) | Reference, every iteration | `Iref = Iref·solvent_mask + lowpass(Iref·M)` | Keeps the micelle's low-resolution shape, which helps alignment, and drops its high-resolution noise; Blush (C) aims for the same treatment. The solvent mask must be protein-only and **disjoint** from M, otherwise the edge is counted twice. **Not saved in `optimiser.star`**: pass it again with `--continue`. |
| C | **Blush**: `--blush` | Reconstruction, every iteration | A denoiser trained on EMDB half maps | The paper describes the denoiser **smoothing, not removing**, micelle density: its training masks keep a 15 Å version of disordered regions. It uses the external-reconstruct path, so it cannot be combined with E. |
| D | Multi-body refinement with a micelle body | Images, every E-step | Subtracts the other bodies' projections on the fly | Bodies under ~100 kDa cannot be aligned. A body with all sigmas at 0 is kept fixed (`keep_fixed_bodies`, `ml_optimiser.cpp` l.4433, 10155): it is never aligned or reconstructed, and only its **consensus** density is subtracted. That equals A followed by a local refinement. |
| E | `--external_reconstruct` with `$RELION_EXTERNAL_RECONSTRUCT_EXECUTABLE` | Reconstruction, every iteration, per half | RELION writes the data and weight terms and runs the program, which returns the map | A hook for prototyping per-iteration micelle models (low-pass, flattening, a parametric toroid) without changing C++. |

**Not available:** a micelle that follows the protein's *refined* pose and is rebuilt from the
current reconstruction, subtracted from the images at every iteration. Multi-body refinement holds a
fixed body at the consensus pose, and particle subtraction uses one fixed reference and pose.
Adding it would mean a `--subtract_mask` option in `relion_refine` that reuses the multi-body
subtraction in `getFourierTransformsAndCtfs`.

## 2. Literature

- **Bai *et al.* 2015**, eLife 4:e11182, [doi:10.7554/eLife.11182](https://doi.org/10.7554/eLife.11182).
  Introduced signal subtraction.
  - The masked consensus map is projected with each particle's consensus orientation and CTF, and
    subtracted. No scaling is described.
  - It was used for **classification with orientations fixed**, and what was subtracted was protein
    plus the amphipol belt.
  - No resolution gain is attributed to the subtraction.
- **Nakane *et al.* 2018**, eLife 7:e36861, [doi:10.7554/eLife.36861](https://doi.org/10.7554/eLife.36861).
  Multi-body refinement.
  - The subtraction is recomputed at every E-step, using the other bodies' previous poses.
  - Bodies under ~100 kDa become hard to align.
  - Membranes are not discussed.
- **Kimanius *et al.* 2024**, Nat Methods 21:1216, [doi:10.1038/s41592-024-02304-8](https://doi.org/10.1038/s41592-024-02304-8).
  Blush regularisation.
  - Ste2 (≈67 kDa ordered, mostly within a detergent micelle) improved from 3.8 to 3.4 Å, and
    PfCRT (EMPIAR-10330) from 3.8 to 3.4 Å.
- **Punjani *et al.* 2020**, Nat Methods 17:1214, [doi:10.1038/s41592-020-00990-8](https://doi.org/10.1038/s41592-020-00990-8).
  Non-uniform refinement.
  - Adaptive regularisation, motivated by micelles and nanodiscs, in cryoSPARC; it is not
    subtraction.
- **Bharadwaj, Veerbeek & Jakobi 2026**, Acta Cryst D 82:848, [doi:10.1107/S205979832600598X](https://doi.org/10.1107/S205979832600598X).
  SURFER.
  - A ChimeraX network that segments membrane-mimic density from unfiltered half maps; an
    alternative way to make the micelle mask.
  - The authors note that its use for particle subtraction has not been tested.
- **Scheres 2016**, Methods Enzymol 579:125, [doi:10.1016/bs.mie.2016.04.012](https://doi.org/10.1016/bs.mie.2016.04.012).
  Practice for masked and subtracted classification in RELION.

**Gap:** none of these reports a gain in alignment or resolution from subtracting a micelle or
nanodisc from the particle images. Treat it as a hypothesis to test.

## 3. Why the micelle matters, and why it is separable

Measured on EMD-32299 (gastric H+,K+-ATPase, about 147 kDa of protein in detergent;
[worked example](../examples/emd-32299)):
- **The deposited map hides the belt.** The deposited map is masked and sharpened, and at the
  author contour the belt does not show.
- **The half maps show it plainly.** Low-passed to 6–20 Å, their average shows a belt 24–49
  solvent sd above bulk solvent, about 44 Å thick (FWHM).
- **The belt carries no high-resolution signal.** The local correlation between the two half
  maps, band-passed 8–4 Å, has a median of 0.95 on the protein but 0.26 in the belt. That makes
  the micelle low-resolution structured signal, which is what every option above treats.

## 4. The masks

`micelle_masks.py` takes two half maps and writes:
- `protein_mask.mrc`: what `relion_particle_subtract --mask` keeps, or the `--solvent_mask` to use
  with `--lowpass_mask`;
- `micelle_mask.mrc`: the region for `--lowpass_mask` or `--solvent_mask2`, or what subtraction
  removes;
- `envelope_mask.mrc`: protein plus micelle, the `--solvent_mask` to use with `--solvent_mask2`.

All three are soft-edged, and the protein and micelle masks are disjoint (P + M ≤ 1). No atomic
model is needed. How the masks are built:
- **Protein is defined by order, not density:** where the half maps correlate locally at 8–4 Å. A
  density threshold (Otsu) separates protein from detergent poorly, because the belt is about
  0.6× the protein's density at 6 Å. On this map it produced a "protein" mask 3–4× the protein's
  volume.
- **Micelle** is the rest of the 15 Å envelope, inside the belt's height band along the membrane
  normal.

Validated on EMD-32299 against PDB 7W4A, which is used only to score (`sweep_masks.py`, 18
settings; defaults `k_env 20, lcc_t 0.6, dilate 2`):

| Score | Value |
|---|---|
| protein atoms inside the protein mask | 99.3% |
| protein atoms touching the micelle mask | 0% |
| micelle mask within 3 Å of protein | 0% |
| measured belt in micelle / protein / neither mask | 76% / 19% / 5% |
| micelle mask inside the belt's height band | 99% |
| protein-mask volume / model volume | 2.5× (a generous keep-mask) |

Both kinds of miss are conservative:
- the belt's inner rim stays inside the padded protein mask, so subtraction keeps it;
- the belt's tails beyond its FWHM fall in neither mask.

## 5. A test plan

Use paired comparisons: the same particles, the same random half-set split (`rlnRandomSubset`),
the same initial model and the same final mask throughout.

1. **Baseline:** refine with a loose protein + micelle envelope mask.
2. A tight protein-only solvent mask (B1).
3. `--solvent_mask envelope --solvent_mask2 micelle` (B2).
4. `--solvent_mask protein --lowpass_mask micelle --lowpass 15` (B3).
5. The baseline with `--blush` (C).
6. Particle subtraction with the protein mask as the keep-mask (A). Then:
   1. re-refine with local searches only;
   2. `--revert` for CTF refinement and polishing;
   3. refine again.
7. Only if 6 beats 4 and 5: per-iteration subtraction, through E or a new `--subtract_mask`.

**Measure:** FSC 0.143 within the same protein mask; local resolution in the membrane domain;
ligand density; RELION's reported angular accuracy; and how closely the per-particle poses agree
with the baseline.

**Practical points:**
- **Rebuild the masks for every refinement.** Use that job's own unfiltered half maps
  (`run_half{1,2}_class001_unfil.mrc`); masks do not transfer between frames or box sizes.
- **Subtraction is only as good as the poses.** Poor poses leave residual micelle and can remove
  protein, so subtract after a converged refinement and re-refine locally.
