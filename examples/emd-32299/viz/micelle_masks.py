#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Diamond Light Source
"""micelle_masks.py -- model-free protein / micelle / envelope masks for RELION micelle handling.

Input: the two half maps of a refinement, in that refinement's frame. No atomic model is needed to
BUILD the masks; if one is given (--model, mmCIF) it is used only to SCORE them.

  1. avg = mean of the half maps; L15 = avg low-passed to 15 A.
  2. ENVELOPE (two passes): a first envelope with the solvent taken from the box's outer shell
     (r > RHO_SOLV of the half box; that shell is ~4x noisier than solvent near the particle, so the
     first envelope is tight); then SOLVENT = inside 0.95 of the half box and > SOLV_GAP_A from the
     first envelope, and ENVELOPE = L15 > solvent median + k_env * solvent sd, largest piece, holes
     filled.
  3. PROTEIN = ORDER, not density. Local correlation of the two half maps, band-passed 8-4 A, in a
     Gaussian window (sigma LCC_SIGMA_A). On EMD-32299 against 7W4A its median is 0.95 on protein,
     0.26 in the micelle belt, 0.07 around the soluble domains, 0.01 in solvent -- whereas DENSITY
     separates protein from detergent poorly (belt/protein ~0.6 at 6 A; a density/Otsu version grew
     a 'protein' mask 3-4x the protein's volume). Core = (lcc > lcc_t) inside the envelope, largest
     piece, holes filled, dilated by `dilate` A.
  4. MEMBRANE AXIS: longest inertia axis of the protein core. BAND: heights where (envelope -
     protein) has the largest cross-section (FWHM of its area profile).
  5. MICELLE = largest piece of (envelope - protein) within the band.
  6. Soft cosine edges (SOFT_A); micelle multiplied by (1 - protein), so the two are DISJOINT and
     protein + micelle <= 1 (RELION's --lowpass_mask ADDS the low-passed region to the
     solvent-masked map, so overlap would count the edge twice).

Outputs (default <data>/viz/masks/):
  protein_mask.mrc   relion_particle_subtract --mask (what is KEPT); or --solvent_mask with --lowpass_mask
  micelle_mask.mrc   --lowpass_mask or --solvent_mask2 region; what particle subtraction removes
  envelope_mask.mrc  protein + micelle: --solvent_mask when the micelle goes in --solvent_mask2
                     (mask2's mean is taken AFTER the main mask is applied, so the main mask must
                     include the micelle or the 'flattened' micelle is just zero)
  masks.json         levels, volumes, axis, band, and (with --model) the scores
    python3 -I viz/micelle_masks.py [--half1 H1 --half2 H2] [--model 7w4a.cif] [--out DIR]
                                    [--k_env 20] [--lcc_t 0.6] [--dilate 3]
"""
import argparse, json, os
import numpy as np
import mrcfile
from scipy import ndimage

K_ENV = 20.0
K_ENV0 = 10.0          # first pass only (outer-shell sd); 20 there leaves ~10 voxels
LCC_T = 0.6
DILATE_A = 2.0
LCC_SIGMA_A = 4.0
SOFT_A = 6.0
RHO_SOLV = 0.80
SOLV_GAP_A = 35.0
NONPROTEIN = ["PCW", "NAG", "8CK", "MG", "NA", "BEF", "HOH"]


def filt(v, apix, lo=None, hi=None):
    """Raised-cosine low-pass at `hi` A, optionally high-pass at `lo` A."""
    n = v.shape[0]
    f = np.fft.fftfreq(n, apix)
    fz, fy, fx = np.meshgrid(f, f, np.fft.rfftfreq(n, apix), indexing="ij", sparse=True)
    s = np.sqrt(fx**2 + fy**2 + fz**2)
    def edge(c):            # 1 below 1/c, 0 above
        f1, f2 = 0.85 / c, 1.15 / c
        w = np.clip((f2 - s) / (f2 - f1), 0, 1)
        return 0.5 - 0.5 * np.cos(np.pi * w)
    w = edge(hi)
    if lo is not None:
        w = w * (1 - edge(lo))
    return np.fft.irfftn(np.fft.rfftn(v) * w, s=v.shape).astype(np.float32)


def largest(b):
    lab, n = ndimage.label(b)
    if n == 0:
        return b
    return lab == (np.argmax(np.bincount(lab.ravel())[1:]) + 1)


def soft(b, apix, width):
    d = ndimage.distance_transform_edt(~b, sampling=apix)
    return np.where(d <= 0, 1.0, np.where(d >= width, 0.0, 0.5 + 0.5 * np.cos(np.pi * d / width))).astype(np.float32)


def write(path, v, apix):
    with mrcfile.new(path, overwrite=True) as m:
        m.set_data(v.astype(np.float32)); m.voxel_size = apix


def main():
    here = os.path.dirname(os.path.abspath(__file__)); data = os.path.dirname(here)
    ap = argparse.ArgumentParser()
    ap.add_argument("--half1", default=f"{data}/emd_32299_half_map_1.map")
    ap.add_argument("--half2", default=f"{data}/emd_32299_half_map_2.map")
    ap.add_argument("--model", default=f"{data}/7w4a.cif", help="mmCIF used ONLY to score; '' to skip")
    ap.add_argument("--out", default=f"{here}/masks")
    ap.add_argument("--k_env", type=float, default=K_ENV)
    ap.add_argument("--lcc_t", type=float, default=LCC_T)
    ap.add_argument("--dilate", type=float, default=DILATE_A)
    a = ap.parse_args()

    h1 = mrcfile.read(a.half1).astype(np.float32); h2 = mrcfile.read(a.half2).astype(np.float32)
    with mrcfile.open(a.half1, header_only=True) as m:
        apix = float(m.voxel_size.x)
    n = h1.shape[0]
    avg = (h1 + h2) / 2
    ax = (np.arange(n) - (n - 1) / 2) * apix
    Z, Y, X = np.meshgrid(ax, ax, ax, indexing="ij", sparse=True)
    r = np.sqrt(X**2 + Y**2 + Z**2)

    # 2. envelope
    L15 = filt(avg, apix, hi=15)
    outer = (r > RHO_SOLV * n * apix / 2) & (r < 0.95 * n * apix / 2)
    env0 = largest(L15 > float(np.median(L15[outer]) + K_ENV0 * L15[outer].std()))
    solv = (r < 0.95 * n * apix / 2) & (ndimage.distance_transform_edt(~env0, sampling=apix) > SOLV_GAP_A)
    env_level = float(np.median(L15[solv]) + a.k_env * L15[solv].std())
    env = ndimage.binary_fill_holes(largest(L15 > env_level))

    # 3. protein from local half-map correlation (order)
    b1, b2 = filt(h1, apix, lo=8, hi=4), filt(h2, apix, lo=8, hi=4)
    sg = LCC_SIGMA_A / apix
    lcc = ndimage.gaussian_filter(b1 * b2, sg) / (np.sqrt(ndimage.gaussian_filter(b1 * b1, sg) * ndimage.gaussian_filter(b2 * b2, sg)) + 1e-12)
    del b1, b2
    core = ndimage.binary_fill_holes(largest((lcc > a.lcc_t) & env))
    prot = ndimage.distance_transform_edt(~core, sampling=apix) <= a.dilate

    # 4-5. axis, band, micelle
    idx = np.argwhere(core) * apix
    cen = idx.mean(0)
    w, vec = np.linalg.eigh(np.cov((idx - cen).T))
    up = vec[:, -1]                                     # zyx
    g = np.arange(n) * apix
    zz, yy, xx = np.meshgrid(g, g, g, indexing="ij", sparse=True)
    h = (zz - cen[0]) * up[0] + (yy - cen[1]) * up[1] + (xx - cen[2]) * up[2]
    rest = env & ~prot
    hb = np.floor(h / 2.0).astype(int)
    hs = np.arange(hb.min(), hb.max() + 1)
    area = np.array([np.count_nonzero(rest & (hb == k)) for k in hs], float)
    ip = int(np.argmax(area)); lo = hi = ip
    while lo > 0 and area[lo - 1] >= area[ip] / 2: lo -= 1
    while hi < len(area) - 1 and area[hi + 1] >= area[ip] / 2: hi += 1
    mic = largest(rest & (hb >= hs[lo]) & (hb <= hs[hi]))

    # 6. soft, disjoint
    P = soft(prot, apix, SOFT_A)
    M = soft(mic, apix, SOFT_A) * (1 - P)
    E = np.maximum(soft(env, apix, SOFT_A), P + M)
    os.makedirs(a.out, exist_ok=True)
    for name, v in (("protein_mask", P), ("micelle_mask", M), ("envelope_mask", E)):
        write(f"{a.out}/{name}.mrc", v, apix)
    write(f"{a.out}/local_cc_8-4A.mrc", lcc, apix)

    rep = dict(source=[os.path.basename(a.half1), os.path.basename(a.half2)], apix=apix, k_env=a.k_env, lcc_t=a.lcc_t, dilate_A=a.dilate,
               lcc_sigma_A=LCC_SIGMA_A, soft_A=SOFT_A, envelope_level_L15=env_level,
               volumes_A3={k: float(b.sum() * apix**3) for k, b in (("envelope", env), ("protein", prot), ("micelle", mic))},
               axis_zyx=up.tolist(), band_A=[float(hs[lo] * 2.0), float((hs[hi] + 1) * 2.0)],
               max_P_plus_M=float((P + M).max()))
    if a.model:
        import sys
        sys.path.insert(0, here)
        from micelle_maps import read_atoms
        xyz, res, ch = read_atoms(a.model)
        def at(v, pts):
            ijk = np.clip(np.rint(pts / apix).astype(int), 0, n - 1)
            return v[ijk[:, 2], ijk[:, 1], ijk[:, 0]]
        pa = xyz[~np.isin(res, NONPROTEIN)]
        occ = np.zeros(h1.shape, bool)
        ijk = np.clip(np.rint(pa / apix).astype(int), 0, n - 1); occ[ijk[:, 2], ijk[:, 1], ijk[:, 0]] = True
        dist = ndimage.distance_transform_edt(~occ, sampling=apix)
        model_vol = float((dist <= 2.5).sum() * apix**3)
        Pb, Mb = P >= 0.5, M >= 0.5
        rep["model_score"] = dict(
            protein_atoms_in_P=float(np.mean(at(P, pa) >= 0.5)),           # sensitivity (want ~1)
            protein_atoms_touching_M=float(np.mean(at(M, pa) > 0.1)),      # want ~0
            P_volume_over_model_volume=float(Pb.sum() * apix**3 / model_vol),
            P_frac_further_than_6A=float(np.mean(dist[Pb] > 6)),           # false 'protein' (want small)
            M_frac_within_3A_of_protein=float(np.mean(dist[Mb] <= 3)) if Mb.any() else float("nan"),
            lipid_PCW_in_M=float(np.mean(at(M, xyz[res == "PCW"]) > 0.1)),
            model_volume_A3=model_vol)
    json.dump(rep, open(f"{a.out}/masks.json", "w"), indent=1)
    print(json.dumps(rep, indent=1))


if __name__ == "__main__":
    main()
