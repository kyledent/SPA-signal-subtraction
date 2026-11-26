#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Diamond Light Source
"""micelle_maps.py -- can the EMD-32299 detergent micelle be made visible by low-pass and threshold?

The deposited primary map is masked (zero outside emd_32299_msk_1: 92% of voxels are exactly 0), so
nothing outside the mask can be recovered from it. This works on the average of the two UNMASKED,
unsharpened half maps.

Measurement (continuous, so it cannot saturate):
  * SHELL: voxels 10-25 A from the nearest 7W4A atom (beyond the protein's own blur at these
    low-passes) and inside the inscribed sphere. SOLVENT: > 35 A from any atom, inside the sphere.
  * PROFILE: the mean of the shell, in 2 A slices along the membrane normal (longest inertia axis of
    7W4A; + = cytoplasm), as (shell mean - solvent median) / solvent sd. A micelle is a band of
    large excess at the membrane (where the bound lipid PCW and revaprazan sit); the soluble domains
    sit near 0.
  * BELT vs THRESHOLD: within the band (FWHM of the profile), the shell is split into 36 azimuthal
    sectors; coverage(k) = fraction of sectors whose mean exceeds solvent median + k * solvent sd.
    Coverage 1 = the belt is closed at that threshold.

Record of what did NOT work (each saturated, reading as a full belt at every height and setting):
thresholded "unexplained density" (> 3 A from atoms) with thresholds from the box corners, from the
half-map-difference noise, or from the solvent; and thresholded model-subtracted residuals. After a
10-20 A low-pass the protein's own density spreads 5-10 A past its atoms, and any residual inside
the protein fills every azimuthal sector of every slice through it.

Control: the deposited map, unfiltered and low-passed to 15 A, judged by the belt-to-protein ratio
(mean over the band shell / mean within 2 A of atoms), which needs no solvent reference; the same
ratio is reported for every half-map average.

Writes (into viz/): maps/avg_unfiltered.mrc, maps/avg_lp<c>.mrc, maps/deposited_lp15.mrc, micelle_profile.tsv,
micelle_belt.tsv, levels.json (render levels), renders/micelle_profiles.png.
    python3 -I viz/micelle_maps.py [DATA_DIR]
"""
import json, os, sys
import numpy as np
import mrcfile
from scipy import ndimage

VIZ = os.path.dirname(os.path.abspath(__file__))
DATA = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(VIZ)
APIX = 0.752
CUTOFFS = [None, 6, 10, 15, 20]           # low-pass cutoff, A (None = unfiltered half-map average)
SHELL_A = (10.0, 25.0)
SOLVENT_A = 35.0
SLICE_A, NSECT = 2.0, 36
KS = [1, 2, 3, 5, 8, 12, 20, 30, 50]      # belt thresholds, x solvent sd
KS_RENDER = [5, 10, 20, 40]               # levels written to levels.json for render_micelle.py


def tag(c):
    return "unfiltered" if c is None else f"lp{c}"


def read_atoms(cif):
    """(xyz, resname, auth chain) of non-water atoms from the _atom_site loop."""
    cols, rows, inloop = [], [], False
    for line in open(cif):
        if line.startswith("_atom_site."):
            cols.append(line.split(".", 1)[1].strip()); inloop = True; continue
        if inloop:
            if line.startswith(("ATOM", "HETATM")):
                rows.append(line.split())
            elif rows and line.startswith("#"):
                break
    c = {k: cols.index(k) for k in ("Cartn_x", "Cartn_y", "Cartn_z", "label_comp_id", "auth_asym_id")}
    rows = [r for r in rows if r[c["label_comp_id"]] != "HOH"]
    xyz = np.array([[float(r[c["Cartn_x"]]), float(r[c["Cartn_y"]]), float(r[c["Cartn_z"]])] for r in rows])
    return xyz, np.array([r[c["label_comp_id"]] for r in rows]), np.array([r[c["auth_asym_id"]] for r in rows])


def lowpass(v, cutoff):
    if cutoff is None:
        return v
    n = v.shape[0]
    f = np.fft.fftfreq(n, APIX)
    fz, fy, fx = np.meshgrid(f, f, np.fft.rfftfreq(n, APIX), indexing="ij", sparse=True)
    s = np.sqrt(fx**2 + fy**2 + fz**2)
    f1, f2 = 0.85 / cutoff, 1.15 / cutoff                    # raised-cosine edge around 1/cutoff
    w = np.clip((f2 - s) / (f2 - f1), 0, 1)
    w = 0.5 - 0.5 * np.cos(np.pi * w)
    return np.fft.irfftn(np.fft.rfftn(v) * w, s=v.shape).astype(np.float32)


def main():
    avg = (mrcfile.read(f"{DATA}/emd_32299_half_map_1.map").astype(np.float32)
           + mrcfile.read(f"{DATA}/emd_32299_half_map_2.map")) / 2
    dmask = mrcfile.read(f"{DATA}/emd_32299_msk_1.map")
    n = avg.shape[0]

    xyz, res, ch = read_atoms(f"{DATA}/7w4a.cif")
    cen = xyz.mean(0)
    w, vec = np.linalg.eigh(np.cov((xyz - cen).T))
    up, perp = vec[:, -1], vec[:, 1]
    if np.dot(xyz[ch == "B"].mean(0) - cen, up) > 0:       # 'up' = away from the beta ectodomain
        up = -up
    third = np.cross(up, perp)
    marks = {r: float(np.mean((xyz[res == r] - cen) @ up)) for r in ("PCW", "8CK")}

    occ = np.zeros(avg.shape, bool)
    ijk = np.clip(np.rint(xyz / APIX).astype(int), 0, n - 1)
    occ[ijk[:, 2], ijk[:, 1], ijk[:, 0]] = True
    dist = ndimage.distance_transform_edt(~occ, sampling=APIX)
    del occ
    ax = np.arange(n) * APIX
    Z, Y, X = np.meshgrid(ax, ax, ax, indexing="ij")
    relx, rely, relz = X - cen[0], Y - cen[1], Z - cen[2]
    hgt = relx * up[0] + rely * up[1] + relz * up[2]
    a = relx * perp[0] + rely * perp[1] + relz * perp[2]
    b = relx * third[0] + rely * third[1] + relz * third[2]
    sect = ((np.arctan2(b, a) + np.pi) / (2 * np.pi) * NSECT).astype(int) % NSECT
    c0 = (n - 1) * APIX / 2
    inner = (X - c0)**2 + (Y - c0)**2 + (Z - c0)**2 < (0.95 * n * APIX / 2)**2
    del X, Y, Z, relx, rely, relz, a, b
    solvent = inner & (dist > SOLVENT_A)
    shell = inner & (dist > SHELL_A[0]) & (dist <= SHELL_A[1])
    hb = np.floor(hgt[shell] / SLICE_A).astype(int)
    sh_sect, sh_h = sect[shell], hgt[shell]
    hbins = np.arange(-40, 41)
    in_mask_shell = dmask[shell] > 0.5

    near = dist <= 2.0
    os.makedirs(f"{VIZ}/maps", exist_ok=True)
    levels, prof_rows, belt_rows, summary = {}, [], [], {}
    # control: the DEPOSITED map (masked + sharpened). Its solvent is exactly 0, so it is judged only
    # by the belt-to-protein ratio, which needs no solvent reference.
    prim = mrcfile.read(f"{DATA}/emd_32299.map").astype(np.float32)
    for c in (None, 15):
        v = lowpass(prim, c)
        t = "deposited_" + tag(c)
        if c is not None:
            with mrcfile.new(f"{VIZ}/maps/{t}.mrc", overwrite=True) as m:
                m.set_data(v); m.voxel_size = APIX
        band_sh = shell & (hgt >= -42) & (hgt < 2)
        ratio = float(v[band_sh].mean() / v[near].mean())
        summary[t] = dict(belt_to_protein_ratio=ratio, protein_mean=float(v[near].mean()))
        if c is not None:      # no solvent sd here: render levels as fractions of the protein mean
            levels[t] = {f"{f}xprot": f * float(v[near].mean()) for f in (0.5, 0.25, 0.1, 0.05)}
        print(f"{t:>20}: belt-to-protein ratio {ratio:.3f}", flush=True)
    del prim
    for c in CUTOFFS:
        t = tag(c)
        v = lowpass(avg, c)
        with mrcfile.new(f"{VIZ}/maps/avg_{t}.mrc", overwrite=True) as m:
            m.set_data(v); m.voxel_size = APIX
        base, sd = float(np.median(v[solvent])), float(v[solvent].std())
        vs = v[shell]
        exc = np.array([(vs[hb == h].mean() - base) / sd if (hb == h).any() else np.nan for h in hbins])
        for h, e in zip(hbins, exc):
            prof_rows.append((t, h * SLICE_A + SLICE_A / 2, e))
        # band = FWHM of the profile around its peak
        ip = int(np.nanargmax(exc)); half = exc[ip] / 2
        lo, hi = ip, ip
        while lo > 0 and exc[lo - 1] >= half: lo -= 1
        while hi < len(exc) - 1 and exc[hi + 1] >= half: hi += 1
        hlo, hhi = hbins[lo] * SLICE_A, (hbins[hi] + 1) * SLICE_A
        inband = (sh_h >= hlo) & (sh_h < hhi)
        smeans = np.array([vs[inband & (sh_sect == s)].mean() for s in range(NSECT)])
        for k in KS:
            belt_rows.append((t, k, float(np.mean(smeans > base + k * sd))))
        frac_mask = float(in_mask_shell[inband].mean())
        ratio = float(vs[inband].mean() / v[near].mean())
        summary[t] = dict(solvent_median=base, solvent_sd=sd, peak_excess_sd=float(exc[ip]),
                          peak_h_A=float(hbins[ip] * SLICE_A + SLICE_A / 2), band_A=[float(hlo), float(hhi)],
                          band_shell_frac_in_deposit_mask=frac_mask,
                          min_sector_excess_sd=float((smeans.min() - base) / sd),
                          belt_to_protein_ratio=ratio)
        levels[t] = {str(k): base + k * sd for k in KS_RENDER}
        print(f"{t:>10}: solvent sd {sd:.5f}  peak excess {exc[ip]:5.1f} sd at h={hbins[ip]*SLICE_A+1:+.0f} A  "
              f"band {hlo:+.0f}..{hhi:+.0f} A  weakest sector {(smeans.min()-base)/sd:5.1f} sd  "
              f"in deposited mask {frac_mask:.2f}  belt/protein {ratio:.3f}", flush=True)

    json.dump(dict(levels=levels, summary=summary, centre=cen.tolist(), up=up.tolist(), perp=perp.tolist(),
                   markers_h=marks, shell_A=SHELL_A, solvent_A=SOLVENT_A), open(f"{VIZ}/levels.json", "w"), indent=1)
    with open(f"{VIZ}/micelle_profile.tsv", "w") as f:
        f.write("map\th_A\tshell_excess_sd\n")
        f.writelines("%s\t%.1f\t%.3f\n" % r for r in prof_rows)
    with open(f"{VIZ}/micelle_belt.tsv", "w") as f:
        f.write("map\tk_sd\tsector_coverage\n")
        f.writelines("%s\t%d\t%.4f\n" % r for r in belt_rows)
    plot(prof_rows, belt_rows, marks, summary)


def plot(prof_rows, belt_rows, marks, summary):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    # ordered cutoffs -> one blue ramp, light (unfiltered) to dark (20 A)
    ramp = ["#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#104281"]
    ink, muted, grid = "#0b0b0b", "#52514e", "#e4e3df"
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": muted, "axes.labelcolor": ink,
                         "xtick.color": muted, "ytick.color": muted})
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10, 4.6), gridspec_kw={"width_ratios": [1.1, 1]})
    for col, c in zip(ramp, CUTOFFS):
        t = tag(c)
        h = [r[1] for r in prof_rows if r[0] == t]; e = [r[2] for r in prof_rows if r[0] == t]
        lab = "unfiltered" if c is None else f"low-pass {c} Å"
        a1.plot(e, h, color=col, lw=2, label=lab)
        k = [r[1] for r in belt_rows if r[0] == t]; cv = [r[2] for r in belt_rows if r[0] == t]
        a2.plot(k, cv, color=col, lw=2, marker="o", ms=4, label=lab)
    for name, hm in marks.items():
        a1.axhline(hm, color=muted, lw=0.8, ls=":")
        a1.text(a1.get_xlim()[1] * 0.98 if False else 1, hm + 1.2, "bound lipid" if name == "PCW" else "revaprazan",
                color=muted, fontsize=8)
    a1.axvline(0, color=muted, lw=0.8)
    a1.set_xlabel("shell density above bulk solvent (solvent sd)")
    a1.set_ylabel("height along membrane normal (Å; + = cytoplasm)")
    a1.set_title("Density 10–25 Å outside the model, by height", color=ink, fontsize=10, loc="left")
    a1.grid(color=grid, lw=0.6); a1.set_axisbelow(True)
    a2.set_xscale("log")
    a2.set_xticks(KS, [str(k) for k in KS])
    a2.set_ylim(-0.03, 1.03)
    a2.set_xlabel("threshold (solvent median + k × solvent sd), k")
    a2.set_ylabel("fraction of 36 sectors above threshold")
    a2.set_title("Is the belt closed at that threshold?", color=ink, fontsize=10, loc="left")
    a2.grid(color=grid, lw=0.6); a2.set_axisbelow(True)
    a2.legend(frameon=False, fontsize=8, loc="lower left")
    for a in (a1, a2):
        for s in ("top", "right"):
            a.spines[s].set_visible(False)
    fig.suptitle("EMD-32299 half-map average: the detergent micelle around the membrane domain",
                 color=ink, fontsize=11, x=0.01, ha="left")
    fig.tight_layout()
    os.makedirs(f"{VIZ}/renders", exist_ok=True)
    fig.savefig(f"{VIZ}/renders/micelle_profiles.png", dpi=150, facecolor="#fcfcfb")


if __name__ == "__main__":
    main()
