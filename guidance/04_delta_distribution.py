#!/usr/bin/env python3
# Per-field and pooled distributions of delta (step vs fiber angle, 0-90 deg), real vs each field's
# own within-frame shuffle; same scoring as 03_guidance.py, which keeps only summaries.
# Field- and donor-level Wilcoxon tests; donor IDs come from alignment_scales.csv.
# Out: results/guidance/delta_per_field.csv, delta_pooled.npz, guidance_score_pooled.csv

import os
import glob
import numpy as np
import pandas as pd
from scipy import stats

# USER INPUTS
ROOT = "."
OUT = f"{ROOT}/results/guidance"
TRACK_DIR = f"{OUT}/traj"

MIN_FRAMES, AGG, N_SAMPLE, N_PERM = 40, 2, 5, 200
N_PERM_POOLED = 10          # permutations whose per-step deltas go into the pooled null
SKIP_FIELDS = ["P170713_field_13"]

rng = np.random.default_rng(0)

# STEP 1. Recompute delta per field, real and shuffled
rows, pooled_real, pooled_null = [], [], []
cos_real, cos_null = [], []

for npz_path in sorted(glob.glob(f"{OUT}/fiber_blocks/*/*_block8.npz")):
    category = os.path.basename(os.path.dirname(npz_path))
    field = os.path.basename(npz_path).replace("_block8.npz", "")
    if field in SKIP_FIELDS:
        continue
    traj_path = f"{TRACK_DIR}/{category}/trajectories_{field}.csv"
    if not os.path.exists(traj_path):
        continue

    st = np.load(npz_path)
    C, S = st["C"], st["S"]
    patch_px = int(st["block"]) * AGG
    n_frames, n_rows, n_cols = C.shape
    r_, c_ = n_rows // AGG, n_cols // AGG
    crop = (slice(None), slice(0, r_ * AGG), slice(0, c_ * AGG))
    C = C[crop].reshape(n_frames, r_, AGG, c_, AGG).sum(axis=(2, 4))
    S = S[crop].reshape(n_frames, r_, AGG, c_, AGG).sum(axis=(2, 4))

    traj = pd.read_csv(traj_path).sort_values(["particle", "frame"])
    traj = traj[traj.groupby("particle").frame.transform("size") >= MIN_FRAMES]
    cell = traj.groupby("particle", sort=False)
    steps = pd.DataFrame({
        "frame_from": traj.frame - cell.frame.diff(),
        "x_from": traj.x - cell.x.diff(), "y_from": traj.y - cell.y.diff(),
        "dx": cell.x.diff(), "dy": cell.y.diff(), "dframe": cell.frame.diff()}).dropna()
    steps = steps[steps.dframe == 1]
    if len(steps) < 100:
        continue

    theta_move = np.arctan2(-steps.dy.values, steps.dx.values)   # y is a row index (down)
    along = np.linspace(0.1, 0.9, N_SAMPLE)[None, :]
    sx = steps.x_from.values[:, None] + steps.dx.values[:, None] * along
    sy = steps.y_from.values[:, None] + steps.dy.values[:, None] * along
    frame = steps.frame_from.values.astype(int)

    bx = np.clip((sx // patch_px).astype(int), 0, c_ - 1)
    by = np.clip((sy // patch_px).astype(int), 0, r_ - 1)
    bf = np.repeat(np.clip(frame, 0, n_frames - 1)[:, None], N_SAMPLE, axis=1)
    theta_fiber = 0.5 * np.arctan2(S[bf, by, bx].sum(axis=1), C[bf, by, bx].sum(axis=1))
    delta = np.degrees(np.arccos(np.clip(np.abs(np.cos(theta_move - theta_fiber)), 0, 1)))

    # null: each step scored against another cell's patch from the same frame. 45 deg alone is
    # not the yardstick: a field-wide fiber axis plus field-wide drift lowers delta on its own.
    idx_of_frame = [np.where(frame == f)[0] for f in np.unique(frame)]
    null_means = np.empty(N_PERM)
    for p in range(N_PERM):
        sh = np.arange(len(steps))
        for idx in idx_of_frame:
            sh[idx] = rng.permutation(idx)
        bxp = np.clip((sx[sh] // patch_px).astype(int), 0, c_ - 1)
        byp = np.clip((sy[sh] // patch_px).astype(int), 0, r_ - 1)
        bfp = np.repeat(np.clip(frame[sh], 0, n_frames - 1)[:, None], N_SAMPLE, axis=1)
        tfp = 0.5 * np.arctan2(S[bfp, byp, bxp].sum(axis=1), C[bfp, byp, bxp].sum(axis=1))
        dnull = np.degrees(np.arccos(np.clip(np.abs(np.cos(theta_move - tfp)), 0, 1)))
        null_means[p] = dnull.mean()
        if p < N_PERM_POOLED:
            pooled_null.append(dnull)
            cos_null.append(np.cos(2 * np.radians(dnull)))
    pooled_real.append(delta)
    cos_real.append(np.cos(2 * np.radians(delta)))

    # one-sided empirical p: how often does a shuffle beat the real pairing?
    p_perm = (np.sum(null_means <= delta.mean()) + 1) / (N_PERM + 1)
    rows.append(dict(field=field, category=category, cafType=category.split("_")[0],
                     n_steps=len(delta), delta_real=delta.mean(),
                     delta_null=null_means.mean(),
                     null_lo=np.percentile(null_means, 2.5),
                     null_hi=np.percentile(null_means, 97.5), p_perm=p_perm))
    print(f"  {field:22s} {len(delta):6d} steps  delta {delta.mean():5.1f} vs null "
          f"{null_means.mean():5.1f}  p={p_perm:.4f}", flush=True)

res = pd.DataFrame(rows)
res.to_csv(f"{OUT}/delta_per_field.csv", index=False, float_format="%.5g")
real = np.concatenate(pooled_real)
null = np.concatenate(pooled_null)
_b = np.arange(0, 91, 3)
cr, cn = np.concatenate(cos_real), np.concatenate(cos_null)
_cb = np.linspace(-1, 1, 41)
np.savez(f"{OUT}/delta_pooled.npz", bins=_b,
         h_real=np.histogram(real, bins=_b, density=True)[0],
         h_null=np.histogram(null, bins=_b, density=True)[0],
         mean_real=real.mean(), mean_null=null.mean(), n_steps=len(real),
         cos_bins=_cb,
         hc_real=np.histogram(cr, bins=_cb, density=True)[0],
         hc_null=np.histogram(cn, bins=_cb, density=True)[0],
         cos_mean_real=cr.mean(), cos_mean_null=cn.mean())
# per-step guidance score (cos 2 delta) histogram as a plain table for R
pd.DataFrame({"bin_lo": _cb[:-1], "bin_hi": _cb[1:], "mid": (_cb[:-1] + _cb[1:]) / 2,
              "real": np.histogram(cr, bins=_cb, density=True)[0] * 2,
              "shuffled": np.histogram(cn, bins=_cb, density=True)[0] * 2}
             ).to_csv(f"{OUT}/guidance_score_pooled.csv", index=False)
print(f"        <cos 2delta>  real {cr.mean():+.4f}   shuffled {cn.mean():+.4f}   "
      f"G_pooled = {cr.mean() - cn.mean():+.4f}")

# STEP 2. Report
w = stats.wilcoxon(res.delta_real, res.delta_null)
t45 = stats.wilcoxon(res.delta_real - 45)
print(f"\n{len(res)} fields, {res.n_steps.sum():,} steps")
print(f"  pooled mean delta   real {real.mean():.2f} deg   shuffled {null.mean():.2f} deg   (45 = chance)")
print(f"  fields below their own null: {int((res.delta_real < res.delta_null).sum())}/{len(res)}"
      f"   permutation p <= 0.05: {int((res.p_perm <= 0.05).sum())}/{len(res)}")
print(f"  paired across fields, real vs null : Wilcoxon p = {w.pvalue:.2g}")
print(f"  fields vs the 45 deg chance line   : Wilcoxon p = {t45.pvalue:.2g}")

# collapse to donors: 54 fields come from only 14 cfECM donors, and the matrix is the donor's
meta = pd.read_csv(f"{OUT}/alignment_scales.csv")[["field", "donor"]]
dn = res.merge(meta, on="field").groupby("donor")[["delta_real", "delta_null"]].mean()
wd = stats.wilcoxon(dn.delta_real, dn.delta_null)
print(f"  collapsed to {len(dn)} donors, real vs null : Wilcoxon p = {wd.pvalue:.4f}"
      f"   ({int((dn.delta_real < dn.delta_null).sum())}/{len(dn)} donors below their null)")

print(f"\nwrote delta_per_field.csv, delta_pooled.csv/.npz, guidance_score_pooled.csv\n"
      f"     figures are drawn by 05_fig_delta.R")
