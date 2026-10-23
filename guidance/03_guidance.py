#!/usr/bin/env python3
# Guidance index G per field: are steps aligned with the fiber under the cell, beyond a shuffled
# control? G = <cos 2 delta>real - <cos 2 delta>shuffled, -1 to 1 (> 0 along the fiber, < 0 across).
# In: results/guidance/fiber_blocks (01), results/guidance/traj (02)
# Out: results/guidance/guidance_per_field.csv, results/guidance/fig_guidance.png

import os
import glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# USER INPUTS
ROOT = "."
TRACK_DIR = f"{ROOT}/results/guidance/traj"      # written by 02_spots_to_trajectories.py
OUT = f"{ROOT}/results/guidance"

UM_PER_PX = 0.64      # micrometers per pixel
MIN_FRAMES = 40       # min frames detected (of 48); short, lost tracks carry the most linking errors
AGG = 2               # 2 x 2 stored 8 px blocks -> 10.24 um fiber patch, about one cell
N_SAMPLE = 5          # points sampled along each step, so a long step sees all the matrix it crossed
N_PERM = 200          # shuffles used to build each field's null

# P170713_field_13 shares 27 of 48 raw frames with field_12 (upstream export bug): not independent
SKIP_FIELDS = ["P170713_field_13"]

rng = np.random.default_rng(0)

# plot colors
C_PRO, C_REST = "#2a78d6", "#eb6834"
INK, INK2, GRID, MUTED = "#0b0b0b", "#52514e", "#e3e2df", "#c3c2b7"

# STEP 1. Per field: build steps, score them against the fiber, compare with shuffles
results = []
all_delta, all_delta_null = [], []

for npz_path in sorted(glob.glob(f"{OUT}/fiber_blocks/*/*_block8.npz")):
    category = os.path.basename(os.path.dirname(npz_path))
    field = os.path.basename(npz_path).replace("_block8.npz", "")
    if field in SKIP_FIELDS:
        continue

    # 1a. fiber field: 01 stored additive block sums C = sum w cos 2theta, S = sum w sin 2theta,
    # so summing AGG x AGG blocks gives the patch exactly
    stored = np.load(npz_path)
    C, S = stored["C"], stored["S"]
    patch_px = int(stored["block"]) * AGG
    n_frames, n_rows, n_cols = C.shape
    rows, cols = n_rows // AGG, n_cols // AGG
    crop = (slice(None), slice(0, rows * AGG), slice(0, cols * AGG))
    C = C[crop].reshape(n_frames, rows, AGG, cols, AGG).sum(axis=(2, 4))
    S = S[crop].reshape(n_frames, rows, AGG, cols, AGG).sum(axis=(2, 4))

    # 1b. trajectories, tracks with at least MIN_FRAMES detections
    traj = pd.read_csv(f"{TRACK_DIR}/{category}/trajectories_{field}.csv")
    traj = traj.sort_values(["particle", "frame"])
    n_detections = traj.groupby("particle").frame.transform("size")
    traj = traj[n_detections >= MIN_FRAMES]

    # consecutive detections of a cell become steps (first row of each track is NaN, dropped)
    cell = traj.groupby("particle", sort=False)
    steps = pd.DataFrame({
        "particle": traj.particle,
        "frame_from": traj.frame - cell.frame.diff(),
        "x_from": traj.x - cell.x.diff(),
        "y_from": traj.y - cell.y.diff(),
        "dx": cell.x.diff(),
        "dy": cell.y.diff(),
        "dframe": cell.frame.diff(),
    }).dropna()
    steps = steps[steps.dframe == 1]          # only moves between adjacent frames
    if len(steps) < 100:
        print(f"  {field}: only {len(steps)} steps, skipped")
        continue

    # 1c. step direction; trajectory y is a row index (down), fiber angles have y up
    theta_move = np.arctan2(-steps.dy.values, steps.dx.values)

    # 1d. fiber orientation under the step, averaged over N_SAMPLE points along it
    along = np.linspace(0.1, 0.9, N_SAMPLE)[None, :]
    sample_x = steps.x_from.values[:, None] + steps.dx.values[:, None] * along
    sample_y = steps.y_from.values[:, None] + steps.dy.values[:, None] * along
    frame = steps.frame_from.values.astype(int)

    bx = np.clip((sample_x // patch_px).astype(int), 0, cols - 1)
    by = np.clip((sample_y // patch_px).astype(int), 0, rows - 1)
    bf = np.repeat(np.clip(frame, 0, n_frames - 1)[:, None], N_SAMPLE, axis=1)
    theta_fiber = 0.5 * np.arctan2(S[bf, by, bx].sum(axis=1), C[bf, by, bx].sum(axis=1))

    A_real = np.cos(2 * (theta_move - theta_fiber)).mean()

    # 1e. shuffled control: each step keeps its direction but is scored against another cell's
    # patch from the same frame, so a field-wide fiber axis plus field-wide drift cancels out of G
    rows_of_frame = [np.where(frame == f)[0] for f in np.unique(frame)]
    null = np.empty(N_PERM)
    for p in range(N_PERM):
        shuffled = np.arange(len(steps))
        for idx in rows_of_frame:
            shuffled[idx] = rng.permutation(idx)
        bx_p = np.clip((sample_x[shuffled] // patch_px).astype(int), 0, cols - 1)
        by_p = np.clip((sample_y[shuffled] // patch_px).astype(int), 0, rows - 1)
        bf_p = np.repeat(np.clip(frame[shuffled], 0, n_frames - 1)[:, None], N_SAMPLE, axis=1)
        theta_fiber_p = 0.5 * np.arctan2(S[bf_p, by_p, bx_p].sum(axis=1),
                                         C[bf_p, by_p, bx_p].sum(axis=1))
        null[p] = np.cos(2 * (theta_move - theta_fiber_p)).mean()

    G = A_real - null.mean()
    # two-sided empirical p (G < 0, moving across fibers, counts too); floor is 2/(N_PERM+1)
    p_along = (np.sum(null >= A_real) + 1) / (N_PERM + 1)
    p_across = (np.sum(null <= A_real) + 1) / (N_PERM + 1)
    p_perm = min(1.0, 2 * min(p_along, p_across))

    # delta in degrees, folded to 0-90 (0 = along, 45 = unrelated, 90 = across)
    delta = np.degrees(np.arccos(np.clip(np.abs(np.cos(theta_move - theta_fiber)), 0, 1)))
    delta_null = np.degrees(np.arccos(np.clip(np.abs(np.cos(theta_move - theta_fiber_p)), 0, 1)))
    all_delta.append(delta)
    all_delta_null.append(delta_null)

    results.append(dict(field=field, category=category, cafType=category.split("_")[0],
                        n_tracks=steps.particle.nunique(), n_steps=len(steps),
                        A_real=A_real, A_null=null.mean(), G=G, p_perm=p_perm,
                        mean_delta_deg=delta.mean(),
                        frac_within_30deg=float((delta < 30).mean())))
    print(f"  {field:22s} {steps.particle.nunique():4d} tracks {len(steps):6d} steps   "
          f"G={G:+.4f}  p={p_perm:.4f}  mean|delta|={delta.mean():.1f} deg", flush=True)

res = pd.DataFrame(results)
res.to_csv(f"{OUT}/guidance_per_field.csv", index=False, float_format="%.5g")

delta = np.concatenate(all_delta)
delta_null = np.concatenate(all_delta_null)
n_along = int(((res.G > 0) & (res.p_perm <= 0.05)).sum())
n_across = int(((res.G < 0) & (res.p_perm <= 0.05)).sum())

print(f"\nSTEP 1  {len(res)} fields, {res.n_tracks.sum():,} tracks, {res.n_steps.sum():,} steps")
print(f"        guided ALONG the fiber:  {n_along} fields")
print(f"        guided ACROSS it:        {n_across} fields")
print(f"        not significant:         {int((res.p_perm > 0.05).sum())} fields")
print(f"        mean |delta| = {delta.mean():.1f} deg   (45.0 = unrelated)")
print(f"        {100*(delta < 30).mean():.1f}% of steps fall within 30 deg of the fiber "
      f"(33.3% expected by chance)")
a = res[res.cafType == "pro"].G
b = res[res.cafType == "rest"].G
print(f"        G: proCAF {a.mean():+.4f} (n={len(a)})   restCAF {b.mean():+.4f} (n={len(b)})")

# STEP 2. Figure
plt.rcParams.update({
    "figure.facecolor": "#fcfcfb", "axes.facecolor": "#fcfcfb",
    "axes.edgecolor": GRID, "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": INK2, "ytick.color": INK2, "font.size": 10,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": .6, "grid.alpha": .8,
})
fig, ax = plt.subplots(1, 2, figsize=(12.5, 5.2))

# panel a: delta distribution, real vs shuffled pairing
bins = np.linspace(0, 90, 19)
centers = bins[:-1] + 2.5
h_real, _ = np.histogram(delta, bins=bins, density=True)
h_null, _ = np.histogram(delta_null, bins=bins, density=True)
ax[0].fill_between(centers, h_null * 90, h_real * 90, where=(h_real >= h_null),
                   color="#cde2fb", alpha=.85, zorder=1)
ax[0].plot(centers, h_real * 90, lw=2.6, color="#1c5cab", zorder=3,
           label="real: each step vs the fiber under it")
ax[0].plot(centers, h_null * 90, lw=2.6, color=MUTED, zorder=2,
           label="shuffled: same steps, another cell's fiber")
ax[0].axhline(1, color=INK, ls="--", lw=1.2, zorder=4)
ax[0].annotate("flat at 1.0 = no relationship", (34, 1.02), fontsize=9, color=INK)
ax[0].set_xticks([0, 15, 30, 45, 60, 75, 90])
ax[0].set_xlim(0, 90)
ax[0].set_xlabel("δ, angle between the step and the fiber under it (°)\n"
                 "0 = along the fiber        45 = unrelated        90 = across it")
ax[0].set_ylabel("relative frequency  (1.0 = flat)")
ax[0].set_title("a   What the guidance index measures", loc="left", fontweight="bold")
ax[0].legend(frameon=False, fontsize=9, loc="upper right")

# panel b: G for every field, ranked
ranked = res.sort_values("G").reset_index(drop=True)
x = np.arange(len(ranked))
for caf, color, label in [("pro", C_PRO, "proCAF matrix"), ("rest", C_REST, "restCAF matrix")]:
    m = (ranked.cafType == caf).values
    ax[1].vlines(x[m], 0, ranked.G[m], color=color, lw=2.2, alpha=.85, zorder=2)
    ax[1].scatter(x[m], ranked.G[m], s=34, color=color, zorder=3, edgecolors="none", label=label)
ns = (ranked.p_perm > 0.05).values
ax[1].scatter(x[ns], ranked.G[ns], s=34, facecolors="#fcfcfb", zorder=4,
              edgecolors=INK2, linewidths=1.2, label="not significant")
ax[1].axhline(0, color=INK, lw=1.2, zorder=1)
ax[1].set_xlim(-1, len(ranked))
ax[1].set_xlabel(f"the {len(ranked)} fields, ranked by G")
ax[1].set_ylabel("guidance index  G = ⟨cos 2δ⟩real − ⟨cos 2δ⟩shuffled")
ax[1].set_title(f"b   {n_along} of {len(ranked)} fields guided along the fiber; "
                f"{n_across} across", loc="left", fontweight="bold")
ax[1].legend(frameon=False, fontsize=9, loc="upper left")

fig.suptitle("Cancer cells migrate along the matrix fiber beneath them",
             fontsize=13, fontweight="bold", x=.012, ha="left", y=.985)
fig.text(.012, .925, f"{len(res)} fields · tracks detected in at least {MIN_FRAMES} of 48 "
                     f"frames ({res.n_steps.sum():,} steps) · 10.2 µm fiber patch",
         fontsize=9.5, color=INK2, ha="left")
fig.tight_layout(rect=[0, 0, 1, .90])
fig.savefig(f"{OUT}/fig_guidance.png", dpi=170)
print(f"\nSTEP 2  wrote fig_guidance.png and guidance_per_field.csv")
