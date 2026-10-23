#!/usr/bin/env python3
# Matrix alignment per field at three scales, from frame 0 of the fiber blocks: global order S_global,
# local order R_local (one 10.2 um patch) and correlation length xi_um; joined to G and turning stats.
# Earlier trackpy-stage script: the thesis uses only the matrix columns (R_local, xi_um, S_global,
# donor, cafType); G comes from 03_guidance.py. Out: results/guidance/alignment_scales.csv, fig4_chain.png

import glob
import os
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

# USER INPUTS
ROOT = "."
OUT = f"{ROOT}/results/guidance"
FRAME = 0                 # the matrix is essentially static over 24 h, so one frame is enough
MAX_LAG_BLOCKS = 120      # 120 blocks x 5.12 um = 614 um
LO_UM, HI_UM = 6.4, 20.9  # credible step band

# STEP 1. Per field: global order S, mean local order R, correlation length xi
rows = []
for npz_path in sorted(glob.glob(f"{OUT}/fiber_blocks/*/*_block8.npz")):
    category = os.path.basename(os.path.dirname(npz_path))
    field = os.path.basename(npz_path).replace("_block8.npz", "")
    z = np.load(npz_path)
    C = z["C"][FRAME].astype(np.float64)
    S = z["S"][FRAME].astype(np.float64)
    W = z["W"][FRAME].astype(np.float64)
    block_um = float(z["block"]) * float(z["um_per_px"])     # 8 px = 5.12 um

    # global order: sum every block's axial vector over the field, 0..1
    S_global = np.hypot(C.sum(), S.sum()) / W.sum()

    # local order inside each stored block, 0..1
    R_block = np.hypot(C, S) / np.maximum(W, 1e-12)
    good = W > np.percentile(W, 10)          # ignore blocks with essentially no signal
    R_local = float(np.average(R_block[good], weights=W[good]))

    # correlation of the axial direction with distance, C(r) = <cos 2(theta(x) - theta(x+r))>;
    # low-signal blocks get zero weight
    norm = np.maximum(np.hypot(C, S), 1e-12)
    ux, uy = C / norm, S / norm
    w = np.where(good, W, 0.0)
    ux, uy = ux * (w > 0), uy * (w > 0)

    lags = np.unique(np.round(np.logspace(0, np.log10(MAX_LAG_BLOCKS), 22)).astype(int))
    corr = []
    for lag in lags:
        # horizontal and vertical offsets pooled
        num = (w[:, :-lag] * w[:, lag:] *
               (ux[:, :-lag] * ux[:, lag:] + uy[:, :-lag] * uy[:, lag:])).sum()
        den = (w[:, :-lag] * w[:, lag:]).sum()
        num += (w[:-lag] * w[lag:] *
                (ux[:-lag] * ux[lag:] + uy[:-lag] * uy[lag:])).sum()
        den += (w[:-lag] * w[lag:]).sum()
        corr.append(num / den)
    corr = np.array(corr)
    r_um = lags * block_um

    # fit C(r) = (1 - c) exp(-r/xi) + c; plateau c is free because any global anisotropy
    # makes C(r) level off near S^2 rather than 0
    try:
        (xi, c), _ = curve_fit(lambda r, xi, c: (1 - c) * np.exp(-r / xi) + c,
                               r_um, corr, p0=[30.0, max(corr[-1], 0.0)],
                               bounds=([1.0, -0.5], [2000.0, 0.99]), maxfev=20000)
    except Exception:
        xi, c = np.nan, np.nan

    rows.append(dict(field=field, category=category,
                     cafType=category.split("_")[0], tumorType=category.split("_")[1],
                     S_global=S_global, R_local=R_local, xi_um=xi, plateau=c,
                     corr_at_50um=float(np.interp(50, r_um, corr)),
                     corr_at_200um=float(np.interp(200, r_um, corr))))
    print(f"  {field:22s} S_global={S_global:.3f}  R_local={R_local:.3f}  "
          f"xi={xi:6.1f} um  plateau={c:+.3f}", flush=True)

align = pd.DataFrame(rows)

# STEP 2. Join G (guidance_per_field.csv) and turning stats (steps_all.csv)
guid = pd.read_csv(f"{OUT}/guidance_per_field.csv")[
    ["field", "G", "A_real", "mean_delta_deg", "n_steps"]]

steps = pd.read_csv(f"{OUT}/steps_all.csv",
                    usecols=["field", "donor", "step_um", "cos_turn", "turn_deg"])
credible = steps[(steps.step_um >= LO_UM) & (steps.step_um < HI_UM)]
# <cos turn> over steps in the credible band: +1 straight, 0 no memory, -1 reversals
migr = credible.groupby("field").agg(cos_turn=("cos_turn", "mean"),
                                     mean_turn_deg=("turn_deg", "mean")).reset_index()
donor = steps[["field", "donor"]].drop_duplicates()

out = align.merge(guid, on="field").merge(migr, on="field").merge(donor, on="field")
out.to_csv(f"{OUT}/alignment_scales.csv", index=False, float_format="%.5g")
print(f"\nwrote {OUT}/alignment_scales.csv")
print(out[["field", "donor", "cafType", "S_global", "R_local", "xi_um",
           "G", "cos_turn"]].round(3).to_string(index=False))


# STEP 3. Figure: matrix alignment -> guidance -> directional migration, one panel per link,
# using the measured matrix properties rather than the CAF label
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

C_PRO, C_REST = "#2a78d6", "#eb6834"
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e3e2df"
plt.rcParams.update({
    "figure.facecolor": "#fcfcfb", "axes.facecolor": "#fcfcfb",
    "axes.edgecolor": GRID, "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": INK2, "ytick.color": INK2, "font.size": 10,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": .6, "grid.alpha": .8,
})

out["log_xi"] = np.log10(out.xi_um)
fig, ax = plt.subplots(1, 3, figsize=(14.5, 4.8))

panels = [
    (0, "R_local", "G",
     "local fiber alignment R of the patch under the cell\n"
     "0 = no dominant direction        1 = perfectly parallel",
     "guidance index G", "a   matrix → guidance", False),
    (1, "log_xi", "G",
     "correlation length ξ (µm, log scale)\nhow far the fiber direction persists",
     "guidance index G", "b   how big the domains are → guidance", True),
    (2, "G", "cos_turn",
     "guidance index G", "migration directionality ⟨cos turn⟩\n"
     "+1 perfectly straight    0 no directional memory",
     "c   guidance → directional migration", False),
]
for i, xcol, ycol, xlabel, ylabel, title, logx in panels:
    for caf, color, label in [("pro", C_PRO, "proCAF matrix"), ("rest", C_REST, "restCAF matrix")]:
        d = out[out.cafType == caf]
        ax[i].scatter(10 ** d[xcol] if logx else d[xcol], d[ycol], s=46, color=color,
                      edgecolors="none", alpha=.85, label=label)
    # trend line is an OLS fit; the statistic is Spearman's rho (assumes only monotonicity)
    xv, yv = out[xcol].values, out[ycol].values
    slope, intercept = np.polyfit(xv, yv, 1)
    xs = np.linspace(xv.min(), xv.max(), 50)
    ax[i].plot(10 ** xs if logx else xs, slope * xs + intercept, color=INK2, lw=1.6, ls="--")
    rho, p = stats.spearmanr(xv, yv)
    ax[i].annotate(f"Spearman ρ = {rho:+.2f}\np = {p:.4f}   n = {len(out)}",
                   (0.04, 0.90), xycoords="axes fraction", fontsize=10,
                   va="top", fontweight="bold", color=INK)
    if logx:
        ax[i].set_xscale("log")
    if ycol == "cos_turn":
        ax[i].axhline(0, color=INK, ls=":", lw=1.1)
    ax[i].set_xlabel(xlabel)
    ax[i].set_ylabel(ylabel)
    ax[i].set_title(title, loc="left", fontweight="bold")
ax[0].legend(frameon=False, fontsize=9, loc="lower right")

fig.suptitle("Matrix alignment vs guidance vs directional migration (field-level correlations)",
             fontsize=13, fontweight="bold", x=.010, ha="left", y=.985)
fig.text(.010, .915, "CAUTION: these are raw Spearman correlations over 54 fields. The 54 fields are only "
         "14 cfECM donors x 4 tumor cell lines, and once that structure is accounted for none of "
         "these\nassociations is significant. Treat the panels as descriptive, not as evidence for the chain.",
         fontsize=8.5, color="#eb6834", ha="left")
fig.tight_layout(rect=[0, 0, 1, .89])
fig.savefig(f"{OUT}/fig4_chain.png", dpi=170)
print(f"STEP 3  wrote fig4_chain.png")
