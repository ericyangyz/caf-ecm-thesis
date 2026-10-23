# Per-sample comparison of the old Keyence 20x alignment from Mathematica (the variances in
# mathematica_variance_test.Rmd) and from UMA/OrientationPy (sd_around_0_deg).
# in:  results/old_keyence_alignment_per_sample.csv
# out: results/old_keyence_uma_vs_mathematica.csv (+ _correlations.csv), results/fig12_uma_vs_mathematica.pdf/.png

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

# USER INPUTS
ROOT      = "."
UMA_FILE  = f"{ROOT}/results/old_keyence_alignment_per_sample.csv"
OUT_CSV   = f"{ROOT}/results/old_keyence_uma_vs_mathematica.csv"
OUT_FIG   = f"{ROOT}/results/fig12_uma_vs_mathematica"

PRO_COLOR, REST_COLOR = "#E74095", "#288D92"   # proCAF, restCAF
CLUSTER_CUT = 20.0     # deg; the 15 samples fall into a group below this and a group above

D2R  = np.pi / 180
UNIF = np.pi**2 / 12          # rad^2, variance of a uniform axial distribution (= 51.96 deg sd)
assert abs(np.sqrt(UNIF) / D2R - 51.96) < 0.01, "the random-matrix ceiling must agree in both units"


# STEP 1. Mathematica variances (rad^2), in the notebook's proSamples / restSamples order
# names use the UMA spelling (notebook: "229", "1207", "0118", "201")
pro_order  = ["015", "229T", "804N", "804T", "815", "908", "1207T"]
rest_order = ["0118N", "201N", "212", "509", "519", "711", "807", "1108"]
pro_var    = [0.657998, 0.618538, 0.51705, 0.772362, 0.407781, 0.583037, 0.589696]
rest_var   = [0.66782, 0.465514, 0.457228, 0.42473, 0.391601, 0.466767, 0.452664, 0.430305]

assert len(pro_order) == len(pro_var) == 7 and len(rest_order) == len(rest_var) == 8
math_map = dict(zip(pro_order + rest_order, pro_var + rest_var))
assert max(math_map.values()) <= UNIF, "a variance above pi^2/12 means these are not rad^2"


# STEP 2. Join onto the UMA table
d = pd.read_csv(UMA_FILE)
d["math_var_rad2"] = d["sample"].map(math_map)

assert len(d) == 15 and d.math_var_rad2.notna().all(), "sample names did not all match"

# Mathematica = unweighted variance (rad^2), UMA = |grad I|^2-weighted sd (deg); both are
# centered on the modal direction, so sqrt + rad->deg puts them on one axis
d["math_sd_deg"]        = np.sqrt(d.math_var_rad2) / D2R
d["math_pct_of_random"] = 100 * d.math_var_rad2 / UNIF
d["uma_pct_of_random"]  = 100 * d.sd_around_0_deg**2 * D2R**2 / UNIF


# STEP 3. Correlations, overall and within groups
# by UMA: nine at 10.2-16.5 deg, four at 22.0-23.5, two above (908 at 36.6, 804T at 51.0);
# the whole-set r mostly reflects that grouping, so within-group r is reported too
tight = d[d.sd_around_0_deg < CLUSTER_CUT]
mid   = d[(d.sd_around_0_deg >= CLUSTER_CUT) & (d.sd_around_0_deg < 24)]
subsets = {
    "all 15":                         d,
    "all 15 except 908":              d[d["sample"] != "908"],
    f"the {len(tight)} aligned ones (<{CLUSTER_CUT:.0f} deg)":  tight,
    f"the {len(mid)} at 22-23.5 deg":                           mid,
}
print(f"{'subset':>32s}  {'n':>2s}  {'Pearson r':>18s}  {'Spearman rho':>18s}   UMA span")
rows = []
for lab, sub in subsets.items():
    if len(sub) < 3:
        continue
    rp = stats.pearsonr(sub.sd_around_0_deg, sub.math_sd_deg)
    rs = stats.spearmanr(sub.sd_around_0_deg, sub.math_sd_deg)
    span = f"{sub.sd_around_0_deg.min():.1f}-{sub.sd_around_0_deg.max():.1f} deg"
    print(f"{lab:>32s}  {len(sub):2d}  {rp[0]:+.3f} (p={rp[1]:.3f})  {rs.statistic:+.3f} (p={rs.pvalue:.3f})   {span}")
    rows.append(dict(subset=lab, n=len(sub), pearson_r=rp[0], pearson_p=rp[1],
                     spearman_rho=rs.statistic, spearman_p=rs.pvalue, uma_span_deg=span))


# STEP 4. Write the joined table
d.sort_values(["group", "sd_around_0_deg"])[
    ["group", "sample", "donor", "sd_around_0_deg", "uma_pct_of_random",
     "math_var_rad2", "math_sd_deg", "math_pct_of_random"]
].to_csv(OUT_CSV, index=False)
pd.DataFrame(rows).to_csv(OUT_CSV.replace(".csv", "_correlations.csv"), index=False)
print(f"\nwritten: {OUT_CSV}")


# STEP 5. Figure (points labeled by sample so color is not the only cue)
fig, ax = plt.subplots(figsize=(7.4, 6.4))

ax.axvline(51.96, color="#999999", ls="--", lw=0.9, zorder=1)
ax.axhline(51.96, color="#999999", ls="--", lw=0.9, zorder=1)
ax.text(51.2, 40.0, "random matrix, 51.96°", rotation=90, ha="right", va="bottom",
        fontsize=8, color="#666666")
ax.text(1.5, 52.3, "random matrix, 51.96°", ha="left", va="bottom", fontsize=8, color="#666666")

# shade the two clusters both methods separate
ax.axvspan(8.5, 18.0, color="#000000", alpha=0.045, zorder=0)
ax.axvspan(21.0, 24.5, color="#000000", alpha=0.045, zorder=0)
ax.text(13.2, 31.0, f"n = {len(tight)}", ha="center", fontsize=8.5, color="#777777")
ax.text(22.75, 31.0, f"n = {len(mid)}", ha="center", fontsize=8.5, color="#777777")

for g, c, m in (("proCAF", PRO_COLOR, "o"), ("restCAF", REST_COLOR, "s")):
    sg = d[d.group == g]
    ax.scatter(sg.sd_around_0_deg, sg.math_sd_deg, s=66, c=c, marker=m,
               edgecolors="white", linewidths=1.2, zorder=3, label=f"{g} (n = {len(sg)})")
    for _, r in sg.iterrows():
        dx, dy = (7, -3.5) if r["sample"] not in ("711", "509", "229T") else (7, 3)
        ax.annotate(r["sample"], (r.sd_around_0_deg, r.math_sd_deg),
                    textcoords="offset points", xytext=(dx, dy),
                    fontsize=8.5, color=c, zorder=4)

# 908 is the one sample the two methods rank in opposite directions
r908 = d[d["sample"] == "908"].iloc[0]
ax.scatter([r908.sd_around_0_deg], [r908.math_sd_deg], s=210, facecolors="none",
           edgecolors="#B0711A", linewidths=1.6, zorder=2)
ax.annotate("908: UMA ranks it 2nd-worst,\nMathematica ranks it mid-pack",
            (r908.sd_around_0_deg, r908.math_sd_deg), xytext=(-6, -34),
            textcoords="offset points", fontsize=8.5, color="#8A5714", ha="right",
            arrowprops=dict(arrowstyle="-", color="#B0711A", lw=0.9))

r_all   = stats.pearsonr(d.sd_around_0_deg, d.math_sd_deg)
r_tight = stats.pearsonr(tight.sd_around_0_deg, tight.math_sd_deg)

ax.set_xlabel("UMA / OrientationPy   sd_around_0_deg   (°, |∇I|²-weighted)")
ax.set_ylabel("Mathematica   √Variance[oCentered]   (°, unweighted)")
ax.set_title("The two analyses agree on the grouping, not within a group",
             fontsize=12, pad=12)
ax.text(0.035, 0.965,
        f"all 15:  r = {r_all[0]:+.2f}  (p < 0.001)\n"
        f"within the {len(tight)} aligned ones:  r = {r_tight[0]:+.2f}  (p = {r_tight[1]:.2f})",
        transform=ax.transAxes, va="top", ha="left", fontsize=9.5,
        bbox=dict(boxstyle="round,pad=0.45", fc="white", ec="#CCCCCC", lw=0.8))
ax.set_xlim(0, 56); ax.set_ylim(30, 56)
ax.legend(loc="lower right", frameon=False, fontsize=9.5)
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(f"{OUT_FIG}.{ext}", dpi=200)
print(f"written: {OUT_FIG}.pdf / .png")
