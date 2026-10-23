#!/usr/bin/env python3
# Supplementary Fig S1: the fibronectin IF image of every cfECM donor in the fiber-orientation analysis
# (Fig 3; 7 proCAF, 8 restCAF), on one page at final size, ordered from most to least aligned matrix.
# In:  data/IF_imaging/{proCAF,restCAF}/*.tif (Fiji exports of Zeiss .czi, channel 1), ECM_CSV below
# Out: results/thesis/suppfig1/suppfig1_if_images.pdf (+ 200-dpi PNG preview), suppfig1_images.csv
import os, re, glob
import numpy as np
import pandas as pd
from PIL import Image
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# USER INPUTS
ROOT    = "."
IF_DIR  = f"{ROOT}/data/IF_imaging"
ECM_CSV = "./results/old_keyence_alignment_per_sample.csv"
OUT_DIR = f"{ROOT}/results/thesis/suppfig1"

# file-name prefix (text before the first '#' or '_') -> sample code in ECM_CSV; each must match once
FILE_SAMPLE = {"015sz": "015", "1207T1": "1207T", "229T1": "229T", "804N": "804N", "804T1": "804T",
               "815": "815", "908": "908", "0118N1": "0118N", "1108": "1108", "201N": "201N",
               "212": "212", "509": "509", "519": "519", "711": "711", "807": "807"}
DONOR_FILL = {"509": "P240509N1", "0118N": "P220118N1"}   # never tracked, so no donor string in the csv

PX_PER_UM_EXPECTED = 8.095163          # 0.12353 um/px; checked against every file's TIFF tags
SCALE_UM   = 100                       # um, scale bar length
N_COL      = 4                         # images per row
PAGE_W     = 6.5 * 72                  # pt, the 1-in frame of a letter page
GAP        = 6                         # pt between images
LABEL_H    = 9                         # pt, room for the label above each image
TITLE_H    = 12                        # pt, block title
BLOCK_GAP  = 8                         # pt between the two blocks
IMG_DPI    = 450                       # resolution of the embedded images at their printed size
PT_LABEL   = 6.5;  PT_TITLE = 7.0          # pt at final size; block titles are bold
BAR_H      = 2.0                       # pt, scale bar thickness
BAR_MARGIN = 5.0                       # pt from the right and bottom edges of the image
PRO, REST  = "#E74095", "#288D92"
FONT       = "Liberation Sans"         # metric-identical to Arial

os.makedirs(OUT_DIR, exist_ok=True)
plt.rcParams.update({"font.family": FONT, "pdf.fonttype": 42, "ps.fonttype": 42})   # 42 = text stays editable

# STEP 1. Every image -> sample code, donor, CAF subtype, fiber spread; check calibration
ecm = pd.read_csv(ECM_CSV, dtype={"sample": str})
assert len(ecm) == 15 and ecm["sample"].is_unique
ecm["donor"] = ecm["donor"].fillna(ecm["sample"].map(DONOR_FILL))
assert ecm["donor"].notna().all() and ecm["donor"].is_unique

rows = []
for f in sorted(glob.glob(f"{IF_DIR}/*/*.tif")):
    name = os.path.basename(f)
    prefix = re.match(r"^([0-9A-Za-z]+?)(#|_)", name).group(1)
    assert prefix in FILE_SAMPLE, name
    im = Image.open(f)
    desc = str(im.tag_v2.get(270, ""))
    xres, yres = float(im.tag_v2.get(282)), float(im.tag_v2.get(283))
    assert "unit=micron" in desc and abs(xres - PX_PER_UM_EXPECTED) < 1e-5 and abs(yres - xres) < 1e-9, name
    assert im.size == (2586, 2586) and im.mode == "RGB", name    # 8-bit RGB, green LUT: signal is in G
    rows.append(dict(file=f, folder=os.path.basename(os.path.dirname(f)), sample=FILE_SAMPLE[prefix],
                     px_per_um=xres, width_px=im.size[0]))
img = pd.DataFrame(rows)
assert len(img) == 15 and img["sample"].is_unique and set(img["sample"]) == set(ecm["sample"])
img = img.merge(ecm[["sample", "group", "donor", "sd_around_0_deg", "perc_within_15deg"]], on="sample")
assert (img["folder"] == img["group"]).all()                          # folder agrees with the csv subtype
img["alias"] = "P" + img["donor"].str[4:]                               # P220815T1 -> P815T1, as on every figure
assert img["alias"].str.len().eq(6).all() and img["alias"].is_unique
img = img.sort_values(["group", "sd_around_0_deg"]).reset_index(drop=True)
print(f"STEP 1  {len(img)} images ({(img.group == 'proCAF').sum()} proCAF, {(img.group == 'restCAF').sum()} restCAF); "
      f"{PX_PER_UM_EXPECTED} px/um = {1 / PX_PER_UM_EXPECTED:.5f} um/px, image width {2586 / PX_PER_UM_EXPECTED:.1f} um")

# STEP 2. Page geometry, all in pt from the top-left corner
IMG = (PAGE_W - (N_COL - 1) * GAP) / N_COL                                # pt, side of one image
n_rows = {g: int(np.ceil((img.group == g).sum() / N_COL)) for g in ["proCAF", "restCAF"]}
block_h = {g: TITLE_H + n_rows[g] * (LABEL_H + IMG) + (n_rows[g] - 1) * GAP for g in n_rows}
PAGE_H = block_h["proCAF"] + BLOCK_GAP + block_h["restCAF"]
bar_len = SCALE_UM * PX_PER_UM_EXPECTED / 2586 * IMG                      # pt
px_out = int(round(IMG / 72 * IMG_DPI))
print(f"STEP 2  page {PAGE_W:.1f} x {PAGE_H:.1f} pt; image {IMG:.1f} pt ({px_out} px embedded, {IMG_DPI} dpi); "
      f"{SCALE_UM}-um bar = {bar_len:.2f} pt")
assert PAGE_H <= 9 * 72

# STEP 3. Draw block titles, labels, images (G channel only) and scale bars
# 908 and 0118N1 lack one strip of the tile scan (rows 795-903, zero signal; 908 marks it in blue)
# brightness is not equalized (exports used different display settings): not comparable between images
fig = plt.figure(figsize=(PAGE_W / 72, PAGE_H / 72))
W, H = PAGE_W, PAGE_H
y_top = 0.0
for g, col in [("proCAF", PRO), ("restCAF", REST)]:
    d = img[img.group == g].reset_index(drop=True)
    fig.text(0, 1 - (y_top + PT_TITLE) / H, f"{g} matrix ({len(d)} cfECM donors)", color=col,
             fontsize=PT_TITLE, fontweight="bold", ha="left", va="baseline")
    for i, r in d.iterrows():
        row, c = divmod(i, N_COL)
        x0 = c * (IMG + GAP)
        y_img = y_top + TITLE_H + row * (LABEL_H + IMG + GAP) + LABEL_H       # top edge of the image
        fig.text(x0 / W, 1 - (y_img - 2) / H, f"{r.alias}  ({r.sd_around_0_deg:.1f}°)",
                 fontsize=PT_LABEL, ha="left", va="baseline")
        a = np.asarray(Image.open(r.file).convert("RGB"))
        green = np.zeros_like(a); green[..., 1] = a[..., 1]                  # the signal only; drops the blue zero marker
        small = np.asarray(Image.fromarray(green).resize((px_out, px_out), Image.LANCZOS))
        ax = fig.add_axes([x0 / W, 1 - (y_img + IMG) / H, IMG / W, IMG / H])
        ax.imshow(small, interpolation="none", extent=(0, IMG, 0, IMG))       # axis units = pt
        ax.set_xlim(0, IMG); ax.set_ylim(0, IMG); ax.axis("off")
        bx = IMG - BAR_MARGIN - bar_len
        ax.add_patch(Rectangle((bx, BAR_MARGIN), bar_len, BAR_H, facecolor="white", edgecolor="none"))
        ax.text(bx + bar_len / 2, BAR_MARGIN + BAR_H + 1.5, f"{SCALE_UM} µm", color="white",
                fontsize=PT_LABEL, ha="center", va="bottom")
    y_top += block_h[g] + BLOCK_GAP

out_pdf = f"{OUT_DIR}/suppfig1_if_images.pdf"
fig.savefig(out_pdf)
fig.savefig(out_pdf.replace(".pdf", "_200dpi.png"), dpi=200, facecolor="white")
plt.close(fig)
print(f"STEP 3  {out_pdf}  ({os.path.getsize(out_pdf) / 1e6:.1f} MB)")

# STEP 4. One row per image, for the legend
img["order_in_block"] = img.groupby("group").cumcount() + 1
img["file"] = img["file"].str.replace(ROOT + "/", "", regex=False)
img[["group", "order_in_block", "alias", "donor", "sample", "sd_around_0_deg", "perc_within_15deg", "px_per_um", "file"]] \
    .round(2).to_csv(f"{OUT_DIR}/suppfig1_images.csv", index=False)
print(f"STEP 4  spread proCAF {img[img.group=='proCAF'].sd_around_0_deg.min():.1f}-{img[img.group=='proCAF'].sd_around_0_deg.max():.1f} deg, "
      f"restCAF {img[img.group=='restCAF'].sd_around_0_deg.min():.1f}-{img[img.group=='restCAF'].sd_around_0_deg.max():.1f} deg")
