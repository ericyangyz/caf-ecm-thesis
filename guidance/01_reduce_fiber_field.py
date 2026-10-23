#!/usr/bin/env python3
# Reduce one field's per-pixel fiber orientation (results/fib_wholefield/, 48 frames) to 8 x 8 px
# block sums, so later steps never re-read the 950 MB per-frame pixel CSVs.
# Out: results/guidance/fiber_blocks/<category>/<field>_block8.npz (C, S, W: float32, (48, H//8, W//8))
# Usage: python guidance/01_reduce_fiber_field.py rest_basal P201015_field_01

import sys
import os
import glob
import time
import numpy as np
import pandas as pd


# USER INPUTS
ROOT = "."
BLOCK = 8            # 8 px = 5.12 um. Aggregate up later; never down.
N_FRAMES = 48

CATEGORY = sys.argv[1]          # e.g. rest_basal
FIELD = sys.argv[2]             # e.g. P201015_field_01

# tracking side spells it PancT6, ECM side PanT6 (deliberate)
ECM_FIELD = FIELD.replace("PancT6", "PanT6")


# STEP 1. Find this field's Pixel folder
found = glob.glob(f"{ROOT}/results/fib_wholefield/{CATEGORY}/{ECM_FIELD}/"
                  f"Alignment_assay_results_*/Pixel")
assert len(found) == 1, f"expected exactly one Pixel dir for {ECM_FIELD}, got {found}"
PIXEL_DIR = found[0]

OUT_DIR = f"{ROOT}/results/guidance/fiber_blocks/{CATEGORY}"
os.makedirs(OUT_DIR, exist_ok=True)
OUT_PATH = f"{OUT_DIR}/{FIELD}_block{BLOCK}.npz"

print(f"STEP 1  {FIELD} ({CATEGORY})  ->  {OUT_PATH}")


# STEP 2. Reduce each of the 48 frames to block sums
started = time.time()
C_frames, S_frames, W_frames = [], [], []
height = width = None

for t in range(N_FRAMES):
    csv = f"{PIXEL_DIR}/{ECM_FIELD}_T_{t:04d}_fibronectin_processed_pixel.csv"
    pixels = pd.read_csv(csv, usecols=["x", "y", "theta_deg", "coherence", "energy"],
                         dtype={"x": np.int32, "y": np.int32, "theta_deg": np.float32,
                                "coherence": np.float64, "energy": np.float32})

    # the reshape below assumes rows are in C order (x fastest, then y); checked once per field
    if height is None:
        width, height = int(pixels.x.max()) + 1, int(pixels.y.max()) + 1
        assert pixels.x.values[1] == 1 and pixels.y.values[1] == 0, "rows are not in C order"
        assert pixels.x.values[width] == 0 and pixels.y.values[width] == 1, "row stride != width"
        assert len(pixels) == width * height, "the pixel grid is not complete"
    theta = pixels.theta_deg.values.reshape(height, width)
    coherence = pixels.coherence.values.reshape(height, width)
    energy = pixels.energy.values.reshape(height, width)

    # weight in [0,1]: raw coherence is heavy-tailed (~0 to 3e12), so use UMA's own
    # log-normalization (clip to [10, 1e8], log, rescale); zero where energy == 0 (no gradient)
    w = np.log(np.clip(coherence, 10.0, 1e8))
    w -= w.min()
    w /= w.max()
    w[energy == 0] = 0.0

    # fibers are axial (theta == theta + 180), so sum the doubled-angle vector instead of angles:
    # C = sum w cos 2theta, S = sum w sin 2theta, W = sum w. Then angle = 0.5 atan2(S, C) and
    # R = sqrt(C^2 + S^2) / W (0..1). Sums are additive, so blocks can be merged exactly later.
    doubled = np.deg2rad(2.0 * theta)
    cos_part = w * np.cos(doubled)
    sin_part = w * np.sin(doubled)

    # sum over BLOCK x BLOCK tiles; the ragged last rows/columns are cropped
    n_rows, n_cols = height // BLOCK, width // BLOCK
    keep = (slice(0, n_rows * BLOCK), slice(0, n_cols * BLOCK))
    C_frames.append(cos_part[keep].reshape(n_rows, BLOCK, n_cols, BLOCK)
                    .sum(axis=(1, 3)).astype(np.float32))
    S_frames.append(sin_part[keep].reshape(n_rows, BLOCK, n_cols, BLOCK)
                    .sum(axis=(1, 3)).astype(np.float32))
    W_frames.append(w[keep].reshape(n_rows, BLOCK, n_cols, BLOCK)
                    .sum(axis=(1, 3)).astype(np.float32))

    print(f"        T_{t:04d}  {time.time()-started:6.0f}s", flush=True)


# STEP 3. Save
np.savez_compressed(OUT_PATH + ".tmp.npz",
                    C=np.stack(C_frames), S=np.stack(S_frames), W=np.stack(W_frames),
                    block=BLOCK, um_per_px=0.64, frames=N_FRAMES,
                    field_h=height, field_w=width)
# write to .tmp first, then rename, so a concurrent reader never sees a partial file
os.replace(OUT_PATH + ".tmp.npz", OUT_PATH)

print(f"STEP 3  saved {OUT_PATH}  "
      f"({os.path.getsize(OUT_PATH)/1e6:.0f} MB, {time.time()-started:.0f}s)")
