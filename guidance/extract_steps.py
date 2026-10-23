#!/usr/bin/env python3
# Per-step and per-track tables from the trackpy trajectories in data/data_tracking/analysis.
# Earlier-stage script: trackpy is a separate tracking from the TrackMate export, so the
# TrackMate-named track features below are not a reproduction of TrackMate's values.
# Out: results/guidance/steps_all.csv (~0.93 M steps), results/guidance/tracks_all.csv (~78 k tracks)

import os
import glob
import zipfile
import xml.etree.ElementTree as ET
import numpy as np
import pandas as pd


# USER INPUTS
ROOT = "."
TRACK_DIR = f"{ROOT}/data/data_tracking/analysis"
OUT_DIR = f"{ROOT}/results/guidance"

UM_PER_PX = 0.64          # micrometers per pixel
MIN_PER_FRAME = 30.0      # minutes between consecutive frames

CATEGORIES = ["pro_basal", "rest_basal", "pro_classical", "rest_classical"]

os.makedirs(OUT_DIR, exist_ok=True)


# STEP 1. Read tracking_metadata.xlsx: which cfECM donor is each field?
# no openpyxl in the env, so read the xlsx (a zip of XML) directly;
# cells with t="s" are indexes into sharedStrings.xml
NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"

with zipfile.ZipFile(f"{TRACK_DIR}/tracking_metadata.xlsx") as z:
    shared = ["".join(t.text or "" for t in si.iter(f"{NS}t"))
              for si in ET.fromstring(z.read("xl/sharedStrings.xml"))]
    sheet = ET.fromstring(z.read("xl/worksheets/sheet1.xml"))

meta_rows = []
for row in sheet.iter(f"{NS}row"):
    vals = []
    for cell in row:
        v = cell.find(f"{NS}v")
        txt = "" if v is None else v.text
        if cell.get("t") == "s":
            txt = shared[int(txt)]
        vals.append(txt)
    meta_rows.append(vals)

meta = pd.DataFrame(meta_rows[1:], columns=meta_rows[0])
meta["field"] = meta.Tumor + "_field_" + meta.Field.astype(int).map("{:02d}".format)
meta = meta.rename(columns={"cfECM": "donor"})[["field", "donor", "cafType", "tumorType"]]

# field number is not a stable donor ID across tumor lines (P201015_field_11 = P220212T1,
# PancT6_field_11 = P240711N1), so always join on donor
print(f"STEP 1  metadata: {len(meta)} field entries, {meta.donor.nunique()} cfECM donors")


# STEP 2. Read every trajectories_*.csv into one long table
pieces = []
for category in CATEGORIES:
    for path in sorted(glob.glob(f"{TRACK_DIR}/{category}/trajectories_*.csv")):
        field = os.path.basename(path)[len("trajectories_"):-len(".csv")]
        t = pd.read_csv(path)                       # columns: frame, x, y (px), particle
        t["field"] = field
        t["category"] = category
        pieces.append(t)
        print(f"        {field:22s} {len(t):6d} detections", flush=True)

traj = pd.concat(pieces, ignore_index=True)
traj["cafType"] = traj.category.str.split("_").str[0]     # 'pro' / 'rest'
traj["tumorType"] = traj.category.str.split("_").str[1]   # 'basal' / 'classical'
traj = traj.merge(meta[["field", "donor"]], on="field", how="left")

# sort so that .diff() below walks each cell forward in time
traj = traj.sort_values(["field", "particle", "frame"]).reset_index(drop=True)
print(f"STEP 2  {len(traj):,} detections, {traj.field.nunique()} fields")


# STEP 3. Consecutive detections -> steps
# each row is the move that ends at frame_to; the first row of every cell is NaN and dropped
cell = traj.groupby(["field", "particle"], sort=False)

steps = pd.DataFrame({
    "field": traj.field, "category": traj.category, "cafType": traj.cafType,
    "tumorType": traj.tumorType, "donor": traj.donor, "particle": traj.particle,
    "frame_to": traj.frame,
    "dframe": cell.frame.diff(),
    "dx_px": cell.x.diff(),
    "dy_px": cell.y.diff(),
    "x_to": traj.x, "y_to": traj.y,
}).dropna(subset=["dx_px"])

steps["frame_from"] = (steps.frame_to - steps.dframe).astype(int)
steps["x_from"] = steps.x_to - steps.dx_px
steps["y_from"] = steps.y_to - steps.dy_px

steps = steps[steps.dframe == 1].copy()

steps["step_um"] = np.hypot(steps.dx_px, steps.dy_px) * UM_PER_PX
steps["speed_um_min"] = steps.step_um / MIN_PER_FRAME

# direction in (-180, 180], 0 = east, +90 = north; image y is down, fiber angles have y up,
# hence atan2(-dy, dx)
steps["step_dir_deg"] = np.degrees(np.arctan2(-steps.dy_px, steps.dx_px))

steps["mid_x_px"] = (steps.x_from + steps.x_to) / 2.0
steps["mid_y_px"] = (steps.y_from + steps.y_to) / 2.0

# turn angle at the start of each step, vs the previous step of the same cell.
# shift(1) is the adjacent step only if tracks have no frame gaps, so check that.
_gaps = traj.groupby(["field", "particle"]).frame.diff().dropna()
assert (_gaps == 1).all(), "some tracks contain frame gaps; shift(1) would span them"
step_of_cell = steps.groupby(["field", "particle"], sort=False)
prev_dir = step_of_cell.step_dir_deg.shift(1)
steps["prev_step_um"] = step_of_cell.step_um.shift(1)

# unsigned turn in [0, 180]: 0 = straight on, 90 = right angle, 180 = reversal
steps["turn_deg"] = ((steps.step_dir_deg - prev_dir + 180.0) % 360.0 - 180.0).abs()

# persistence: +1 straight, 0 no memory, -1 reversals
steps["cos_turn"] = np.cos(np.radians(steps.turn_deg))
# collinearity: +1 for straight and reversal, -1 for right angles (motion confined to a 1-D track)
steps["cos2_turn"] = np.cos(2 * np.radians(steps.turn_deg))

print(f"STEP 3  {len(steps):,} steps, {steps.turn_deg.notna().sum():,} with a turn angle")


# STEP 4. Per-track features, named after TrackMate's (TrackMate name in the column comments)
tracks = steps.groupby(["field", "category", "cafType", "tumorType", "donor", "particle"],
                       sort=False).agg(
    n_steps=("step_um", "size"),
    total_distance_um=("step_um", "sum"),             # TOTAL_DISTANCE_TRAVELED
    median_step_um=("step_um", "median"),
    mean_speed_um_min=("speed_um_min", "mean"),       # TRACK_MEAN_SPEED
    mean_turn_deg=("turn_deg", "mean"),
    x_start=("x_from", "first"), y_start=("y_from", "first"),
    x_end=("x_to", "last"), y_end=("y_to", "last"),
).reset_index()

tracks["n_detections"] = tracks.n_steps + 1
tracks["duration_min"] = tracks.n_steps * MIN_PER_FRAME                     # TRACK_DURATION
tracks["displacement_um"] = np.hypot(tracks.x_end - tracks.x_start,         # TRACK_DISPLACEMENT
                                     tracks.y_end - tracks.y_start) * UM_PER_PX

# LINEARITY_OF_FORWARD_PROGRESSION = net displacement / path length, 0..1. For a random walk it
# decays like 1/sqrt(n_steps), so it is not comparable between tracks of different length.
tracks["linearity"] = tracks.displacement_um / tracks.total_distance_um

# MEAN_DIRECTIONAL_CHANGE_RATE: mean turn (rad) / 30 min, rad/min, 0..pi/30. 0.0524 = random
# directions (mean turn 90 deg), 0.0698 = pure localization noise (120 deg), so jitter raises it.
tracks["mean_dcr_rad_min"] = np.radians(tracks.mean_turn_deg) / MIN_PER_FRAME

tracks = tracks.drop(columns=["x_start", "y_start", "x_end", "y_end"])
print(f"STEP 4  {len(tracks):,} tracks")


# STEP 5. Write both tables
keep = ["field", "category", "cafType", "tumorType", "donor", "particle",
        "frame_from", "frame_to",
        # raw positions in px: x_from, y_from at frame_from; x_to, y_to at frame_to
        "x_from", "y_from", "x_to", "y_to",
        "mid_x_px", "mid_y_px", "dx_px", "dy_px",
        "step_um", "speed_um_min", "step_dir_deg", "turn_deg", "prev_step_um",
        "cos_turn", "cos2_turn"]
steps[keep].to_csv(f"{OUT_DIR}/steps_all.csv", index=False, float_format="%.5g")
tracks.to_csv(f"{OUT_DIR}/tracks_all.csv", index=False, float_format="%.6g")

print(f"STEP 5  wrote steps_all.csv and tracks_all.csv to {OUT_DIR}")
print(f"        fields {steps.field.nunique()}   donors {steps.donor.nunique()}")
