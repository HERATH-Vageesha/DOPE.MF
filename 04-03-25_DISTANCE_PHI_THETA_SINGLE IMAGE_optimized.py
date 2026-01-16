##################################################################################################
###########################   SMLM IMAGE ANALYSIS (FULL + FIXED + ALL PLOTS)   #####################
##################################################################################################

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from scipy.spatial import cKDTree
from scipy.stats import pearsonr
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist


# =============================================================================
# USER THRESHOLDS
# =============================================================================
lower_threshold_c1 = 0
upper_threshold_c1 = 40
lower_threshold_c2 = 0
upper_threshold_c2 = 40

x_lower = 0
x_upper = 80000
y_lower = 0
y_upper = 80000

intensity_lower_c1 = 0
intensity_upper_c1 = 100000
intensity_lower_c2 = 0
intensity_upper_c2 = 100000


# =============================================================================
# PARAMETERS
# =============================================================================
RADIUS_NM = 232.0
TRACK_LINK_NM = 400.0
ROD_LENGTH_NM = 120.0

ONE_TO_ONE_PER_FRAME = True
DROP_AMBIGUOUS_TRACKS = True

C1_COLOR = "green"
C2_COLOR = "red"


# =============================================================================
# COLUMN NAMES
# =============================================================================
ID_COL = "id"
FRAME_COL = "frame"
XCOL = "x [nm]"
YCOL = "y [nm]"
UNCOL = "uncertainty_xy [nm]"
ICOL = "intensity [photon]"


# =============================================================================
# LOAD DATA
# =============================================================================
cols = [ID_COL, FRAME_COL, XCOL, YCOL, UNCOL, ICOL]
df_c1 = pd.read_csv("TIRF560_imageregperformed.csv", usecols=cols)
df_c2 = pd.read_csv("TIRF560_imageregperformed.csv", usecols=cols)


# =============================================================================
# FILTERING
# =============================================================================
def apply_filters(df, u_lo, u_hi, x_lo, x_hi, y_lo, y_hi, i_lo, i_hi):
    df = df.dropna().copy()
    return df[
        (df[UNCOL] >= u_lo) & (df[UNCOL] <= u_hi) &
        (df[ICOL] >= i_lo) & (df[ICOL] <= i_hi) &
        (df[XCOL] >= x_lo) & (df[XCOL] <= x_hi) &
        (df[YCOL] >= y_lo) & (df[YCOL] <= y_hi)
    ].reset_index(drop=True)


df_c1 = apply_filters(df_c1, lower_threshold_c1, upper_threshold_c1,
                      x_lower, x_upper, y_lower, y_upper,
                      intensity_lower_c1, intensity_upper_c1)

df_c2 = apply_filters(df_c2, lower_threshold_c2, upper_threshold_c2,
                      x_lower, x_upper, y_lower, y_upper,
                      intensity_lower_c2, intensity_upper_c2)


# =============================================================================
# INITIAL SCATTER (UNCHANGED)
# =============================================================================
plt.figure(figsize=(8, 6))
plt.scatter(df_c1[XCOL], df_c1[YCOL], color=C1_COLOR, alpha=0.3, label="TIRF 560")
plt.scatter(df_c2[XCOL], df_c2[YCOL], color=C2_COLOR, alpha=0.3, label="TIRF 647")
plt.xlabel("X Position (nm)")
plt.ylabel("Y Position (nm)")
plt.title("Scatter Plot of C1 and C2 Channels")
plt.legend()
plt.grid(True)
plt.gca().invert_yaxis()
plt.savefig("scatter_c1_c2.png", dpi=300, bbox_inches="tight")
plt.show()


# =============================================================================
# SAME-CHANNEL CROWDING + OPPOSITE CHANNEL REMOVAL
# =============================================================================
def remove_ambiguous_triplets_framewise(df1, df2, r):
    rem1 = np.zeros(len(df1), dtype=bool)
    rem2 = np.zeros(len(df2), dtype=bool)

    frames = set(df1[FRAME_COL]).union(df2[FRAME_COL])
    for fr in frames:
        idx1 = df1.index[df1[FRAME_COL] == fr].to_numpy()
        idx2 = df2.index[df2[FRAME_COL] == fr].to_numpy()
        if len(idx1) == 0 or len(idx2) == 0:
            continue

        xy1 = df1.loc[idx1, [XCOL, YCOL]].to_numpy()
        xy2 = df2.loc[idx2, [XCOL, YCOL]].to_numpy()

        t1, t2 = cKDTree(xy1), cKDTree(xy2)

        for i, j in t1.query_pairs(r):
            if t2.query_ball_point(xy1[i], r) or t2.query_ball_point(xy1[j], r):
                rem1[idx1[[i, j]]] = True
                rem2[idx2[t2.query_ball_point(xy1[i], r)]] = True
                rem2[idx2[t2.query_ball_point(xy1[j], r)]] = True

        for i, j in t2.query_pairs(r):
            if t1.query_ball_point(xy2[i], r) or t1.query_ball_point(xy2[j], r):
                rem2[idx2[[i, j]]] = True
                rem1[idx1[t1.query_ball_point(xy2[i], r)]] = True
                rem1[idx1[t1.query_ball_point(xy2[j], r)]] = True

    return df1.loc[~rem1].reset_index(drop=True), df2.loc[~rem2].reset_index(drop=True)


df_c1, df_c2 = remove_ambiguous_triplets_framewise(df_c1, df_c2, RADIUS_NM)


# =============================================================================
# ONE-TO-ONE PAIRING PER FRAME
# =============================================================================
rows = []

for fr in sorted(set(df_c1[FRAME_COL]).intersection(df_c2[FRAME_COL])):
    g1 = df_c1[df_c1[FRAME_COL] == fr].reset_index(drop=True)
    g2 = df_c2[df_c2[FRAME_COL] == fr].reset_index(drop=True)

    xy1 = g1[[XCOL, YCOL]].to_numpy()
    xy2 = g2[[XCOL, YCOL]].to_numpy()
    if len(xy1) == 0 or len(xy2) == 0:
        continue

    D = cdist(xy1, xy2)
    D[D > RADIUS_NM] = 1e9
    r, c = linear_sum_assignment(D)

    for i, j in zip(r, c):
        if D[i, j] > RADIUS_NM:
            continue

        dx = xy2[j, 0] - xy1[i, 0]
        dy = xy2[j, 1] - xy1[i, 1]
        dist = np.hypot(dx, dy)

        rows.append([
            g1.loc[i, ID_COL], fr, xy1[i, 0], xy1[i, 1],
            g2.loc[j, ID_COL], fr, xy2[j, 0], xy2[j, 1],
            dist, g1.loc[i, UNCOL], g2.loc[j, UNCOL]
        ])

distance_df = pd.DataFrame(rows, columns=[
    "C1 id", "C1 Frame", "C1 X (nm)", "C1 Y (nm)",
    "C2 id", "C2 Frame", "C2 X (nm)", "C2 Y (nm)",
    "Distance (nm)", "C1 Uncertainty (nm)", "C2 Uncertainty (nm)"
])


# =============================================================================
# ANGLES + MIDPOINTS (θ AS NaN WHEN INVALID)
# =============================================================================
dx = distance_df["C2 X (nm)"] - distance_df["C1 X (nm)"]
dy = distance_df["C2 Y (nm)"] - distance_df["C1 Y (nm)"]

distance_df["Φ (degrees)"] = (np.degrees(np.arctan2(dy, dx)) + 360) % 360

theta = np.full(len(distance_df), np.nan)
valid = distance_df["Distance (nm)"] <= ROD_LENGTH_NM
theta[valid] = np.degrees(np.arccos(distance_df.loc[valid, "Distance (nm)"] / ROD_LENGTH_NM))
distance_df["θ (degrees)"] = theta

distance_df["mid_x"] = (distance_df["C1 X (nm)"] + distance_df["C2 X (nm)"]) / 2
distance_df["mid_y"] = (distance_df["C1 Y (nm)"] + distance_df["C2 Y (nm)"]) / 2


# =============================================================================
# TRACKING
# =============================================================================
distance_df["Track ID"] = np.nan
next_id = 1
prev = {}

for fr in sorted(distance_df["C1 Frame"].unique()):
    mask = distance_df["C1 Frame"] == fr
    pts = distance_df.loc[mask, ["mid_x", "mid_y"]].to_numpy()

    if fr == distance_df["C1 Frame"].min():
        ids = np.arange(next_id, next_id + len(pts))
        distance_df.loc[mask, "Track ID"] = ids
        next_id += len(pts)
        for i, idx in enumerate(distance_df[mask].index):
            prev[int(ids[i])] = pts[i]
        continue

    prev_ids = list(prev.keys())
    prev_pts = np.vstack([prev[k] for k in prev_ids])
    cost = cdist(prev_pts, pts)
    r, c = linear_sum_assignment(cost)

    assigned = np.full(len(pts), np.nan)
    for i, j in zip(r, c):
        if cost[i, j] < TRACK_LINK_NM:
            assigned[j] = prev_ids[i]

    for i in range(len(pts)):
        if np.isnan(assigned[i]):
            assigned[i] = next_id
            next_id += 1

    distance_df.loc[mask, "Track ID"] = assigned
    for i, idx in enumerate(distance_df[mask].index):
        prev[int(assigned[i])] = pts[i]


# =============================================================================
# MASKS FOR ZERO DISTANCE
# =============================================================================
zero_dist_mask = distance_df["Distance (nm)"] == 0
nonzero_dist_mask = ~zero_dist_mask


# =============================================================================
# ALL ORIGINAL PLOTS (UNCHANGED, WITH ZERO-DISTANCE FIX)
# =============================================================================

# Distance histogram
plt.figure(figsize=(10, 6))
plt.hist(distance_df["Distance (nm)"], bins=30, edgecolor="black")
plt.xlabel("Distance (nm)")
plt.ylabel("Frequency")
plt.title("End-to-End Distance Distribution (nm)")
plt.show()

# Φ vs Distance
plt.figure(figsize=(8, 6))
plt.scatter(distance_df["Φ (degrees)"], distance_df["Distance (nm)"], alpha=0.5)
plt.xlabel("Φ (degrees)")
plt.ylabel("Distance (nm)")
plt.grid(True)
plt.show()

# θ-colored midpoints
plt.figure(figsize=(10, 8))
plt.scatter(distance_df["mid_x"], distance_df["mid_y"],
            c=distance_df["θ (degrees)"], cmap="viridis", vmin=0, vmax=90)
plt.colorbar(label="θ (degrees)")
plt.gca().invert_yaxis()
plt.show()

# Φ-colored midpoints
plt.figure(figsize=(10, 8))
plt.scatter(distance_df["mid_x"], distance_df["mid_y"],
            c=distance_df["Φ (degrees)"], cmap="plasma", vmin=0, vmax=360)
plt.colorbar(label="Φ (degrees)")
plt.gca().invert_yaxis()
plt.show()

# Φ arrows (NO arrows for distance == 0)
plt.figure(figsize=(10, 8))
plt.scatter(distance_df["mid_x"], distance_df["mid_y"],
            c=distance_df["Φ (degrees)"], cmap="plasma", vmin=0, vmax=360)

phi = np.radians(distance_df.loc[nonzero_dist_mask, "Φ (degrees)"])
plt.quiver(distance_df.loc[nonzero_dist_mask, "mid_x"],
           distance_df.loc[nonzero_dist_mask, "mid_y"],
           np.cos(phi)*1800, np.sin(phi)*1800)

for _, r in distance_df.loc[zero_dist_mask].iterrows():
    plt.text(r["mid_x"], r["mid_y"], "NaN", ha="center")

plt.colorbar(label="Φ (degrees)")
plt.gca().invert_yaxis()
plt.show()

# θ-colored with Φ arrows
plt.figure(figsize=(10, 8))
plt.scatter(distance_df["mid_x"], distance_df["mid_y"],
            c=distance_df["θ (degrees)"], cmap="viridis", vmin=0, vmax=90)

plt.quiver(distance_df.loc[nonzero_dist_mask, "mid_x"],
           distance_df.loc[nonzero_dist_mask, "mid_y"],
           np.cos(phi)*1800, np.sin(phi)*1800)

for _, r in distance_df.loc[zero_dist_mask].iterrows():
    plt.text(r["mid_x"], r["mid_y"], "NaN", ha="center")

plt.colorbar(label="θ (degrees)")
plt.gca().invert_yaxis()
plt.show()


