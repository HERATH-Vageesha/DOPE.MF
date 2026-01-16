#################################################################################################################################
#################################   SMLM IMAGE ANALYSIS TEMPLATE (FULL PIPELINE)   ##############################################
#################################################################################################################################

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from datetime import datetime
from zoneinfo import ZoneInfo

from scipy.spatial import cKDTree
from scipy.stats import pearsonr
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist

from matplotlib.ticker import MaxNLocator


# =================================================================================================
# USER-ADJUSTABLE THRESHOLDS (edit these as needed)
# =================================================================================================
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

# =================================================================================================
# COLUMN NAMES (edit only if your CSV columns differ)
# =================================================================================================
ID_COL = "id"
FRAME_COL = "frame"
XCOL = "x [nm]"
YCOL = "y [nm]"
UNCERTAINTY_COL = "uncertainty_xy [nm]"
INTENSITY_COL = "intensity [photon]"   # set to None if your file does not have intensity

# =================================================================================================
# PARAMETERS
# =================================================================================================
RADIUS_NM = 232.0              # C1–C2 pairing radius
TRACK_LINK_NM = 400.0          # linking threshold for tracking midpoints across frames

# Colors (enforced consistently)
C1_COLOR = "green"
C2_COLOR = "red"

# Timestamped outputs (NY time)
ts = datetime.now(ZoneInfo("America/New_York")).strftime("%Y%m%d_%H%M%S")
OUT_XLSX = f"End-to-end distance_{ts}.xlsx"
TRACKED_XLSX = f"Tracked_Dipoles_{ts}.xlsx"


# =================================================================================================
# LOAD + FILTER
#   - Loads the needed columns once
#   - Applies uncertainty, (optional) intensity, and XY window thresholds
# =================================================================================================
def load_and_filter(
    csv_path: str,
    lower_unc: float, upper_unc: float,
    x_lo: float, x_hi: float,
    y_lo: float, y_hi: float,
    intensity_lo: float, intensity_hi: float,
    id_col: str = ID_COL,
    frame_col: str = FRAME_COL,
    xcol: str = XCOL,
    ycol: str = YCOL,
    ucol: str = UNCERTAINTY_COL,
    icol: str | None = INTENSITY_COL,
) -> pd.DataFrame:

    usecols = [id_col, frame_col, xcol, ycol, ucol]
    if icol is not None:
        usecols.append(icol)

    df = pd.read_csv(csv_path, usecols=usecols).dropna(subset=usecols).reset_index(drop=True)

    m = (
        (df[ucol] >= lower_unc) & (df[ucol] <= upper_unc) &
        (df[xcol] >= x_lo) & (df[xcol] <= x_hi) &
        (df[ycol] >= y_lo) & (df[ycol] <= y_hi)
    )
    if icol is not None:
        m = m & (df[icol] >= intensity_lo) & (df[icol] <= intensity_hi)

    return df.loc[m].reset_index(drop=True)


# -------------------------------------------------------------------------------------------------
# Load CSV files at the top (as requested)
# -------------------------------------------------------------------------------------------------
df_c1 = load_and_filter(
    "TIRF560_imageregperformed.csv",
    lower_unc=lower_threshold_c1, upper_unc=upper_threshold_c1,
    x_lo=x_lower, x_hi=x_upper, y_lo=y_lower, y_hi=y_upper,
    intensity_lo=intensity_lower_c1, intensity_hi=intensity_upper_c1
)

df_c2 = load_and_filter(
    "TIRF647_imageregperformed.csv",
    lower_unc=lower_threshold_c2, upper_unc=upper_threshold_c2,
    x_lo=x_lower, x_hi=x_upper, y_lo=y_lower, y_hi=y_upper,
    intensity_lo=intensity_lower_c2, intensity_hi=intensity_upper_c2
)

print(f"C1 after thresholds: {len(df_c1)}")
print(f"C2 after thresholds: {len(df_c2)}")


# =================================================================================================
# HIGH-POPULATION AMBIGUITY DELETION (your rule)
#   If a same-channel close pair (<= RADIUS_NM) exists AND either member is close to opposite channel
#   (<= RADIUS_NM), remove BOTH same-channel puncta AND the opposite-channel puncta within RADIUS_NM.
# =================================================================================================
def remove_ambiguous_triplets(df_c1: pd.DataFrame, df_c2: pd.DataFrame, r_nm: float = 232.0) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    c1_xy = df_c1[[XCOL, YCOL]].to_numpy(dtype=float)
    c2_xy = df_c2[[XCOL, YCOL]].to_numpy(dtype=float)

    remove_c1 = np.zeros(len(df_c1), dtype=bool)
    remove_c2 = np.zeros(len(df_c2), dtype=bool)

    if len(c1_xy) == 0 or len(c2_xy) == 0:
        report = {
            "r_nm": r_nm,
            "initial_c1": len(df_c1), "initial_c2": len(df_c2),
            "same_channel_pairs_checked_c1": 0, "same_channel_pairs_checked_c2": 0,
            "removed_c1": 0, "removed_c2": 0,
            "final_c1": len(df_c1), "final_c2": len(df_c2),
        }
        return df_c1, df_c2, report

    c1_tree = cKDTree(c1_xy)
    c2_tree = cKDTree(c2_xy)

    def process_same_channel_pairs(same_tree, same_xy, other_tree, remove_same, remove_other) -> int:
        pairs = same_tree.query_pairs(r=r_nm)
        if not pairs:
            return 0

        for (i, j) in pairs:
            neigh_i = other_tree.query_ball_point(same_xy[i], r=r_nm)
            neigh_j = other_tree.query_ball_point(same_xy[j], r=r_nm)

            if neigh_i or neigh_j:
                remove_same[i] = True
                remove_same[j] = True
                if neigh_i:
                    remove_other[neigh_i] = True
                if neigh_j:
                    remove_other[neigh_j] = True

        return len(pairs)

    n_c1_pairs = process_same_channel_pairs(c1_tree, c1_xy, c2_tree, remove_c1, remove_c2)
    n_c2_pairs = process_same_channel_pairs(c2_tree, c2_xy, c1_tree, remove_c2, remove_c1)

    df_c1_filt = df_c1.loc[~remove_c1].reset_index(drop=True)
    df_c2_filt = df_c2.loc[~remove_c2].reset_index(drop=True)

    report = {
        "r_nm": r_nm,
        "initial_c1": len(df_c1), "initial_c2": len(df_c2),
        "same_channel_pairs_checked_c1": n_c1_pairs, "same_channel_pairs_checked_c2": n_c2_pairs,
        "removed_c1": int(remove_c1.sum()), "removed_c2": int(remove_c2.sum()),
        "final_c1": len(df_c1_filt), "final_c2": len(df_c2_filt),
    }
    return df_c1_filt, df_c2_filt, report


# -------------------------------------------------------------------------------------------------
# APPLY AMBIGUITY DELETION (THIS IS THE HIGH-POPULATION REMOVAL STEP)
# -------------------------------------------------------------------------------------------------
df_c1, df_c2, deletion_report = remove_ambiguous_triplets(df_c1, df_c2, r_nm=RADIUS_NM)
print("High-population ambiguity deletion report:", deletion_report)


# =================================================================================================
# QC SCATTER PLOT (C1 green, C2 red)
# =================================================================================================
plt.figure(figsize=(8, 6))
plt.scatter(df_c1[XCOL], df_c1[YCOL], color=C1_COLOR, alpha=0.2, label="C1 (TIRF 560)")
plt.scatter(df_c2[XCOL], df_c2[YCOL], color=C2_COLOR, alpha=0.2, label="C2 (TIRF 647)")
plt.xlabel("X Position (nm)")
plt.ylabel("Y Position (nm)")
plt.title("Scatter Plot: CHANNELS BEFORE IMAGE REGISTRATION")
plt.legend()
plt.grid(True)
plt.savefig("scatter_plot.png", dpi=300, bbox_inches="tight")
plt.show()


# =================================================================================================
# PART A: FRAME-AGNOSTIC PAIRING (your template output)
#   - Computes all C1–C2 pairs within RADIUS_NM (across all frames)
#   - Writes a timestamped Excel output: OUT_XLSX
# =================================================================================================
c1_xy = df_c1[[XCOL, YCOL]].to_numpy(dtype=float)
c2_xy = df_c2[[XCOL, YCOL]].to_numpy(dtype=float)

if len(c1_xy) == 0 or len(c2_xy) == 0:
    distance_df = pd.DataFrame(columns=[
        "C1 X (nm)", "C1 Y (nm)", "C2 X (nm)", "C2 Y (nm)",
        "Distance (nm)", "C1 Uncertainty (nm)", "C2 Uncertainty (nm)",
        "Dipole Angle (degrees)"
    ])
else:
    c2_tree_global = cKDTree(c2_xy)
    neighbors = c2_tree_global.query_ball_point(c1_xy, r=RADIUS_NM)

    c1_idx = np.repeat(np.arange(len(neighbors)), [len(n) for n in neighbors])
    if c1_idx.size == 0:
        distance_df = pd.DataFrame(columns=[
            "C1 X (nm)", "C1 Y (nm)", "C2 X (nm)", "C2 Y (nm)",
            "Distance (nm)", "C1 Uncertainty (nm)", "C2 Uncertainty (nm)",
            "Dipole Angle (degrees)"
        ])
    else:
        c2_idx = np.concatenate([np.asarray(n, dtype=int) for n in neighbors])
        c1_sel = c1_xy[c1_idx]
        c2_sel = c2_xy[c2_idx]

        dxy = c2_sel - c1_sel
        dist = np.sqrt((dxy ** 2).sum(axis=1))
        angles = (np.degrees(np.arctan2(dxy[:, 1], dxy[:, 0])) + 360) % 360

        distance_df = pd.DataFrame({
            "C1 X (nm)": c1_sel[:, 0],
            "C1 Y (nm)": c1_sel[:, 1],
            "C2 X (nm)": c2_sel[:, 0],
            "C2 Y (nm)": c2_sel[:, 1],
            "Distance (nm)": dist,
            "C1 Uncertainty (nm)": df_c1.loc[c1_idx, UNCERTAINTY_COL].to_numpy(),
            "C2 Uncertainty (nm)": df_c2.loc[c2_idx, UNCERTAINTY_COL].to_numpy(),
            "Dipole Angle (degrees)": angles,
        })

distance_df.to_excel(OUT_XLSX, index=False)
print(f"Saved output: {OUT_XLSX}")


# -------------------------------------------------------------------------------------------------
# End-to-end distance histogram (NO GRIDLINES + LESS CLUTTERED Y AXIS)
# -------------------------------------------------------------------------------------------------
plt.figure(figsize=(10, 6))
plt.hist(distance_df["Distance (nm)"], bins=30, color="blue", alpha=0.7, edgecolor="black")
plt.xlabel("Distance (nm)")
plt.ylabel("Frequency")
plt.title("END-TO-END DISTANCE DISTRIBUTION (nm)")

# Less clutter: cap number of major ticks
ax = plt.gca()
ax.yaxis.set_major_locator(MaxNLocator(nbins=6, integer=True))

# No gridlines here by request
plt.savefig("Distance_Histogram.png", dpi=300, bbox_inches="tight")
plt.show()

# Polar histogram of azimuth angles (template)
plt.figure(figsize=(8, 8))
ax = plt.subplot(111, projection="polar")
angles_rad = np.radians(distance_df["Dipole Angle (degrees)"].to_numpy())
ax.hist(angles_rad, bins=30, color="pink", alpha=0.3, edgecolor="black")
ax.set_theta_zero_location("E")
ax.set_theta_direction(1)
ax.set_title("Dipole Angle Distribution (Degrees)")
plt.show()


# =================================================================================================
# PART B: FRAME-AWARE PAIRING + Φ/θ + MIDPOINTS + TRACKING + TRACKED OUTPUT + ADDITIONAL PLOTS
#   - Pairs C1↔C2 ONLY within the same frame
#   - Computes Φ and θ (as in your second script)
#   - Tracks dipole midpoints across frames using Hungarian assignment
#   - Saves TRACKED_XLSX
# =================================================================================================
# Verify required columns exist (id + frame)
required_cols = {ID_COL, FRAME_COL, XCOL, YCOL, UNCERTAINTY_COL}
missing_c1 = required_cols - set(df_c1.columns)
missing_c2 = required_cols - set(df_c2.columns)
if missing_c1 or missing_c2:
    raise ValueError(
        f"Missing required columns for frame-aware pairing/tracking.\n"
        f"C1 missing: {sorted(missing_c1)}\n"
        f"C2 missing: {sorted(missing_c2)}"
    )

# Build per-frame C2 trees
grouped_c2 = {frame: group.reset_index(drop=True) for frame, group in df_c2.groupby(FRAME_COL)}
tree_dict = {
    frame: cKDTree(group[[XCOL, YCOL]].to_numpy(dtype=float))
    for frame, group in grouped_c2.items()
}

results = []
for _, row in df_c1.iterrows():
    frame_val = row[FRAME_COL]
    if frame_val not in tree_dict:
        continue

    idxs = tree_dict[frame_val].query_ball_point([row[XCOL], row[YCOL]], RADIUS_NM)
    if not idxs:
        continue

    group = grouped_c2[frame_val]
    for j in idxs:
        other_row = group.iloc[j]
        dx = other_row[XCOL] - row[XCOL]
        dy = other_row[YCOL] - row[YCOL]
        dist = float(np.sqrt(dx * dx + dy * dy))

        results.append([
            row[ID_COL], row[FRAME_COL], row[XCOL], row[YCOL],
            other_row[ID_COL], other_row[FRAME_COL], other_row[XCOL], other_row[YCOL],
            dist, row[UNCERTAINTY_COL], other_row[UNCERTAINTY_COL]
        ])

cols = [
    "C1 id", "C1 Frame", "C1 X (nm)", "C1 Y (nm)",
    "C2 id", "C2 Frame", "C2 X (nm)", "C2 Y (nm)",
    "Distance (nm)", "C1 Uncertainty (nm)", "C2 Uncertainty (nm)"
]
distance_df_tracked = pd.DataFrame(results, columns=cols)

if len(distance_df_tracked) == 0:
    print("No frame-matched dipoles found (after filters). Tracking and tracked plots skipped.")
else:
    # Φ
    dxy_x = distance_df_tracked["C2 X (nm)"].to_numpy() - distance_df_tracked["C1 X (nm)"].to_numpy()
    dxy_y = distance_df_tracked["C2 Y (nm)"].to_numpy() - distance_df_tracked["C1 Y (nm)"].to_numpy()
    phi_deg = (np.degrees(np.arctan2(dxy_y, dxy_x)) + 360) % 360
    distance_df_tracked["Φ (degrees)"] = phi_deg

    # θ (your provided model)
    ratio = np.clip(distance_df_tracked["Distance (nm)"].to_numpy() / 120.0, -1.0, 1.0)
    theta_deg = np.degrees(np.arccos(ratio))
    distance_df_tracked["θ (degrees)"] = theta_deg

    # Remove θ == 0 (as your second script)
    distance_df_tracked = distance_df_tracked[distance_df_tracked["θ (degrees)"] != 0].reset_index(drop=True)

    # Midpoints
    distance_df_tracked["mid_x"] = (distance_df_tracked["C1 X (nm)"] + distance_df_tracked["C2 X (nm)"]) / 2.0
    distance_df_tracked["mid_y"] = (distance_df_tracked["C1 Y (nm)"] + distance_df_tracked["C2 Y (nm)"]) / 2.0

    # TRACKING across frames
    distance_df_tracked = distance_df_tracked.sort_values(by="C1 Frame").reset_index(drop=True)
    distance_df_tracked["Track ID"] = np.nan

    unique_frames = sorted(distance_df_tracked["C1 Frame"].unique())
    next_track_id = 1
    tracks_prev = {}  # track_id -> last midpoint position

    for frame in unique_frames:
        mask = distance_df_tracked["C1 Frame"] == frame
        current = distance_df_tracked.loc[mask].copy()
        curr_pos = current[["mid_x", "mid_y"]].to_numpy()

        if frame == unique_frames[0]:
            n = len(current)
            if n > 0:
                new_ids = np.arange(next_track_id, next_track_id + n)
                distance_df_tracked.loc[mask, "Track ID"] = new_ids
                next_track_id += n
                for k, idx in enumerate(current.index):
                    tracks_prev[int(new_ids[k])] = np.array([
                        distance_df_tracked.loc[idx, "mid_x"],
                        distance_df_tracked.loc[idx, "mid_y"]
                    ])
            continue

        if len(tracks_prev) == 0 or len(curr_pos) == 0:
            n = len(curr_pos)
            if n > 0:
                new_ids = np.arange(next_track_id, next_track_id + n)
                distance_df_tracked.loc[mask, "Track ID"] = new_ids
                next_track_id += n
                for k, idx in enumerate(current.index):
                    tracks_prev[int(new_ids[k])] = np.array([
                        distance_df_tracked.loc[idx, "mid_x"],
                        distance_df_tracked.loc[idx, "mid_y"]
                    ])
            continue

        prev_ids = list(tracks_prev.keys())
        prev_pos = np.vstack([tracks_prev[t] for t in prev_ids])

        cost = cdist(prev_pos, curr_pos)
        r_ind, c_ind = linear_sum_assignment(cost)

        assigned = np.full(len(curr_pos), np.nan)
        for r, c in zip(r_ind, c_ind):
            if cost[r, c] < TRACK_LINK_NM:
                assigned[c] = prev_ids[r]

        for i in range(len(curr_pos)):
            if np.isnan(assigned[i]):
                assigned[i] = next_track_id
                next_track_id += 1

        distance_df_tracked.loc[mask, "Track ID"] = assigned

        current_ids = distance_df_tracked.loc[mask, "Track ID"].to_numpy().astype(int)
        for i, idx in enumerate(current.index):
            tracks_prev[int(current_ids[i])] = np.array([
                distance_df_tracked.loc[idx, "mid_x"],
                distance_df_tracked.loc[idx, "mid_y"]
            ])

    # SAVE TRACKED EXCEL (timestamped)
    with pd.ExcelWriter(TRACKED_XLSX, engine="xlsxwriter") as writer:
        distance_df_tracked.to_excel(writer, sheet_name="Master", index=False)
        for track_id, group in distance_df_tracked.groupby("Track ID"):
            group_sorted = group.sort_values("C1 Frame")
            track_df = group_sorted[[
                "C1 Frame", "Distance (nm)", "Φ (degrees)", "θ (degrees)",
                "C1 X (nm)", "C1 Y (nm)", "C2 X (nm)", "C2 Y (nm)"
            ]]
            track_df.to_excel(writer, sheet_name=f"Track_{int(track_id)}", index=False)

    print(f"Saved tracked output: {TRACKED_XLSX}")

    # GLOBAL ANALYSIS PLOTS (tracked set)
    plt.figure(figsize=(10, 6))
    plt.hist(distance_df_tracked["Distance (nm)"], bins=30, color="blue", alpha=0.7, edgecolor="black")
    plt.xlabel("Distance (nm)")
    plt.ylabel("Frequency")
    plt.title("End-to-End Distance Distribution (nm) - All Tracks")
    plt.savefig("Distance_Histogram_Filtered.png", dpi=300, bbox_inches="tight")
    plt.show()

    plt.figure(figsize=(8, 6))
    plt.scatter(distance_df_tracked["Φ (degrees)"], distance_df_tracked["Distance (nm)"], color="purple", alpha=0.5)
    plt.xlabel("Φ Angle (degrees)")
    plt.ylabel("End-to-End Distance (nm)")
    plt.title("Scatter Plot: Distance vs. Φ Angle (All Tracks)")
    plt.grid(True)
    plt.show()

    if len(distance_df_tracked) > 1:
        corr_phi, p_val_phi = pearsonr(distance_df_tracked["Φ (degrees)"], distance_df_tracked["Distance (nm)"])
        print(f"Φ Pearson Correlation: {corr_phi:.4f}, P-Value: {p_val_phi:.4f}")

    # Midpoints colored by θ
    plt.figure(figsize=(10, 8))
    sc = plt.scatter(
        distance_df_tracked["mid_x"], distance_df_tracked["mid_y"],
        c=distance_df_tracked["θ (degrees)"],
        cmap="viridis", vmin=0, vmax=90,
        s=20, alpha=0.9,
        edgecolors="k", linewidths=0.1
    )
    cbar = plt.colorbar(sc)
    cbar.set_label("θ (degrees)", rotation=270, labelpad=15)
    plt.xlabel("X Position (nm)")
    plt.ylabel("Y Position (nm)")
    plt.title("Dipole Midpoints Colored by θ Angle")
    plt.grid(True)
    plt.gca().invert_yaxis()
    plt.savefig("midpoint_theta_colormap.png", dpi=300, bbox_inches="tight")
    plt.show()

    # Midpoints colored by Φ
    plt.figure(figsize=(10, 8))
    sc_phi = plt.scatter(
        distance_df_tracked["mid_x"], distance_df_tracked["mid_y"],
        c=distance_df_tracked["Φ (degrees)"],
        cmap="plasma", vmin=0, vmax=360,
        s=20, alpha=0.9,
        edgecolors="k", linewidths=0.1
    )
    cbar_phi = plt.colorbar(sc_phi)
    cbar_phi.set_label("Φ (degrees)", rotation=270, labelpad=15)
    plt.xlabel("X Position (nm)")
    plt.ylabel("Y Position (nm)")
    plt.title("Dipole Midpoints Colored by Φ Angle")
    plt.grid(True)
    plt.gca().invert_yaxis()
    plt.savefig("midpoint_phi_colormap.png", dpi=300, bbox_inches="tight")
    plt.show()

    # Midpoints colored by Φ with arrows
    plt.figure(figsize=(10, 8))
    sc_phi2 = plt.scatter(
        distance_df_tracked["mid_x"], distance_df_tracked["mid_y"],
        c=distance_df_tracked["Φ (degrees)"],
        cmap="plasma", vmin=0, vmax=360,
        s=20, alpha=0.9,
        edgecolors="k", linewidths=0.1
    )

    arrow_length = 1800.0
    phi_rad = np.radians(distance_df_tracked["Φ (degrees)"].to_numpy())
    dx = np.cos(phi_rad) * arrow_length
    dy = np.sin(phi_rad) * arrow_length

    plt.quiver(
        distance_df_tracked["mid_x"], distance_df_tracked["mid_y"],
        dx, dy,
        angles="xy", scale_units="xy", scale=1,
        color="black", width=0.0025, alpha=0.7
    )

    cbar_phi2 = plt.colorbar(sc_phi2)
    cbar_phi2.set_label("Φ (degrees)", rotation=270, labelpad=15)
    plt.xlabel("X Position (nm)")
    plt.ylabel("Y Position (nm)")
    plt.title("Dipole Midpoints Colored by Φ Angle with Arrows")
    plt.grid(True)
    plt.gca().invert_yaxis()
    plt.savefig("midpoint_phi_colormap_arrows.png", dpi=300, bbox_inches="tight")
    plt.show()

    # Midpoints colored by θ with Φ arrows
    plt.figure(figsize=(10, 8))
    sc_theta = plt.scatter(
        distance_df_tracked["mid_x"], distance_df_tracked["mid_y"],
        c=distance_df_tracked["θ (degrees)"],
        cmap="viridis", vmin=0, vmax=90,
        s=20, alpha=0.9,
        edgecolors="k", linewidths=0.1
    )

    plt.quiver(
        distance_df_tracked["mid_x"], distance_df_tracked["mid_y"],
        dx, dy,
        angles="xy", scale_units="xy", scale=1,
        color="black", width=0.0025, alpha=0.7
    )

    cbar_theta = plt.colorbar(sc_theta)
    cbar_theta.set_label("θ (degrees)", rotation=270, labelpad=15)
    plt.xlabel("X Position (nm)")
    plt.ylabel("Y Position (nm)")
    plt.title("Dipole Midpoints Colored by θ Angle with Φ Arrows")
    plt.grid(True)
    plt.gca().invert_yaxis()
    plt.savefig("midpoint_theta_colormap_arrows.png", dpi=300, bbox_inches="tight")
    plt.show()
