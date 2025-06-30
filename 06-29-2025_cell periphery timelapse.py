import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
from scipy.stats import pearsonr
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist
import seaborn as sns

# =============================================================================
# SET THRESHOLDS FOR FILTERING
# =============================================================================

lower_threshold_c1 = 0
upper_threshold_c1 = 10
lower_threshold_c2 = 0
upper_threshold_c2 = 10
x_lower = 0
x_upper = 80000
y_lower = 0
y_upper = 80000
intensity_lower_c1 = 1000
intensity_upper_c1 = 100000
intensity_lower_c2 = 800
intensity_upper_c2 = 100000

# =============================================================================
# LOAD CSV DATA
# =============================================================================

columns_to_use = ["id", "frame", "x [nm]", "y [nm]", "uncertainty_xy [nm]", "intensity [photon]"]
df_c1 = pd.read_csv("C1-Cy3b_ThunderStorm_ROI_ImageReg Performed.csv", usecols=columns_to_use)
df_c2 = pd.read_csv("C2-ATTO_ThunderStorm_ROI_ImageReg Performed.csv", usecols=columns_to_use)

# =============================================================================
# FILTER DATA
# =============================================================================

df_c1 = df_c1[(df_c1["uncertainty_xy [nm]"] >= lower_threshold_c1) &
              (df_c1["uncertainty_xy [nm]"] <= upper_threshold_c1)]
df_c2 = df_c2[(df_c2["uncertainty_xy [nm]"] >= lower_threshold_c2) &
              (df_c2["uncertainty_xy [nm]"] <= upper_threshold_c2)]
df_c1 = df_c1[(df_c1["intensity [photon]"] >= intensity_lower_c1) &
              (df_c1["intensity [photon]"] <= intensity_upper_c1)]
df_c2 = df_c2[(df_c2["intensity [photon]"] >= intensity_lower_c2) &
              (df_c2["intensity [photon]"] <= intensity_upper_c2)]
df_c1 = df_c1[(df_c1["x [nm]"] >= x_lower) & (df_c1["x [nm]"] <= x_upper) &
              (df_c1["y [nm]"] >= y_lower) & (df_c1["y [nm]"] <= y_upper)]
df_c2 = df_c2[(df_c2["x [nm]"] >= x_lower) & (df_c2["x [nm]"] <= x_upper) &
              (df_c2["y [nm]"] >= y_lower) & (df_c2["y [nm]"] <= y_upper)]

# =============================================================================
# SCATTER PLOT WITH TIRF 560/647 CHANNELS
# =============================================================================

plt.figure(figsize=(8,6))
plt.scatter(df_c1["x [nm]"], df_c1["y [nm]"], color='green', alpha=0.3, label='TIRF 560')
plt.scatter(df_c2["x [nm]"], df_c2["y [nm]"], color='red', alpha=0.3, label='TIRF 647')
plt.xlabel("X Position (nm)")
plt.ylabel("Y Position (nm)")
plt.title("Scatter Plot of C1 and C2 Channels")
plt.legend()
plt.grid(True)
plt.gca().invert_yaxis()
plt.savefig("scatter_c1_c2.png", dpi=300, bbox_inches="tight")
plt.show()

# =============================================================================
# DISTANCE CALCULATIONS OF DIPOLES FOR TIRF
# =============================================================================

radius_threshold = 232
results = []
grouped_c2 = {frame: group for frame, group in df_c2.groupby("frame")}
tree_dict = {frame: cKDTree(group[["x [nm]", "y [nm]"]].values) for frame, group in grouped_c2.items()}

for i, row in df_c1.iterrows():
    frame_val = row["frame"]
    if frame_val in tree_dict:
        indices = tree_dict[frame_val].query_ball_point([row["x [nm]"], row["y [nm]"]], radius_threshold)
        group = grouped_c2[frame_val]
        for idx in indices:
            other_row = group.iloc[idx]
            dist = np.sqrt((other_row["x [nm]"] - row["x [nm]"])**2 +
                           (other_row["y [nm]"] - row["y [nm]"])**2)
            results.append([
                row["id"], row["frame"], row["x [nm]"], row["y [nm]"],
                other_row["id"], other_row["frame"], other_row["x [nm]"], other_row["y [nm]"],
                dist, row["uncertainty_xy [nm]"], other_row["uncertainty_xy [nm]"]
            ])

cols = ["C1 id", "C1 Frame", "C1 X (nm)", "C1 Y (nm)",
        "C2 id", "C2 Frame", "C2 X (nm)", "C2 Y (nm)",
        "Distance (nm)", "C1 Uncertainty (nm)", "C2 Uncertainty (nm)"]
distance_df = pd.DataFrame(results, columns=cols)

# =============================================================================
# ANGLES AND MIDPOINT
# =============================================================================

angles_phi_rad = np.arctan2(distance_df["C2 Y (nm)"] - distance_df["C1 Y (nm)"],
                            distance_df["C2 X (nm)"] - distance_df["C1 X (nm)"])
distance_df["\u03a6 (degrees)"] = (np.degrees(angles_phi_rad) + 360) % 360
distance_ratio = np.clip(distance_df["Distance (nm)"] / 120, -1, 1)
theta_rad = np.arccos(distance_ratio)
distance_df["\u03b8 (degrees)"] = np.degrees(theta_rad)
distance_df["mid_x"] = (distance_df["C1 X (nm)"] + distance_df["C2 X (nm)"]) / 2
distance_df["mid_y"] = (distance_df["C1 Y (nm)"] + distance_df["C2 Y (nm)"]) / 2

# =============================================================================
# TRACKING PARTICLES ACROSS FRAMES (with 3-frame gap allowed)
# =============================================================================

distance_df = distance_df.sort_values(by="C1 Frame").reset_index(drop=True)
distance_df["Track ID"] = np.nan
unique_frames = sorted(distance_df["C1 Frame"].unique())
max_frame_gap = 3
track_id_counter = 1

tracks_prev = {}

for frame in unique_frames:
    current_mask = distance_df["C1 Frame"] == frame
    current_data = distance_df[current_mask].copy()
    current_positions = current_data[["mid_x", "mid_y"]].to_numpy()

    assigned_ids = np.full(len(current_data), np.nan)

    valid_prev_tracks = {tid: (pos, last_frame) for tid, (pos, last_frame) in tracks_prev.items()
                         if (frame - last_frame) <= max_frame_gap}

    if valid_prev_tracks and len(current_positions) > 0:
        prev_track_ids = list(valid_prev_tracks.keys())
        prev_positions = np.array([valid_prev_tracks[tid][0] for tid in prev_track_ids])
        cost_matrix = cdist(prev_positions, current_positions)
        row_ind, col_ind = linear_sum_assignment(cost_matrix)

        for r, c in zip(row_ind, col_ind):
            if cost_matrix[r, c] < 400:
                assigned_ids[c] = prev_track_ids[r]

    for i in range(len(current_positions)):
        if np.isnan(assigned_ids[i]):
            assigned_ids[i] = track_id_counter
            track_id_counter += 1

    distance_df.loc[current_mask, "Track ID"] = assigned_ids
    current_ids = distance_df.loc[current_mask, "Track ID"].values.astype(int)
    for i, idx in enumerate(current_data.index):
        tracks_prev[int(current_ids[i])] = (
            distance_df.loc[idx, ["mid_x", "mid_y"]].to_numpy(),
            frame
        )

# =============================================================================
# SAVE TRACKED DIPOLES TO EXCEL
# =============================================================================

with pd.ExcelWriter("Tracked_Dipoles.xlsx", engine="xlsxwriter") as writer:
    distance_df.to_excel(writer, sheet_name="Master", index=False)
    for track_id, group in distance_df.groupby("Track ID"):
        group_sorted = group.sort_values("C1 Frame")
        track_df = group_sorted[["C1 Frame", "Distance (nm)", "\u03a6 (degrees)", "\u03b8 (degrees)",
                                 "C1 X (nm)", "C1 Y (nm)", "C2 X (nm)", "C2 Y (nm)"]]
        sheet_name = f"Track_{int(track_id)}"
        track_df.to_excel(writer, sheet_name=sheet_name, index=False)

# =============================================================================
# GLOBAL ANALYSIS PLOTS
# =============================================================================

plt.figure(figsize=(10, 6))
plt.hist(distance_df["Distance (nm)"], bins=30, color='blue', alpha=0.7, edgecolor='black')
plt.xlabel("Distance (nm)")
plt.ylabel("Frequency")
plt.title("End-to-End Distance Distribution (nm) - All Tracks")
plt.savefig("Distance_Histogram_Filtered.png", dpi=300, bbox_inches="tight")
plt.show()

plt.figure(figsize=(8,6))
plt.scatter(distance_df["\u03a6 (degrees)"], distance_df["Distance (nm)"], color='purple', alpha=0.5)
plt.xlabel("\u03a6 Angle (degrees)")
plt.ylabel("End-to-End Distance (nm)")
plt.title("Scatter Plot: Distance vs. \u03a6 Angle (All Tracks)")
plt.grid(True)
corr_phi, p_val_phi = pearsonr(distance_df["\u03a6 (degrees)"], distance_df["Distance (nm)"])
plt.show()
print(f"\u03a6 Pearson Correlation: {corr_phi:.4f}, P-Value: {p_val_phi:.4f}")

# =============================================================================
# DIPLOE MIDPOINT COLOR-CODED BY THETA ANGLE
# =============================================================================

plt.figure(figsize=(10, 8))
sc = plt.scatter(
    distance_df["mid_x"],
    distance_df["mid_y"],
    c=distance_df["\u03b8 (degrees)"],
    cmap="viridis",
    vmin=0,
    vmax=90,
    s=20,
    alpha=0.9,
    edgecolors="k",
    linewidths=0.1
)
cbar = plt.colorbar(sc)
cbar.set_label("\u03b8 (degrees)", rotation=270, labelpad=15)
plt.xlabel("X Position (nm)")
plt.ylabel("Y Position (nm)")
plt.title("Dipole Midpoints Colored by \u03b8 Angle")
plt.grid(True)
plt.gca().invert_yaxis()
plt.savefig("midpoint_theta_colormap.png", dpi=300, bbox_inches="tight")
plt.show()

# =============================================================================
# DIPLOE MIDPOINT COLOR-CODED BY PHI ANGLE WITH ARROWS
# =============================================================================

plt.figure(figsize=(10, 8))
sc_phi = plt.scatter(
    distance_df["mid_x"],
    distance_df["mid_y"],
    c=distance_df["\u03a6 (degrees)"],
    cmap="plasma",
    vmin=0,
    vmax=360,
    s=20,
    alpha=0.9,
    edgecolors="k",
    linewidths=0.1
)
arrow_length = 1800
phi_rad = np.radians(distance_df["\u03a6 (degrees)"])
dx = np.cos(phi_rad) * arrow_length
dy = np.sin(phi_rad) * arrow_length
plt.quiver(
    distance_df["mid_x"],
    distance_df["mid_y"],
    dx,
    dy,
    angles='xy',
    scale_units='xy',
    scale=1,
    color='black',
    width=0.0025,
    headwidth=300,
    headlength=400,
    headaxislength=350,
    alpha=0.7
)
cbar_phi = plt.colorbar(sc_phi)
cbar_phi.set_label("\u03a6 (degrees)", rotation=270, labelpad=15)
plt.xlabel("X Position (nm)")
plt.ylabel("Y Position (nm)")
plt.title("Dipole Midpoints Colored by \u03a6 Angle with Arrows")
plt.grid(True)
plt.gca().invert_yaxis()
plt.savefig("midpoint_phi_colormap_arrows.png", dpi=300, bbox_inches="tight")
plt.show()
