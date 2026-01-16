import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable

# =============================================================================
# INPUT (EDIT THIS ONE LINE AS NEEDED)
# =============================================================================
# Load the Excel file
filename = 'End-to-end distance_20260115_204819.xlsx'
data = pd.read_excel(filename)

# =============================================================================
# COMPUTE VECTORS / MIDPOINTS / DISTANCES / ANGLES
# =============================================================================
data["dx"] = data["C2 X (nm)"] - data["C1 X (nm)"]
data["dy"] = data["C2 Y (nm)"] - data["C1 Y (nm)"]

# Distance
data["Distance (nm)"] = np.sqrt(data["dx"] ** 2 + data["dy"] ** 2)

# Angle for coloring + normalization and for arrow rotation
# NOTE: For text rotation, using the signed angle in degrees keeps the orientation correct.
data["angle_signed (deg)"] = np.degrees(np.arctan2(data["dy"], data["dx"]))

# Optional: also store 0–360 representation (useful if you want it later)
data["angle_0to360 (deg)"] = (data["angle_signed (deg)"] + 360) % 360

# Midpoints (for arrow placement)
data["x_mid"] = data["C1 X (nm)"] + data["dx"] / 2
data["y_mid"] = data["C1 Y (nm)"] + data["dy"] / 2

# =============================================================================
# COLORMAP: RED -> GREEN (short -> long or long -> short depending on vmin/vmax)
# =============================================================================
cmap = LinearSegmentedColormap.from_list("red_green", ["red", "green"])
norm = Normalize(vmin=data["Distance (nm)"].min(), vmax=data["Distance (nm)"].max())

# =============================================================================
# VECTOR PLOT 1: ARROWS + DISTANCE/ANGLE ANNOTATIONS (from your first script)
# =============================================================================
fig, ax = plt.subplots(figsize=(16, 12))

for _, row in data.iterrows():
    arrow_color = cmap(norm(row["Distance (nm)"]))

    # Arrow rendered as text (matches your original style)
    ax.text(
        row["x_mid"], row["y_mid"],
        " ---> ",
        color=arrow_color,
        rotation=row["angle_signed (deg)"],     # signed angle gives correct rotation
        rotation_mode="anchor",
        ha="center", va="center",
        fontsize=16
    )

    # Annotation below arrow: distance + angle
    annotation = f"{row['Distance (nm)']:.1f} nm\n{row['angle_signed (deg)']:.1f}°"
    ax.text(
        row["x_mid"], row["y_mid"] - 400,
        annotation,
        color="black",
        ha="center", va="center",
        fontsize=10
    )

ax.set_xlabel("X (nm)")
ax.set_ylabel("Y (nm)")
ax.set_title("Vector Plot with Distance and Angle Annotations")
ax.grid(True)

# Fixed axis limits / aspect (as in your scripts)
ax.set_xlim(0, 80000)
ax.set_ylim(0, 80000)
ax.set_aspect("equal")
ax.set_autoscale_on(False)

# Colorbar
sm = ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
fig.colorbar(sm, label="Distance (nm)", ax=ax)

plt.show()

# =============================================================================
# VECTOR PLOT 2: CLEANER "DIPOLE DIRECTIONS" VIEW (from your second script)
#   - Only arrows, colored by distance, no extra annotation text
# =============================================================================
fig2, ax2 = plt.subplots(figsize=(12, 10))

for _, row in data.iterrows():
    ax2.text(
        row["x_mid"], row["y_mid"],
        "--->",
        color=cmap(norm(row["Distance (nm)"])),
        rotation=row["angle_signed (deg)"],  # signed angle aligns with arrow direction
        ha="center", va="center",
        fontsize=10
    )

ax2.set(
    xlabel="X (nm)",
    ylabel="Y (nm)",
    title="Dipole Directions (Colored by Distance)",
    xlim=(0, 80000),
    ylim=(0, 80000),
    aspect="equal",
)
ax2.grid(True)

fig2.colorbar(ScalarMappable(norm=norm, cmap=cmap), ax=ax2, label="Distance (nm)")
plt.show()
