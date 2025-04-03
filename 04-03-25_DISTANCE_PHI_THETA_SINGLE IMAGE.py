#################################################################################################################################
################################# LOADING DATA, FILTERING BY UNCERTAINTY, INTENSITY & COORDINATES, AND CALCULATE DISTANCE AND ANGLE VALUES  ############
#################################################################################################################################

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree

# UNCERTAINTY FILTERATION
lower_threshold = 0    # Set your desired lower threshold for uncertainty
upper_threshold = 10   # Set your desired upper threshold for uncertainty

# ROI FILTERATION (X and Y coordinate cutoff)
x_lower = 00000
x_upper = 80000
y_lower = 0000
y_upper = 80000

# INTENSITY FILTERATION
intensity_lower = 0       # Set your desired lower threshold for intensity [photon]
intensity_upper = 100000  # Set your desired upper threshold for intensity [photon]

# Load CSV files (including the intensity column)
columns_to_use = ["x [nm]", "y [nm]", "uncertainty_xy [nm]", "intensity [photon]"]
df_c1 = pd.read_csv("Channel-01.csv", usecols=columns_to_use)  # MUST BE THE CY3B CHANNEL
df_c2 = pd.read_csv("Channel-02.csv", usecols=columns_to_use)  # MUST BE THE ATTO647N CHANNEL

# Filter data based on uncertainty values for both channels
df_c1 = df_c1[(df_c1["uncertainty_xy [nm]"] >= lower_threshold) & (df_c1["uncertainty_xy [nm]"] <= upper_threshold)]
df_c2 = df_c2[(df_c2["uncertainty_xy [nm]"] >= lower_threshold) & (df_c2["uncertainty_xy [nm]"] <= upper_threshold)]

# Filter data based on intensity values for both channels
df_c1 = df_c1[(df_c1["intensity [photon]"] >= intensity_lower) & (df_c1["intensity [photon]"] <= intensity_upper)]
df_c2 = df_c2[(df_c2["intensity [photon]"] >= intensity_lower) & (df_c2["intensity [photon]"] <= intensity_upper)]

# Filter data based on X and Y coordinates for both channels
df_c1 = df_c1[(df_c1["x [nm]"] >= x_lower) & (df_c1["x [nm]"] <= x_upper) &
              (df_c1["y [nm]"] >= y_lower) & (df_c1["y [nm]"] <= y_upper)]
df_c2 = df_c2[(df_c2["x [nm]"] >= x_lower) & (df_c2["x [nm]"] <= x_upper) &
              (df_c2["y [nm]"] >= y_lower) & (df_c2["y [nm]"] <= y_upper)]

# Scatter plot before image registration (using filtered data)
plt.figure(figsize=(8, 6))
plt.scatter(df_c1["x [nm]"], df_c1["y [nm]"], color='green', alpha=0.2, label='TIRF 560')
plt.scatter(df_c2["x [nm]"], df_c2["y [nm]"], color='red', alpha=0.2, label='TIRF 647')
plt.xlabel("X Position (nm)")
plt.ylabel("Y Position (nm)")
plt.title("Scatter Plot: CHANNELS BEFORE IMAGE REGISTRATION\n(Filtered by Uncertainty, Intensity & Coordinates)")
plt.legend()
plt.grid(True)
plt.savefig("scatter_plot_filtered.png", dpi=300, bbox_inches="tight")
plt.show()

# Calculate distances within 232 nm radius
radius_threshold = 232
c2_tree = cKDTree(df_c2[["x [nm]", "y [nm]"]])
results = []
for i, row in df_c1.iterrows():
    indices = c2_tree.query_ball_point([row["x [nm]"], row["y [nm]"]], radius_threshold)
    for idx in indices:
        dist = np.sqrt((df_c2.iloc[idx]["x [nm]"] - row["x [nm]"])**2 + 
                       (df_c2.iloc[idx]["y [nm]"] - row["y [nm]"])**2)
        results.append([row["x [nm]"], row["y [nm]"],
                        df_c2.iloc[idx]["x [nm]"], df_c2.iloc[idx]["y [nm]"],
                        dist, row["uncertainty_xy [nm]"], df_c2.iloc[idx]["uncertainty_xy [nm]"]])

# Save distance data
columns = ["C1 X (nm)", "C1 Y (nm)", "C2 X (nm)", "C2 Y (nm)", 
           "Distance (nm)", "C1 Uncertainty (nm)", "C2 Uncertainty (nm)"]
distance_df = pd.DataFrame(results, columns=columns)
distance_df.to_excel("End-to-End Distance_PUNCTA_Filtered.xlsx", index=False)

# Plot histogram of distances
plt.figure(figsize=(10, 6))
plt.hist(distance_df["Distance (nm)"], bins=30, color='blue', alpha=0.7, edgecolor='black')
plt.xlabel("Distance (nm)")
plt.ylabel("Frequency")
plt.title("END-TO-END DISTANCE DISTRIBUTION (nm)")
plt.yticks(range(0, int(max(plt.hist(distance_df["Distance (nm)"], bins=30)[0])) + 1))
plt.savefig("Distance_Histogram_Filtered.png", dpi=300, bbox_inches="tight")
plt.show()

# -----------------------------------------------------------------------------
# Compute angles and update the Excel file with both PHI and THETA angles
# -----------------------------------------------------------------------------

# Compute PHI angle (dipole angle) and denote with "Φ" symbol
angles_phi_rad = np.arctan2(distance_df["C2 Y (nm)"] - distance_df["C1 Y (nm)"],
                             distance_df["C2 X (nm)"] - distance_df["C1 X (nm)"])
distance_df["Φ (degrees)"] = (np.degrees(angles_phi_rad) + 360) % 360

# Compute THETA angle (denoted by "θ") using cosine relation: cos(θ) = Distance / 120 nm
# Clip the ratio to the valid range for arccos
distance_ratio = np.clip(distance_df["Distance (nm)"] / 120, -1, 1)
theta_rad = np.arccos(distance_ratio)
distance_df["θ (degrees)"] = np.degrees(theta_rad)

# Save updated distance data with both angles
distance_df.to_excel("End-to-End Distance_PUNCTA_Filtered.xlsx", index=False)

# -----------------------------------------------------------------------------
# Plot polar histogram for PHI angle (Φ)
plt.figure(figsize=(8, 8))
ax_phi = plt.subplot(111, projection='polar')
phi_rad = np.radians(distance_df["Φ (degrees)"])
ax_phi.hist(phi_rad, bins=30, color='pink', alpha=0.3, edgecolor='black')
ax_phi.set_theta_zero_location('E')
ax_phi.set_theta_direction(1)
ax_phi.set_yticks(range(0, int(max(ax_phi.hist(phi_rad, bins=30)[0])) + 1))
ax_phi.set_title("Φ Angle Distribution (Degrees)")
plt.show()

# -----------------------------------------------------------------------------
# Plot polar histogram for THETA angle (θ)
plt.figure(figsize=(8, 8))
ax_theta = plt.subplot(111, projection='polar')
theta_rad_plot = np.radians(distance_df["θ (degrees)"])
ax_theta.hist(theta_rad_plot, bins=30, color='lightblue', alpha=0.3, edgecolor='black')
ax_theta.set_theta_zero_location('E')
ax_theta.set_theta_direction(1)
ax_theta.set_yticks(range(0, int(max(ax_theta.hist(theta_rad_plot, bins=30)[0])) + 1))
ax_theta.set_title("θ Angle Distribution (Degrees)")
plt.show()


################################################################################################################################
############## DISTANCE/ANGLE ANALYSIS : DISTANCE VS ANGLE & DISTANCE VS UNCERTAINTY & CHROMATIC ABBERATION DETECTION
###############################################################################################################################
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import pearsonr

# Load the saved Excel file (with updated angles)
file_name = "End-to-End Distance_PUNCTA_Filtered.xlsx"
distance_df = pd.read_excel(file_name)

# Scatter Plot: Distance vs. Φ Angle (PHI)
plt.figure(figsize=(8, 6))
plt.scatter(distance_df["Φ (degrees)"], distance_df["Distance (nm)"], color='purple', alpha=0.5)
plt.xlabel("Φ Angle (degrees)")
plt.ylabel("End-to-End Distance (nm)")
plt.title("Scatter Plot: Distance vs. Φ Angle")
plt.grid(True)
correlation_phi, p_value_phi = pearsonr(distance_df["Φ (degrees)"], distance_df["Distance (nm)"])
plt.show()
print(f"Φ Pearson Correlation Coefficient: {correlation_phi:.4f}, P-Value: {p_value_phi:.4f}")

# Scatter Plot: End-to-End Distance vs. Uncertainty
plt.figure(figsize=(8, 6))
plt.scatter(distance_df["Distance (nm)"], distance_df["C1 Uncertainty (nm)"], color='green', alpha=0.5, label="C1 Uncertainty")
plt.scatter(distance_df["Distance (nm)"], distance_df["C2 Uncertainty (nm)"], color='red', alpha=0.5, label="C2 Uncertainty")
plt.xlabel("End-to-End Distance (nm)")
plt.ylabel("Uncertainty (nm)")
plt.title("End-to-End Distance vs. Uncertainty (C1 & C2)")
plt.legend()
plt.grid(True)
correlation_c1, p_value_c1 = pearsonr(distance_df["Distance (nm)"], distance_df["C1 Uncertainty (nm)"])
correlation_c2, p_value_c2 = pearsonr(distance_df["Distance (nm)"], distance_df["C2 Uncertainty (nm)"])
plt.show()
print(f"C1 Pearson Correlation: {correlation_c1:.4f}, P-Value: {p_value_c1:.4f}")
print(f"C2 Pearson Correlation: {correlation_c2:.4f}, P-Value: {p_value_c2:.4f}")

# Scatter Plot: Distance from Center vs. Dipole Length
roi_size, pixel_size_nm = 512, 160
image_center_x, image_center_y = (roi_size / 2) * pixel_size_nm, (roi_size / 2) * pixel_size_nm
dipole_center_x = (distance_df["C1 X (nm)"] + distance_df["C2 X (nm)"]) / 2
dipole_center_y = (distance_df["C1 Y (nm)"] + distance_df["C2 Y (nm)"]) / 2
distance_from_center = np.sqrt((dipole_center_x - image_center_x) ** 2 + (dipole_center_y - image_center_y) ** 2)

plt.figure(figsize=(8, 6))
plt.scatter(distance_from_center, distance_df["Distance (nm)"], color='blue', alpha=0.5)
plt.xlabel("Distance from Image Center (nm)")
plt.ylabel("End-to-End Distance (nm)")
plt.title("Chromatic Aberration Analysis: Distance from Center vs. Dipole Length")
plt.grid(True)
plt.show()


##########################################################################################
#########################  DIPOLE WITH DIRECTIONS  ########################################
##########################################################################################
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load the saved Excel file (with updated angles)
file_name = "End-to-End Distance_PUNCTA_Filtered.xlsx"
distance_df = pd.read_excel(file_name)

# Extract C1 and C2 coordinates
x_c1, y_c1 = distance_df["C1 X (nm)"], distance_df["C1 Y (nm)"]
x_c2, y_c2 = distance_df["C2 X (nm)"], distance_df["C2 Y (nm)"]

# Create scatter plot with symbolic arrows indicating directionality
plt.figure(figsize=(10, 8))
plt.scatter(x_c1, y_c1, color='green', alpha=0.7, label="TIRF 560")
plt.scatter(x_c2, y_c2, color='red', alpha=0.7, label="TIRF 647")

# Draw symbolic arrows from C1 to C2 and add annotations (with both Φ and θ)
for i in range(len(distance_df)):
    # Draw the dashed arrow line
    plt.plot([x_c1[i], x_c2[i]], [y_c1[i], y_c2[i]], 'k--', alpha=0.7)  # Dashed line for trajectory
    
    # Draw the arrow text as in the original code
    plt.text((x_c1[i] + x_c2[i]) / 2, (y_c1[i] + y_c2[i]) / 2, "--->", fontsize=12, color='black', 
             ha='center', va='center', fontweight='bold',
             rotation=np.degrees(np.arctan2(y_c2[i] - y_c1[i], x_c2[i] - x_c1[i])))
    
    # Calculate additional annotation: midpoint, distance, PHI and THETA angles
    mid_x = (x_c1[i] + x_c2[i]) / 2
    mid_y = (y_c1[i] + y_c2[i]) / 2
    phi_angle = distance_df["Φ (degrees)"].iloc[i]
    theta_angle = distance_df["θ (degrees)"].iloc[i]
    distance_nm = distance_df["Distance (nm)"].iloc[i]
    
    # Create annotation text with distance, Φ and θ angle values
    annotation_text = f"{distance_nm:.1f} nm, Φ={phi_angle:.1f}°, θ={theta_angle:.1f}°"
    
    # Place the annotation slightly offset from the midpoint (adjust offset as needed)
    plt.text(mid_x, mid_y + 10, annotation_text, fontsize=8, color='black', 
             ha='center', va='center', fontweight='bold')

# Labels, title, and legend
plt.xlabel("X Position (nm)")
plt.ylabel("Y Position (nm)")
plt.title("Scatter Plot with Symbolic Directional Arrows: C1 to C2")
plt.legend()
plt.grid(True)

# Show the plot
plt.show()
