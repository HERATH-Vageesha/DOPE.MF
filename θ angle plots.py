import pandas as pd
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d

# === Step 1: Load the Excel file ===
excel_path = "Trackwise_θ_degrees_vs_Frame.xlsx"  # Excel file in same folder
theta_df = pd.read_excel(excel_path)

# === Step 2: Prepare frame and track columns ===
frame_numbers = theta_df['C1 Frame']
track_columns = theta_df.columns[1:]  # Skip 'C1 Frame' column

# === Step 3: Compute variance for each track ===
track_variances = theta_df[track_columns].var().sort_values()

# === Step 4: Select 10 most representative tracks (median variance) ===
midpoint = len(track_variances) // 2
representative_tracks = track_variances.iloc[midpoint - 5: midpoint + 5].index.tolist()

# === Step 5: Plot representative tracks with smoothing ===
fig, ax = plt.subplots(figsize=(12, 6))

for track in representative_tracks:
    y_values = theta_df[track]

    # Smooth with Gaussian filter after filling missing values
    smoothed = gaussian_filter1d(
        y_values.fillna(method='ffill').fillna(method='bfill'),
        sigma=1
    )

    # Plot with small circle markers and solid lines
    ax.plot(frame_numbers, smoothed, marker='o', markersize=3, linestyle='-',
            linewidth=1.0, alpha=0.9, label=track)

# === Step 6: Format plot ===
ax.set_title("Representative Theta Angle (θ)", fontsize=14)
ax.set_xlabel("C1 Frame", fontsize=12)
ax.set_ylabel("θ (degrees)", fontsize=12)
ax.set_ylim(40, 90)  # crop Y-axis to 40–90°
ax.grid(True, linestyle='--', alpha=0.5)

# Add legend to identify tracks
ax.legend(loc='center left', bbox_to_anchor=(1.02, 0.5), title="Tracks")

plt.tight_layout()
plt.show()
