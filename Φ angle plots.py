import pandas as pd
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d

# === Step 1: Load the Excel file ===
excel_path = "Trackwise_Φ_degrees_vs_Frame.xlsx"  # Excel file in same folder
phi_df = pd.read_excel(excel_path)
# Extract frame and track columns
frame_numbers_phi = phi_df['C1 Frame']
track_columns_phi = phi_df.columns[1:]

# Compute variance and select 10 representative tracks (middle variance)
track_variances_phi = phi_df[track_columns_phi].var().sort_values()
midpoint_phi = len(track_variances_phi) // 2
representative_tracks_phi = track_variances_phi.iloc[midpoint_phi - 5: midpoint_phi + 5].index.tolist()

# Plot with Gaussian smoothing
fig, ax = plt.subplots(figsize=(12, 6))

for track in representative_tracks_phi:
    y_values = phi_df[track]
    smoothed = gaussian_filter1d(
        y_values.fillna(method='ffill').fillna(method='bfill'),
        sigma=1
    )
    ax.plot(frame_numbers_phi, smoothed, marker='o', markersize=3, linestyle='-', linewidth=1.0, alpha=0.9, label=track)

ax.set_title("Representative Phi Angle (Φ) Dynamics", fontsize=14)
ax.set_xlabel("C1 Frame", fontsize=12)
ax.set_ylabel("Φ (degrees)", fontsize=12)
ax.set_ylim(0, 360)
ax.grid(True, linestyle='--', alpha=0.5)
ax.legend(loc='center left', bbox_to_anchor=(1.02, 0.5), title="Tracks")

plt.tight_layout()
plt.show()