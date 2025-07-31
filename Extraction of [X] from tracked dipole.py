import pandas as pd

# === ✅ SET THIS LINE ONLY TO EXTRACT A DIFFERENT COLUMN ===
value_column = "θ (degrees)"  # <-- Change this to e.g., "Φ (degrees)", "θ (degrees)", etc.

# Load the Excel file (assumes it is in the same folder as the script)
excel_path = "BETA Sample_timelapse_analysis.xlsx"
xls = pd.ExcelFile(excel_path)

# Get all sheet names and exclude the "Master" sheet
sheet_names = xls.sheet_names
non_master_sheets = [sheet for sheet in sheet_names if sheet.lower() != "master"]

# Create a DataFrame with frame numbers from 1 to 50
frames_df = pd.DataFrame({'C1 Frame': range(1, 51)})

# Initialize result DataFrame
result_df = frames_df.copy()

# Loop through each track sheet
for sheet in non_master_sheets:
    # Read the worksheet
    df = xls.parse(sheet)

    # Check if the desired column exists in the current sheet
    if "C1 Frame" in df.columns and value_column in df.columns:
        # Extract relevant columns and drop rows with missing values
        track_values = df[['C1 Frame', value_column]].dropna()

        # Rename the value column to the sheet name (Track name)
        track_values = track_values.rename(columns={value_column: sheet})

        # Merge with the result DataFrame on C1 Frame
        result_df = result_df.merge(track_values, on='C1 Frame', how='left')

# Save the final DataFrame to an Excel file
output_filename = f"Trackwise_{value_column.replace(' ', '_').replace('(', '').replace(')', '')}_vs_Frame.xlsx"
result_df.to_excel(output_filename, index=False)
