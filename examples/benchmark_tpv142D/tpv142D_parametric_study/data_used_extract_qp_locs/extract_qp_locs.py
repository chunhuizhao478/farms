import pandas as pd

## ------------main fault----------------##

# Read the CSV file
input_filename = "tpv142D_tria_csv_main_fault_0010.csv"  # Replace with your input CSV file name
df = pd.read_csv(input_filename)

# Extract the columns "x", "y", "z"
extracted_df = df[['x', 'y', 'z']]

# Save the extracted columns to a new CSV file
output_filename = "main_fault_qp_coords.csv"  # Replace with your desired output file name
extracted_df.to_csv(output_filename, index=False)

print(f"Extracted columns saved to {output_filename}")

## ------------branch fault----------------##

# Read the CSV file
input_filename = "tpv142D_tria_csv_branch_fault_0010.csv"  # Replace with your input CSV file name
df = pd.read_csv(input_filename)

# Extract the columns "x", "y", "z"
extracted_df = df[['x', 'y', 'z']]

# Save the extracted columns to a new CSV file
output_filename = "branch_fault_qp_coords.csv"  # Replace with your desired output file name
extracted_df.to_csv(output_filename, index=False)

print(f"Extracted columns saved to {output_filename}")