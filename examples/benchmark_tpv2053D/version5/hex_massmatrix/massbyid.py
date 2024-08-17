# Python script to process the input file
input_file = 'massbyidpreprocess.txt'
output_file = 'massbyidpostprocess.txt'

# Initialize a list to store the non-zero values from the second column
non_zero_values = []

# Open the input file for reading
with open(input_file, 'r') as file:
    lines = file.readlines()
    start_index = 0
    end_index = 0

    # Find the start and end indices of the relevant data section
    for i, line in enumerate(lines):
        if '#	Value' in line:
            start_index = i + 1
        if 'Solve Converged!' in line:
            end_index = i
            break

    # Extract the lines and filter nonzero values from the second column
    for i in range(start_index, end_index):
        parts = lines[i].split()
        if len(parts) == 2:  # Ensure there are exactly two columns
            value = int(parts[1])
            if value != 0:  # Keep only non-zero values
                non_zero_values.append(str(value))

# Save the non-zero values to the output file
with open(output_file, 'w') as file:
    file.write('\n'.join(non_zero_values))

print(f'Filtered non-zero values from the second column saved to {output_file}')