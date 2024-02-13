def calculate_average(filepath):
    with open(filepath, 'r') as file:
        numbers = [int(line.strip()) for line in file if line.strip()]
    
    if numbers:
        average = sum(numbers) / len(numbers)
        return average
    else:
        return None

# Filepath
filepath = '/home/projects/MAAG/LAI_sim/read_simulation/pestis.insize'

# Calculate and print the average
avg = calculate_average(filepath)
if avg is not None:
    print(f"The average number is: {avg}")
else:
    print("No numbers found in the file.")


