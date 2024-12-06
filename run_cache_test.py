import os
import time
import re
import subprocess
import csv
import matplotlib.pyplot as plt
import pandas as pd
from data_gen import generate_test_strings, generate_genomic_strings

# Define paths and parameters
base_length = 10  # Starting length of the test strings
num_changes = 8  # Number of insertions/deletions/substitutions
num_tests = 6  # Number of test string pairs to generate - 1
num_pairs = 1
penalty_file = "penalty.csv"  # Path to penalty function CSV

lengths = [(base_length * (2 **i)) for i in range(0, num_tests)]
changes = [((num_changes * 2 * i) + 10) for i in range(0, num_tests)]

input_files = []
for i, length in enumerate(lengths):
    for file_number in range(1, num_pairs + 1):
        dna_string = generate_test_strings(length, changes[i], 1)
        input_file = f"test_length_{length}_file_{file_number}.txt"
        with open(input_file, 'w') as file:
            file.write(f"{dna_string[0][0]}\n{dna_string[0][1]}")
        input_files.append(input_file)


# Define methods to run
alignment_methods = [
    {
        "name": "Hirschberg",
        "command": lambda input_file: f"./build/hirschberg {input_file} {penalty_file}"
    },
    {
        "name": "Affine_0",
        "command": lambda input_file, mode="affine": f"./build/alignment {input_file} {mode} 0 {penalty_file} 5 20"
    },
    {
        "name": "Affine_1",
        "command": lambda input_file, mode="affine": f"./build/alignment {input_file} {mode} 1 {penalty_file} 5 20"
    },
    {
        "name": "Affine_2",
        "command": lambda input_file, mode="affine": f"./build/alignment {input_file} {mode} 2 {penalty_file} 5 20"
    },
    {
        "name": "Global_0",
        "command": lambda input_file, mode="global": f"./build/alignment {input_file} {mode} 0 {penalty_file}"
    },
    {
        "name": "Global_1",
        "command": lambda input_file, mode="global": f"./build/alignment {input_file} {mode} 1 {penalty_file}"
    },
    {
        "name": "Global_2",
        "command": lambda input_file, mode="global": f"./build/alignment {input_file} {mode} 2 {penalty_file}"
    },
    {
        "name": "Local_0",
        "command": lambda input_file, mode="local": f"./build/alignment {input_file} {mode} 0 {penalty_file}"
    },
    {
        "name": "Local_1",
        "command": lambda input_file, mode="local": f"./build/alignment {input_file} {mode} 1 {penalty_file}"
    },
    {
        "name": "Local_2",
        "command": lambda input_file, mode="local": f"./build/alignment {input_file} {mode} 2 {penalty_file}"
    },
    # {
    #     "name": "Four Russians",
    #     "command": lambda input_file: f"./build/four_russians {input_file} {3}"
    # }
]

cachegrind_dir = "cachegrind_outputs"
os.makedirs(cachegrind_dir, exist_ok=True)

# Run each method on each input file and record cache misses
results = []
for method in alignment_methods:
    for input_file in input_files:
        # Call the lambda function to get the actual command string
        command_string = method["command"](input_file)
        
        # Specify the Cachegrind output file
        output_file = os.path.join(cachegrind_dir, f"cachegrind.out.{method['name']}")
        
        # Full Cachegrind command
        command = f"valgrind --tool=cachegrind --cache-sim=yes --cachegrind-out-file={output_file} {command_string}"
        
        try:
            # Run the command
            process = subprocess.run(command, shell=True, check=True, stderr=subprocess.PIPE, text=True)
            valgrind_output = process.stderr
        except subprocess.CalledProcessError as e:
            print(f"Error running {method['name']} on {input_file}: {e}")
            continue

        # Parse cachegrind summary using cg_annotate
        annotate_command = f"cg_annotate {output_file}"
        annotate_output = subprocess.check_output(annotate_command, shell=True).decode()

        # print(annotate_output)

        L1_misses = 0
        for line in valgrind_output.splitlines():
            match = re.search(r"D1  misses:\s+([\d,]+)", line)
            if match:
                L1_misses = int(match.group(1).replace(",", ""))
        print("L1 Misses: ",L1_misses)
        break
        
        results.append({
            "method": method["name"],
            "input_file": input_file,
            "L1_misses": L1_misses
        })
        print(f"{method['name']} on {input_file} completed")

# Save results to CSV
csv_file = "cachegrind_comparison_parallel.csv"
with open(csv_file, 'w', newline='') as csvfile:
    fieldnames = ["method", "input_file", "L1_misses"]
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()
    for result in results:
        writer.writerow(result)

print(f"Cachegrind comparison saved to {csv_file}")
