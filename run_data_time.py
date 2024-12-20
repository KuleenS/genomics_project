import os
import time
import subprocess
import csv
import matplotlib.pyplot as plt
import pandas as pd
from data_gen import generate_test_strings, generate_genomic_strings

# Define paths and parameters
base_length = 500  # Starting length of the test strings
num_changes = 10  # Number of insertions/deletions/substitutions
num_tests = 6  # Number of test string pairs to generate - 1
num_pairs = 5 # Number of points to average
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
    {
        "name": "Four Russians",
        "command": lambda input_file: f"./build/four_russians {input_file} {3}"
    }
]

# Run each method on each input file and record runtime
results = []
for method in alignment_methods:
    for input_file in input_files:
        start_time = time.time()
        command = method["command"](input_file)
        try:
            subprocess.run(command, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error running {method['name']} on {input_file}: {e}")
            continue
        runtime = time.time() - start_time
        results.append({
            "method": method["name"],
            "input_file": input_file,
            "runtime": runtime
        })
        print(f"{method['name']} on {input_file} completed in {runtime:.4f} seconds")

# Save results to CSV
output_csv = "runtime_results_not_parallel.csv" # or runtime_results_parallel.csv
with open(output_csv, 'w', newline='') as csvfile:
    fieldnames = ["method", "input_file", "runtime"]
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()
    for result in results:
        writer.writerow(result)

print(f"Runtime results saved to {output_csv}")
