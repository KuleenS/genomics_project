import os
import time
import subprocess
import csv
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# Load results from CSV
data_parallel = pd.read_csv('runtime_results_parallel.csv', skiprows=1, names=['test_type','test_length_file', 'time'])
data_not_parallel = pd.read_csv('runtime_results_not_parallel.csv', skiprows=1, names=['test_type','test_length_file', 'time'])


methods = data_parallel['test_type'].unique()
print(methods)

# Plot runtimes for each method
for method in methods:
    fig = plt.figure(figsize=(10, 6))
    method_data_parallel = data_parallel[data_parallel['test_type'] == method]
    method_data_not_parallel = data_not_parallel[data_not_parallel['test_type'] == method]
    # Extract input lengths from file names
    method_data_parallel['test_length_file'] = method_data_parallel['test_length_file'].str.extract(r'(\d+)').astype(int)
    method_data_not_parallel['test_length_file'] = method_data_not_parallel['test_length_file'].str.extract(r'(\d+)').astype(int)
    
    # Sort by input length
    method_data_parallel = method_data_parallel.sort_values(by='test_length_file')
    method_data_not_parallel = method_data_not_parallel.sort_values(by='test_length_file')
    plt.plot(method_data_parallel['test_length_file'], method_data_parallel['time'], 'ro', label='parallel')
    plt.plot(method_data_not_parallel['test_length_file'], method_data_not_parallel['time'], 'bo', label='not parallel')
    plt.xlabel("Number Bases")
    plt.ylabel("Runtime (s)")
    plt.xticks(np.linspace(method_data_parallel['test_length_file'].min(), method_data_parallel['test_length_file'].max(), 5).astype(int))
    plt.title(f"Runtime Comparison of {method}")
    plt.legend()
    plt.savefig(f'runtime_plot_{method}.png', format='png', dpi=300)
