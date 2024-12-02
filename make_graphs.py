import csv
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import matplotlib.colors as mcolors

# Load results from CSV
data_parallel = pd.read_csv('runtime_results_parallel.csv', skiprows=1, names=['test_type','test_length_file', 'time'])
data_not_parallel = pd.read_csv('runtime_results_not_parallel.csv', skiprows=1, names=['test_type','test_length_file', 'time'])

method_types = ['Hirschberg','Local', 'Global', 'Affine']
colors = list(mcolors.TABLEAU_COLORS.values())

for method in method_types:
    fig = plt.figure(figsize=(10, 6))
    if method != 'Hirschberg':
        unique_modes = {0, 1, 2}
        for idx, mode in enumerate(unique_modes):
            color_base = colors[idx % len(colors)]  # Rotate through colors if there are more modes than colors
        
            # Filter parallel and non-parallel data for the current mode
            method_data_parallel = data_parallel[data_parallel['test_type'] == f"{method}_{mode}"]
            method_data_not_parallel = data_not_parallel[data_not_parallel['test_type'] == f"{method}_{mode}"]

            # Extract input lengths from file names
            method_data_parallel['test_length_file'] = method_data_parallel['test_length_file'].str.extract(r'(\d+)').astype(int)
            method_data_not_parallel['test_length_file'] = method_data_not_parallel['test_length_file'].str.extract(r'(\d+)').astype(int)

            method_data_parallel = method_data_parallel.sort_values(by='test_length_file')
            method_data_not_parallel = method_data_not_parallel.sort_values(by='test_length_file')

            plt.plot(method_data_parallel['test_length_file'], method_data_parallel['time'], color=color_base, label=f'{method}_{mode} parallel', linestyle='-')
            plt.plot(method_data_not_parallel['test_length_file'], method_data_not_parallel['time'], color=color_base, label=f'{method}_{mode} serial', linestyle='--')
        plt.xlabel("Number Bases")
        plt.ylabel("Runtime (s)")
        plt.xticks(np.linspace(1000, 5000, 5).astype(int))
        plt.title(f"Runtime Comparison of {method}")
        plt.legend()
        plt.savefig(f'runtime_plot_{method}.png', format='png', dpi=300)
    else:
        color_base = colors[3 % len(colors)]
        fig = plt.figure(figsize=(10, 6))
        method_data_parallel = data_parallel[data_parallel['test_type'] == method]
        method_data_not_parallel = data_not_parallel[data_not_parallel['test_type'] == method]
        # Extract input lengths from file names
        method_data_parallel['test_length_file'] = method_data_parallel['test_length_file'].str.extract(r'(\d+)').astype(int)
        method_data_not_parallel['test_length_file'] = method_data_not_parallel['test_length_file'].str.extract(r'(\d+)').astype(int)
    
        # Sort by input length
        method_data_parallel = method_data_parallel.sort_values(by='test_length_file')
        method_data_not_parallel = method_data_not_parallel.sort_values(by='test_length_file')
        plt.plot(method_data_parallel['test_length_file'], method_data_parallel['time'], color=color_base, label=f'{method} parallel', linestyle='-')
        plt.plot(method_data_not_parallel['test_length_file'], method_data_not_parallel['time'], color=color_base, label=f'{method} serial', linestyle='--')
    
        plt.xlabel("Number Bases")
        plt.ylabel("Runtime (s)")
        plt.xticks(np.linspace(1000, 5000, 5).astype(int))
        plt.title(f"Runtime Comparison of {method}")
        plt.legend()
        plt.savefig(f'runtime_plot_{method}.png', format='png', dpi=300)
    