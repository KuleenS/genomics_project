import csv
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import numpy as np
import matplotlib.colors as mcolors
from scipy.optimize import curve_fit


mpl.rcParams.update({'font.size': 14})

def fitting_equation(x,a,b,c):
    y = a*(x**2) + b*x + c
    return y

# Load results from CSV
data_parallel = pd.read_csv('runtime_results_parallel.csv', skiprows=1, names=['test_type','test_length_file', 'time'])
data_not_parallel = pd.read_csv('runtime_results_not_parallel.csv', skiprows=1, names=['test_type','test_length_file', 'time'])

method_types = ['Hirschberg','Local', 'Global', 'Affine','Four Russians'] # 
colors = list(mcolors.TABLEAU_COLORS.values())

names = {
    0: 'Vertical',
    1: 'Horizontal',
    2: 'Antidiagonal',
}

results = {}

for method in method_types:
    fig = plt.figure(figsize=(10, 6))
    if method != 'Hirschberg' and method != 'Four Russians':
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

            popt_parallel, pcov_parallel = curve_fit(fitting_equation, method_data_parallel['test_length_file'], method_data_parallel['time'], maxfev=20000)
            popt_nonparallel, pcov_nonparallel = curve_fit(fitting_equation, method_data_not_parallel['test_length_file'], method_data_not_parallel['time'],  maxfev=20000)
            
            results[f'{method}{mode}'] = (popt_parallel,popt_nonparallel)
            length = np.linspace(0, 16000, 100)
            parallel_fit = fitting_equation(length, *popt_parallel)
            nonparallel_fit = fitting_equation(length, *popt_nonparallel)

            plt.plot(method_data_parallel['test_length_file'], method_data_parallel['time'], 'o', color=color_base)
            plt.plot(length, parallel_fit, linestyle='-', color=color_base, label=f'{method} {names[mode]} parallel')

            plt.plot(method_data_not_parallel['test_length_file'], method_data_not_parallel['time'],'o', color=color_base)
            plt.plot(length, nonparallel_fit, linestyle='--', color=color_base, label=f'{method} {names[mode]} serial')

        plt.xlabel("Number Bases")
        plt.ylabel("Runtime (s)")
        plt.legend()
        plt.savefig(f'runtime_plot_{method}.png', format='png', dpi=300, bbox_inches = 'tight')
    elif method == 'Four Russians':
        color_base = colors[3 % len(colors)]
        fig = plt.figure(figsize=(10, 6))
        method_data_parallel = data_parallel[data_parallel['test_type'] == method]
        #method_data_not_parallel = data_not_parallel[data_not_parallel['test_type'] == method]
        # Extract input lengths from file names
        method_data_parallel['test_length_file'] = method_data_parallel['test_length_file'].str.extract(r'(\d+)').astype(int)
        #method_data_not_parallel['test_length_file'] = method_data_not_parallel['test_length_file'].str.extract(r'(\d+)').astype(int)
    
        # Sort by input length
        popt_parallel, pcov_parallel = curve_fit(fitting_equation, method_data_parallel['test_length_file'], method_data_parallel['time'], maxfev=20000)
        #popt_nonparallel, pcov_nonparallel = curve_fit(fitting_equation, method_data_not_parallel['test_length_file'], method_data_not_parallel['time'],  maxfev=20000)
            
        results[method] = (popt_parallel)

        length = np.linspace(0, 16000, 100)
        parallel_fit = fitting_equation(length, *popt_parallel)

        plt.plot(method_data_parallel['test_length_file'], method_data_parallel['time'], 'o', color=color_base)
        plt.plot(length, parallel_fit, linestyle='-', color=color_base, label=f'{method}')

        plt.xlabel("Number Bases")
        plt.ylabel("Runtime (s)")
        plt.legend()
        plt.savefig(f'runtime_plot_{method}.png', format='png', dpi=300, bbox_inches = 'tight')
    else:
        color_base = colors[3 % len(colors)]
        fig = plt.figure(figsize=(10, 6))
        method_data_parallel = data_parallel[data_parallel['test_type'] == method]
        method_data_not_parallel = data_not_parallel[data_not_parallel['test_type'] == method]
        # Extract input lengths from file names
        method_data_parallel['test_length_file'] = method_data_parallel['test_length_file'].str.extract(r'(\d+)').astype(int)
        method_data_not_parallel['test_length_file'] = method_data_not_parallel['test_length_file'].str.extract(r'(\d+)').astype(int)
    
        # Sort by input length
        popt_parallel, pcov_parallel = curve_fit(fitting_equation, method_data_parallel['test_length_file'], method_data_parallel['time'], maxfev=20000)
        popt_nonparallel, pcov_nonparallel = curve_fit(fitting_equation, method_data_not_parallel['test_length_file'], method_data_not_parallel['time'],  maxfev=20000)
            
        results[method] = (popt_parallel,popt_nonparallel)

        length = np.linspace(0, 16000, 100)
        parallel_fit = fitting_equation(length, *popt_parallel)
        nonparallel_fit = fitting_equation(length, *popt_nonparallel)

        plt.plot(method_data_parallel['test_length_file'], method_data_parallel['time'], 'o', color=color_base)
        plt.plot(length, parallel_fit, linestyle='-', color=color_base, label=f'{method} serial')

        plt.plot(method_data_not_parallel['test_length_file'], method_data_not_parallel['time'],'o', color=color_base)
        plt.plot(length, nonparallel_fit, linestyle='--', color=color_base, label=f'{method} parallel')

        plt.xlabel("Number Bases")
        plt.ylabel("Runtime (s)")
        plt.legend()
        plt.savefig(f'runtime_plot_{method}.png', format='png', dpi=300, bbox_inches = 'tight')
    


df = pd.DataFrame(results)

df.to_csv("output_results_fourrussians.csv", index=False)
    