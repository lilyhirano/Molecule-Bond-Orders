import h5py
import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
import re


label_map = {
    "H2 (0-1)": "H2 : H-H",
    "HF (0-1)": "HF : H-F",
    "H2O (0-1)": "H2O : O-H (1)",
    "H2O (0-2)": "H2O : O-H (2)",
    "HO (0-1)": "HO : O-H (Radical)",
    "HCN (0-1)": "HCN : C#N",
    "HCN (0-2)": "HCN : C-H",
    "C2H4 (0-1)": "C2H4 : C=C",
    "C2H4 (0-2)": "C2H4 : C-H (1)",
    "C2H4 (0-3)": "C2H4 : C-H (2)",
    "C2H4 (1-4)": "C2H4 : C-H (3)",
    "C2H4 (1-5)": "C2H4 : C-H (4)",
}

all_mols = ["H2", "HF", "H2O", "HO", "HCN", "C2H4"]
methods = ["wiberg", "mayer", "mulliken"]
raw_results = []

for mol in all_mols:
    h5_path = f"project_output/final_proj_hdf5/{mol}.hdf5"
    csv_path = f"tests/test_checks/{mol}.csv"
    
    if not os.path.exists(h5_path) or not os.path.exists(csv_path): continue
    adj_matrix = pd.read_csv(csv_path, header=None).values
    
    with h5py.File(h5_path, 'r') as f:
        n = adj_matrix.shape[0]
        for i in range(n):
            for j in range(i + 1, n):
                if adj_matrix[i, j] > 0.5:
                    vals = [f[f"{m}_bond_order_matrix"][i, j] for m in methods]
                    raw_key = f"{mol} ({i}-{j})"
                    full_label = label_map.get(raw_key, raw_key)
                    
                    raw_results.append({
                        "FullLabel": full_label,
                        "Mean": np.mean(vals),
                        "StdDev": np.std(vals)
                    })

df = pd.DataFrame(raw_results)

def clean_bond_label(label):
    return re.sub(r'\s*\(\d+\)', '', label)

df['CleanLabel'] = df['FullLabel'].apply(clean_bond_label)

df_unique = df.drop_duplicates(subset=['CleanLabel']).copy()

plt.figure(figsize=(12, 7))
plt.errorbar(df_unique['CleanLabel'], df_unique['Mean'], yerr=df_unique['StdDev'], 
             fmt='o', color='#3498db', capsize=6, markersize=8, elinewidth=2, markeredgecolor='black')

plt.xticks(rotation=45, ha='right', fontsize=10, fontweight='bold')
plt.ylabel('Bond Order (Mean +/- StdDev)', fontsize=12, fontweight='bold')
plt.title('Bond Order Analysis Across All Test Molecules', fontsize=14, fontweight='bold', pad=20)
plt.grid(axis='y', linestyle='--', alpha=0.7)
plt.tight_layout()

plt.savefig('project_output/bonds_comparison_plot.png')
print("Combined plot saved to project_output/bonds_comparison_plot.png")
