import pandas as pd
import numpy as np

URL = "https://raw.githubusercontent.com/statOmics/pda/data/quantification/cancer/proteinGroups.txt"

print("Downloading dataset...")
print("Source: StatOmics/pda — Squamous cell carcinoma vs Control (LC-MS/MS)")
df_raw = pd.read_csv(URL, sep='\t', low_memory=False)
print(f"Raw data: {len(df_raw)} proteins, {len(df_raw.columns)} columns")

intensity_cols = [c for c in df_raw.columns if c.startswith('Intensity ')]
sample_names   = [c.replace('Intensity ', '') for c in intensity_cols]
print(f"Samples ({len(sample_names)}): {sample_names}")

groups  = ['Disease' if s.startswith('OR') else 'Control' for s in sample_names]
batches = ['Batch1']*6 + ['Batch2']*6 + ['Batch3']*6

metadata = pd.DataFrame({
    'sample': sample_names,
    'group' : groups,
    'batch' : batches
}).set_index('sample')

print(f"\nGroup distribution:")
print(metadata['group'].value_counts().to_string())
print(f"\nBatch distribution:")
print(metadata['batch'].value_counts().to_string())

protein_matrix = df_raw[intensity_cols].copy()
protein_matrix.columns = sample_names
protein_matrix.index   = df_raw['Majority protein IDs'].fillna(df_raw['Protein IDs'])
protein_matrix = protein_matrix.replace(0, np.nan)
protein_matrix = np.log2(protein_matrix)
protein_matrix = protein_matrix.replace([np.inf, -np.inf], np.nan)

if 'Potential contaminant' in df_raw.columns:
    contam = df_raw['Potential contaminant'].fillna('').str.contains(r'\+')
    protein_matrix = protein_matrix[~contam.values]

print(f"\nFinal matrix: {protein_matrix.shape[0]} proteins x {protein_matrix.shape[1]} samples")
print(f"Overall missing value rate: {protein_matrix.isna().mean().mean()*100:.1f}%")

protein_matrix.to_csv('data/intensity_matrix.csv')
metadata.to_csv('data/metadata.csv')

print("\nStep 1 COMPLETE")
print("  data/intensity_matrix.csv")
print("  data/metadata.csv")
