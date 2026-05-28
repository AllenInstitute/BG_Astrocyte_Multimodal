import os
from pathlib import Path
import glob
import json
import numpy as np
import anndata as ad
import pandas as pd
import matplotlib.pyplot as plt
import sys

sys.path.append("/code")

# ── Subgroup params ────────────────────────────────────────────────────────────

subgroup_cluster_dict = {
    'GM STR': ['h-7', 'h-14', 'h-31'],
    'GM exSTR': [
        'h-145',
        'h-146',
        'h-149',
        'h-150',
        'h-152',
        'h-159',
        'h-160',
        'h-450',
        'h-452'
    ],
    'WM': ['h-28', 'h-227', 'h-230', 'h-231', 'h-472', 'h-479', 'h-485'],
    'Mixed': ['h-143', 'h-232']
}

# replace h- with Human-
subgroup_cluster_dict = {k: [c.replace('h-', 'Human-') for c in v] for k, v in subgroup_cluster_dict.items()}

# invert the dictionary to get cluster to subgroup mapping
cluster_subgroup_dict = {c: k for k, v in subgroup_cluster_dict.items() for c in v}

# ── Load data ──────────────────────────────────────────────────────────────────

adata_spatial = ad.read_h5ad('/data/cirro_files_hmba_human_20250922/hmba_release_format_human_cirro_h5ad_20250922.h5ad')

adata_spatial_astro = adata_spatial[adata_spatial.obs['Group'] == 'Astrocyte', :].copy()

human_mmc_output_dir = "/data/mapping_output/mapping_output"
extended_results_dir = "/data/mapping_output/extended_results"

h5ad_files = glob.glob(os.path.join(human_mmc_output_dir, "**/*.h5ad"), recursive=True)
dropped_cols_files = glob.glob(os.path.join(human_mmc_output_dir, "**/*dropped_cols.csv"), recursive=True)
extended_json_files = glob.glob(os.path.join(extended_results_dir, "**/*.json"), recursive=True)

# ── Helper functions ───────────────────────────────────────────────────────────

def extract_runner_ups(level_dict, n=5):
    """Extract up to n runner-ups from a taxonomy level dict."""
    runner_up_assignments = level_dict.get('runner_up_assignment', [])[:n]
    runner_up_correlations = level_dict.get('runner_up_correlation', [])[:n]
    runner_up_probabilities = level_dict.get('runner_up_probability', [])[:n]

    result = {}
    for i in range(n):
        result[f'runner_up_{i+1}'] = runner_up_assignments[i] if i < len(runner_up_assignments) else None
        result[f'runner_up_{i+1}_correlation'] = runner_up_correlations[i] if i < len(runner_up_correlations) else None
        result[f'runner_up_{i+1}_probability'] = runner_up_probabilities[i] if i < len(runner_up_probabilities) else None
    return result

def parse_cell_result(cell, specimen_name):
    """Parse a single cell result dict into a flat record."""
    row = {
        'cell_id': cell['cell_id'],
        'specimen_name': specimen_name,
    }

    levels = ['Neighborhood', 'Class', 'Subclass', 'Group', 'cluster_id']
    for level in levels:
        level_dict = cell.get(level, {})
        row[f'{level}'] = level_dict.get('assignment')
        row[f'{level}_bootstrapping_probability'] = level_dict.get('bootstrapping_probability')
        row[f'{level}_avg_correlation'] = level_dict.get('avg_correlation')
        # Add runner-ups only for cluster_id and Subclass levels (most granular)
        runner_ups = extract_runner_ups(level_dict, n=5)
        for k, v in runner_ups.items():
            row[f'{level}_{k}'] = v

    return row

# ── Parse extended JSON results ────────────────────────────────────────────────

rows = []
for json_file in extended_json_files:
    # Extract specimen name from filename
    specimen_name = os.path.basename(json_file).replace('_mmc_extended_results.json', '').replace('.json', '')

    with open(json_file, 'r') as f:
        extended_results = json.load(f)

    for cell in extended_results['results']:
        subclass_assignment = cell.get('Subclass', {}).get('assignment', '')
        if subclass_assignment == 'Astrocyte':
            row = parse_cell_result(cell, specimen_name)
            rows.append(row)

astrocyte_runner_up_df = pd.DataFrame(rows).set_index('cell_id')

# ── Merge extended results into AnnData ────────────────────────────────────────

adata_spatial_astro_extended = adata_spatial_astro.copy()

adata_spatial_astro_extended.obs = adata_spatial_astro_extended.obs.merge(
    astrocyte_runner_up_df,
    left_index=True,
    right_index=True,
    how='left',
    suffixes=('', '_extended')
)

# ── Fix None values and write out ─────────────────────────────────────────────

for col in adata_spatial_astro_extended.obs.columns:
    if adata_spatial_astro_extended.obs[col].dtype == object:
        adata_spatial_astro_extended.obs[col] = (
            adata_spatial_astro_extended.obs[col].fillna("None").astype(str)
        )

for cluster in adata_spatial_astro_extended.obs['cluster_id_extended'].unique():
    if cluster not in cluster_subgroup_dict:
        cluster_subgroup_dict[cluster] = 'Non-Astrocyte'

all_subGroups_ordered = ['GM STR', 'GM exSTR', 'WM', 'Mixed', 'Non-Astrocyte']

adata_spatial_astro_extended.obs['SubGroup_extended'] = adata_spatial_astro_extended.obs['cluster_id_extended'].map(cluster_subgroup_dict)
adata_spatial_astro_extended.obs['SubGroup_extended'] = adata_spatial_astro_extended.obs['SubGroup_extended'].astype('category')
adata_spatial_astro_extended.obs['SubGroup_extended'] = adata_spatial_astro_extended.obs['SubGroup_extended'].cat.set_categories(all_subGroups_ordered, ordered=True)
for runnerup in range(1, 6):
    adata_spatial_astro_extended.obs[f'SubGroup_runner_up_{runnerup}'] = adata_spatial_astro_extended.obs[f'cluster_id_runner_up_{runnerup}'].map(cluster_subgroup_dict)
    adata_spatial_astro_extended.obs[f'SubGroup_runner_up_{runnerup}'] = adata_spatial_astro_extended.obs[f'SubGroup_runner_up_{runnerup}'].fillna('Non-Astrocyte')
    adata_spatial_astro_extended.obs[f'SubGroup_runner_up_{runnerup}'] = adata_spatial_astro_extended.obs[f'SubGroup_runner_up_{runnerup}'].astype('category')
    adata_spatial_astro_extended.obs[f'SubGroup_runner_up_{runnerup}'] = adata_spatial_astro_extended.obs[f'SubGroup_runner_up_{runnerup}'].cat.set_categories(all_subGroups_ordered, ordered=True)

adata_spatial_astro_extended.write_h5ad('/scratch/astrocytes/human_spatial_astrocytes_extended_mmc_results_test.h5ad')