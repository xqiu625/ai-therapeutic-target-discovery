#!/usr/bin/env python
"""
01 - Data Preparation Script
Sepsis Target Discovery Pipeline

This script loads the extracted Seurat data and prepares it for downstream analysis.
"""

import scanpy as sc
import pandas as pd
import numpy as np
import scipy.io
import scipy.sparse as sp
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
import seaborn as sns
import gzip
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Settings
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=100, facecolor='white', frameon=False)

# Paths - Import from config or use defaults
try:
    from config import PROCESSED_DIR, FIGURES_DIR
except ImportError:
    BASE_DIR = Path(__file__).parent.parent.resolve()
    PROCESSED_DIR = BASE_DIR / 'data' / 'processed'
    FIGURES_DIR = BASE_DIR / 'figures'

PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

print("="*60)
print("01 - DATA PREPARATION")
print("="*60)
print(f"Scanpy version: {sc.__version__}")
print(f"Data directory: {PROCESSED_DIR}")

# ============================================================
# 1. Load Data
# ============================================================
print("\n" + "="*60)
print("1. LOADING DATA")
print("="*60)

# Load count matrix
print("\nLoading count matrix...")
counts_file = PROCESSED_DIR / 'counts.mtx.gz'
with gzip.open(counts_file, 'rb') as f:
    counts = scipy.io.mmread(f)
counts = counts.T.tocsr()  # Transpose: genes x cells -> cells x genes
print(f"Count matrix shape: {counts.shape} (cells x genes)")

# Load gene names
print("Loading gene names...")
genes = pd.read_csv(PROCESSED_DIR / 'genes.csv')
print(f"Loaded {len(genes)} genes")

# Load cell metadata
print("Loading cell metadata...")
metadata = pd.read_csv(PROCESSED_DIR / 'metadata.csv')
print(f"Loaded metadata for {len(metadata)} cells")
print(f"Columns: {metadata.columns.tolist()}")

# Create AnnData
print("\nCreating AnnData object...")
adata = sc.AnnData(
    X=counts,
    obs=metadata.set_index('barcode') if 'barcode' in metadata.columns else metadata,
    var=genes.set_index('gene_name') if 'gene_name' in genes.columns else genes
)
adata.var_names = adata.var_names.astype(str)
adata.var_names_make_unique()
print(f"AnnData created: {adata.n_obs} cells x {adata.n_vars} genes")

# Load embeddings if available
print("\nLoading embeddings...")
pca_file = PROCESSED_DIR / 'pca.csv'
if pca_file.exists():
    pca = pd.read_csv(pca_file, index_col=0)
    adata.obsm['X_pca'] = pca.values
    print(f"Loaded PCA: {pca.shape}")

umap_file = PROCESSED_DIR / 'umap.csv'
if umap_file.exists():
    umap = pd.read_csv(umap_file, index_col=0)
    adata.obsm['X_umap'] = umap.values
    print(f"Loaded UMAP: {umap.shape}")

tsne_file = PROCESSED_DIR / 'tsne.csv'
if tsne_file.exists():
    tsne = pd.read_csv(tsne_file, index_col=0)
    adata.obsm['X_tsne'] = tsne.values
    print(f"Loaded tSNE: {tsne.shape}")

# ============================================================
# 2. Inspect Metadata
# ============================================================
print("\n" + "="*60)
print("2. METADATA INSPECTION")
print("="*60)

for col in adata.obs.columns:
    n_unique = adata.obs[col].nunique()
    print(f"\n{col}: {n_unique} unique values")
    if n_unique <= 10:
        print(adata.obs[col].value_counts().to_string())

# ============================================================
# 3. Standardize Metadata
# ============================================================
print("\n" + "="*60)
print("3. STANDARDIZING METADATA")
print("="*60)

# Map condition to outcome
# From metadata: condition has HC (healthy), and likely sepsis survivors/non-survivors
# time column has timepoints
# cell_type already exists

# Check what's in condition column
print("\nCondition column:")
print(adata.obs['condition'].value_counts())

# Create outcome column based on condition
# Assuming condition contains: HC (healthy control), S (survivor), NS (non-survivor)
# or similar patterns
condition_values = adata.obs['condition'].unique()
print(f"\nUnique conditions: {condition_values}")

# Create outcome mapping
outcome_map = {}
for val in condition_values:
    val_lower = str(val).lower()
    if 'hc' in val_lower or 'healthy' in val_lower or 'control' in val_lower:
        outcome_map[val] = 'Healthy'
    elif 'ns' in val_lower or 'non' in val_lower:
        outcome_map[val] = 'Non-survivor'
    elif 's' in val_lower or 'survivor' in val_lower:
        outcome_map[val] = 'Survivor'
    else:
        outcome_map[val] = val  # Keep original

print(f"Outcome mapping: {outcome_map}")
adata.obs['outcome'] = adata.obs['condition'].map(outcome_map)
print("\nOutcome distribution:")
print(adata.obs['outcome'].value_counts())

# Create sepsis/control condition
adata.obs['disease'] = adata.obs['outcome'].apply(
    lambda x: 'Control' if x == 'Healthy' else 'Sepsis'
)
print("\nDisease distribution:")
print(adata.obs['disease'].value_counts())

# Standardize timepoint if exists
if 'time' in adata.obs.columns:
    adata.obs['timepoint'] = adata.obs['time']
    print("\nTimepoint distribution:")
    print(adata.obs['timepoint'].value_counts())

# Use sample_id as patient_id
if 'sample_id' in adata.obs.columns:
    adata.obs['patient_id'] = adata.obs['sample_id']

# ============================================================
# 4. Quality Control
# ============================================================
print("\n" + "="*60)
print("4. QUALITY CONTROL")
print("="*60)

# Store raw counts
adata.layers['counts'] = adata.X.copy()
print("Stored raw counts in layers['counts']")

# Identify MT and ribo genes
adata.var['mt'] = adata.var_names.str.startswith('MT-')
adata.var['ribo'] = adata.var_names.str.startswith(('RPS', 'RPL'))
print(f"Mitochondrial genes: {adata.var['mt'].sum()}")
print(f"Ribosomal genes: {adata.var['ribo'].sum()}")

# Calculate QC metrics
sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=['mt', 'ribo'],
    percent_top=None,
    log1p=False,
    inplace=True
)

print(f"\nQC Summary:")
print(f"  n_genes median: {adata.obs['n_genes_by_counts'].median():.0f}")
print(f"  total_counts median: {adata.obs['total_counts'].median():.0f}")
print(f"  pct_mt median: {adata.obs['pct_counts_mt'].median():.1f}%")

# Plot QC
fig, axes = plt.subplots(1, 4, figsize=(16, 4))
sc.pl.violin(adata, ['n_genes_by_counts'], jitter=0.4, ax=axes[0], show=False)
sc.pl.violin(adata, ['total_counts'], jitter=0.4, ax=axes[1], show=False)
sc.pl.violin(adata, ['pct_counts_mt'], jitter=0.4, ax=axes[2], show=False)
sc.pl.violin(adata, ['pct_counts_ribo'], jitter=0.4, ax=axes[3], show=False)
plt.tight_layout()
plt.savefig(FIGURES_DIR / 'qc_metrics_violin.png', dpi=150, bbox_inches='tight')
plt.close()
print("Saved QC violin plot")

# ============================================================
# 5. Normalization
# ============================================================
print("\n" + "="*60)
print("5. NORMALIZATION")
print("="*60)

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.layers['normalized'] = adata.X.copy()
print("Normalization complete")

# ============================================================
# 6. Highly Variable Genes
# ============================================================
print("\n" + "="*60)
print("6. HIGHLY VARIABLE GENES")
print("="*60)

sc.pp.highly_variable_genes(
    adata,
    n_top_genes=3000,
    subset=False,
    flavor='seurat_v3',
    layer='counts'
)
print(f"Highly variable genes: {adata.var['highly_variable'].sum()}")

sc.pl.highly_variable_genes(adata, show=False)
plt.savefig(FIGURES_DIR / 'hvg_selection.png', dpi=150, bbox_inches='tight')
plt.close()

# ============================================================
# 7. Dimensionality Reduction
# ============================================================
print("\n" + "="*60)
print("7. DIMENSIONALITY REDUCTION")
print("="*60)

if 'X_pca' not in adata.obsm:
    print("Computing PCA...")
    adata_hvg = adata[:, adata.var['highly_variable']].copy()
    sc.pp.scale(adata_hvg, max_value=10)
    sc.tl.pca(adata_hvg, n_comps=50)
    adata.obsm['X_pca'] = adata_hvg.obsm['X_pca']
    adata.uns['pca'] = adata_hvg.uns['pca']
else:
    print(f"Using existing PCA: {adata.obsm['X_pca'].shape}")

if 'X_umap' not in adata.obsm:
    print("Computing UMAP...")
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
    sc.tl.umap(adata)
else:
    print(f"Using existing UMAP: {adata.obsm['X_umap'].shape}")
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)

# Plot UMAP
fig, axes = plt.subplots(2, 2, figsize=(12, 12))
sc.pl.umap(adata, color='outcome', ax=axes[0, 0], show=False, title='Outcome')
sc.pl.umap(adata, color='cell_type', ax=axes[0, 1], show=False, title='Cell Type')
if 'timepoint' in adata.obs.columns:
    sc.pl.umap(adata, color='timepoint', ax=axes[1, 0], show=False, title='Timepoint')
sc.pl.umap(adata, color='patient_id', ax=axes[1, 1], show=False, title='Patient')
plt.tight_layout()
plt.savefig(FIGURES_DIR / 'umap_metadata.png', dpi=150, bbox_inches='tight')
plt.close()
print("Saved UMAP plot")

# ============================================================
# 8. Differential Expression
# ============================================================
print("\n" + "="*60)
print("8. DIFFERENTIAL EXPRESSION")
print("="*60)

# Subset to sepsis samples
adata_sepsis = adata[adata.obs['disease'] == 'Sepsis'].copy()
print(f"Sepsis samples: {adata_sepsis.n_obs} cells")
print(adata_sepsis.obs['outcome'].value_counts())

# Run DE if we have both survivor and non-survivor
if adata_sepsis.obs['outcome'].nunique() >= 2:
    sc.tl.rank_genes_groups(
        adata_sepsis,
        groupby='outcome',
        method='wilcoxon',
        use_raw=False
    )
    adata.uns['rank_genes_groups_outcome'] = adata_sepsis.uns['rank_genes_groups']

    # Get results
    de_results = sc.get.rank_genes_groups_df(adata_sepsis, group='Survivor')
    print("\nTop 20 DE genes (Survivor vs Non-survivor):")
    print(de_results.head(20).to_string())

    # Save DE results
    de_results.to_csv(PROCESSED_DIR / 'de_results_outcome.csv', index=False)
    print("\nSaved DE results")
else:
    print("Cannot run DE: need both Survivor and Non-survivor groups")
    de_results = None

# ============================================================
# 9. CD52 Analysis
# ============================================================
print("\n" + "="*60)
print("9. CD52 ANALYSIS - Proof of Concept Target")
print("="*60)

if 'CD52' in adata.var_names:
    print("CD52 found in gene list!")

    if de_results is not None:
        cd52_de = de_results[de_results['names'] == 'CD52']
        if not cd52_de.empty:
            print("\nCD52 DE results (Survivor vs Non-survivor):")
            print(cd52_de.to_string(index=False))

    # Plot CD52
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    sc.pl.umap(adata, color='CD52', ax=axes[0], show=False, title='CD52 Expression')
    sc.pl.violin(adata, 'CD52', groupby='outcome', ax=axes[1], show=False)
    axes[1].set_title('CD52 by Outcome')
    sc.pl.violin(adata, 'CD52', groupby='cell_type', rotation=45, ax=axes[2], show=False)
    axes[2].set_title('CD52 by Cell Type')
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / 'cd52_expression_overview.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Saved CD52 plots")
else:
    print("WARNING: CD52 not found in gene names!")

# ============================================================
# 10. Save
# ============================================================
print("\n" + "="*60)
print("10. SAVING DATA")
print("="*60)

# Summary
print(f"\nTotal cells: {adata.n_obs}")
print(f"Total genes: {adata.n_vars}")
print(f"Highly variable genes: {adata.var['highly_variable'].sum()}")
print(f"\nCell types:\n{adata.obs['cell_type'].value_counts().to_string()}")
print(f"\nOutcome:\n{adata.obs['outcome'].value_counts().to_string()}")
print(f"\nLayers: {list(adata.layers.keys())}")
print(f"Obsm: {list(adata.obsm.keys())}")

# Save
output_path = PROCESSED_DIR / 'adata_processed.h5ad'
adata.write(output_path)
print(f"\nSaved processed data to: {output_path}")

print("\n" + "="*60)
print("DATA PREPARATION COMPLETE!")
print("="*60)
