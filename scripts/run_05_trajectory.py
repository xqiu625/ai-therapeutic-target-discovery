#!/usr/bin/env python
"""
05 - Trajectory Analysis
Sepsis Target Discovery Pipeline

Performs trajectory/pseudotime analysis to capture:
- Disease progression from T0 to T6
- Gene expression dynamics along trajectory
- CD52 dynamics along trajectory
"""

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy.stats import spearmanr
import warnings
warnings.filterwarnings('ignore')

# Paths
# Import paths from config
try:
    from config import BASE_DIR, PROCESSED_DIR, TABLES_DIR, FIGURES_DIR, EXTERNAL_DIR
except ImportError:
    # Fallback for standalone execution
    BASE_DIR = Path(__file__).parent.parent.resolve()
    PROCESSED_DIR = BASE_DIR / 'data' / 'processed'
    TABLES_DIR = BASE_DIR / 'results' / 'tables'
    FIGURES_DIR = BASE_DIR / 'figures'
    EXTERNAL_DIR = BASE_DIR / 'data' / 'external'
PROCESSED_DIR = BASE_DIR / 'data' / 'processed'
RESULTS_DIR = BASE_DIR / 'results'
TABLES_DIR = RESULTS_DIR / 'tables'
FIGURES_DIR = BASE_DIR / 'figures'

TABLES_DIR.mkdir(parents=True, exist_ok=True)

print("="*60)
print("05 - TRAJECTORY ANALYSIS")
print("="*60)

# ============================================================
# 1. Load Data
# ============================================================
print("\n" + "="*60)
print("1. LOADING DATA")
print("="*60)

# Try to load in order of preference
for fname in ['adata_grn.h5ad', 'adata_embeddings.h5ad', 'adata_processed.h5ad']:
    adata_file = PROCESSED_DIR / fname
    if adata_file.exists():
        break

adata = sc.read_h5ad(adata_file)
print(f"Loaded: {adata.n_obs} cells, {adata.n_vars} genes")
print(f"From: {adata_file}")

# ============================================================
# 2. Subset to Sepsis Samples
# ============================================================
print("\n" + "="*60)
print("2. SUBSETTING TO SEPSIS")
print("="*60)

adata_sepsis = adata[adata.obs['disease'] == 'Sepsis'].copy()
print(f"Sepsis cells: {adata_sepsis.n_obs}")
print(f"\nOutcome distribution:")
print(adata_sepsis.obs['outcome'].value_counts())

# Check timepoint column
timepoint_col = 'timepoint' if 'timepoint' in adata_sepsis.obs.columns else 'time'
if timepoint_col in adata_sepsis.obs.columns:
    print(f"\nTimepoint distribution:")
    print(adata_sepsis.obs[timepoint_col].value_counts())
else:
    print("\nNo timepoint column found")
    timepoint_col = None

# ============================================================
# 3. Compute Neighbors (if not present)
# ============================================================
print("\n" + "="*60)
print("3. COMPUTING NEIGHBORS")
print("="*60)

if 'neighbors' not in adata_sepsis.uns:
    print("Computing neighbors...")
    sc.pp.neighbors(adata_sepsis, n_neighbors=15, n_pcs=30)
else:
    print("Using existing neighbors")

# ============================================================
# 4. Diffusion Map
# ============================================================
print("\n" + "="*60)
print("4. DIFFUSION MAP")
print("="*60)

sc.tl.diffmap(adata_sepsis, n_comps=15)
print("Diffusion map computed")

# Visualize
fig, axes = plt.subplots(1, 3, figsize=(15, 5))
sc.pl.diffmap(adata_sepsis, color='outcome', components=['1,2'], ax=axes[0], show=False, title='Outcome')
sc.pl.diffmap(adata_sepsis, color='cell_type', components=['1,2'], ax=axes[1], show=False, title='Cell Type')
if timepoint_col:
    sc.pl.diffmap(adata_sepsis, color=timepoint_col, components=['1,2'], ax=axes[2], show=False, title='Timepoint')
plt.tight_layout()
plt.savefig(FIGURES_DIR / 'diffmap_components.png', dpi=150, bbox_inches='tight')
plt.close()
print("Saved diffusion map plot")

# ============================================================
# 5. Diffusion Pseudotime
# ============================================================
print("\n" + "="*60)
print("5. DIFFUSION PSEUDOTIME")
print("="*60)

# Find root cell (earliest timepoint or random if no timepoint)
dc1 = adata_sepsis.obsm['X_diffmap'][:, 0]

if timepoint_col and 'T0' in adata_sepsis.obs[timepoint_col].values:
    t0_mask = adata_sepsis.obs[timepoint_col] == 'T0'
    t0_dc1 = dc1[t0_mask]
    root_idx = np.where(t0_mask)[0][np.argmin(t0_dc1)]
    print(f"Root cell from T0: {adata_sepsis.obs_names[root_idx]}")
else:
    # Use cell with minimum diffusion component 1
    root_idx = np.argmin(dc1)
    print(f"Root cell (min DC1): {adata_sepsis.obs_names[root_idx]}")

adata_sepsis.uns['iroot'] = root_idx

# Compute DPT
sc.tl.dpt(adata_sepsis, n_branchings=0)
pseudotime = adata_sepsis.obs['dpt_pseudotime'].values
print(f"Pseudotime range: {pseudotime.min():.3f} - {pseudotime.max():.3f}")

# Visualize
fig, axes = plt.subplots(1, 3, figsize=(15, 5))
sc.pl.umap(adata_sepsis, color='dpt_pseudotime', ax=axes[0], show=False, title='Pseudotime')
sc.pl.umap(adata_sepsis, color='outcome', ax=axes[1], show=False, title='Outcome')
if timepoint_col:
    sc.pl.umap(adata_sepsis, color=timepoint_col, ax=axes[2], show=False, title='Timepoint')
plt.tight_layout()
plt.savefig(FIGURES_DIR / 'fig5_trajectory_pseudotime.png', dpi=150, bbox_inches='tight')
plt.close()
print("Saved pseudotime plot")

# ============================================================
# 6. Gene Dynamics Along Trajectory
# ============================================================
print("\n" + "="*60)
print("6. GENE TRAJECTORY CORRELATIONS")
print("="*60)

gene_names = adata_sepsis.var_names.tolist()
X = adata_sepsis.X
if hasattr(X, 'toarray'):
    X = X.toarray()

correlations = []
for i, gene in enumerate(gene_names):
    expr = X[:, i]
    if expr.std() > 0:
        corr, pval = spearmanr(pseudotime, expr)
        correlations.append({
            'gene': gene,
            'correlation': corr,
            'pvalue': pval,
            'abs_correlation': abs(corr)
        })
    if i % 5000 == 0:
        print(f"  Processed {i}/{len(gene_names)} genes")

trajectory_genes = pd.DataFrame(correlations)
trajectory_genes = trajectory_genes.sort_values('abs_correlation', ascending=False)

print(f"\nTop 20 genes correlated with pseudotime:")
print(trajectory_genes.head(20).to_string())

# Check CD52
cd52_traj = trajectory_genes[trajectory_genes['gene'] == 'CD52']
if not cd52_traj.empty:
    cd52_idx = trajectory_genes[trajectory_genes['gene'] == 'CD52'].index[0]
    cd52_rank = (trajectory_genes['abs_correlation'] > trajectory_genes.loc[cd52_idx, 'abs_correlation']).sum() + 1
    cd52_corr = cd52_traj['correlation'].values[0]
    print(f"\n*** CD52 trajectory rank: {cd52_rank} (correlation: {cd52_corr:.4f}) ***")
else:
    cd52_rank = None
    print("\nCD52 not found")

# ============================================================
# 7. CD52 Trajectory Plot
# ============================================================
print("\n" + "="*60)
print("7. CD52 TRAJECTORY VISUALIZATION")
print("="*60)

if 'CD52' in adata_sepsis.var_names:
    cd52_expr = adata_sepsis[:, 'CD52'].X
    if hasattr(cd52_expr, 'toarray'):
        cd52_expr = cd52_expr.toarray().flatten()
    else:
        cd52_expr = np.array(cd52_expr).flatten()

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Scatter by outcome
    for outcome, color in [('Survivor', 'blue'), ('Non-survivor', 'red')]:
        mask = adata_sepsis.obs['outcome'] == outcome
        axes[0].scatter(pseudotime[mask], cd52_expr[mask], c=color, alpha=0.2, s=5, label=outcome)
    axes[0].set_xlabel('Pseudotime')
    axes[0].set_ylabel('CD52 Expression')
    axes[0].set_title('CD52 Along Trajectory')
    axes[0].legend()

    # Smoothed trajectory by outcome
    for outcome, color in [('Survivor', 'blue'), ('Non-survivor', 'red')]:
        mask = adata_sepsis.obs['outcome'] == outcome
        pt_sub = pseudotime[mask]
        expr_sub = cd52_expr[mask]

        bins = np.linspace(0, 1, 20)
        bin_means = []
        bin_centers = []
        for j in range(len(bins)-1):
            bin_mask = (pt_sub >= bins[j]) & (pt_sub < bins[j+1])
            if bin_mask.sum() > 10:
                bin_means.append(expr_sub[bin_mask].mean())
                bin_centers.append((bins[j] + bins[j+1]) / 2)

        axes[1].plot(bin_centers, bin_means, '-o', color=color, label=outcome, linewidth=2)

    axes[1].set_xlabel('Pseudotime')
    axes[1].set_ylabel('Mean CD52 Expression')
    axes[1].set_title('CD52 Trajectory by Outcome')
    axes[1].legend()

    plt.tight_layout()
    plt.savefig(FIGURES_DIR / 'cd52_trajectory.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Saved CD52 trajectory plot")

# ============================================================
# 8. Pseudotime by Outcome
# ============================================================
print("\n" + "="*60)
print("8. PSEUDOTIME DISTRIBUTIONS")
print("="*60)

fig, ax = plt.subplots(figsize=(8, 5))
for outcome, color in [('Survivor', 'blue'), ('Non-survivor', 'red')]:
    mask = adata_sepsis.obs['outcome'] == outcome
    ax.hist(pseudotime[mask], bins=30, alpha=0.5, label=outcome, density=True, color=color)
ax.set_xlabel('Pseudotime')
ax.set_ylabel('Density')
ax.set_title('Pseudotime Distribution by Outcome')
ax.legend()
plt.tight_layout()
plt.savefig(FIGURES_DIR / 'pseudotime_by_outcome.png', dpi=150, bbox_inches='tight')
plt.close()

# Stats
print("\nMean pseudotime by outcome:")
print(adata_sepsis.obs.groupby('outcome')['dpt_pseudotime'].mean())

# ============================================================
# 9. Save Results
# ============================================================
print("\n" + "="*60)
print("9. SAVING RESULTS")
print("="*60)

# Save trajectory correlations
trajectory_genes.to_csv(TABLES_DIR / 'trajectory_correlations.csv', index=False)
print("Saved trajectory correlations")

# Add pseudotime to main adata
adata.obs['dpt_pseudotime'] = np.nan
adata.obs.loc[adata_sepsis.obs_names, 'dpt_pseudotime'] = adata_sepsis.obs['dpt_pseudotime']

# Save
adata.write(PROCESSED_DIR / 'adata_trajectory.h5ad')
print("Saved adata_trajectory.h5ad")

# ============================================================
# Summary
# ============================================================
print("\n" + "="*60)
print("TRAJECTORY ANALYSIS SUMMARY")
print("="*60)
print(f"Cells analyzed: {adata_sepsis.n_obs}")
print(f"Pseudotime range: {pseudotime.min():.3f} - {pseudotime.max():.3f}")
print(f"\nTop 10 trajectory-correlated genes:")
for i, (_, row) in enumerate(trajectory_genes.head(10).iterrows()):
    marker = " *** CD52 ***" if row['gene'] == 'CD52' else ""
    print(f"  {i+1}. {row['gene']}: r={row['correlation']:.3f}{marker}")

if cd52_rank:
    print(f"\nCD52 trajectory rank: {cd52_rank}")

print("\n" + "="*60)
print("TRAJECTORY ANALYSIS COMPLETE!")
print("="*60)
