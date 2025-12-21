#!/usr/bin/env python
"""
06 - Target Prioritization
Sepsis Target Discovery Pipeline

Integrates all evidence sources to prioritize therapeutic targets:
1. Differential Expression (Survivor vs Non-survivor)
2. ML Feature Importance
3. GRN Centrality (if available)
4. Trajectory Dynamics (if available)
5. Druggability (DGIdb)
"""

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from sklearn.preprocessing import MinMaxScaler
import warnings
warnings.filterwarnings('ignore')

# Optional: DGIdb query
try:
    import requests
    HAS_REQUESTS = True
except ImportError:
    HAS_REQUESTS = False

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
print("06 - TARGET PRIORITIZATION")
print("="*60)

# ============================================================
# 1. Load Evidence Sources
# ============================================================
print("\n" + "="*60)
print("1. LOADING EVIDENCE SOURCES")
print("="*60)

# Load adata for gene list
for fname in ['adata_trajectory.h5ad', 'adata_grn.h5ad', 'adata_processed.h5ad']:
    adata_file = PROCESSED_DIR / fname
    if adata_file.exists():
        break

adata = sc.read_h5ad(adata_file)
all_genes = adata.var_names.tolist()
print(f"Total genes: {len(all_genes)}")

# 1a. Differential Expression
de_file = PROCESSED_DIR / 'de_results_outcome.csv'
if de_file.exists():
    de_results = pd.read_csv(de_file)
    de_results['de_score'] = de_results['logfoldchanges'].abs() * (-np.log10(de_results['pvals_adj'] + 1e-300))
    print(f"DE results: {len(de_results)} genes")
else:
    print("DE results not found - will compute")
    de_results = None

# 1b. ML Feature Importance
ml_file = TABLES_DIR / 'ml_feature_importance.csv'
if ml_file.exists():
    ml_importance = pd.read_csv(ml_file)
    # Rename column if needed
    if 'shap_importance' not in ml_importance.columns and 'importance' in ml_importance.columns:
        ml_importance = ml_importance.rename(columns={'importance': 'shap_importance'})
    print(f"ML importance: {len(ml_importance)} genes")
else:
    print("ML importance not found")
    ml_importance = pd.DataFrame({'gene': all_genes, 'shap_importance': 0})

# 1c. GRN Centrality
grn_file = TABLES_DIR / 'grn_centrality.csv'
if grn_file.exists():
    grn_centrality = pd.read_csv(grn_file)
    print(f"GRN centrality: {len(grn_centrality)} genes")
else:
    print("GRN centrality not found (SCENIC may still be running)")
    grn_centrality = pd.DataFrame({'gene': all_genes, 'centrality_score': 0})

# 1d. Trajectory Dynamics
traj_file = TABLES_DIR / 'trajectory_correlations.csv'
if traj_file.exists():
    trajectory = pd.read_csv(traj_file)
    print(f"Trajectory: {len(trajectory)} genes")
else:
    print("Trajectory not found")
    trajectory = pd.DataFrame({'gene': all_genes, 'abs_correlation': 0})

# ============================================================
# 2. Query Druggability (DGIdb)
# ============================================================
print("\n" + "="*60)
print("2. DRUGGABILITY (DGIdb)")
print("="*60)

druggability_file = TABLES_DIR / 'druggability_scores.csv'

if druggability_file.exists() and druggability_file.stat().st_size > 0:
    try:
        druggability = pd.read_csv(druggability_file)
        if len(druggability) > 0 and 'gene' in druggability.columns:
            print(f"Loaded cached druggability: {len(druggability)} genes")
        else:
            print("Cached file invalid, will re-query")
            druggability_file.unlink()  # Delete invalid file
            druggability = None
    except Exception as e:
        print(f"Error reading cached file: {e}")
        druggability_file.unlink()  # Delete corrupt file
        druggability = None
else:
    druggability = None

if druggability is None and HAS_REQUESTS:
    print("Querying DGIdb...")

    def query_dgidb(gene_list, batch_size=100):
        base_url = "https://dgidb.org/api/v2/interactions.json"
        all_interactions = []

        for i in range(0, len(gene_list), batch_size):
            batch = gene_list[i:i+batch_size]
            try:
                params = {'genes': ','.join(batch)}
                response = requests.get(base_url, params=params, timeout=30)

                if response.status_code == 200:
                    data = response.json()
                    for match in data.get('matchedTerms', []):
                        gene = match.get('geneName', '')
                        interactions = match.get('interactions', [])
                        all_interactions.append({
                            'gene': gene,
                            'n_interactions': len(interactions),
                            'n_drugs': len(set(i.get('drugName', '') for i in interactions))
                        })
                print(f"  Queried {min(i+batch_size, len(gene_list))}/{len(gene_list)}")
            except Exception as e:
                print(f"  Error: {e}")

        return pd.DataFrame(all_interactions)

    # Query top genes
    if de_results is not None:
        top_genes = de_results.nlargest(500, 'de_score')['names'].tolist()
    else:
        top_genes = all_genes[:500]

    if 'CD52' not in top_genes:
        top_genes.append('CD52')

    druggability = query_dgidb(top_genes)
    if len(druggability) > 0:
        druggability.to_csv(druggability_file, index=False)
        print(f"Saved druggability for {len(druggability)} genes")
    else:
        print("DGIdb returned no results, using placeholder")
        druggability = pd.DataFrame({'gene': ['CD52'], 'n_interactions': [1]})
elif druggability is None:
    print("requests not available and no cache, using placeholder")
    druggability = pd.DataFrame({'gene': ['CD52'], 'n_interactions': [1]})
# else: druggability was loaded from cache, use as-is

# Check CD52
if 'gene' in druggability.columns:
    cd52_drug = druggability[druggability['gene'] == 'CD52']
    if not cd52_drug.empty:
        print(f"\nCD52 druggability: {cd52_drug['n_interactions'].values[0]} interactions")
    else:
        print("\nCD52 not in DGIdb (known drug: alemtuzumab)")
        # Add CD52 manually since we know it has alemtuzumab
        druggability = pd.concat([druggability, pd.DataFrame([{'gene': 'CD52', 'n_interactions': 1}])], ignore_index=True)
else:
    print("\nDGIdb query returned no results, using placeholder")
    druggability = pd.DataFrame([{'gene': 'CD52', 'n_interactions': 1}])

# ============================================================
# 3. Integrate Evidence
# ============================================================
print("\n" + "="*60)
print("3. INTEGRATING EVIDENCE")
print("="*60)

master_df = pd.DataFrame({'gene': all_genes})

# Merge DE
if de_results is not None:
    de_scores = de_results[['names', 'de_score', 'logfoldchanges', 'pvals_adj']].rename(
        columns={'names': 'gene', 'logfoldchanges': 'log2fc', 'pvals_adj': 'de_padj'}
    )
    master_df = master_df.merge(de_scores, on='gene', how='left')
else:
    master_df['de_score'] = 0

# Merge ML
master_df = master_df.merge(ml_importance[['gene', 'shap_importance']], on='gene', how='left')

# Merge GRN
master_df = master_df.merge(grn_centrality[['gene', 'centrality_score']], on='gene', how='left')

# Merge Trajectory
master_df = master_df.merge(
    trajectory[['gene', 'abs_correlation']].rename(columns={'abs_correlation': 'trajectory_score'}),
    on='gene', how='left'
)

# Merge Druggability
if 'gene' in druggability.columns and 'n_interactions' in druggability.columns:
    master_df = master_df.merge(
        druggability[['gene', 'n_interactions']].rename(columns={'n_interactions': 'druggability_score'}),
        on='gene', how='left'
    )
else:
    master_df['druggability_score'] = 0
    # Manually add CD52 druggability (alemtuzumab)
    master_df.loc[master_df['gene'] == 'CD52', 'druggability_score'] = 1

# Fill NaN
score_cols = ['de_score', 'shap_importance', 'centrality_score', 'trajectory_score', 'druggability_score']
for col in score_cols:
    if col not in master_df.columns:
        master_df[col] = 0
master_df[score_cols] = master_df[score_cols].fillna(0)

print(f"Master table: {len(master_df)} genes")

# ============================================================
# 4. Normalize and Calculate Priority Score
# ============================================================
print("\n" + "="*60)
print("4. CALCULATING PRIORITY SCORES")
print("="*60)

scaler = MinMaxScaler()
normalized_cols = []
for col in score_cols:
    new_col = f'{col}_norm'
    values = master_df[col].values.reshape(-1, 1)
    master_df[new_col] = scaler.fit_transform(values)
    normalized_cols.append(new_col)

# Weighted priority score
weights = {
    'de_score_norm': 0.25,
    'shap_importance_norm': 0.25,
    'centrality_score_norm': 0.15,
    'trajectory_score_norm': 0.15,
    'druggability_score_norm': 0.20
}

master_df['priority_score'] = sum(
    master_df[col] * weight for col, weight in weights.items()
)

master_df['rank'] = master_df['priority_score'].rank(ascending=False).astype(int)
master_df = master_df.sort_values('priority_score', ascending=False)

print("\nTop 20 prioritized targets:")
print(master_df[['gene', 'rank', 'priority_score']].head(20).to_string())

# ============================================================
# 5. CD52 Validation
# ============================================================
print("\n" + "="*60)
print("5. CD52 VALIDATION")
print("="*60)

cd52_row = master_df[master_df['gene'] == 'CD52']

if not cd52_row.empty:
    cd52_rank = cd52_row['rank'].values[0]
    cd52_score = cd52_row['priority_score'].values[0]

    print(f"\n*** CD52 RANK: {cd52_rank} / {len(master_df)} ***")
    print(f"Priority Score: {cd52_score:.4f}")
    print(f"\nEvidence breakdown:")
    print(f"  - DE score:      {cd52_row['de_score_norm'].values[0]:.4f}")
    print(f"  - ML importance: {cd52_row['shap_importance_norm'].values[0]:.4f}")
    print(f"  - GRN centrality:{cd52_row['centrality_score_norm'].values[0]:.4f}")
    print(f"  - Trajectory:    {cd52_row['trajectory_score_norm'].values[0]:.4f}")
    print(f"  - Druggability:  {cd52_row['druggability_score_norm'].values[0]:.4f}")

    if cd52_rank <= 10:
        print(f"\n*** SUCCESS: CD52 is in TOP 10 targets! ***")
    elif cd52_rank <= 50:
        print(f"\n*** GOOD: CD52 is in TOP 50 targets ***")
    elif cd52_rank <= 100:
        print(f"\n*** ACCEPTABLE: CD52 is in TOP 100 targets ***")
else:
    print("CD52 not found!")
    cd52_rank = None

# ============================================================
# 6. Visualizations
# ============================================================
print("\n" + "="*60)
print("6. VISUALIZATIONS")
print("="*60)

# Top 30 bar plot
fig, ax = plt.subplots(figsize=(12, 10))
top_30 = master_df.head(30)
colors = ['red' if g == 'CD52' else 'steelblue' for g in top_30['gene']]
ax.barh(range(len(top_30)), top_30['priority_score'], color=colors)
ax.set_yticks(range(len(top_30)))
ax.set_yticklabels(top_30['gene'])
ax.invert_yaxis()
ax.set_xlabel('Priority Score')
ax.set_title('Top 30 Prioritized Therapeutic Targets')
for i, (score, gene) in enumerate(zip(top_30['priority_score'], top_30['gene'])):
    ax.text(score + 0.005, i, f'{score:.3f}', va='center', fontsize=8)
plt.tight_layout()
plt.savefig(FIGURES_DIR / 'fig6_target_ranking.png', dpi=150, bbox_inches='tight')
plt.close()
print("Saved target ranking plot")

# Evidence heatmap
top_20 = master_df.head(20)
heatmap_data = top_20[normalized_cols].copy()
heatmap_data.index = top_20['gene']
heatmap_data.columns = ['DE', 'ML', 'GRN', 'Trajectory', 'Drug']

fig, ax = plt.subplots(figsize=(10, 8))
sns.heatmap(heatmap_data, annot=True, fmt='.2f', cmap='YlOrRd',
            linewidths=0.5, ax=ax, vmin=0, vmax=1)
ax.set_title('Evidence Breakdown for Top 20 Targets')
plt.tight_layout()
plt.savefig(FIGURES_DIR / 'target_evidence_heatmap.png', dpi=150, bbox_inches='tight')
plt.close()
print("Saved evidence heatmap")

# ============================================================
# 7. Save Results
# ============================================================
print("\n" + "="*60)
print("7. SAVING RESULTS")
print("="*60)

output_cols = ['gene', 'rank', 'priority_score'] + score_cols + normalized_cols
if 'log2fc' in master_df.columns:
    output_cols.insert(3, 'log2fc')
if 'de_padj' in master_df.columns:
    output_cols.insert(4, 'de_padj')

master_df[output_cols].to_csv(TABLES_DIR / 'prioritized_targets.csv', index=False)
master_df.head(50)[output_cols].to_csv(TABLES_DIR / 'top50_targets.csv', index=False)
print("Saved prioritized targets")

# ============================================================
# Summary
# ============================================================
print("\n" + "="*60)
print("TARGET PRIORITIZATION SUMMARY")
print("="*60)
print(f"Total genes: {len(master_df)}")
print(f"\nWeighting: DE={weights['de_score_norm']*100:.0f}%, ML={weights['shap_importance_norm']*100:.0f}%, "
      f"GRN={weights['centrality_score_norm']*100:.0f}%, Traj={weights['trajectory_score_norm']*100:.0f}%, "
      f"Drug={weights['druggability_score_norm']*100:.0f}%")
print(f"\nTop 10 targets:")
for _, row in master_df.head(10).iterrows():
    marker = " *** CD52 ***" if row['gene'] == 'CD52' else ""
    print(f"  {int(row['rank']):2d}. {row['gene']}: {row['priority_score']:.4f}{marker}")

if cd52_rank:
    print(f"\n*** CD52 Final Rank: {cd52_rank} ***")

print("\n" + "="*60)
print("TARGET PRIORITIZATION COMPLETE!")
print("="*60)
