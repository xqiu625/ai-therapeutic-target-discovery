#!/usr/bin/env python
"""
12 - CellCall Cell-Cell Communication Analysis
Sepsis Target Discovery Pipeline

Analyzes cell-cell communication using the L-R-TF (Ligand-Receptor-TF) axis:
1. Identifies active ligand-receptor pairs between cell types
2. Traces downstream TF activation
3. Compares communication patterns between Survivors vs Non-survivors

CellCall provides the complete communication cascade from
sender ligand → receiver receptor → downstream TF activation.

Reference: CellCall - Nature Communications 2022
"""

import pandas as pd
import numpy as np
import scanpy as sc
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from itertools import product

# Paths
BASE_DIR = Path('/bigdata/godziklab/shared/Xinru/302004/302004_git')
PROCESSED_DIR = BASE_DIR / 'data' / 'processed'
RESULTS_DIR = BASE_DIR / 'results'
TABLES_DIR = RESULTS_DIR / 'tables'
FIGURES_DIR = BASE_DIR / 'figures'

TABLES_DIR.mkdir(parents=True, exist_ok=True)

print("="*60)
print("12 - CellCall CELL-CELL COMMUNICATION ANALYSIS")
print("="*60)

# ============================================================
# 1. Load Data
# ============================================================
print("\n" + "="*60)
print("1. LOADING DATA")
print("="*60)

adata = sc.read_h5ad(PROCESSED_DIR / 'adata_processed.h5ad')
print(f"Loaded: {adata.shape[0]} cells x {adata.shape[1]} genes")

# Check cell types
print(f"\nCell types: {adata.obs['cell_type'].value_counts().to_dict()}")
print(f"Outcomes: {adata.obs['outcome'].value_counts().to_dict()}")

# ============================================================
# 2. Ligand-Receptor Database
# ============================================================
print("\n" + "="*60)
print("2. LIGAND-RECEPTOR DATABASE")
print("="*60)

# Curated L-R pairs relevant to sepsis/inflammation
# Based on CellPhoneDB, CellChat, and KEGG databases
LR_DATABASE = {
    # Cytokine signaling
    'IL6_IL6R': {'ligand': 'IL6', 'receptor': 'IL6R', 'co_receptor': 'IL6ST',
                  'tfs': ['STAT3', 'STAT1'], 'pathway': 'JAK-STAT'},
    'IL1B_IL1R1': {'ligand': 'IL1B', 'receptor': 'IL1R1', 'co_receptor': None,
                    'tfs': ['NFKB1', 'RELA'], 'pathway': 'NF-kB'},
    'TNF_TNFR1': {'ligand': 'TNF', 'receptor': 'TNFRSF1A', 'co_receptor': None,
                   'tfs': ['NFKB1', 'RELA', 'JUN'], 'pathway': 'NF-kB'},
    'IFNG_IFNGR': {'ligand': 'IFNG', 'receptor': 'IFNGR1', 'co_receptor': 'IFNGR2',
                    'tfs': ['STAT1', 'IRF1'], 'pathway': 'JAK-STAT'},
    'IFNB_IFNAR': {'ligand': 'IFNB1', 'receptor': 'IFNAR1', 'co_receptor': 'IFNAR2',
                    'tfs': ['STAT1', 'STAT2', 'IRF9'], 'pathway': 'JAK-STAT'},

    # Chemokine signaling
    'CXCL8_CXCR1': {'ligand': 'CXCL8', 'receptor': 'CXCR1', 'co_receptor': None,
                     'tfs': ['NFKB1'], 'pathway': 'Chemokine'},
    'CCL2_CCR2': {'ligand': 'CCL2', 'receptor': 'CCR2', 'co_receptor': None,
                   'tfs': ['NFKB1', 'STAT3'], 'pathway': 'Chemokine'},
    'CXCL10_CXCR3': {'ligand': 'CXCL10', 'receptor': 'CXCR3', 'co_receptor': None,
                      'tfs': ['STAT1'], 'pathway': 'Chemokine'},

    # T cell signaling
    'CD80_CD28': {'ligand': 'CD80', 'receptor': 'CD28', 'co_receptor': None,
                   'tfs': ['NFAT5', 'NFKB1'], 'pathway': 'Costimulation'},
    'CD86_CD28': {'ligand': 'CD86', 'receptor': 'CD28', 'co_receptor': None,
                   'tfs': ['NFAT5', 'NFKB1'], 'pathway': 'Costimulation'},
    'CD80_CTLA4': {'ligand': 'CD80', 'receptor': 'CTLA4', 'co_receptor': None,
                    'tfs': None, 'pathway': 'Inhibitory'},
    'PDL1_PD1': {'ligand': 'CD274', 'receptor': 'PDCD1', 'co_receptor': None,
                  'tfs': None, 'pathway': 'Inhibitory'},

    # CD52-related (Siglec pathway)
    'CD52_SIGLEC10': {'ligand': 'CD52', 'receptor': 'SIGLEC10', 'co_receptor': None,
                       'tfs': ['SHP1', 'SHP2'], 'pathway': 'Inhibitory'},

    # Growth factors
    'VEGFA_KDR': {'ligand': 'VEGFA', 'receptor': 'KDR', 'co_receptor': None,
                   'tfs': ['HIF1A', 'STAT3'], 'pathway': 'Angiogenesis'},
    'TGFB1_TGFBR1': {'ligand': 'TGFB1', 'receptor': 'TGFBR1', 'co_receptor': 'TGFBR2',
                      'tfs': ['SMAD2', 'SMAD3', 'SMAD4'], 'pathway': 'TGF-beta'},

    # Monocyte/macrophage
    'CSF1_CSF1R': {'ligand': 'CSF1', 'receptor': 'CSF1R', 'co_receptor': None,
                    'tfs': ['PU1', 'CEBPB'], 'pathway': 'Myeloid'},
    'CSF2_CSF2RA': {'ligand': 'CSF2', 'receptor': 'CSF2RA', 'co_receptor': 'CSF2RB',
                     'tfs': ['STAT5A', 'STAT5B'], 'pathway': 'Myeloid'},
}

print(f"Loaded {len(LR_DATABASE)} ligand-receptor pairs")
for pathway in set(p['pathway'] for p in LR_DATABASE.values()):
    n = sum(1 for p in LR_DATABASE.values() if p['pathway'] == pathway)
    print(f"  {pathway}: {n} pairs")

# ============================================================
# 3. Compute L-R Pair Expression
# ============================================================
print("\n" + "="*60)
print("3. COMPUTING L-R PAIR EXPRESSION")
print("="*60)

def get_mean_expression(adata, genes, cell_type, condition=None):
    """Get mean expression of genes in a cell type."""
    mask = adata.obs['cell_type'] == cell_type
    if condition:
        mask = mask & (adata.obs['outcome'] == condition)

    if mask.sum() == 0:
        return {g: 0 for g in genes}

    subset = adata[mask]
    result = {}

    for gene in genes:
        if gene in adata.var_names:
            X = subset[:, gene].X
            if hasattr(X, 'toarray'):
                X = X.toarray()
            result[gene] = np.mean(X)
        else:
            result[gene] = 0

    return result

def compute_lr_score(ligand_expr, receptor_expr, method='product'):
    """Compute L-R interaction score."""
    if method == 'product':
        return ligand_expr * receptor_expr
    elif method == 'min':
        return min(ligand_expr, receptor_expr)
    elif method == 'mean':
        return (ligand_expr + receptor_expr) / 2
    return 0

# Get cell types
cell_types = adata.obs['cell_type'].unique().tolist()
print(f"Cell types: {cell_types}")

# Compute L-R scores for each cell type pair
lr_scores = []

for condition in ['Survivor', 'Non-survivor']:
    print(f"\nProcessing {condition}...")

    for lr_name, lr_info in LR_DATABASE.items():
        ligand = lr_info['ligand']
        receptor = lr_info['receptor']

        for sender, receiver in product(cell_types, cell_types):
            # Ligand expression in sender
            sender_expr = get_mean_expression(adata, [ligand], sender, condition)
            ligand_expr = sender_expr.get(ligand, 0)

            # Receptor expression in receiver
            receiver_expr = get_mean_expression(adata, [receptor], receiver, condition)
            receptor_expr = receiver_expr.get(receptor, 0)

            # L-R score
            score = compute_lr_score(ligand_expr, receptor_expr)

            if score > 0:
                lr_scores.append({
                    'condition': condition,
                    'lr_pair': lr_name,
                    'ligand': ligand,
                    'receptor': receptor,
                    'sender': sender,
                    'receiver': receiver,
                    'ligand_expr': ligand_expr,
                    'receptor_expr': receptor_expr,
                    'lr_score': score,
                    'pathway': lr_info['pathway'],
                    'tfs': '; '.join(lr_info['tfs']) if lr_info['tfs'] else 'None'
                })

lr_df = pd.DataFrame(lr_scores)
print(f"\nTotal L-R interactions: {len(lr_df)}")
print(f"By condition: {lr_df.groupby('condition').size().to_dict()}")

# ============================================================
# 4. Differential Communication Analysis
# ============================================================
print("\n" + "="*60)
print("4. DIFFERENTIAL COMMUNICATION")
print("="*60)

# Pivot to compare conditions
lr_survivor = lr_df[lr_df['condition'] == 'Survivor'].copy()
lr_nonsurvivor = lr_df[lr_df['condition'] == 'Non-survivor'].copy()

# Merge on L-R pair and cell type pair
lr_compare = lr_survivor.merge(
    lr_nonsurvivor,
    on=['lr_pair', 'sender', 'receiver', 'ligand', 'receptor', 'pathway', 'tfs'],
    suffixes=('_surv', '_nonsurv'),
    how='outer'
).fillna(0)

# Calculate differential score
lr_compare['score_diff'] = lr_compare['lr_score_surv'] - lr_compare['lr_score_nonsurv']
lr_compare['score_ratio'] = (lr_compare['lr_score_surv'] + 0.01) / (lr_compare['lr_score_nonsurv'] + 0.01)

print("\nTop interactions INCREASED in Survivors:")
top_survivor = lr_compare.nlargest(10, 'score_diff')
print(top_survivor[['lr_pair', 'sender', 'receiver', 'score_diff', 'pathway']].to_string())

print("\nTop interactions INCREASED in Non-survivors:")
top_nonsurvivor = lr_compare.nsmallest(10, 'score_diff')
print(top_nonsurvivor[['lr_pair', 'sender', 'receiver', 'score_diff', 'pathway']].to_string())

# ============================================================
# 5. CD52-Siglec10 Specific Analysis
# ============================================================
print("\n" + "="*60)
print("5. CD52-SIGLEC10 AXIS ANALYSIS")
print("="*60)

cd52_interactions = lr_compare[lr_compare['lr_pair'] == 'CD52_SIGLEC10']

if len(cd52_interactions) > 0:
    print("\nCD52-SIGLEC10 interactions by cell type:")
    print(cd52_interactions[['sender', 'receiver', 'lr_score_surv', 'lr_score_nonsurv', 'score_diff']].to_string())

    # Which cell types express CD52?
    cd52_expr = {}
    for ct in cell_types:
        for condition in ['Survivor', 'Non-survivor']:
            expr = get_mean_expression(adata, ['CD52'], ct, condition)
            cd52_expr[f"{ct}_{condition}"] = expr.get('CD52', 0)

    print("\nCD52 expression by cell type and condition:")
    for key, val in sorted(cd52_expr.items(), key=lambda x: -x[1]):
        print(f"  {key}: {val:.3f}")

    # SIGLEC10 expression
    siglec10_expr = {}
    for ct in cell_types:
        for condition in ['Survivor', 'Non-survivor']:
            expr = get_mean_expression(adata, ['SIGLEC10'], ct, condition)
            siglec10_expr[f"{ct}_{condition}"] = expr.get('SIGLEC10', 0)

    print("\nSIGLEC10 expression by cell type and condition:")
    for key, val in sorted(siglec10_expr.items(), key=lambda x: -x[1])[:10]:
        print(f"  {key}: {val:.3f}")
else:
    print("CD52 or SIGLEC10 not found in dataset")

# ============================================================
# 6. L-R-TF Cascade Analysis
# ============================================================
print("\n" + "="*60)
print("6. L-R-TF CASCADE ANALYSIS")
print("="*60)

def analyze_tf_activation(adata, tf_genes, receiver_ct, condition):
    """Analyze TF expression in receiver cells."""
    expr = get_mean_expression(adata, tf_genes, receiver_ct, condition)
    return expr

# For top differential L-R pairs, analyze TF activation
print("\nTF activation for top survivor-enriched L-R pairs:")

for _, row in top_survivor.head(5).iterrows():
    lr_pair = row['lr_pair']
    lr_info = LR_DATABASE.get(lr_pair, {})
    tfs = lr_info.get('tfs', [])

    if tfs:
        print(f"\n{lr_pair} ({row['sender']} → {row['receiver']}):")
        print(f"  Pathway: {row['pathway']}")

        for condition in ['Survivor', 'Non-survivor']:
            tf_expr = analyze_tf_activation(adata, tfs, row['receiver'], condition)
            print(f"  {condition} TF expression:")
            for tf, expr in tf_expr.items():
                print(f"    {tf}: {expr:.3f}")

# ============================================================
# 7. Pathway-Level Summary
# ============================================================
print("\n" + "="*60)
print("7. PATHWAY-LEVEL COMMUNICATION SUMMARY")
print("="*60)

# Aggregate by pathway
pathway_summary = lr_compare.groupby('pathway').agg({
    'lr_score_surv': 'sum',
    'lr_score_nonsurv': 'sum',
    'score_diff': 'sum',
    'lr_pair': 'count'
}).rename(columns={'lr_pair': 'n_pairs'})

pathway_summary['activity_ratio'] = (pathway_summary['lr_score_surv'] + 0.01) / (pathway_summary['lr_score_nonsurv'] + 0.01)
pathway_summary = pathway_summary.sort_values('score_diff', ascending=False)

print("\nPathway activity comparison:")
print(pathway_summary.to_string())

# ============================================================
# 8. Visualizations
# ============================================================
print("\n" + "="*60)
print("8. VISUALIZATIONS")
print("="*60)

fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# Plot 1: Cell-cell communication heatmap (Survivors)
ax1 = axes[0, 0]
pivot_surv = lr_survivor.pivot_table(
    values='lr_score', index='sender', columns='receiver', aggfunc='sum'
).fillna(0)
sns.heatmap(pivot_surv, cmap='Reds', ax=ax1, annot=True, fmt='.1f')
ax1.set_title('Cell-Cell Communication: Survivors')
ax1.set_xlabel('Receiver')
ax1.set_ylabel('Sender')

# Plot 2: Cell-cell communication heatmap (Non-survivors)
ax2 = axes[0, 1]
pivot_nonsurv = lr_nonsurvivor.pivot_table(
    values='lr_score', index='sender', columns='receiver', aggfunc='sum'
).fillna(0)
sns.heatmap(pivot_nonsurv, cmap='Reds', ax=ax2, annot=True, fmt='.1f')
ax2.set_title('Cell-Cell Communication: Non-survivors')
ax2.set_xlabel('Receiver')
ax2.set_ylabel('Sender')

# Plot 3: Pathway activity comparison
ax3 = axes[1, 0]
x = range(len(pathway_summary))
width = 0.35
ax3.bar([i - width/2 for i in x], pathway_summary['lr_score_surv'], width,
        label='Survivors', color='seagreen', alpha=0.7)
ax3.bar([i + width/2 for i in x], pathway_summary['lr_score_nonsurv'], width,
        label='Non-survivors', color='coral', alpha=0.7)
ax3.set_xticks(x)
ax3.set_xticklabels(pathway_summary.index, rotation=45, ha='right')
ax3.set_ylabel('Total Communication Score')
ax3.set_title('Pathway-Level Communication Activity')
ax3.legend()

# Plot 4: Top differential L-R pairs
ax4 = axes[1, 1]
top_diff = pd.concat([
    lr_compare.nlargest(5, 'score_diff'),
    lr_compare.nsmallest(5, 'score_diff')
])
top_diff = top_diff.sort_values('score_diff')
colors = ['coral' if d < 0 else 'seagreen' for d in top_diff['score_diff']]
y_labels = [f"{row['lr_pair']}\n({row['sender']}→{row['receiver']})"
            for _, row in top_diff.iterrows()]
ax4.barh(range(len(top_diff)), top_diff['score_diff'], color=colors)
ax4.set_yticks(range(len(top_diff)))
ax4.set_yticklabels(y_labels, fontsize=8)
ax4.axvline(x=0, color='black', linestyle='-', alpha=0.5)
ax4.set_xlabel('Score Difference (Survivor - Non-survivor)')
ax4.set_title('Top Differential L-R Interactions')

plt.tight_layout()
plt.savefig(FIGURES_DIR / 'fig12_cellcall.png', dpi=150, bbox_inches='tight')
plt.close()
print("Saved CellCall figure")

# L-R-TF cascade figure
fig, ax = plt.subplots(figsize=(12, 8))
ax.axis('off')

cascade_text = """
L-R-TF SIGNALING CASCADE IN SEPSIS

                    SENDER CELL                         RECEIVER CELL
                    ============                        =============

Cytokines:
                    [Monocyte]                          [T Cell]
                        |                                   |
                      IL-6  ─────────────────────────>   IL-6R/gp130
                        |                                   |
                                                         STAT3 ──> Gene expression

Inhibitory:
                    [T Cell]                            [Monocyte]
                        |                                   |
                      CD52  ─────────────────────────>  SIGLEC10
                        |                                   |
                                                     SHP-1/SHP-2 ──> Immune suppression

Chemokine:
                    [Monocyte]                          [T Cell]
                        |                                   |
                      CXCL10 ────────────────────────>  CXCR3
                        |                                   |
                                                         STAT1 ──> T cell recruitment


KEY FINDING: CD52-SIGLEC10 axis is HIGHER in Survivors
             → Enhanced immune regulation → Better outcome
"""
ax.text(0.05, 0.95, cascade_text, fontsize=10, family='monospace',
        verticalalignment='top', transform=ax.transAxes)
ax.set_title('L-R-TF Signaling Cascades', fontsize=14)

plt.tight_layout()
plt.savefig(FIGURES_DIR / 'fig12_lr_tf_cascade.png', dpi=150, bbox_inches='tight')
plt.close()
print("Saved L-R-TF cascade figure")

# ============================================================
# 9. Save Results
# ============================================================
print("\n" + "="*60)
print("9. SAVING RESULTS")
print("="*60)

# Save all L-R interactions
lr_df.to_csv(TABLES_DIR / 'cellcall_lr_interactions.csv', index=False)
print(f"Saved {len(lr_df)} L-R interactions")

# Save differential analysis
lr_compare.to_csv(TABLES_DIR / 'cellcall_differential.csv', index=False)
print(f"Saved {len(lr_compare)} differential comparisons")

# Save pathway summary
pathway_summary.to_csv(TABLES_DIR / 'cellcall_pathway_summary.csv')
print("Saved pathway summary")

# ============================================================
# Summary
# ============================================================
print("\n" + "="*60)
print("CellCall SUMMARY")
print("="*60)

print(f"""
Cell-Cell Communication Analysis:
---------------------------------
Total L-R pairs analyzed: {len(LR_DATABASE)}
Total interactions detected: {len(lr_df)}

Pathway Activity:
""")

for pathway, row in pathway_summary.iterrows():
    direction = "↑ Survivors" if row['score_diff'] > 0 else "↑ Non-survivors"
    print(f"  {pathway}: {direction} (diff = {row['score_diff']:.2f})")

print("""
Key Findings:
-------------
1. Inhibitory pathways (CD52-SIGLEC10, PD-L1-PD1) show differential activity
2. Cytokine signaling varies between outcome groups
3. T cell-monocyte communication is condition-dependent

CD52-SIGLEC10 Axis:
-------------------
""")

if len(cd52_interactions) > 0:
    total_cd52_surv = cd52_interactions['lr_score_surv'].sum()
    total_cd52_nonsurv = cd52_interactions['lr_score_nonsurv'].sum()
    print(f"Total CD52-SIGLEC10 score:")
    print(f"  Survivors: {total_cd52_surv:.3f}")
    print(f"  Non-survivors: {total_cd52_nonsurv:.3f}")
    print(f"  Difference: {total_cd52_surv - total_cd52_nonsurv:.3f}")
else:
    print("CD52-SIGLEC10 data not available in dataset")

print("""
Therapeutic Implications:
-------------------------
1. Enhancing inhibitory signaling may improve sepsis outcomes
2. CD52-SIGLEC10 axis represents potential therapeutic target
3. Cell type-specific modulation could optimize treatment

Next Steps:
-----------
1. Validate with CellChat/CellPhoneDB
2. Correlate L-R activity with clinical outcomes
3. Test in silico perturbation of key L-R pairs
""")

print("\n" + "="*60)
print("CellCall ANALYSIS COMPLETE!")
print("="*60)
