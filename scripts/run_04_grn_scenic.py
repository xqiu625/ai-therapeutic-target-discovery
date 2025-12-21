#!/usr/bin/env python
"""
04 - Gene Regulatory Network Analysis (SCENIC)
Sepsis Target Discovery Pipeline

Performs GRN inference using pySCENIC:
1. GRNBoost2 - TF-target inference
2. cisTarget - Motif-based pruning
3. AUCell - Regulon activity scoring

Note: Run on HPCC with pyscenic_env
"""

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import warnings
import pickle
warnings.filterwarnings('ignore')

# SCENIC imports
try:
    from pyscenic.utils import modules_from_adjacencies
    from pyscenic.prune import prune2df, df2regulons
    from pyscenic.aucell import aucell
    from arboreto.algo import grnboost2
    from ctxcore.rnkdb import FeatherRankingDatabase
    SCENIC_AVAILABLE = True
    print("pySCENIC imported successfully")
except ImportError as e:
    print(f"pySCENIC not available: {e}")
    SCENIC_AVAILABLE = False

import networkx as nx
from scipy.stats import mannwhitneyu
from sklearn.preprocessing import MinMaxScaler

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
EXTERNAL_DIR = BASE_DIR / 'data' / 'external'
GRN_DIR = BASE_DIR / 'data' / 'grn'
RESULTS_DIR = BASE_DIR / 'results'
TABLES_DIR = RESULTS_DIR / 'tables'
FIGURES_DIR = BASE_DIR / 'figures'

GRN_DIR.mkdir(parents=True, exist_ok=True)
TABLES_DIR.mkdir(parents=True, exist_ok=True)

# Number of workers for parallel processing
N_WORKERS = 16

print("="*60)
print("04 - GRN/SCENIC ANALYSIS")
print("="*60)

# ============================================================
# 1. Load Data
# ============================================================
print("\n" + "="*60)
print("1. LOADING DATA")
print("="*60)

# Try adata_embeddings first, fall back to adata_processed
adata_file = PROCESSED_DIR / 'adata_embeddings.h5ad'
if not adata_file.exists():
    adata_file = PROCESSED_DIR / 'adata_processed.h5ad'

adata = sc.read_h5ad(adata_file)
print(f"Loaded: {adata.n_obs} cells, {adata.n_vars} genes")
print(f"From: {adata_file}")

# ============================================================
# 2. Load TF List
# ============================================================
print("\n" + "="*60)
print("2. LOADING TRANSCRIPTION FACTORS")
print("="*60)

tf_file = EXTERNAL_DIR / 'hs_hgnc_tfs.txt'

# Try multiple sources for TF list
TF_URLS = [
    "https://resources.aertslab.org/cistarget/tf_lists/allTFs_hg38.txt",
    "https://raw.githubusercontent.com/aertslab/pySCENIC/master/resources/hs_hgnc_tfs.txt",
]

# Download if missing or empty
if not tf_file.exists() or tf_file.stat().st_size == 0:
    import urllib.request
    for url in TF_URLS:
        print(f"Trying: {url}...")
        try:
            urllib.request.urlretrieve(url, tf_file)
            if tf_file.stat().st_size > 0:
                print(f"Downloaded successfully!")
                break
        except Exception as e:
            print(f"Failed: {e}")

# Load TF list
if tf_file.exists() and tf_file.stat().st_size > 0:
    tf_names = pd.read_csv(tf_file, header=None)[0].tolist()
    tf_names = [tf for tf in tf_names if tf in adata.var_names]
    print(f"TFs in dataset: {len(tf_names)} / {len(pd.read_csv(tf_file, header=None))} total")
else:
    # Fallback: use common TF gene name patterns
    print("Using fallback: detecting TFs by gene name patterns...")
    # Common TF families
    tf_patterns = ['ZNF', 'KLF', 'SOX', 'FOX', 'GATA', 'IRF', 'STAT', 'NFE', 'NFK',
                   'ETV', 'ELF', 'ETS', 'PAX', 'HOX', 'POU', 'RUNX', 'CEBP', 'ATF',
                   'JUN', 'FOS', 'MAF', 'BACH', 'NRF', 'RELA', 'RELB', 'REL', 'MYC',
                   'MAX', 'MYB', 'E2F', 'TP53', 'TP63', 'TP73', 'NFAT', 'TFAP', 'TCF',
                   'LEF', 'SMAD', 'NOTCH', 'HIF', 'ARNT', 'EPAS', 'SP1', 'SP2', 'SP3']
    tf_names = [g for g in adata.var_names if any(g.startswith(p) for p in tf_patterns)]
    print(f"Found {len(tf_names)} potential TFs by pattern matching")

# ============================================================
# 3. Prepare Expression Matrix
# ============================================================
print("\n" + "="*60)
print("3. PREPARING EXPRESSION MATRIX")
print("="*60)

# Use raw counts
if 'counts' in adata.layers:
    expression_matrix = adata.layers['counts'].copy()
    print("Using counts layer")
else:
    expression_matrix = adata.X.copy()
    print("Using X matrix (assuming counts)")

if hasattr(expression_matrix, 'toarray'):
    expression_matrix = expression_matrix.toarray()

# Create DataFrame
ex_matrix = pd.DataFrame(
    expression_matrix,
    index=adata.obs_names,
    columns=adata.var_names
)

print(f"Expression matrix: {ex_matrix.shape}")

# ============================================================
# 4. GRNBoost2 - TF-Target Inference
# ============================================================
print("\n" + "="*60)
print("4. GRNBoost2 - TF-TARGET INFERENCE")
print("="*60)

adj_file = GRN_DIR / 'adjacencies.csv'

if adj_file.exists():
    print(f"Loading cached adjacencies from {adj_file}")
    adjacencies = pd.read_csv(adj_file)
elif SCENIC_AVAILABLE and len(tf_names) > 0:
    print("Running GRNBoost2... (this may take hours)")
    print(f"Using {N_WORKERS} workers")

    adjacencies = grnboost2(
        expression_data=ex_matrix,
        tf_names=tf_names,
        verbose=True,
        seed=42
    )

    adjacencies.to_csv(adj_file, index=False)
    print(f"Saved adjacencies to {adj_file}")
else:
    print("Cannot run GRNBoost2: SCENIC not available or no TFs found")
    adjacencies = pd.DataFrame()

if not adjacencies.empty:
    print(f"Adjacencies: {len(adjacencies)} TF-target pairs")
    print(f"Unique TFs: {adjacencies['TF'].nunique()}")
    print(f"\nTop 10 TF-target pairs by importance:")
    print(adjacencies.nlargest(10, 'importance').to_string())

# ============================================================
# 5. cisTarget - Motif Pruning
# ============================================================
print("\n" + "="*60)
print("5. CISTARGET - MOTIF PRUNING")
print("="*60)

db_files = list(EXTERNAL_DIR.glob('*.feather'))
print(f"Found {len(db_files)} database files:")
for f in db_files:
    print(f"  - {f.name}")

regulons_file = GRN_DIR / 'regulons.pkl'

if regulons_file.exists():
    print(f"Loading cached regulons from {regulons_file}")
    with open(regulons_file, 'rb') as f:
        regulons = pickle.load(f)
elif SCENIC_AVAILABLE and len(db_files) > 0 and not adjacencies.empty:
    print("Running cisTarget pruning...")

    dbs = [FeatherRankingDatabase(fname=str(f)) for f in db_files]

    motif_files = list(EXTERNAL_DIR.glob('*motif*.tbl'))
    if motif_files:
        motif_annotations = str(motif_files[0])
        print(f"Using motif file: {motif_annotations}")
    else:
        print("Motif annotation file not found!")
        motif_annotations = None

    modules = list(modules_from_adjacencies(adjacencies, ex_matrix))
    print(f"Created {len(modules)} modules")

    if motif_annotations:
        df = prune2df(dbs, modules, motif_annotations)
        regulons = df2regulons(df)

        with open(regulons_file, 'wb') as f:
            pickle.dump(regulons, f)
        print(f"Saved {len(regulons)} regulons")
    else:
        regulons = []
else:
    print("Cannot run cisTarget: missing dependencies or databases")
    print("Please download from: https://resources.aertslab.org/cistarget/databases/")
    regulons = []

if regulons:
    print(f"Total regulons: {len(regulons)}")
    print(f"\nTop 10 regulons by size:")
    regulon_sizes = [(r.name, len(r.genes)) for r in regulons]
    regulon_sizes = sorted(regulon_sizes, key=lambda x: x[1], reverse=True)[:10]
    for name, size in regulon_sizes:
        print(f"  {name}: {size} target genes")

# ============================================================
# 6. AUCell - Regulon Activity Scoring
# ============================================================
print("\n" + "="*60)
print("6. AUCELL - REGULON ACTIVITY")
print("="*60)

auc_file = GRN_DIR / 'auc_matrix.csv'

if auc_file.exists():
    print(f"Loading cached AUCell scores from {auc_file}")
    auc_mtx = pd.read_csv(auc_file, index_col=0)
elif SCENIC_AVAILABLE and regulons:
    print(f"Running AUCell with {N_WORKERS} workers...")

    auc_mtx = aucell(ex_matrix, regulons, num_workers=N_WORKERS)

    auc_mtx.to_csv(auc_file)
    print(f"Saved AUCell scores to {auc_file}")
else:
    print("Cannot run AUCell: missing dependencies")
    auc_mtx = pd.DataFrame()

if not auc_mtx.empty:
    print(f"AUCell matrix: {auc_mtx.shape}")
    adata.obsm['X_aucell'] = auc_mtx.values

    for regulon_name in auc_mtx.columns:
        adata.obs[f'AUC_{regulon_name}'] = auc_mtx[regulon_name].values

# ============================================================
# 7. Differential Regulon Activity
# ============================================================
print("\n" + "="*60)
print("7. DIFFERENTIAL REGULON ACTIVITY")
print("="*60)

if not auc_mtx.empty and 'outcome' in adata.obs.columns:
    adata_sepsis = adata[adata.obs['disease'] == 'Sepsis'].copy()

    if 'Survivor' in adata_sepsis.obs['outcome'].values and 'Non-survivor' in adata_sepsis.obs['outcome'].values:
        survivor_mask = adata_sepsis.obs['outcome'] == 'Survivor'

        results = []
        for regulon in auc_mtx.columns:
            col = f'AUC_{regulon}'
            if col in adata_sepsis.obs.columns:
                surv_vals = adata_sepsis.obs.loc[survivor_mask, col]
                nonsurv_vals = adata_sepsis.obs.loc[~survivor_mask, col]

                stat, pval = mannwhitneyu(surv_vals, nonsurv_vals, alternative='two-sided')
                fc = surv_vals.mean() / (nonsurv_vals.mean() + 1e-10)

                results.append({
                    'regulon': regulon,
                    'mean_survivor': surv_vals.mean(),
                    'mean_nonsurvivor': nonsurv_vals.mean(),
                    'fold_change': fc,
                    'pvalue': pval
                })

        diff_regulons = pd.DataFrame(results)
        diff_regulons['padj'] = diff_regulons['pvalue'] * len(diff_regulons)
        diff_regulons = diff_regulons.sort_values('pvalue')

        print("Top 20 differential regulons (Survivor vs Non-survivor):")
        print(diff_regulons.head(20).to_string())

        diff_regulons.to_csv(TABLES_DIR / 'differential_regulons.csv', index=False)
        print("\nSaved differential regulons")
    else:
        print("Need both Survivor and Non-survivor for comparison")
        diff_regulons = pd.DataFrame()
else:
    print("No AUCell matrix or outcome column")
    diff_regulons = pd.DataFrame()

# ============================================================
# 8. Network Centrality Analysis
# ============================================================
print("\n" + "="*60)
print("8. NETWORK CENTRALITY ANALYSIS")
print("="*60)

if not adjacencies.empty:
    top_adj = adjacencies.nlargest(5000, 'importance')

    G = nx.from_pandas_edgelist(
        top_adj,
        source='TF',
        target='target',
        edge_attr='importance',
        create_using=nx.DiGraph()
    )

    print(f"Network: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")

    # Calculate centrality
    pagerank = nx.pagerank(G)
    in_degree = dict(G.in_degree())
    out_degree = dict(G.out_degree())
    betweenness = nx.betweenness_centrality(G)

    all_genes = set(G.nodes())
    centrality_df = pd.DataFrame({
        'gene': list(all_genes),
        'pagerank': [pagerank.get(g, 0) for g in all_genes],
        'in_degree': [in_degree.get(g, 0) for g in all_genes],
        'out_degree': [out_degree.get(g, 0) for g in all_genes],
        'betweenness': [betweenness.get(g, 0) for g in all_genes]
    })

    scaler = MinMaxScaler()
    centrality_df['centrality_score'] = scaler.fit_transform(
        centrality_df[['pagerank', 'in_degree', 'betweenness']]
    ).mean(axis=1)

    centrality_df = centrality_df.sort_values('centrality_score', ascending=False)

    print("\nTop 20 genes by network centrality:")
    print(centrality_df.head(20).to_string())

    centrality_df.to_csv(TABLES_DIR / 'grn_centrality.csv', index=False)
    print("\nSaved centrality scores")

    # Check CD52
    cd52_centrality = centrality_df[centrality_df['gene'] == 'CD52']
    if not cd52_centrality.empty:
        print(f"\n*** CD52 Network Metrics ***")
        print(cd52_centrality.to_string())

        cd52_regulators = adjacencies[adjacencies['target'] == 'CD52'].sort_values('importance', ascending=False)
        print(f"\nTop TFs regulating CD52:")
        print(cd52_regulators.head(10).to_string())
    else:
        print("\nCD52 not in top network edges")
else:
    centrality_df = pd.DataFrame()

# ============================================================
# 9. Network Visualization
# ============================================================
print("\n" + "="*60)
print("9. NETWORK VISUALIZATION")
print("="*60)

if not adjacencies.empty:
    top_tfs = centrality_df[centrality_df['out_degree'] > 0].head(20)['gene'].tolist()
    sub_adj = adjacencies[adjacencies['TF'].isin(top_tfs)].nlargest(200, 'importance')

    plt.figure(figsize=(14, 14))

    sub_G = nx.from_pandas_edgelist(
        sub_adj,
        source='TF',
        target='target',
        edge_attr='importance',
        create_using=nx.DiGraph()
    )

    pos = nx.spring_layout(sub_G, k=2, iterations=50, seed=42)

    tfs_in_net = set(sub_adj['TF'])
    node_colors = []
    node_sizes = []
    for node in sub_G.nodes():
        if node == 'CD52':
            node_colors.append('red')
            node_sizes.append(500)
        elif node in tfs_in_net:
            node_colors.append('orange')
            node_sizes.append(300)
        else:
            node_colors.append('lightblue')
            node_sizes.append(100)

    nx.draw(sub_G, pos, node_color=node_colors, node_size=node_sizes,
            with_labels=True, font_size=8, arrows=True, alpha=0.7)

    plt.title('Top TF-Target Network (Top 20 TFs)')
    plt.savefig(FIGURES_DIR / 'fig4_grn_network.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Saved GRN network plot")

# ============================================================
# 10. Save Results
# ============================================================
print("\n" + "="*60)
print("10. SAVING RESULTS")
print("="*60)

adata.write(PROCESSED_DIR / 'adata_grn.h5ad')
print(f"Saved adata_grn.h5ad")

# ============================================================
# Summary
# ============================================================
print("\n" + "="*60)
print("GRN ANALYSIS SUMMARY")
print("="*60)
if not adjacencies.empty:
    print(f"TF-target pairs: {len(adjacencies)}")
    print(f"Unique TFs: {adjacencies['TF'].nunique()}")
if regulons:
    print(f"Regulons: {len(regulons)}")
if not auc_mtx.empty:
    print(f"AUCell matrix: {auc_mtx.shape}")
if not centrality_df.empty:
    print(f"\nTop 5 central genes:")
    for _, row in centrality_df.head(5).iterrows():
        marker = " *** CD52 ***" if row['gene'] == 'CD52' else ""
        print(f"  {row['gene']}: centrality={row['centrality_score']:.4f}{marker}")

print("\n" + "="*60)
print("GRN ANALYSIS COMPLETE!")
print("="*60)
