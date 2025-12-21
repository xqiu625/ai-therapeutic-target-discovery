#!/usr/bin/env python
"""
08 - CMap/LINCS Drug Signature Query
Sepsis Target Discovery Pipeline

Queries Connectivity Map (CMap) L1000 to find drugs that:
1. REVERSE the sepsis poor-outcome (non-survivor) signature
2. MIMIC the sepsis good-outcome (survivor) signature

This identifies drug repurposing candidates for sepsis treatment.
"""

import pandas as pd
import numpy as np
import requests
import json
import time
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
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
print("08 - CMap/LINCS DRUG QUERY")
print("="*60)

# ============================================================
# 1. Load DE Results
# ============================================================
print("\n" + "="*60)
print("1. LOADING DE RESULTS")
print("="*60)

de_file = PROCESSED_DIR / 'de_results_outcome.csv'
de_results = pd.read_csv(de_file)
print(f"Loaded {len(de_results)} DE genes")
print(de_results.head())

# Separate upregulated and downregulated genes
# In survivors vs non-survivors:
# - Positive logFC = higher in survivors (good outcome genes)
# - Negative logFC = higher in non-survivors (bad outcome genes)

# Filter significant genes
sig_de = de_results[de_results['pvals_adj'] < 0.05].copy()
print(f"\nSignificant genes (padj < 0.05): {len(sig_de)}")

# Get top upregulated (higher in survivors) and downregulated genes
up_genes = sig_de[sig_de['logfoldchanges'] > 0.5].nlargest(150, 'logfoldchanges')['names'].tolist()
down_genes = sig_de[sig_de['logfoldchanges'] < -0.5].nsmallest(150, 'logfoldchanges')['names'].tolist()

print(f"Upregulated in survivors (top 150): {len(up_genes)}")
print(f"Downregulated in survivors (top 150): {len(down_genes)}")

# Check CD52
cd52_de = de_results[de_results['names'] == 'CD52']
if not cd52_de.empty:
    cd52_lfc = cd52_de['logfoldchanges'].values[0]
    print(f"\nCD52 log2FC: {cd52_lfc:.3f} ({'UP in survivors' if cd52_lfc > 0 else 'DOWN in survivors'})")

# ============================================================
# 2. L1000FWD API Query
# ============================================================
print("\n" + "="*60)
print("2. QUERYING L1000FWD API")
print("="*60)

L1000FWD_URL = "https://maayanlab.cloud/L1000FWD"

def query_l1000fwd(up_genes, down_genes, description="sepsis_signature"):
    """
    Query L1000FWD for drugs that reverse/mimic a gene signature.
    L1000FWD uses a two-step process: submit query, then retrieve results.
    """
    payload = {
        "up_genes": up_genes[:100],  # API limit
        "down_genes": down_genes[:100],
    }

    try:
        print(f"Sending query with {len(payload['up_genes'])} up and {len(payload['down_genes'])} down genes...")

        # Step 1: Submit signature search
        response = requests.post(f"{L1000FWD_URL}/sig_search", json=payload, timeout=60)

        if response.status_code == 200:
            result = response.json()
            print(f"Response keys: {list(result.keys())}")

            # Check if we got a search_id (async result)
            if 'search_id' in result:
                search_id = result['search_id']
                print(f"Got search_id: {search_id}, retrieving results...")
                time.sleep(2)  # Wait for processing

                # Step 2: Retrieve results
                result_response = requests.get(f"{L1000FWD_URL}/result/{search_id}", timeout=60)
                if result_response.status_code == 200:
                    result = result_response.json()
                    print(f"Result keys: {list(result.keys())}")

            # Print sample of result structure for debugging
            for key in list(result.keys())[:5]:
                val = result[key]
                if isinstance(val, list) and len(val) > 0:
                    print(f"  {key}: list with {len(val)} items, first item keys: {list(val[0].keys()) if isinstance(val[0], dict) else type(val[0])}")
                elif isinstance(val, dict):
                    print(f"  {key}: dict with keys {list(val.keys())[:5]}")
                else:
                    print(f"  {key}: {type(val).__name__}")

            return result
        else:
            print(f"Error: {response.status_code}")
            print(response.text[:500])
            return None
    except Exception as e:
        print(f"Request failed: {e}")
        return None

# Query 1: Find drugs that MIMIC survivor signature (therapeutic)
print("\n--- Query 1: Drugs mimicking SURVIVOR signature ---")
print("(These drugs push cells toward good-outcome state)")
mimic_result = query_l1000fwd(up_genes, down_genes, "mimic_survivor_signature")

# Query 2: Find drugs that REVERSE non-survivor signature
# Swap up/down to find drugs that do the opposite of poor outcome
print("\n--- Query 2: Drugs REVERSING NON-SURVIVOR signature ---")
print("(These drugs counter the bad-outcome gene program)")
reverse_result = query_l1000fwd(down_genes, up_genes, "reverse_nonsurvivor_signature")

# ============================================================
# 3. Process Results
# ============================================================
print("\n" + "="*60)
print("3. PROCESSING RESULTS")
print("="*60)

def process_l1000_results(result, query_type):
    """Extract drug information from L1000FWD results."""
    if result is None:
        return pd.DataFrame()

    drugs = []

    # Try multiple possible response formats
    possible_keys = ['similar', 'signatures', 'topMeta', 'result', 'results']

    for key in possible_keys:
        if key in result and isinstance(result[key], list):
            print(f"  Processing '{key}' with {len(result[key])} entries...")
            for sig in result[key][:100]:
                if isinstance(sig, dict):
                    # Try multiple field name variations
                    pert_name = (sig.get('pert_name') or sig.get('pert_iname') or
                                sig.get('drug_name') or sig.get('name') or
                                sig.get('Drug') or sig.get('perturbagen') or '')
                    drugs.append({
                        'query_type': query_type,
                        'sig_id': sig.get('sig_id', sig.get('id', '')),
                        'pert_name': pert_name,
                        'pert_type': sig.get('pert_type', sig.get('type', '')),
                        'cell_line': sig.get('cell_id', sig.get('cell_line', sig.get('cell', ''))),
                        'similarity_score': sig.get('score', sig.get('similarity', sig.get('combined_score', sig.get('zscore', 0)))),
                        'pert_dose': sig.get('pert_dose', sig.get('dose', '')),
                        'pert_time': sig.get('pert_time', sig.get('time', ''))
                    })

    # If no drugs found, try treating entire result as list
    if len(drugs) == 0 and isinstance(result, list):
        print(f"  Processing top-level list with {len(result)} entries...")
        for sig in result[:100]:
            if isinstance(sig, dict):
                pert_name = (sig.get('pert_name') or sig.get('pert_iname') or
                            sig.get('drug_name') or sig.get('name') or '')
                if pert_name:
                    drugs.append({
                        'query_type': query_type,
                        'sig_id': sig.get('sig_id', ''),
                        'pert_name': pert_name,
                        'pert_type': sig.get('pert_type', ''),
                        'cell_line': sig.get('cell_id', ''),
                        'similarity_score': sig.get('score', 0),
                        'pert_dose': sig.get('pert_dose', ''),
                        'pert_time': sig.get('pert_time', '')
                    })

    return pd.DataFrame(drugs)

mimic_drugs = process_l1000_results(mimic_result, 'mimic_survivor')
reverse_drugs = process_l1000_results(reverse_result, 'reverse_nonsurvivor')

print(f"\nMimic survivor results: {len(mimic_drugs)}")
print(f"Reverse non-survivor results: {len(reverse_drugs)}")

# Combine results
all_drugs = pd.concat([mimic_drugs, reverse_drugs], ignore_index=True)

if len(all_drugs) > 0:
    # Filter to small molecules (drugs)
    drug_results = all_drugs[all_drugs['pert_type'].isin(['trt_cp', 'ctl_vehicle', ''])].copy()

    # Remove empty names
    drug_results = drug_results[drug_results['pert_name'].str.len() > 0]

    # Aggregate by drug name (take best score)
    drug_summary = drug_results.groupby('pert_name').agg({
        'similarity_score': 'max',
        'query_type': 'first',
        'cell_line': lambda x: ', '.join(x.unique()[:3]),
        'pert_dose': 'first'
    }).reset_index()

    drug_summary = drug_summary.sort_values('similarity_score', ascending=False)

    print(f"\nUnique drugs found: {len(drug_summary)}")
    print("\nTop 20 drug candidates:")
    print(drug_summary.head(20).to_string())
else:
    drug_summary = pd.DataFrame()
    print("\nNo drugs found from L1000FWD query")

# ============================================================
# 4. Alternative: SigCom LINCS Query
# ============================================================
print("\n" + "="*60)
print("4. SIGCOM LINCS QUERY (Alternative)")
print("="*60)

SIGCOM_URL = "https://maayanlab.cloud/sigcom-lincs/api/v1/enrichment/overlap"

def query_sigcom(up_genes, down_genes, top_n=50):
    """Query SigCom LINCS for drug enrichment."""
    payload = {
        "up": up_genes[:100],
        "down": down_genes[:100],
        "filterTerm": "CP",  # Filter to compounds
        "offset": 0,
        "limit": top_n
    }

    try:
        response = requests.post(SIGCOM_URL, json=payload, timeout=60)
        if response.status_code == 200:
            return response.json()
        else:
            print(f"SigCom error: {response.status_code}")
            return None
    except Exception as e:
        print(f"SigCom request failed: {e}")
        return None

# Query SigCom
print("Querying SigCom LINCS...")
sigcom_result = query_sigcom(up_genes, down_genes)

sigcom_drugs = []
if sigcom_result and 'results' in sigcom_result:
    for r in sigcom_result['results'][:50]:
        sigcom_drugs.append({
            'pert_name': r.get('pert', ''),
            'cell_line': r.get('cell', ''),
            'score': r.get('score', 0),
            'pvalue': r.get('pvalue', 1),
            'direction': r.get('direction', '')
        })

sigcom_df = pd.DataFrame(sigcom_drugs)
if len(sigcom_df) > 0:
    print(f"\nSigCom drugs found: {len(sigcom_df)}")
    print(sigcom_df.head(10).to_string())

# ============================================================
# 4b. Alternative: Enrichr Drug Signatures
# ============================================================
print("\n" + "="*60)
print("4b. ENRICHR DRUG SIGNATURES (Alternative)")
print("="*60)

ENRICHR_URL = "https://maayanlab.cloud/Enrichr"

def query_enrichr_drugs(gene_list, libraries=['DSigDB', 'Drug_Perturbations_from_GEO_down',
                                               'Drug_Perturbations_from_GEO_up',
                                               'LINCS_L1000_Chem_Pert_Consensus_Sigs']):
    """Query Enrichr for drug-related gene set enrichment."""
    all_results = []

    try:
        # Step 1: Add gene list
        add_url = f"{ENRICHR_URL}/addList"
        payload = {
            'list': (None, '\n'.join(gene_list)),
            'description': (None, 'sepsis_signature')
        }
        response = requests.post(add_url, files=payload, timeout=30)
        if response.status_code != 200:
            print(f"Enrichr addList error: {response.status_code}")
            return pd.DataFrame()

        result = response.json()
        user_list_id = result.get('userListId')
        print(f"Enrichr list ID: {user_list_id}")

        # Step 2: Query each library
        for library in libraries:
            enrich_url = f"{ENRICHR_URL}/enrich"
            params = {
                'userListId': user_list_id,
                'backgroundType': library
            }
            response = requests.get(enrich_url, params=params, timeout=30)

            if response.status_code == 200:
                data = response.json()
                if library in data:
                    for term in data[library][:20]:
                        # Enrichr format: [rank, term, pval, zscore, combined_score, genes, adj_pval, ...]
                        all_results.append({
                            'library': library,
                            'term': term[1] if len(term) > 1 else '',
                            'pvalue': term[2] if len(term) > 2 else 1,
                            'zscore': term[3] if len(term) > 3 else 0,
                            'combined_score': term[4] if len(term) > 4 else 0,
                            'genes': ';'.join(term[5]) if len(term) > 5 else '',
                            'adj_pvalue': term[6] if len(term) > 6 else 1
                        })
                    print(f"  {library}: {len(data[library])} results")

    except Exception as e:
        print(f"Enrichr query failed: {e}")

    return pd.DataFrame(all_results)

# Query Enrichr with up-regulated genes (survivor signature)
print("\nQuerying Enrichr with survivor UP-regulated genes...")
enrichr_up = query_enrichr_drugs(up_genes[:100])

print("\nQuerying Enrichr with survivor DOWN-regulated genes...")
enrichr_down = query_enrichr_drugs(down_genes[:100])

enrichr_drugs = pd.concat([enrichr_up, enrichr_down], ignore_index=True)
if len(enrichr_drugs) > 0:
    enrichr_drugs = enrichr_drugs.sort_values('combined_score', ascending=False)
    print(f"\nEnrichr drug terms found: {len(enrichr_drugs)}")
    print("\nTop 15 Enrichr drug hits:")
    print(enrichr_drugs.head(15)[['library', 'term', 'combined_score', 'pvalue']].to_string())

# ============================================================
# 5. Known Drug Annotations
# ============================================================
print("\n" + "="*60)
print("5. ANNOTATING KNOWN DRUGS")
print("="*60)

# FDA-approved drugs relevant to sepsis/inflammation
known_drugs = {
    'dexamethasone': 'Corticosteroid, anti-inflammatory, used in sepsis',
    'hydrocortisone': 'Corticosteroid, used in septic shock',
    'methylprednisolone': 'Corticosteroid, anti-inflammatory',
    'ibuprofen': 'NSAID, anti-inflammatory',
    'aspirin': 'NSAID, anti-inflammatory',
    'tocilizumab': 'IL-6 inhibitor, used in cytokine storm',
    'anakinra': 'IL-1 receptor antagonist',
    'baricitinib': 'JAK inhibitor, anti-inflammatory',
    'ruxolitinib': 'JAK inhibitor',
    'tofacitinib': 'JAK inhibitor',
    'sirolimus': 'mTOR inhibitor, immunomodulator',
    'everolimus': 'mTOR inhibitor',
    'cyclosporine': 'Calcineurin inhibitor, immunosuppressant',
    'tacrolimus': 'Calcineurin inhibitor',
    'metformin': 'Diabetes drug with anti-inflammatory effects',
    'statins': 'Anti-inflammatory pleiotropic effects',
    'vitamin-d': 'Immunomodulator',
    'n-acetylcysteine': 'Antioxidant, used in sepsis trials',
    'pentoxifylline': 'TNF inhibitor, hemorheologic',
    'interferon': 'Immunomodulator',
}

# Check if any known drugs appear in results
if len(drug_summary) > 0:
    drug_summary['known_annotation'] = drug_summary['pert_name'].str.lower().map(
        lambda x: next((v for k, v in known_drugs.items() if k in x.lower()), '')
    )

    known_hits = drug_summary[drug_summary['known_annotation'] != '']
    if len(known_hits) > 0:
        print("\nKnown drugs in results:")
        print(known_hits[['pert_name', 'similarity_score', 'known_annotation']].to_string())

# ============================================================
# 6. CD52-Related Drug Search
# ============================================================
print("\n" + "="*60)
print("6. CD52-RELATED DRUGS")
print("="*60)

# Drugs known to affect CD52 or related pathways
cd52_related = {
    'alemtuzumab': 'Anti-CD52 antibody (OPPOSITE of what we want)',
    'siglec': 'CD52 binds Siglec-10',
    'gpi-anchor': 'CD52 is GPI-anchored',
}

print("""
CD52 Therapeutic Context:
------------------------
- CD52 is HIGHER in survivors (protective)
- Soluble CD52 is anti-inflammatory
- CD52 binds Siglec-10 to suppress T cells

Therapeutic Strategy:
- DO NOT use anti-CD52 (alemtuzumab) - this depletes CD52+ cells
- LOOK FOR drugs that:
  1. Upregulate CD52 expression
  2. Mimic soluble CD52 effects
  3. Activate Siglec-10 pathway
  4. Have similar anti-inflammatory mechanisms
""")

# ============================================================
# 7. Visualizations
# ============================================================
print("\n" + "="*60)
print("7. VISUALIZATIONS")
print("="*60)

if len(drug_summary) > 0:
    # Top drugs bar plot
    fig, ax = plt.subplots(figsize=(12, 10))
    top_drugs = drug_summary.head(30)

    colors = ['green' if a != '' else 'steelblue'
              for a in top_drugs.get('known_annotation', [''] * len(top_drugs))]

    ax.barh(range(len(top_drugs)), top_drugs['similarity_score'], color=colors)
    ax.set_yticks(range(len(top_drugs)))
    ax.set_yticklabels(top_drugs['pert_name'])
    ax.invert_yaxis()
    ax.set_xlabel('CMap Similarity Score')
    ax.set_title('Top 30 Drug Candidates from CMap Query\n(Green = Known anti-inflammatory)')

    plt.tight_layout()
    plt.savefig(FIGURES_DIR / 'fig8_cmap_drugs.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Saved CMap drug plot")

# Enrichr results plot
if len(enrichr_drugs) > 0:
    fig, ax = plt.subplots(figsize=(14, 10))
    top_enrichr = enrichr_drugs.head(25)

    # Color by library
    lib_colors = {'DSigDB': 'steelblue', 'Drug_Perturbations_from_GEO_down': 'coral',
                  'Drug_Perturbations_from_GEO_up': 'seagreen',
                  'LINCS_L1000_Chem_Pert_Consensus_Sigs': 'purple'}
    colors = [lib_colors.get(lib, 'gray') for lib in top_enrichr['library']]

    ax.barh(range(len(top_enrichr)), top_enrichr['combined_score'], color=colors)
    ax.set_yticks(range(len(top_enrichr)))

    # Truncate long term names
    labels = [t[:50] + '...' if len(t) > 50 else t for t in top_enrichr['term']]
    ax.set_yticklabels(labels, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel('Enrichr Combined Score')
    ax.set_title('Top 25 Drug Signatures from Enrichr\n(Blue=DSigDB, Coral=GEO_down, Green=GEO_up, Purple=LINCS)')

    plt.tight_layout()
    plt.savefig(FIGURES_DIR / 'fig8_enrichr_drugs.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Saved Enrichr drug plot")

# ============================================================
# 8. Save Results
# ============================================================
print("\n" + "="*60)
print("8. SAVING RESULTS")
print("="*60)

# Save all results
if len(drug_summary) > 0:
    drug_summary.to_csv(TABLES_DIR / 'cmap_drug_candidates.csv', index=False)
    print(f"Saved {len(drug_summary)} L1000FWD drug candidates")

if len(sigcom_df) > 0:
    sigcom_df.to_csv(TABLES_DIR / 'sigcom_drug_candidates.csv', index=False)
    print(f"Saved {len(sigcom_df)} SigCom results")

if len(enrichr_drugs) > 0:
    enrichr_drugs.to_csv(TABLES_DIR / 'enrichr_drug_candidates.csv', index=False)
    print(f"Saved {len(enrichr_drugs)} Enrichr drug results")

# Save query signature for reproducibility
query_sig = pd.DataFrame({
    'gene': up_genes + down_genes,
    'direction': ['up'] * len(up_genes) + ['down'] * len(down_genes)
})
query_sig.to_csv(TABLES_DIR / 'cmap_query_signature.csv', index=False)
print("Saved query signature")

# ============================================================
# Summary
# ============================================================
print("\n" + "="*60)
print("CMAP QUERY SUMMARY")
print("="*60)
print(f"""
Query Parameters:
- Up-regulated genes (survivors): {len(up_genes)}
- Down-regulated genes (survivors): {len(down_genes)}

Results:
- L1000FWD drugs: {len(drug_summary) if len(drug_summary) > 0 else 'No results'}
- SigCom drugs: {len(sigcom_df)}
- Enrichr drug terms: {len(enrichr_drugs)}

Top Drug Candidates:
""")

if len(drug_summary) > 0:
    print("From L1000FWD:")
    for i, (_, row) in enumerate(drug_summary.head(10).iterrows()):
        print(f"  {i+1}. {row['pert_name']}: score={row['similarity_score']:.3f}")
elif len(sigcom_df) > 0:
    print("From SigCom:")
    for i, (_, row) in enumerate(sigcom_df.head(10).iterrows()):
        print(f"  {i+1}. {row['pert_name']}: score={row['score']:.3f}")

if len(enrichr_drugs) > 0:
    print("\nFrom Enrichr (drug signature databases):")
    for i, (_, row) in enumerate(enrichr_drugs.head(10).iterrows()):
        print(f"  {i+1}. {row['term']}: score={row['combined_score']:.1f} (p={row['pvalue']:.2e})")

print(f"""
Next Steps:
1. Validate top candidates in literature
2. Check for sepsis/inflammation clinical trials
3. Cross-reference with DGIdb druggability
4. Prioritize FDA-approved drugs for repurposing

Key CD52 Insight:
- Look for drugs that ENHANCE immune regulation
- Avoid drugs that deplete lymphocytes
- Focus on anti-inflammatory mechanisms matching CD52 pathway
""")

print("\n" + "="*60)
print("CMAP QUERY COMPLETE!")
print("="*60)
