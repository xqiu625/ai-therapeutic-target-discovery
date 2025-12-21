# AI-Enabled Therapeutic Target Discovery Pipeline

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

An end-to-end computational pipeline for identifying and prioritizing therapeutic targets from single-cell RNA-seq data, integrating machine learning, gene regulatory networks, trajectory analysis, and drug repurposing.

## Overview

This pipeline demonstrates **AI-enabled target discovery** by combining multiple analytical approaches to identify high-confidence therapeutic targets:

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                    AI-Enabled Target Discovery Pipeline                      │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│   scRNA-seq ──► Preprocessing ──► Multi-Modal Analysis ──► Target Ranking  │
│                                          │                                  │
│                          ┌───────────────┼───────────────┐                  │
│                          ▼               ▼               ▼                  │
│                    ┌─────────┐    ┌───────────┐   ┌────────────┐           │
│                    │   ML    │    │    GRN    │   │ Trajectory │           │
│                    │ Models  │    │ (SCENIC)  │   │ (Pseudotime)│          │
│                    └────┬────┘    └─────┬─────┘   └──────┬─────┘           │
│                         │               │                │                  │
│                         └───────────────┼────────────────┘                  │
│                                         ▼                                   │
│                              ┌─────────────────────┐                        │
│                              │ Evidence Integration │                       │
│                              │  + Druggability      │                       │
│                              └──────────┬──────────┘                        │
│                                         ▼                                   │
│                              ┌─────────────────────┐                        │
│                              │  Drug Repurposing   │                        │
│                              │  (CMap/Enrichr)     │                        │
│                              └─────────────────────┘                        │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
```

## Case Study: Sepsis Therapeutic Target Discovery

**Objective:** Identify therapeutic targets that differentiate survivors from non-survivors in sepsis patients.

### Key Results

| Analysis | Method | Key Finding |
|----------|--------|-------------|
| **ML Prediction** | RF, GradientBoosting, LightGBM | AUC = 0.99 for outcome prediction |
| **Differential Expression** | Wilcoxon rank-sum | 6,014 significant genes (padj < 0.05) |
| **Target Prioritization** | Multi-evidence integration | Top target ranked #3/33,538 genes |
| **Drug Discovery** | Enrichr + CMap L1000 | 160 drug candidates identified |

### Validated Target: CD52

The pipeline successfully identified **CD52** as a high-priority therapeutic target:

- **ML Importance**: Ranked #12/3,000 highly variable genes
- **Fold Change**: 2.93x higher in survivors vs non-survivors (p < 0.001)
- **Final Priority Rank**: #3 out of 33,538 genes
- **Druggability**: Known drug exists (alemtuzumab), but therapeutic strategy requires *agonism* not antagonism

## Pipeline Components

### 1. Data Preparation (`run_01_data_prep.py`)
- Converts Seurat objects to AnnData format
- Quality control and normalization
- Differential expression analysis

### 2. Cell Type Annotation (`run_02_cell_annotation.py`)
- Automated cell type identification
- Marker gene validation
- UMAP visualization

### 3. ML Outcome Prediction (`run_03_ml_prediction.py`)
- **Models**: Random Forest, Gradient Boosting, LightGBM
- **Feature Selection**: Top 3,000 highly variable genes
- **Validation**: Patient-stratified train/test split (prevents data leakage)
- **Interpretation**: SHAP values for feature importance

### 4. Gene Regulatory Networks (`run_04_grn_scenic.py`)
- **GRNBoost2**: TF-target inference
- **cisTarget**: Motif-based regulon pruning
- **AUCell**: Regulon activity scoring
- **Network Analysis**: Centrality metrics for target prioritization

### 5. Trajectory Analysis (`run_05_trajectory.py`)
- Diffusion pseudotime computation
- Gene dynamics along disease progression
- Temporal expression patterns

### 6. Target Prioritization (`run_06_prioritization.py`)
Integrates multiple evidence sources with weighted scoring:

| Evidence Source | Weight | Description |
|-----------------|--------|-------------|
| Differential Expression | 25% | Statistical significance + effect size |
| ML Feature Importance | 25% | SHAP-based predictive power |
| GRN Centrality | 15% | Network hub importance |
| Trajectory Dynamics | 15% | Temporal expression patterns |
| Druggability (DGIdb) | 20% | Existing drug-gene interactions |

### 7. Target Deep Dive (`run_07_cd52_deep_dive.py`)
- Expression patterns across cell types
- Temporal dynamics (T0 vs T6)
- Upstream regulator identification
- Correlation network analysis

### 8. Drug Repurposing (`run_08_cmap_query.py`)
- **L1000FWD**: Connectivity Map signature matching
- **Enrichr**: Drug signature database queries
  - DSigDB
  - Drug_Perturbations_from_GEO
  - LINCS_L1000_Chem_Pert_Consensus_Sigs

## Installation

### Environment Setup

```bash
# Clone repository
git clone https://github.com/xqiu625/ai-therapeutic-target-discovery.git
cd ai-therapeutic-target-discovery

# Create conda environment
conda env create -f environment.yml
conda activate target-discovery

# Or install dependencies manually
pip install scanpy pandas numpy scikit-learn matplotlib seaborn
pip install pyscenic  # For GRN analysis
```

### Required Data

Place your data in the `data/` directory:
```
data/
├── raw/           # Raw count matrices
├── processed/     # Processed AnnData objects
└── external/      # External databases (TF lists, etc.)
```

## Usage

### Quick Start

```bash
# 1. Data preparation
python scripts/run_01_data_prep.py

# 2. ML prediction
python scripts/run_03_ml_prediction.py

# 3. Target prioritization
python scripts/run_06_prioritization.py

# 4. Drug repurposing
python scripts/run_08_cmap_query.py
```

### HPC Execution (SLURM)

```bash
# Submit individual jobs
sbatch scripts/hpcc/submit_ml_prediction.sh
sbatch scripts/hpcc/submit_grn_scenic.sh
sbatch scripts/hpcc/submit_trajectory.sh

# Or run the full pipeline
sbatch scripts/hpcc/submit_full_pipeline.sh
```

## Output Structure

```
results/
├── tables/
│   ├── prioritized_targets.csv      # Full ranked target list
│   ├── top50_targets.csv            # Top 50 candidates
│   ├── ml_feature_importance.csv    # ML-based gene rankings
│   ├── grn_centrality.csv           # Network centrality scores
│   ├── trajectory_correlations.csv  # Pseudotime correlations
│   └── enrichr_drug_candidates.csv  # Drug repurposing results
└── figures/
    ├── fig3_ml_performance.png      # ROC curves, confusion matrix
    ├── fig5_trajectory.png          # Pseudotime analysis
    └── fig6_target_ranking.png      # Priority score visualization

figures/
├── target_evidence_heatmap.png      # Multi-evidence summary
└── fig8_enrichr_drugs.png           # Drug candidates
```

## Key Dependencies

| Package | Version | Purpose |
|---------|---------|---------|
| scanpy | ≥1.9.0 | Single-cell analysis |
| pyscenic | ≥0.12.0 | Gene regulatory networks |
| scikit-learn | ≥1.0.0 | Machine learning models |
| pandas | ≥1.3.0 | Data manipulation |
| matplotlib/seaborn | - | Visualization |

## Customization

### Modify Prioritization Weights

Edit `run_06_prioritization.py`:
```python
weights = {
    'de_score_norm': 0.25,        # Differential expression
    'shap_importance_norm': 0.25,  # ML importance
    'centrality_score_norm': 0.15, # GRN centrality
    'trajectory_score_norm': 0.15, # Trajectory dynamics
    'druggability_score_norm': 0.20 # Drug-gene interactions
}
```

### Add Custom Evidence Sources

The pipeline is modular - add new evidence by:
1. Creating a scoring script that outputs `gene, score` CSV
2. Merging into the master dataframe in `run_06_prioritization.py`
3. Adding to the weighted scoring formula

## Citation

If you use this pipeline, please cite:

```bibtex
@software{qiu2024target,
  author = {Qiu, Xinru},
  title = {AI-Enabled Therapeutic Target Discovery Pipeline},
  year = {2024},
  url = {https://github.com/xqiu625/ai-therapeutic-target-discovery}
}
```

## License

MIT License - see [LICENSE](LICENSE) for details.

## Contact

**Xinru Qiu, PhD**
AI-enabled Target Discovery × Mechanistic Biology

- Website: [xqiu625.github.io](https://xqiu625.github.io)
- LinkedIn: [linkedin.com/in/xinru-qiu](https://linkedin.com/in/xinru-qiu)
- Email: xinru.qiu@email.com

---

*This pipeline demonstrates end-to-end therapeutic target discovery from single-cell transcriptomics, integrating machine learning, network biology, and drug repurposing approaches.*
