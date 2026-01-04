#!/usr/bin/env python3
"""
Script 19: CCTL (Causal-Contrastive Target Learning) Analysis

This script applies our novel CCTL algorithm to the sepsis single-cell data
for drug target discovery. CCTL combines causal inference with contrastive
learning to identify genes that causally influence survival outcome.

Key Features:
1. Learns causal DAG structure from observational data
2. Uses causal structure to guide positive/negative sampling
3. Learns disentangled representations (causal vs spurious)
4. Provides interpretable drug target rankings

Author: Xinru Qiu
Date: 2024
"""

import os
import sys
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime

# Add parent directory to path for CCTL import
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Import CCTL components
from cctl import CCTLModel, CCTLTrainer, CausalInterpreter, CausalDAGLearner

# Configuration
BASE_DIR = "/bigdata/godziklab/shared/Xinru/302004/302004_git"
DATA_DIR = os.path.join(BASE_DIR, "data/processed")
RESULTS_DIR = os.path.join(BASE_DIR, "results/advanced_analysis")
FIGURES_DIR = os.path.join(BASE_DIR, "figures/cctl")

# CCTL hyperparameters
TEMPERATURE = 0.07
LAMBDA_CAUSAL = 0.1
LAMBDA_DISENTANGLE = 0.1
LAMBDA_HARD = 0.5
EPOCHS = 50  # Reduced from 100
BATCH_SIZE = 512  # Increased from 256
LEARNING_RATE = 1e-3
MAX_CELLS = 5000  # Subsample for faster training
N_TOP_GENES = 500  # Reduced from 2000 for faster DAG learning


def load_data():
    """Load preprocessed sepsis single-cell data."""
    print("=" * 60)
    print("Loading Data")
    print("=" * 60)

    # Load AnnData
    adata_path = os.path.join(DATA_DIR, "adata_processed.h5ad")
    if not os.path.exists(adata_path):
        # Try alternative path
        adata_path = os.path.join(DATA_DIR, "sepsis_monocytes.h5ad")

    if os.path.exists(adata_path):
        adata = sc.read_h5ad(adata_path)
        print(f"Loaded AnnData: {adata.shape[0]} cells x {adata.shape[1]} genes")
    else:
        # Do NOT use synthetic data - require real data
        raise FileNotFoundError(
            f"Real data file not found at {adata_path}. "
            "Please ensure adata_processed.h5ad or sepsis_monocytes.h5ad exists. "
            "Synthetic/simulated data is NOT allowed for generating conclusions."
        )

    return adata


def prepare_data(adata, n_top_genes=N_TOP_GENES):
    """Prepare data for CCTL training."""
    print("\nPreparing data...")

    # Subsample cells if needed
    if adata.n_obs > MAX_CELLS:
        print(f"Subsampling from {adata.n_obs} to {MAX_CELLS} cells...")
        np.random.seed(42)
        indices = np.random.choice(adata.n_obs, size=MAX_CELLS, replace=False)
        adata = adata[indices].copy()
        print(f"Subsampled to {adata.n_obs} cells")

    # Get expression matrix
    if hasattr(adata.X, 'toarray'):
        X = adata.X.toarray()
    else:
        X = np.array(adata.X)

    # Get outcome
    if 'survival_7d' in adata.obs.columns:
        y = adata.obs['survival_7d'].values
    elif 'Mortality' in adata.obs.columns:
        y = adata.obs['Mortality'].values
    elif 'outcome' in adata.obs.columns:
        y = adata.obs['outcome'].values
    else:
        # Create binary outcome based on some criterion
        print("Warning: No outcome column found, creating from median expression")
        y = (X.mean(axis=1) > np.median(X.mean(axis=1))).astype(int)

    # Convert to binary integers
    # Handle string outcomes like "Survivor"/"Non-survivor"/"Healthy"
    y_original = y.copy() if hasattr(y, 'copy') else np.array(y)
    return_mask = None  # Will be set if filtering is needed

    if hasattr(y, 'dtype') and (y.dtype == object or str(y.dtype) == 'category'):
        # Convert to string array
        y_str = np.array([str(v).lower().strip() for v in y])
        print(f"Unique outcome values: {np.unique(y_str)}")

        # Create binary labels: Survivor=1, Non-survivor=0, exclude others
        y_binary = np.full(len(y_str), -1, dtype=int)  # -1 = exclude

        for i, val in enumerate(y_str):
            if 'non-survivor' in val or 'nonsurvivor' in val:
                y_binary[i] = 0  # Non-survivor = 0 (bad outcome)
            elif 'survivor' in val:
                y_binary[i] = 1  # Survivor = 1 (good outcome)
            elif 'healthy' in val:
                y_binary[i] = -1  # Exclude healthy controls from survival analysis
            else:
                y_binary[i] = -1  # Exclude unknown

        # Filter to only valid outcomes (sepsis patients with known outcome)
        valid_mask = y_binary >= 0
        n_excluded = np.sum(~valid_mask)
        if n_excluded > 0:
            print(f"Excluding {n_excluded} samples without valid survival outcome (e.g., Healthy controls)")

        y = y_binary

        # Return mask for filtering X as well
        return_mask = valid_mask
    else:
        y = np.array(y).astype(int)
        return_mask = np.ones(len(y), dtype=bool)

    print(f"Outcome distribution: 0 (Non-survivor): {np.sum(y[y>=0]==0)}, 1 (Survivor): {np.sum(y[y>=0]==1)}")

    # Apply mask to filter to valid samples only
    if 'return_mask' in dir() and return_mask is not None:
        X = X[return_mask]
        y = y[return_mask]
        print(f"After filtering: {X.shape[0]} cells")

    # Select top variable genes
    gene_var = X.var(axis=0)
    top_gene_idx = np.argsort(gene_var)[::-1][:n_top_genes]

    X = X[:, top_gene_idx]
    gene_names = list(adata.var_names[top_gene_idx])

    # Standardize
    X = (X - X.mean(axis=0)) / (X.std(axis=0) + 1e-6)

    print(f"Prepared data: {X.shape[0]} cells x {X.shape[1]} genes")
    print(f"Final outcome distribution: 0={np.sum(y==0)}, 1={np.sum(y==1)}")

    return X, y, gene_names


def train_cctl(X, y, gene_names):
    """Train CCTL model."""
    print("\n" + "=" * 60)
    print("Training CCTL Model")
    print("=" * 60)

    # Split data
    n_samples = len(y)
    indices = np.random.permutation(n_samples)
    train_idx = indices[:int(0.8 * n_samples)]
    val_idx = indices[int(0.8 * n_samples):]

    X_train, y_train = X[train_idx], y[train_idx]
    X_val, y_val = X[val_idx], y[val_idx]

    print(f"Training set: {len(train_idx)} samples")
    print(f"Validation set: {len(val_idx)} samples")

    # Create trainer
    trainer = CCTLTrainer(
        temperature=TEMPERATURE,
        lambda_causal=LAMBDA_CAUSAL,
        lambda_disentangle=LAMBDA_DISENTANGLE,
        lambda_hard=LAMBDA_HARD,
        learning_rate=LEARNING_RATE,
        patience=15
    )

    # Train
    trainer.fit(
        X_train, y_train, gene_names,
        X_val=X_val, y_val=y_val,
        epochs=EPOCHS,
        batch_size=BATCH_SIZE,
        verbose=True
    )

    return trainer, X_train, y_train, X_val, y_val


def interpret_results(trainer, X, y, gene_names):
    """Interpret CCTL results."""
    print("\n" + "=" * 60)
    print("Interpreting Results")
    print("=" * 60)

    # Create interpreter
    interpreter = CausalInterpreter(trainer, n_bootstrap=100)

    # Compute causal importance
    print("\nComputing causal importance scores...")
    importance = interpreter.compute_causal_importance(X, y, method='combined')

    print("\nTop 20 Causal Genes:")
    print("-" * 40)
    for i, (gene, score) in enumerate(list(importance.items())[:20]):
        is_causal = gene in trainer.dag_learner.get_causal_effects()
        causal_marker = "*" if is_causal else " "
        print(f"{i+1:2d}. {gene:15s}  {score:.4f} {causal_marker}")

    # Rank drug targets
    print("\nRanking drug targets...")
    rankings = interpreter.rank_drug_targets(X, y, top_k=20)

    print("\nTop 20 Drug Target Candidates:")
    print("-" * 70)
    print(f"{'Rank':<6}{'Gene':<15}{'Causal':<10}{'Effect':<12}{'Score':<10}")
    print("-" * 70)
    for r in rankings:
        print(f"{r['rank']:<6}{r['gene']:<15}{r['causal_importance']:.4f}    "
              f"{r['intervention_effect']:+.4f}     {r['target_score']:.4f}")

    # Simulate interventions for top targets
    print("\nIntervention Simulations:")
    print("-" * 50)
    for r in rankings[:5]:
        gene = r['gene']
        result = interpreter.simulate_intervention(
            X, y, gene, intervention_value=0.0
        )
        effect = result['causal_effect']
        ci = result['confidence_interval']
        sig = "***" if result['significant'] else ""
        print(f"  Targeting {gene:12s}: {effect:+.4f} (95% CI: [{ci[0]:.4f}, {ci[1]:.4f}]) {sig}")

    # Generate report
    print("\nGenerating interpretation report...")
    report = interpreter.generate_report(X, y)

    return interpreter, importance, rankings, report


def visualize_results(trainer, interpreter, X, y, gene_names, rankings):
    """Create visualizations."""
    print("\n" + "=" * 60)
    print("Creating Visualizations")
    print("=" * 60)

    os.makedirs(FIGURES_DIR, exist_ok=True)

    # 1. Training history
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))

    axes[0].plot(trainer.history['train_loss'], label='Train')
    if trainer.history['val_loss']:
        axes[0].plot(trainer.history['val_loss'], label='Validation')
    axes[0].set_xlabel('Epoch')
    axes[0].set_ylabel('Loss')
    axes[0].set_title('CCTL Training History')
    axes[0].legend()

    # 2. Causal importance distribution
    importance_values = list(interpreter.causal_importance.values())
    axes[1].hist(importance_values, bins=50, edgecolor='black', alpha=0.7)
    axes[1].axvline(np.percentile(importance_values, 95), color='r',
                    linestyle='--', label='95th percentile')
    axes[1].set_xlabel('Causal Importance')
    axes[1].set_ylabel('Count')
    axes[1].set_title('Distribution of Causal Importance Scores')
    axes[1].legend()

    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, 'cctl_training.png'), dpi=300, bbox_inches='tight')
    plt.close()

    # 3. Embedding visualization with UMAP
    print("Computing UMAP of embeddings...")
    embeddings = trainer.transform(X)
    causal_embeddings = trainer.get_causal_embeddings(X)

    # Use sklearn UMAP if available
    try:
        import umap
        reducer = umap.UMAP(n_components=2, random_state=42)
        umap_full = reducer.fit_transform(embeddings)
        umap_causal = reducer.fit_transform(causal_embeddings)

        fig, axes = plt.subplots(1, 2, figsize=(14, 5))

        # Full embedding
        scatter1 = axes[0].scatter(umap_full[:, 0], umap_full[:, 1],
                                   c=y, cmap='coolwarm', alpha=0.5, s=5)
        axes[0].set_title('Full CCTL Embedding')
        axes[0].set_xlabel('UMAP 1')
        axes[0].set_ylabel('UMAP 2')
        plt.colorbar(scatter1, ax=axes[0], label='Outcome')

        # Causal embedding
        scatter2 = axes[1].scatter(umap_causal[:, 0], umap_causal[:, 1],
                                   c=y, cmap='coolwarm', alpha=0.5, s=5)
        axes[1].set_title('Causal Subspace Embedding')
        axes[1].set_xlabel('UMAP 1')
        axes[1].set_ylabel('UMAP 2')
        plt.colorbar(scatter2, ax=axes[1], label='Outcome')

        plt.tight_layout()
        plt.savefig(os.path.join(FIGURES_DIR, 'cctl_embeddings.png'), dpi=300, bbox_inches='tight')
        plt.close()
    except ImportError:
        print("UMAP not available, using PCA instead")
        from sklearn.decomposition import PCA
        pca = PCA(n_components=2)
        pca_full = pca.fit_transform(embeddings)

        fig, ax = plt.subplots(figsize=(8, 6))
        scatter = ax.scatter(pca_full[:, 0], pca_full[:, 1],
                            c=y, cmap='coolwarm', alpha=0.5, s=5)
        ax.set_title('CCTL Embedding (PCA)')
        ax.set_xlabel('PC 1')
        ax.set_ylabel('PC 2')
        plt.colorbar(scatter, label='Outcome')
        plt.savefig(os.path.join(FIGURES_DIR, 'cctl_embeddings.png'), dpi=300, bbox_inches='tight')
        plt.close()

    # 4. Top targets bar plot
    fig, ax = plt.subplots(figsize=(12, 6))

    target_names = [r['gene'] for r in rankings[:15]]
    target_scores = [r['target_score'] for r in rankings[:15]]
    target_effects = [r['intervention_effect'] for r in rankings[:15]]

    colors = ['green' if e < 0 else 'red' for e in target_effects]

    bars = ax.barh(range(len(target_names)), target_scores, color=colors, alpha=0.7)
    ax.set_yticks(range(len(target_names)))
    ax.set_yticklabels(target_names)
    ax.set_xlabel('Target Score')
    ax.set_title('Top 15 Drug Target Candidates (CCTL)')
    ax.invert_yaxis()

    # Add effect size annotations
    for i, (bar, effect) in enumerate(zip(bars, target_effects)):
        ax.text(bar.get_width() + 0.01, bar.get_y() + bar.get_height()/2,
                f'{effect:+.3f}', va='center', fontsize=9)

    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, 'cctl_targets.png'), dpi=300, bbox_inches='tight')
    plt.close()

    # 5. Causal DAG visualization (simplified)
    fig, ax = plt.subplots(figsize=(10, 8))

    causal_effects = trainer.dag_learner.get_causal_effects()
    causal_parents = list(causal_effects.keys())[:10]
    effects = [causal_effects[g] for g in causal_parents]

    # Simple bar plot of causal effects
    colors = ['green' if e < 0 else 'red' for e in effects]
    ax.barh(range(len(causal_parents)), effects, color=colors, alpha=0.7)
    ax.axvline(0, color='black', linestyle='-', linewidth=0.5)
    ax.set_yticks(range(len(causal_parents)))
    ax.set_yticklabels(causal_parents)
    ax.set_xlabel('Causal Effect on Outcome')
    ax.set_title('Direct Causal Effects (from DAG)')
    ax.invert_yaxis()

    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, 'cctl_causal_effects.png'), dpi=300, bbox_inches='tight')
    plt.close()

    print(f"Figures saved to {FIGURES_DIR}")


def save_results(trainer, interpreter, rankings, report):
    """Save analysis results."""
    print("\n" + "=" * 60)
    print("Saving Results")
    print("=" * 60)

    os.makedirs(RESULTS_DIR, exist_ok=True)

    # Save rankings
    rankings_df = pd.DataFrame(rankings)
    rankings_df.to_csv(os.path.join(RESULTS_DIR, 'cctl_target_rankings.csv'), index=False)
    print(f"Saved rankings to {RESULTS_DIR}/cctl_target_rankings.csv")

    # Save causal importance
    importance_df = pd.DataFrame([
        {'gene': gene, 'causal_importance': score}
        for gene, score in interpreter.causal_importance.items()
    ])
    importance_df.to_csv(os.path.join(RESULTS_DIR, 'cctl_causal_importance.csv'), index=False)
    print(f"Saved importance scores to {RESULTS_DIR}/cctl_causal_importance.csv")

    # Save causal effects from DAG
    effects_df = pd.DataFrame([
        {'gene': gene, 'causal_effect': effect}
        for gene, effect in trainer.dag_learner.get_causal_effects().items()
    ])
    effects_df.to_csv(os.path.join(RESULTS_DIR, 'cctl_dag_effects.csv'), index=False)
    print(f"Saved DAG effects to {RESULTS_DIR}/cctl_dag_effects.csv")

    # Save report
    with open(os.path.join(RESULTS_DIR, 'cctl_report.txt'), 'w') as f:
        f.write(report)
    print(f"Saved report to {RESULTS_DIR}/cctl_report.txt")

    # Save model
    trainer.save(os.path.join(RESULTS_DIR, 'cctl_model.pt'))
    print(f"Saved model to {RESULTS_DIR}/cctl_model.pt")


def main():
    """Main analysis pipeline."""
    print("\n" + "=" * 70)
    print("CCTL: Causal-Contrastive Target Learning Analysis")
    print("=" * 70)
    print(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

    # Set random seed
    np.random.seed(42)

    try:
        # Load and prepare data
        adata = load_data()
        X, y, gene_names = prepare_data(adata, n_top_genes=N_TOP_GENES)

        # Train CCTL
        trainer, X_train, y_train, X_val, y_val = train_cctl(X, y, gene_names)

        # Interpret results
        interpreter, importance, rankings, report = interpret_results(
            trainer, X, y, gene_names
        )

        # Visualize
        visualize_results(trainer, interpreter, X, y, gene_names, rankings)

        # Save results
        save_results(trainer, interpreter, rankings, report)

        print("\n" + "=" * 70)
        print("CCTL Analysis Complete!")
        print("=" * 70)
        print(f"End time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

        # Print summary
        print("\nSummary:")
        print(f"  - Identified {len(trainer.dag_learner.get_causal_effects())} causal genes")
        print(f"  - Top target: {rankings[0]['gene']} (score: {rankings[0]['target_score']:.4f})")
        print(f"  - Best intervention effect: {rankings[0]['intervention_effect']:.4f}")

    except Exception as e:
        print(f"\nError: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
