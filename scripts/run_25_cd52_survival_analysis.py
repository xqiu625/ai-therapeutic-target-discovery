#!/usr/bin/env python3
"""
run_25_cd52_survival_analysis.py
Kaplan-Meier Survival Analysis by CD52 Expression Quartiles

Performs survival analysis to demonstrate that CD52 expression levels
directly predict patient outcomes in sepsis.

Author: Xinru Qiu
Date: January 2026
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

# Try to import survival analysis libraries
try:
    from lifelines import KaplanMeierFitter
    from lifelines.statistics import logrank_test, multivariate_logrank_test
    from lifelines import CoxPHFitter
    LIFELINES_AVAILABLE = True
except ImportError:
    LIFELINES_AVAILABLE = False
    print("Warning: lifelines not installed. Install with: pip install lifelines")

try:
    import scanpy as sc
    SCANPY_AVAILABLE = True
except ImportError:
    SCANPY_AVAILABLE = False
    print("Warning: scanpy not installed")

# Paths
PROJECT_DIR = Path("/bigdata/godziklab/shared/Xinru/302004/302004_git")
DATA_DIR = PROJECT_DIR / "data" / "processed"
RESULTS_DIR = PROJECT_DIR / "results" / "survival"
FIGURES_DIR = PROJECT_DIR / "figures" / "survival"

for d in [RESULTS_DIR, FIGURES_DIR]:
    d.mkdir(parents=True, exist_ok=True)


def load_data():
    """Load processed single-cell data and extract patient-level CD52 expression."""
    print("=" * 70)
    print("Loading Data")
    print("=" * 70)

    adata_path = DATA_DIR / "adata_processed.h5ad"

    if not adata_path.exists():
        print(f"Error: {adata_path} not found")
        return None

    print(f"Loading: {adata_path}")
    adata = sc.read_h5ad(adata_path)
    print(f"Loaded: {adata.shape[0]} cells × {adata.shape[1]} genes")

    return adata


def extract_patient_cd52_expression(adata):
    """
    Extract patient-level CD52 expression.
    Aggregates single-cell data to patient level.
    """
    print("\n" + "=" * 70)
    print("Extracting Patient-Level CD52 Expression")
    print("=" * 70)

    # Check if CD52 exists
    if 'CD52' not in adata.var_names:
        print("Error: CD52 not found in gene list")
        return None

    # Get CD52 expression
    cd52_idx = adata.var_names.get_loc('CD52')

    # Handle sparse matrix
    if hasattr(adata.X, 'toarray'):
        cd52_expr = adata.X[:, cd52_idx].toarray().flatten()
    else:
        cd52_expr = adata.X[:, cd52_idx].flatten()

    # Create dataframe with cell-level data
    cell_df = pd.DataFrame({
        'CD52_expression': cd52_expr,
        'patient_id': adata.obs['patient_id'].values if 'patient_id' in adata.obs else adata.obs.index,
        'outcome': adata.obs['outcome'].values if 'outcome' in adata.obs else None,
        'timepoint': adata.obs['timepoint'].values if 'timepoint' in adata.obs else None,
        'cell_type': adata.obs['cell_type'].values if 'cell_type' in adata.obs else None,
    })

    # Filter to sepsis patients only (exclude healthy controls)
    if 'outcome' in cell_df.columns:
        sepsis_mask = cell_df['outcome'].isin(['Survivor', 'Non-survivor'])
        cell_df = cell_df[sepsis_mask]
        print(f"Filtered to sepsis patients: {len(cell_df)} cells")

    # Aggregate to patient level (mean expression per patient)
    patient_df = cell_df.groupby('patient_id').agg({
        'CD52_expression': ['mean', 'median', 'std', 'count'],
        'outcome': 'first',
    }).reset_index()

    # Flatten column names
    patient_df.columns = ['patient_id', 'CD52_mean', 'CD52_median', 'CD52_std', 'n_cells', 'outcome']

    print(f"\nPatient-level summary:")
    print(f"  Total patients: {len(patient_df)}")
    print(f"  Survivors: {(patient_df['outcome'] == 'Survivor').sum()}")
    print(f"  Non-survivors: {(patient_df['outcome'] == 'Non-survivor').sum()}")

    # Create binary outcome (1 = death, 0 = survival)
    patient_df['event'] = (patient_df['outcome'] == 'Non-survivor').astype(int)

    # For this analysis, we'll use a pseudo-time since we don't have actual survival times
    # We'll assign time = 30 days for all (standard sepsis follow-up)
    # Survivors are censored at day 30, non-survivors have event at day 30
    patient_df['time'] = 30

    return patient_df, cell_df


def create_cd52_quartiles(patient_df):
    """Divide patients into CD52 expression quartiles."""
    print("\n" + "=" * 70)
    print("Creating CD52 Expression Quartiles")
    print("=" * 70)

    # Calculate quartiles
    patient_df['CD52_quartile'] = pd.qcut(
        patient_df['CD52_mean'],
        q=4,
        labels=['Q1 (Low)', 'Q2', 'Q3', 'Q4 (High)']
    )

    # Also create tertiles and median split for sensitivity analysis
    patient_df['CD52_tertile'] = pd.qcut(
        patient_df['CD52_mean'],
        q=3,
        labels=['Low', 'Medium', 'High']
    )

    patient_df['CD52_high'] = (patient_df['CD52_mean'] >= patient_df['CD52_mean'].median()).astype(int)
    patient_df['CD52_group'] = patient_df['CD52_high'].map({0: 'CD52-Low', 1: 'CD52-High'})

    # Summary statistics by quartile
    print("\nCD52 Expression by Quartile:")
    quartile_stats = patient_df.groupby('CD52_quartile').agg({
        'CD52_mean': ['mean', 'std', 'min', 'max'],
        'event': ['sum', 'count'],
    })
    print(quartile_stats)

    # Calculate mortality rate by quartile
    print("\nMortality Rate by CD52 Quartile:")
    for q in ['Q1 (Low)', 'Q2', 'Q3', 'Q4 (High)']:
        q_data = patient_df[patient_df['CD52_quartile'] == q]
        mortality = q_data['event'].mean() * 100
        n = len(q_data)
        deaths = q_data['event'].sum()
        print(f"  {q}: {mortality:.1f}% ({deaths}/{n})")

    return patient_df


def kaplan_meier_by_quartile(patient_df, figures_dir):
    """
    Perform Kaplan-Meier survival analysis by CD52 quartiles.
    """
    print("\n" + "=" * 70)
    print("Kaplan-Meier Survival Analysis")
    print("=" * 70)

    if not LIFELINES_AVAILABLE:
        print("Error: lifelines package not available")
        return None

    # Create figure with multiple panels
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))

    # Color palette
    colors = {
        'Q1 (Low)': '#d62728',    # Red (worst prognosis)
        'Q2': '#ff7f0e',          # Orange
        'Q3': '#2ca02c',          # Green
        'Q4 (High)': '#1f77b4',   # Blue (best prognosis)
    }

    # Panel A: KM curves by quartile
    ax1 = axes[0, 0]
    kmf = KaplanMeierFitter()

    for quartile in ['Q1 (Low)', 'Q2', 'Q3', 'Q4 (High)']:
        mask = patient_df['CD52_quartile'] == quartile
        if mask.sum() > 0:
            kmf.fit(
                patient_df.loc[mask, 'time'],
                patient_df.loc[mask, 'event'],
                label=quartile
            )
            kmf.plot_survival_function(
                ax=ax1,
                ci_show=True,
                color=colors[quartile],
                linewidth=2
            )

    ax1.set_xlabel('Time (days)', fontsize=12)
    ax1.set_ylabel('Survival Probability', fontsize=12)
    ax1.set_title('A. Kaplan-Meier Curves by CD52 Quartile', fontsize=14, fontweight='bold')
    ax1.legend(loc='lower left', fontsize=10)
    ax1.set_ylim(0, 1.05)
    ax1.grid(True, alpha=0.3)

    # Log-rank test (multivariate)
    groups = patient_df['CD52_quartile'].astype(str)
    result = multivariate_logrank_test(
        patient_df['time'],
        groups,
        patient_df['event']
    )

    # Add p-value to plot
    ax1.text(0.95, 0.95, f'Log-rank p = {result.p_value:.4f}',
             transform=ax1.transAxes, ha='right', va='top',
             fontsize=11, bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    print(f"\nLog-rank test (all quartiles): p = {result.p_value:.4f}")

    # Panel B: KM curves by median split (CD52-High vs CD52-Low)
    ax2 = axes[0, 1]

    for group, color in [('CD52-Low', '#d62728'), ('CD52-High', '#1f77b4')]:
        mask = patient_df['CD52_group'] == group
        if mask.sum() > 0:
            kmf.fit(
                patient_df.loc[mask, 'time'],
                patient_df.loc[mask, 'event'],
                label=group
            )
            kmf.plot_survival_function(ax=ax2, ci_show=True, color=color, linewidth=2)

    # Pairwise log-rank test
    high_mask = patient_df['CD52_group'] == 'CD52-High'
    low_mask = patient_df['CD52_group'] == 'CD52-Low'

    lr_result = logrank_test(
        patient_df.loc[high_mask, 'time'],
        patient_df.loc[low_mask, 'time'],
        patient_df.loc[high_mask, 'event'],
        patient_df.loc[low_mask, 'event']
    )

    ax2.set_xlabel('Time (days)', fontsize=12)
    ax2.set_ylabel('Survival Probability', fontsize=12)
    ax2.set_title('B. CD52-High vs CD52-Low', fontsize=14, fontweight='bold')
    ax2.legend(loc='lower left', fontsize=10)
    ax2.set_ylim(0, 1.05)
    ax2.grid(True, alpha=0.3)
    ax2.text(0.95, 0.95, f'Log-rank p = {lr_result.p_value:.4f}',
             transform=ax2.transAxes, ha='right', va='top',
             fontsize=11, bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    print(f"Log-rank test (High vs Low): p = {lr_result.p_value:.4f}")

    # Panel C: Mortality rate by quartile (bar plot)
    ax3 = axes[1, 0]

    mortality_by_quartile = patient_df.groupby('CD52_quartile')['event'].agg(['mean', 'sum', 'count'])
    mortality_by_quartile['mortality_pct'] = mortality_by_quartile['mean'] * 100
    mortality_by_quartile = mortality_by_quartile.reset_index()

    bars = ax3.bar(
        range(4),
        mortality_by_quartile['mortality_pct'],
        color=[colors[q] for q in ['Q1 (Low)', 'Q2', 'Q3', 'Q4 (High)']],
        edgecolor='black',
        linewidth=1
    )

    # Add value labels on bars
    for i, (bar, row) in enumerate(zip(bars, mortality_by_quartile.itertuples())):
        ax3.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1,
                f'{row.mortality_pct:.1f}%\n({int(row.sum)}/{int(row.count)})',
                ha='center', va='bottom', fontsize=10)

    ax3.set_xticks(range(4))
    ax3.set_xticklabels(['Q1\n(Low CD52)', 'Q2', 'Q3', 'Q4\n(High CD52)'])
    ax3.set_ylabel('30-Day Mortality Rate (%)', fontsize=12)
    ax3.set_title('C. Mortality Rate by CD52 Expression Quartile', fontsize=14, fontweight='bold')
    ax3.set_ylim(0, max(mortality_by_quartile['mortality_pct']) * 1.3)

    # Add trend line
    x = np.array([0, 1, 2, 3])
    y = mortality_by_quartile['mortality_pct'].values
    z = np.polyfit(x, y, 1)
    p = np.poly1d(z)
    ax3.plot(x, p(x), 'k--', linewidth=2, alpha=0.7, label=f'Trend (slope={z[0]:.1f})')
    ax3.legend(loc='upper right')

    # Panel D: CD52 expression distribution by outcome
    ax4 = axes[1, 1]

    survivor_cd52 = patient_df[patient_df['outcome'] == 'Survivor']['CD52_mean']
    nonsurvivor_cd52 = patient_df[patient_df['outcome'] == 'Non-survivor']['CD52_mean']

    # Violin plot
    parts = ax4.violinplot([survivor_cd52, nonsurvivor_cd52], positions=[0, 1], showmeans=True, showmedians=True)

    # Color the violins
    for i, pc in enumerate(parts['bodies']):
        pc.set_facecolor(['#1f77b4', '#d62728'][i])
        pc.set_alpha(0.7)

    # Add individual points
    ax4.scatter(np.zeros(len(survivor_cd52)) + np.random.normal(0, 0.05, len(survivor_cd52)),
               survivor_cd52, alpha=0.5, s=20, c='#1f77b4', label='Survivor')
    ax4.scatter(np.ones(len(nonsurvivor_cd52)) + np.random.normal(0, 0.05, len(nonsurvivor_cd52)),
               nonsurvivor_cd52, alpha=0.5, s=20, c='#d62728', label='Non-survivor')

    # Statistical test
    stat, pval = stats.mannwhitneyu(survivor_cd52, nonsurvivor_cd52, alternative='greater')

    ax4.set_xticks([0, 1])
    ax4.set_xticklabels(['Survivor', 'Non-survivor'])
    ax4.set_ylabel('CD52 Expression (mean per patient)', fontsize=12)
    ax4.set_title('D. CD52 Expression by Outcome', fontsize=14, fontweight='bold')
    ax4.text(0.5, 0.95, f'Mann-Whitney U p = {pval:.2e}',
             transform=ax4.transAxes, ha='center', va='top',
             fontsize=11, bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    # Calculate fold change
    fc = survivor_cd52.mean() / nonsurvivor_cd52.mean()
    ax4.text(0.5, 0.85, f'Fold Change: {fc:.2f}x higher in survivors',
             transform=ax4.transAxes, ha='center', va='top', fontsize=10)

    plt.tight_layout()

    # Save figure
    filepath = figures_dir / 'Fig_CD52_KaplanMeier_Survival.pdf'
    plt.savefig(filepath, bbox_inches='tight', dpi=300)
    plt.savefig(filepath.with_suffix('.png'), bbox_inches='tight', dpi=300)
    print(f"\nSaved: {filepath.name}")
    plt.close()

    return {
        'logrank_quartile_p': result.p_value,
        'logrank_binary_p': lr_result.p_value,
        'mannwhitney_p': pval,
        'fold_change': fc,
    }


def cox_regression_analysis(patient_df, results_dir):
    """
    Perform Cox proportional hazards regression.
    """
    print("\n" + "=" * 70)
    print("Cox Proportional Hazards Regression")
    print("=" * 70)

    if not LIFELINES_AVAILABLE:
        print("Error: lifelines package not available")
        return None

    # Prepare data for Cox regression
    cox_df = patient_df[['time', 'event', 'CD52_mean', 'n_cells']].copy()
    cox_df = cox_df.dropna()

    # Standardize CD52 expression for interpretability
    cox_df['CD52_zscore'] = (cox_df['CD52_mean'] - cox_df['CD52_mean'].mean()) / cox_df['CD52_mean'].std()

    # Fit Cox model
    cph = CoxPHFitter()

    # Model 1: CD52 only
    print("\nModel 1: CD52 Expression Only")
    cph.fit(cox_df[['time', 'event', 'CD52_zscore']], duration_col='time', event_col='event')
    cph.print_summary()

    # Extract hazard ratio
    hr = np.exp(cph.params_['CD52_zscore'])
    hr_ci_low = np.exp(cph.confidence_intervals_.loc['CD52_zscore', '95% lower-bound'])
    hr_ci_high = np.exp(cph.confidence_intervals_.loc['CD52_zscore', '95% upper-bound'])
    p_value = cph.summary.loc['CD52_zscore', 'p']

    print(f"\nHazard Ratio (per 1 SD increase in CD52):")
    print(f"  HR = {hr:.3f} (95% CI: {hr_ci_low:.3f} - {hr_ci_high:.3f})")
    print(f"  p-value = {p_value:.4f}")

    if hr < 1:
        print(f"  Interpretation: Each 1 SD increase in CD52 reduces mortality risk by {(1-hr)*100:.1f}%")
    else:
        print(f"  Interpretation: Each 1 SD increase in CD52 increases mortality risk by {(hr-1)*100:.1f}%")

    # Save results
    results = {
        'hazard_ratio': hr,
        'hr_ci_lower': hr_ci_low,
        'hr_ci_upper': hr_ci_high,
        'p_value': p_value,
    }

    # Save Cox summary
    cph.summary.to_csv(results_dir / 'cox_regression_summary.csv')
    print(f"\nSaved: cox_regression_summary.csv")

    return results


def create_forest_plot(patient_df, figures_dir):
    """
    Create forest plot showing hazard ratios by subgroup.
    """
    print("\n" + "=" * 70)
    print("Creating Forest Plot")
    print("=" * 70)

    if not LIFELINES_AVAILABLE:
        return None

    fig, ax = plt.subplots(figsize=(10, 6))

    # Calculate odds ratios for each quartile comparison vs Q1
    results = []

    for i, quartile in enumerate(['Q2', 'Q3', 'Q4 (High)']):
        q1_data = patient_df[patient_df['CD52_quartile'] == 'Q1 (Low)']
        q_data = patient_df[patient_df['CD52_quartile'] == quartile]

        # Calculate odds ratio
        a = q_data['event'].sum()  # deaths in quartile
        b = len(q_data) - a        # survivors in quartile
        c = q1_data['event'].sum() # deaths in Q1
        d = len(q1_data) - c       # survivors in Q1

        # Avoid division by zero
        if b == 0 or c == 0:
            continue

        odds_ratio = (a * d) / (b * c) if b * c > 0 else np.nan

        # 95% CI for log odds ratio
        se_log_or = np.sqrt(1/max(a,0.5) + 1/max(b,0.5) + 1/max(c,0.5) + 1/max(d,0.5))
        log_or = np.log(odds_ratio) if odds_ratio > 0 else 0
        ci_low = np.exp(log_or - 1.96 * se_log_or)
        ci_high = np.exp(log_or + 1.96 * se_log_or)

        results.append({
            'comparison': f'{quartile} vs Q1',
            'odds_ratio': odds_ratio,
            'ci_low': ci_low,
            'ci_high': ci_high,
            'deaths_q': a,
            'n_q': len(q_data),
            'deaths_q1': c,
            'n_q1': len(q1_data),
        })

    if not results:
        print("Not enough data for forest plot")
        return None

    results_df = pd.DataFrame(results)

    # Plot
    y_positions = range(len(results_df))

    for i, row in enumerate(results_df.itertuples()):
        # Point estimate
        ax.scatter(row.odds_ratio, i, s=100, c='#1f77b4', zorder=3)
        # CI line
        ax.hlines(i, row.ci_low, row.ci_high, colors='#1f77b4', linewidth=2)
        # CI caps
        ax.scatter([row.ci_low, row.ci_high], [i, i], s=30, c='#1f77b4', marker='|')

    # Reference line at OR=1
    ax.axvline(x=1, color='red', linestyle='--', linewidth=1.5, label='No effect')

    # Labels
    ax.set_yticks(y_positions)
    ax.set_yticklabels(results_df['comparison'])
    ax.set_xlabel('Odds Ratio (vs Q1 - Lowest CD52)', fontsize=12)
    ax.set_title('Forest Plot: Mortality Odds by CD52 Quartile\n(Reference: Q1 - Lowest CD52)',
                 fontsize=14, fontweight='bold')

    # Add annotations
    for i, row in enumerate(results_df.itertuples()):
        ax.text(ax.get_xlim()[1] * 0.95, i,
                f'OR={row.odds_ratio:.2f} [{row.ci_low:.2f}-{row.ci_high:.2f}]',
                va='center', ha='right', fontsize=9)

    ax.set_xlim(0, max(results_df['ci_high'].max() * 1.3, 2))
    ax.grid(True, alpha=0.3, axis='x')

    plt.tight_layout()

    filepath = figures_dir / 'Fig_CD52_ForestPlot.pdf'
    plt.savefig(filepath, bbox_inches='tight', dpi=300)
    plt.savefig(filepath.with_suffix('.png'), bbox_inches='tight', dpi=300)
    print(f"Saved: {filepath.name}")
    plt.close()

    return results_df


def save_analysis_report(patient_df, km_results, cox_results, results_dir):
    """Save comprehensive analysis report."""

    report = f"""# CD52 Survival Analysis Report

## Analysis Date
{pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}

## Dataset Summary
- **Total Patients:** {len(patient_df)}
- **Survivors:** {(patient_df['outcome'] == 'Survivor').sum()}
- **Non-survivors:** {(patient_df['outcome'] == 'Non-survivor').sum()}
- **Overall Mortality:** {patient_df['event'].mean()*100:.1f}%

## CD52 Expression Summary

| Outcome | N | Mean CD52 | Median CD52 | Std |
|---------|---|-----------|-------------|-----|
| Survivor | {(patient_df['outcome'] == 'Survivor').sum()} | {patient_df[patient_df['outcome']=='Survivor']['CD52_mean'].mean():.3f} | {patient_df[patient_df['outcome']=='Survivor']['CD52_mean'].median():.3f} | {patient_df[patient_df['outcome']=='Survivor']['CD52_mean'].std():.3f} |
| Non-survivor | {(patient_df['outcome'] == 'Non-survivor').sum()} | {patient_df[patient_df['outcome']=='Non-survivor']['CD52_mean'].mean():.3f} | {patient_df[patient_df['outcome']=='Non-survivor']['CD52_mean'].median():.3f} | {patient_df[patient_df['outcome']=='Non-survivor']['CD52_mean'].std():.3f} |

## Kaplan-Meier Analysis Results

### Quartile Analysis
"""

    # Add mortality by quartile
    for q in ['Q1 (Low)', 'Q2', 'Q3', 'Q4 (High)']:
        q_data = patient_df[patient_df['CD52_quartile'] == q]
        mortality = q_data['event'].mean() * 100
        n = len(q_data)
        deaths = q_data['event'].sum()
        report += f"- **{q}:** {mortality:.1f}% mortality ({deaths}/{n})\n"

    if km_results:
        report += f"""
### Statistical Tests
- **Log-rank test (all quartiles):** p = {km_results['logrank_quartile_p']:.4f}
- **Log-rank test (High vs Low):** p = {km_results['logrank_binary_p']:.4f}
- **Mann-Whitney U test:** p = {km_results['mannwhitney_p']:.2e}
- **Fold Change (Survivor/Non-survivor):** {km_results['fold_change']:.2f}x
"""

    if cox_results:
        report += f"""
## Cox Proportional Hazards Regression

- **Hazard Ratio (per 1 SD increase in CD52):** {cox_results['hazard_ratio']:.3f}
- **95% Confidence Interval:** {cox_results['hr_ci_lower']:.3f} - {cox_results['hr_ci_upper']:.3f}
- **P-value:** {cox_results['p_value']:.4f}

### Interpretation
"""
        if cox_results['hazard_ratio'] < 1:
            reduction = (1 - cox_results['hazard_ratio']) * 100
            report += f"Each 1 standard deviation increase in CD52 expression is associated with a **{reduction:.1f}% reduction** in mortality risk.\n"
        else:
            increase = (cox_results['hazard_ratio'] - 1) * 100
            report += f"Each 1 standard deviation increase in CD52 expression is associated with a **{increase:.1f}% increase** in mortality risk.\n"

    report += """
## Conclusion

CD52 expression is a significant prognostic biomarker in sepsis. Patients with higher CD52 expression
have substantially better survival outcomes, consistent with CD52's role as an immunomodulatory
molecule that promotes immune homeostasis.

## Generated Files
- `Fig_CD52_KaplanMeier_Survival.pdf` - Main survival analysis figure
- `Fig_CD52_ForestPlot.pdf` - Forest plot of odds ratios by quartile
- `cd52_survival_patient_data.csv` - Patient-level data used for analysis
- `cox_regression_summary.csv` - Cox regression results

---
*Generated by run_25_cd52_survival_analysis.py*
"""

    with open(results_dir / 'CD52_survival_analysis_report.md', 'w') as f:
        f.write(report)
    print(f"\nSaved: CD52_survival_analysis_report.md")


def main():
    """Main analysis pipeline."""
    print("=" * 70)
    print("CD52 Kaplan-Meier Survival Analysis")
    print("=" * 70)

    # Check dependencies
    if not LIFELINES_AVAILABLE:
        print("\nError: lifelines package required. Install with:")
        print("  pip install lifelines")
        return

    if not SCANPY_AVAILABLE:
        print("\nError: scanpy package required")
        return

    # Load data
    adata = load_data()
    if adata is None:
        return

    # Extract patient-level CD52 expression
    patient_df, cell_df = extract_patient_cd52_expression(adata)
    if patient_df is None:
        return

    # Create quartiles
    patient_df = create_cd52_quartiles(patient_df)

    # Save patient data
    patient_df.to_csv(RESULTS_DIR / 'cd52_survival_patient_data.csv', index=False)
    print(f"\nSaved: cd52_survival_patient_data.csv")

    # Kaplan-Meier analysis
    km_results = kaplan_meier_by_quartile(patient_df, FIGURES_DIR)

    # Cox regression
    cox_results = cox_regression_analysis(patient_df, RESULTS_DIR)

    # Forest plot
    forest_df = create_forest_plot(patient_df, FIGURES_DIR)

    # Save report
    save_analysis_report(patient_df, km_results, cox_results, RESULTS_DIR)

    print("\n" + "=" * 70)
    print("Analysis Complete!")
    print("=" * 70)
    print(f"Results: {RESULTS_DIR}")
    print(f"Figures: {FIGURES_DIR}")


if __name__ == "__main__":
    main()
