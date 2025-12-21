#!/usr/bin/env python
"""
03 - Outcome Prediction with Machine Learning
Sepsis Target Discovery Pipeline

Trains ML classifiers to predict sepsis outcome and identifies important features.
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
warnings.filterwarnings('ignore')

# ML libraries
from sklearn.model_selection import StratifiedKFold, cross_val_score, train_test_split
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.metrics import roc_auc_score, roc_curve, classification_report
from sklearn.preprocessing import StandardScaler, LabelEncoder
import joblib

# Optional libraries - graceful fallback
try:
    import xgboost as xgb
    HAS_XGB = True
except ImportError:
    print("xgboost not available, using sklearn GradientBoosting")
    HAS_XGB = False

try:
    import lightgbm as lgb
    HAS_LGB = True
except ImportError:
    print("lightgbm not available, skipping")
    HAS_LGB = False

try:
    import shap
    HAS_SHAP = True
except ImportError:
    print("shap not available, using feature_importances_ instead")
    HAS_SHAP = False

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
FIGURES_DIR = BASE_DIR / 'figures'
MODELS_DIR = RESULTS_DIR / 'models'
TABLES_DIR = RESULTS_DIR / 'tables'

MODELS_DIR.mkdir(parents=True, exist_ok=True)
TABLES_DIR.mkdir(parents=True, exist_ok=True)

print("="*60)
print("03 - ML OUTCOME PREDICTION")
print("="*60)

# ============================================================
# 1. Load Data
# ============================================================
print("\n" + "="*60)
print("1. LOADING DATA")
print("="*60)

adata = sc.read_h5ad(PROCESSED_DIR / 'adata_processed.h5ad')
print(f"Loaded: {adata.n_obs} cells, {adata.n_vars} genes")

# Filter to sepsis samples only
adata_sepsis = adata[adata.obs['disease'] == 'Sepsis'].copy()
print(f"Sepsis samples: {adata_sepsis.n_obs} cells")
print(adata_sepsis.obs['outcome'].value_counts())

# ============================================================
# 2. Prepare Features
# ============================================================
print("\n" + "="*60)
print("2. PREPARING FEATURES")
print("="*60)

# Use HVGs as features
hvg_mask = adata_sepsis.var['highly_variable']
X_cells = adata_sepsis[:, hvg_mask].X
if hasattr(X_cells, 'toarray'):
    X_cells = X_cells.toarray()

hvg_names = adata_sepsis.var_names[hvg_mask].tolist()
print(f"Features (HVGs): {X_cells.shape[1]}")
print(f"Cells: {X_cells.shape[0]}")

# Prepare labels
le = LabelEncoder()
y_cells = le.fit_transform(adata_sepsis.obs['outcome'])
print(f"Classes: {le.classes_}")
print(f"Label distribution: {dict(zip(le.classes_, [sum(y_cells==i) for i in range(len(le.classes_))]))}")

# ============================================================
# 3. Train/Test Split (Patient-based)
# ============================================================
print("\n" + "="*60)
print("3. TRAIN/TEST SPLIT")
print("="*60)

# Patient-based splitting to avoid data leakage
patients = adata_sepsis.obs['patient_id'].unique()
print(f"Total patients: {len(patients)}")

# Split patients
train_patients, test_patients = train_test_split(patients, test_size=0.3, random_state=42)

train_mask = adata_sepsis.obs['patient_id'].isin(train_patients)
test_mask = adata_sepsis.obs['patient_id'].isin(test_patients)

X_train, X_test = X_cells[train_mask], X_cells[test_mask]
y_train, y_test = y_cells[train_mask], y_cells[test_mask]

print(f"Train: {X_train.shape[0]} cells from {len(train_patients)} patients")
print(f"Test: {X_test.shape[0]} cells from {len(test_patients)} patients")

# Scale features
scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled = scaler.transform(X_test)

# ============================================================
# 4. Train Classifiers
# ============================================================
print("\n" + "="*60)
print("4. TRAINING CLASSIFIERS")
print("="*60)

# Store models and predictions
models = {}
predictions = {}

# Random Forest
print("\nTraining Random Forest...")
rf = RandomForestClassifier(
    n_estimators=200,
    max_depth=10,
    min_samples_split=5,
    class_weight='balanced',
    random_state=42,
    n_jobs=-1
)
rf.fit(X_train_scaled, y_train)
y_pred_rf = rf.predict_proba(X_test_scaled)[:, 1]
auc_rf = roc_auc_score(y_test, y_pred_rf)
print(f"Random Forest AUC: {auc_rf:.4f}")
models['rf'] = rf
predictions['Random Forest'] = (y_pred_rf, auc_rf)

# XGBoost or GradientBoosting fallback
if HAS_XGB:
    print("\nTraining XGBoost...")
    scale_pos_weight = sum(y_train==0) / max(sum(y_train==1), 1)
    xgb_clf = xgb.XGBClassifier(
        n_estimators=200,
        max_depth=6,
        learning_rate=0.1,
        scale_pos_weight=scale_pos_weight,
        random_state=42,
        n_jobs=-1,
        eval_metric='logloss'
    )
    xgb_clf.fit(X_train_scaled, y_train)
    y_pred_xgb = xgb_clf.predict_proba(X_test_scaled)[:, 1]
    auc_xgb = roc_auc_score(y_test, y_pred_xgb)
    print(f"XGBoost AUC: {auc_xgb:.4f}")
    models['xgb'] = xgb_clf
    predictions['XGBoost'] = (y_pred_xgb, auc_xgb)
else:
    print("\nTraining GradientBoosting (sklearn)...")
    gb_clf = GradientBoostingClassifier(
        n_estimators=200,
        max_depth=6,
        learning_rate=0.1,
        random_state=42
    )
    gb_clf.fit(X_train_scaled, y_train)
    y_pred_gb = gb_clf.predict_proba(X_test_scaled)[:, 1]
    auc_gb = roc_auc_score(y_test, y_pred_gb)
    print(f"GradientBoosting AUC: {auc_gb:.4f}")
    models['gb'] = gb_clf
    predictions['GradientBoosting'] = (y_pred_gb, auc_gb)

# LightGBM (optional)
if HAS_LGB:
    print("\nTraining LightGBM...")
    lgb_clf = lgb.LGBMClassifier(
        n_estimators=200,
        max_depth=6,
        learning_rate=0.1,
        class_weight='balanced',
        random_state=42,
        n_jobs=-1,
        verbose=-1
    )
    lgb_clf.fit(X_train_scaled, y_train)
    y_pred_lgb = lgb_clf.predict_proba(X_test_scaled)[:, 1]
    auc_lgb = roc_auc_score(y_test, y_pred_lgb)
    print(f"LightGBM AUC: {auc_lgb:.4f}")
    models['lgb'] = lgb_clf
    predictions['LightGBM'] = (y_pred_lgb, auc_lgb)

# ============================================================
# 5. ROC Curves
# ============================================================
print("\n" + "="*60)
print("5. ROC CURVES")
print("="*60)

fig, ax = plt.subplots(figsize=(8, 8))

for name, (y_pred, auc) in predictions.items():
    fpr, tpr, _ = roc_curve(y_test, y_pred)
    ax.plot(fpr, tpr, label=f'{name} (AUC = {auc:.3f})', linewidth=2)

ax.plot([0, 1], [0, 1], 'k--', label='Random')
ax.set_xlabel('False Positive Rate', fontsize=12)
ax.set_ylabel('True Positive Rate', fontsize=12)
ax.set_title('ROC Curves: Sepsis Outcome Prediction', fontsize=14)
ax.legend(loc='lower right', fontsize=10)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(FIGURES_DIR / 'fig3_outcome_prediction.png', dpi=150, bbox_inches='tight')
plt.close()
print("Saved ROC curve")

# ============================================================
# 6. Feature Importance (SHAP or feature_importances_)
# ============================================================
print("\n" + "="*60)
print("6. FEATURE IMPORTANCE")
print("="*60)

# Use best available model for feature importance
best_model = models.get('xgb') or models.get('gb') or models['rf']
best_model_name = 'XGBoost' if 'xgb' in models else ('GradientBoosting' if 'gb' in models else 'RandomForest')

if HAS_SHAP and (HAS_XGB or 'gb' in models):
    print(f"Computing SHAP values for {best_model_name}...")

    # Sample for speed
    n_shap = min(2000, len(X_test_scaled))
    shap_idx = np.random.choice(len(X_test_scaled), n_shap, replace=False)
    X_shap = X_test_scaled[shap_idx]

    explainer = shap.TreeExplainer(best_model)
    shap_values = explainer.shap_values(X_shap)
    print(f"SHAP values shape: {shap_values.shape}")

    # Mean absolute SHAP values
    mean_shap = np.abs(shap_values).mean(axis=0)
    feature_importance = pd.DataFrame({
        'gene': hvg_names,
        'importance': mean_shap
    }).sort_values('importance', ascending=False)

    # SHAP summary plot
    plt.figure(figsize=(10, 12))
    shap.summary_plot(shap_values, X_shap, feature_names=hvg_names, max_display=30, show=False)
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / 'shap_summary.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Saved SHAP summary plot")
else:
    print(f"Using feature_importances_ from {best_model_name}...")
    feature_importance = pd.DataFrame({
        'gene': hvg_names,
        'importance': best_model.feature_importances_
    }).sort_values('importance', ascending=False)

print("\nTop 20 genes by importance:")
print(feature_importance.head(20).to_string())

# Check CD52
cd52_row = feature_importance[feature_importance['gene'] == 'CD52']
if not cd52_row.empty:
    cd52_idx = feature_importance[feature_importance['gene'] == 'CD52'].index[0]
    cd52_rank = (feature_importance['importance'] > feature_importance.loc[cd52_idx, 'importance']).sum() + 1
    cd52_imp = cd52_row['importance'].values[0]
    print(f"\n*** CD52 rank: {cd52_rank} (importance: {cd52_imp:.4f}) ***")
else:
    cd52_rank = None
    print("\nCD52 not in HVGs")

# Bar plot (top 30)
plt.figure(figsize=(10, 10))
top_30 = feature_importance.head(30)
colors = ['red' if g == 'CD52' else 'steelblue' for g in top_30['gene']]
plt.barh(range(len(top_30)), top_30['importance'].values, color=colors)
plt.yticks(range(len(top_30)), top_30['gene'].values)
plt.xlabel('Feature Importance')
plt.ylabel('Gene')
plt.title(f'Top 30 Predictive Genes ({best_model_name})')
plt.gca().invert_yaxis()
plt.tight_layout()
plt.savefig(FIGURES_DIR / 'feature_importance_bar.png', dpi=150, bbox_inches='tight')
plt.close()
print("Saved feature importance bar plot")

# ============================================================
# 7. Save Results
# ============================================================
print("\n" + "="*60)
print("7. SAVING RESULTS")
print("="*60)

# Save models
for name, model in models.items():
    joblib.dump(model, MODELS_DIR / f'{name}_classifier.joblib')
joblib.dump(scaler, MODELS_DIR / 'feature_scaler.joblib')
print(f"Saved {len(models)} models")

# Save feature importance
feature_importance.to_csv(TABLES_DIR / 'ml_feature_importance.csv', index=False)
print("Saved feature importance")

# Save performance
performance = pd.DataFrame([
    {'model': name, 'auc': auc}
    for name, (_, auc) in predictions.items()
])
performance.to_csv(TABLES_DIR / 'ml_model_performance.csv', index=False)
print("Saved model performance")
print(performance)

# ============================================================
# Summary
# ============================================================
print("\n" + "="*60)
print("ML PREDICTION SUMMARY")
print("="*60)
print(f"Best model: {performance.loc[performance['auc'].idxmax(), 'model']} (AUC = {performance['auc'].max():.4f})")
print(f"\nTop 10 predictive genes:")
for i, (_, row) in enumerate(feature_importance.head(10).iterrows()):
    marker = " *** CD52 ***" if row['gene'] == 'CD52' else ""
    print(f"  {i+1}. {row['gene']}: {row['importance']:.4f}{marker}")

if cd52_rank:
    print(f"\nCD52 rank: {cd52_rank}")
else:
    print("\nNote: CD52 not in HVGs, check full gene importance")

print("\n" + "="*60)
print("ML PREDICTION COMPLETE!")
print("="*60)
