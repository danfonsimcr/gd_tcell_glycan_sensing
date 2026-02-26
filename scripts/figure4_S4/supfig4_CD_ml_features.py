#!/usr/bin/env python3
"""
Supplementary Figure 4C–D — ML Feature Importance & CRISPR Correlation

  Sup Fig 4C: ML feature importance (normalised 0–1) per model
              (Adaboost, Random Forest, GB Trees, XGBoost, Logistic Regression)
  Sup Fig 4D: Sum of normalised ML feature importances vs CRISPR β-score
              differential effect

--- Input data ----------------------------------------------------------------

Script looks in demo/figure4_S4/. Same expression format as figure4_ml_classifier.py; CRISPR file: Gene, diff.

Usage:
  conda activate gdt201
  python supfig4_CD_ml_features.py
"""

import os
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import (
    RandomForestClassifier, GradientBoostingClassifier, AdaBoostClassifier,
)
from sklearn.tree import DecisionTreeClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.preprocessing import StandardScaler
from sklearn.exceptions import ConvergenceWarning
from xgboost import XGBClassifier

warnings.filterwarnings("ignore", category=ConvergenceWarning)

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.join(SCRIPT_DIR, "..", "..")
DATA_DIR = os.path.join(REPO_ROOT, "demo", "figure4_S4")

CLASSIFIERS = {
    "Logistic Regression": LogisticRegression(random_state=0),
    "XGBoost": XGBClassifier(random_state=0, eval_metric="logloss"),
    "Decision Tree": DecisionTreeClassifier(random_state=0),
    "Random Forest": RandomForestClassifier(random_state=0),
    "GB Trees": GradientBoostingClassifier(random_state=0),
    "Adaboost": AdaBoostClassifier(random_state=0, algorithm="SAMME"),
    "Naive Bayes": GaussianNB(),
}

FEATURE_IMPORTANCE_MODELS = [
    "Adaboost", "Random Forest", "GB Trees", "XGBoost", "Logistic Regression"
]


def train_all_classifiers(X, y, test_size=0.25, random_state=0):
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    X_train, X_test, y_train, y_test = train_test_split(
        X_scaled, y, stratify=y, test_size=test_size, random_state=random_state
    )
    fitted = {}
    for name, clf in CLASSIFIERS.items():
        clf.fit(X_train, y_train)
        fitted[name] = clf
    return fitted, scaler


def _get_importance(clf, feature_names):
    if hasattr(clf, "feature_importances_"):
        raw = clf.feature_importances_
    elif hasattr(clf, "coef_"):
        raw = np.abs(clf.coef_).flatten()
    else:
        return None
    return pd.Series(raw, index=feature_names)


def _normalise_01(series):
    mn, mx = series.min(), series.max()
    return (series - mn) / (mx - mn) if mx > mn else series * 0.0


def get_all_normalised_importances(fitted_models, feature_names):
    result = {}
    for name in FEATURE_IMPORTANCE_MODELS:
        clf = fitted_models.get(name)
        if clf is None:
            continue
        imp = _get_importance(clf, feature_names)
        if imp is not None:
            result[name] = _normalise_01(imp)
    return result


# =============================================================================
# Supplementary Figure 4C — Feature importance per model (normalised 0–1)
# =============================================================================

def plot_feature_importance(fitted_models, feature_names):
    norm_imps = get_all_normalised_importances(fitted_models, feature_names)
    n_models = len(norm_imps)
    if n_models == 0:
        return None

    fig, axes = plt.subplots(1, n_models, figsize=(5 * n_models, 6), sharey=True)
    if n_models == 1:
        axes = [axes]

    for ax, (name, imp) in zip(axes, norm_imps.items()):
        imp.sort_values().plot.barh(ax=ax, color="#003B7F")
        ax.set_title(name, fontsize=11, fontweight="bold")
        ax.set_xlabel("Normalised Importance (0–1)")
        ax.set_xlim(0, 1.05)

    plt.suptitle("Feature Importance per Model (Normalised 0–1)",
                 fontsize=14, fontweight="bold", y=1.02)
    plt.tight_layout()
    plt.show()
    return norm_imps


# =============================================================================
# Supplementary Figure 4D — Sum of normalised features vs CRISPR β-score
# =============================================================================

def plot_sum_features_vs_beta(crispr_df, norm_imps_dict):
    imp_df = pd.DataFrame(norm_imps_dict)
    imp_df["sum_norm_importance"] = imp_df.sum(axis=1)

    merged = crispr_df.merge(
        imp_df["sum_norm_importance"],
        left_on="Gene", right_index=True, how="inner"
    )
    r, p = pearsonr(merged["diff"], merged["sum_norm_importance"])

    fig, ax = plt.subplots(figsize=(7, 7))
    ax.scatter(merged["diff"], merged["sum_norm_importance"],
               s=50, alpha=0.7, edgecolors="k")
    for _, row in merged.iterrows():
        ax.annotate(row["Gene"], (row["diff"], row["sum_norm_importance"]),
                    fontsize=7, alpha=0.8)

    ax.set_xlabel("CRISPR β-score Differential Effect", fontsize=11)
    ax.set_ylabel("Sum of Normalised ML Feature Importances", fontsize=11)
    ax.set_title(f"CRISPR β-score vs ML Importance (r = {r:.2f}, p = {p:.2e})",
                 fontsize=12, fontweight="bold")
    plt.tight_layout()
    plt.show()


# =============================================================================

if __name__ == "__main__":
    ccle_data_path = os.path.join(DATA_DIR, "ccle_gadeta_expression.csv")
    if not os.path.exists(ccle_data_path):
        ccle_data_path = os.path.join(DATA_DIR, "ccle_gadeta_expression_demo.csv")
    if not os.path.exists(ccle_data_path):
        raise FileNotFoundError(f"Input data not found in {DATA_DIR}")
    if "_demo.csv" in ccle_data_path:
        print("Using demo data:", ccle_data_path)

    df = pd.read_csv(ccle_data_path)
    gene_cols = [c for c in df.columns if c not in ("cell_line", "resp_status")]

    y = df["resp_status"]
    X = df[gene_cols]

    fitted, scaler = train_all_classifiers(X, y)

    # Sup Fig 4C — Feature importance
    norm_imps = plot_feature_importance(fitted, gene_cols)

    # Sup Fig 4D — Sum of normalised features vs CRISPR β-score
    crispr_path = os.path.join(DATA_DIR, "pos_genes_gdt_vs_nc_1_1_mle_combined.csv")
    if not os.path.exists(crispr_path):
        crispr_path = os.path.join(DATA_DIR, "pos_genes_gdt_vs_nc_1_1_mle_combined_demo.csv")
    if os.path.exists(crispr_path) and norm_imps:
        crispr_df = pd.read_csv(crispr_path)
        crispr_df.columns = ["Gene", "diff"]
        plot_sum_features_vs_beta(crispr_df, norm_imps)
