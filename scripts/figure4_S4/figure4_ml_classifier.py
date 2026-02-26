#!/usr/bin/env python3
"""
Figure 4A–D — ML Classifier for GDT201 Reactivity

  Fig 4A: Distribution of CCLE cell lines by # models classifying as responder (like plot_pred_with_indication_distribution)
  Fig 4B: Among CCLE with ≥6 models: count plot by indication (like plot_distribution on df_ccle_filtered_6)
  Fig 4C: Distribution of Normal cell lines by # models classifying as responder
  Fig 4D: Among Normal with ≥3 models: count plot by sub_anatomical_region1 (like plot_distribution on df_norm_filtered_2)

--- Input data (tidy) ---------------------------------------------------------
  One table: rows = cell lines, columns = cell_line + gene names + resp_status.
  resp_status: 0 = non-responder, 1 = responder. Optional metadata: cell_line → indication.

  cell_line   B3GLCT   FUT4   ...   TM9SF3   resp_status
  786O        12.16    14.14   ...  132.29    0
  AGS         8.12     37.18   ...  316.98    1

  Data: demo/figure4_S4/. Files: ccle_gadeta_expression.csv (or _demo.csv), optional ccle_indications.csv, normal_gadeta_expression.csv, normal_indications.csv.

--- Feature engineering -------------------------------------------------------
  Features = gene columns only; target = resp_status. StandardScaler fit on training
  data, transform train and any prediction set (no other feature construction).

--- Model building ------------------------------------------------------------
  Seven classifiers: Logistic Regression, XGBoost, Decision Tree, Random Forest,
  GB Trees, AdaBoost, Naive Bayes. Stratified train/test split (75/25, random_state=0);
  fit on X_train, evaluate on X_test (accuracy, ROC-AUC, recall, precision, F1).
  This script uses tuned hyperparameters via train_tuned_classifiers() (see below).

--- Hyperparameter tuning -----------------------------------------------------
  train_tuned_classifiers() below: GridSearchCV with StratifiedKFold(n_splits=5),
  scoring='accuracy'. Per-model param_grids (LogisticRegression: C, penalty;
  RF/GB/XGB/Adaboost: n_estimators, max_depth, learning_rate, etc.). Best
  estimator per model refit on full training fold and returned for prediction.

--- Outputs -------------------------------------------------------------------
  Per–cell line y_pred / y_prob from each model; consensus = majority vote.
  Plots: 4A/4C full distribution by # models; 4B count by indication (CCLE); 4D count by sub_anatomical_region1 (normal). No NA drop for count plots.

Usage:
  conda activate gdt201
  python figure4_ml_classifier.py
"""

import os
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.model_selection import train_test_split, StratifiedKFold, GridSearchCV
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    roc_auc_score, recall_score, precision_score,
    accuracy_score, f1_score,
)
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


def _input_path(filename):
    # Prefer non-demo name, then demo name (both in demo/figure4_S4/)
    main = os.path.join(DATA_DIR, filename)
    if os.path.exists(main):
        return main
    demo_name = filename.replace(".csv", "_demo.csv") if "demo" not in filename else filename
    demo = os.path.join(DATA_DIR, demo_name)
    if os.path.exists(demo):
        return demo
    return main  # let caller handle missing

CLASSIFIERS = {
    "Logistic Regression": LogisticRegression(random_state=0),
    "XGBoost": XGBClassifier(random_state=0, eval_metric="logloss"),
    "Decision Tree": DecisionTreeClassifier(random_state=0),
    "Random Forest": RandomForestClassifier(random_state=0),
    "GB Trees": GradientBoostingClassifier(random_state=0),
    "Adaboost": AdaBoostClassifier(random_state=0, algorithm="SAMME"),
    "Naive Bayes": GaussianNB(),
}

PREDICTION_COLUMNS = [
    "logistic_regression_y_pred", "xgboost_y_pred", "decision_tree_y_pred",
    "gb_trees_y_pred", "random_forest_y_pred", "adaboost_y_pred",
    "naive_bayes_y_pred",
]

# Hyperparameter tuning: param_grids per model (from gadeta_classifier_ccle.ipynb)
TUNED_PARAM_GRIDS = {
    "Logistic Regression": {"C": [0.001, 0.01, 0.1, 1, 10, 100], "penalty": ["l1", "l2"], "solver": ["saga"]},
    "XGBoost": {"learning_rate": [0.01, 0.1, 0.2], "max_depth": [3, 5, 7], "n_estimators": [50, 100, 200]},
    "Decision Tree": {"max_depth": [3, 5, 7], "min_samples_split": [2, 5, 10]},
    "Random Forest": {"n_estimators": [50, 100], "max_depth": [None, 10, 20]},
    "GB Trees": {"learning_rate": [0.01, 0.1, 0.2], "n_estimators": [100, 200], "max_depth": [3, 5, 7]},
    "Adaboost": {"n_estimators": [50, 100, 200], "learning_rate": [0.01, 0.1, 1.0]},
    "Naive Bayes": {},  # no grid
}


def train_all_classifiers(X, y, test_size=0.25, random_state=0):
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, stratify=y, test_size=test_size, random_state=random_state
    )
    scaler = StandardScaler()
    X_train = scaler.fit_transform(X_train)
    X_test = scaler.transform(X_test)

    metrics_rows = []
    fitted = {}
    for name, clf in CLASSIFIERS.items():
        clf.fit(X_train, y_train)
        y_pred = clf.predict(X_test)
        y_prob = clf.predict_proba(X_test)[:, 1]

        metrics_rows.append({
            "model_name": name,
            "accuracy": accuracy_score(y_test, y_pred),
            "roc_auc": roc_auc_score(y_test, y_prob),
            "recall": recall_score(y_test, y_pred),
            "precision": precision_score(y_test, y_pred),
            "f1": f1_score(y_test, y_pred),
        })
        fitted[name] = clf

    metrics_df = pd.DataFrame(metrics_rows).set_index("model_name").round(3)
    return fitted, metrics_df, scaler


def train_tuned_classifiers(X, y, test_size=0.25, random_state=0, cv_splits=5):
    """Train with hyperparameter tuning via GridSearchCV + StratifiedKFold."""
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, stratify=y, test_size=test_size, random_state=random_state
    )
    scaler = StandardScaler()
    X_train = scaler.fit_transform(X_train)
    X_test = scaler.transform(X_test)
    skf = StratifiedKFold(n_splits=cv_splits, shuffle=True, random_state=random_state)
    fitted = {}
    metrics_rows = []
    for name, base_clf in CLASSIFIERS.items():
        param_grid = TUNED_PARAM_GRIDS.get(name, {})
        if param_grid:
            search = GridSearchCV(
                base_clf, param_grid, cv=skf, scoring="accuracy", refit=True
            )
            search.fit(X_train, y_train)
            clf = search.best_estimator_
        else:
            base_clf.fit(X_train, y_train)
            clf = base_clf
        fitted[name] = clf
        y_pred = clf.predict(X_test)
        y_prob = clf.predict_proba(X_test)[:, 1]
        metrics_rows.append({
            "model_name": name,
            "accuracy": accuracy_score(y_test, y_pred),
            "roc_auc": roc_auc_score(y_test, y_prob),
            "recall": recall_score(y_test, y_pred),
            "precision": precision_score(y_test, y_pred),
            "f1": f1_score(y_test, y_pred),
        })
    metrics_df = pd.DataFrame(metrics_rows).set_index("model_name").round(3)
    return fitted, metrics_df, scaler


def predict_on_dataset(fitted_models, scaler, df, gene_cols):
    X_new = scaler.transform(df[gene_cols])
    result = df.copy()
    for name, clf in fitted_models.items():
        col_base = name.lower().replace(" ", "_")
        result[f"{col_base}_y_pred"] = clf.predict(X_new)
        result[f"{col_base}_y_prob"] = clf.predict_proba(X_new)[:, 1]
    return result


def add_consensus(df, pred_cols=None):
    if pred_cols is None:
        pred_cols = [c for c in df.columns if c.endswith("_y_pred")]
    votes = df[pred_cols].sum(axis=1)
    n_models = len(pred_cols)
    df["consensus"] = np.where(
        votes > n_models / 2, 1,
        np.where(votes < n_models / 2, 0, -1)
    )
    return df


# -- Figure 4A / 4C (notebook: plot_pred_with_indication_distribution) -------------

def plot_responder_model_count(df, prediction_columns, title):
    """Distribution of cell lines by # models that classified as responder. No NA drop."""
    df = df.copy()
    df["count_ones"] = df[prediction_columns].sum(axis=1)
    counts = df["count_ones"].value_counts().sort_index()

    fig, ax = plt.subplots(figsize=(10, 6))
    sns.barplot(x=counts.index, y=counts.values, hue=counts.index, palette="coolwarm", ax=ax)

    total = counts.sum()
    for idx, val in enumerate(counts):
        ax.text(idx, val, f"{(val / total) * 100:.2f}%\n({val})", ha="center", va="bottom", fontsize=9)

    ax.set_title(f"Distribution of {title} Cell Lines Responder by Different Models", fontsize=13, fontweight="bold")
    ax.set_xlabel("Number of Models that Classified one cell line as Responder")
    ax.set_ylabel("Number of Cell Lines")
    ax.set_ylim(0, total / 2)
    ax.get_legend().remove()
    plt.tight_layout()
    plt.show()


# -- Figure 4B (notebook: df_ccle_filtered_6, plot_distribution by indication) ------
# -- Figure 4D (notebook: df_norm_filtered_2, plot_distribution by var) -------------
# No NA drop: include all rows; missing indication appears as a bar if present.

def plot_distribution_by_var(df, var, title, palette="coolwarm"):
    """Count plot by grouping variable (e.g. indication). Matches notebook plot_distribution. No NA drop."""
    fig, ax = plt.subplots(figsize=(12, 6))
    sns.countplot(x=var, data=df, palette=palette, ax=ax)
    total = len(df)
    for p in ax.patches:
        count = p.get_height()
        ax.annotate(f"{(count / total) * 100:.2f}%\n({int(count)})",
                    (p.get_x() + p.get_width() / 2.0, p.get_height()),
                    ha="center", va="bottom", fontsize=9)
    ax.set_title(f"Distribution of {var}", fontsize=13, fontweight="bold")
    ax.set_xlabel(var)
    ax.set_ylabel("Count")
    ax.set_ylim(0, total / 2 if total else 1)
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.show()


def plot_responder_rate_by_indication(df, indication_col, title):
    """Optional: responder rate per indication. Drops NA only for groupby (so rate is per known indication)."""
    if indication_col not in df.columns or df[indication_col].isna().all():
        return
    if "consensus" not in df.columns:
        df = add_consensus(df)
    # Drop NA only for rate calculation so we don't get a NaN row in the table
    rates = (
        df.dropna(subset=[indication_col])
        .groupby(indication_col)["consensus"]
        .apply(lambda s: (s == 1).mean())
        .sort_values(ascending=False)
        .reset_index()
    )
    rates.columns = [indication_col, "responder_rate"]
    fig, ax = plt.subplots(figsize=(14, 7))
    sns.barplot(data=rates, x=indication_col, y="responder_rate", color="#de2d26", ax=ax)
    for idx, row in rates.iterrows():
        ax.text(idx, row["responder_rate"] + 0.02, f"{row['responder_rate']:.0%}", ha="center", fontsize=8)
    ax.set_title(f"Responder Rate by {title}", fontsize=13, fontweight="bold")
    ax.set_xlabel(title)
    ax.set_ylabel("Fraction Classified as Responder")
    ax.set_ylim(0, 1.1)
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.show()


# =============================================================================

if __name__ == "__main__":
    # Assume tidy input: cell_line + gene columns + resp_status (see docstring).
    ccle_data_path = _input_path("ccle_gadeta_expression.csv")
    if not os.path.exists(ccle_data_path):
        raise FileNotFoundError(f"Input not found: {ccle_data_path} (nor in {DATA_DIR})")
    if "_demo.csv" in ccle_data_path:
        print("Using demo data:", ccle_data_path)

    df = pd.read_csv(ccle_data_path)
    gene_cols = [c for c in df.columns if c not in ("cell_line", "resp_status")]
    X = df[gene_cols]
    y = df["resp_status"]

    # Train with hyperparameter tuning (GridSearchCV + StratifiedKFold)
    fitted, metrics, scaler = train_tuned_classifiers(X, y)

    df_pred = predict_on_dataset(fitted, scaler, df, gene_cols)
    df_pred = add_consensus(df_pred)

    # Fig 4A — CCLE: distribution by # models that classified as responder (notebook: plot_pred_with_indication_distribution)
    df_pred["count_ones"] = df_pred[PREDICTION_COLUMNS].sum(axis=1)
    plot_responder_model_count(df_pred, PREDICTION_COLUMNS, "CCLE")

    # Fig 4B — Among CCLE with ≥6 models: count plot by indication (column "indication" only)
    indication_path = _input_path("ccle_indications.csv")
    if os.path.exists(indication_path):
        df_indications = pd.read_csv(indication_path)
        df_ccle_ind = df_pred.merge(df_indications, on="cell_line", how="left")
        df_ccle_filtered_6 = df_ccle_ind[df_ccle_ind["count_ones"] >= 6]
        if len(df_ccle_filtered_6) > 0 and "indication" in df_ccle_ind.columns:
            plot_distribution_by_var(df_ccle_filtered_6, "indication", title="CCLE (≥6 models as responder)")
        elif len(df_ccle_filtered_6) > 0:
            plot_distribution_by_var(df_ccle_filtered_6.assign(grp="≥6 models"), "grp", title="CCLE (≥6 models)")
    else:
        df_ccle_filtered_6 = df_pred[df_pred["count_ones"] >= 6]
        if len(df_ccle_filtered_6) > 0:
            plot_distribution_by_var(df_ccle_filtered_6.assign(grp="≥6 models"), "grp", title="CCLE (≥6 models)")

    # Fig 4C — Normal: distribution by # models that classified as responder
    normal_path = _input_path("normal_gadeta_expression.csv")
    if os.path.exists(normal_path):
        df_normal = pd.read_csv(normal_path)
        df_normal_pred = predict_on_dataset(fitted, scaler, df_normal, gene_cols)
        df_normal_pred = add_consensus(df_normal_pred)
        df_normal_pred["count_ones"] = df_normal_pred[PREDICTION_COLUMNS].sum(axis=1)
        plot_responder_model_count(df_normal_pred, PREDICTION_COLUMNS, "Normal")

        # Fig 4D — Among Normal with ≥3 models: count plot by sub_anatomical_region1 only
        normal_meta_path = _input_path("normal_indications.csv")  # must have cell_line, sub_anatomical_region1
        df_norm_filtered_3 = df_normal_pred[df_normal_pred["count_ones"] >= 3]
        if len(df_norm_filtered_3) > 0 and os.path.exists(normal_meta_path):
            df_norm_meta = pd.read_csv(normal_meta_path)
            df_norm_merged = df_norm_filtered_3.merge(df_norm_meta, on="cell_line", how="left")
            if "sub_anatomical_region1" in df_norm_merged.columns:
                plot_distribution_by_var(df_norm_merged, "sub_anatomical_region1", title="Normal (≥3 models as responder)")
            else:
                plot_distribution_by_var(df_norm_filtered_3.assign(grp="≥3 models"), "grp", title="Normal (≥3 models)")
        elif len(df_norm_filtered_3) > 0:
            plot_distribution_by_var(df_norm_filtered_3.assign(grp="≥3 models"), "grp", title="Normal (≥3 models)")

    # Optional: responder rate by cancer type (uses column "indication")
    if os.path.exists(indication_path):
        df_ind = pd.read_csv(indication_path)
        if "indication" in df_ind.columns:
            df_merged = df_pred.merge(df_ind, on="cell_line", how="left")
            plot_responder_rate_by_indication(df_merged, "indication", "Cancer Type")
