# Surrogate Model Building and Feature Ranking

Surrogate-model training and Shapley (SHAP) feature-ranking analysis on the
prostate cancer (PCa) dataset. A Neural Network, SVM, and Random Forest are
trained as surrogate classifiers, and SHAP is used to rank input features by
their contribution to the predicted output.

- **Models:** Neural Network (Keras), SVM (RBF), Random Forest
- **Input:** all three cohorts (CNT, BR, TF) — boolean node states, ERK/AKT
  concentrations, and PTEN/EGF analog values (22 input features)
- **Output:** binary class, derived by thresholding the NCG value at **0.15**

The workflow has two parts:

1. **[Part 1 — Surrogate model training](#part-1--surrogate-model-training):**
   train and evaluate the NN, SVM, and RF surrogate classifiers
   (`NN_SVM_RF.ipynb`).
2. **[Part 2 — Feature ranking with Shapley](#part-2--feature-ranking-with-shapley):**
   run SHAP on the trained models to rank input features
   (`feature_ranking.py`).

## Repository contents

| File | Description |
|------|-------------|
| `NN_SVM_RF.ipynb` | **Part 1** — trains the surrogate models (NN, SVM, RF) and sweeps NCG thresholds; reports balanced accuracy and produces the accuracy plots. |
| `feature_ranking.py` | **Part 2** — trains NN/SVM/RF at threshold 0.15 and runs the SHAP analysis, saving one summary plot per model. |
| `boolean_input_output_default_run_Testo0.csv` | Input CSV 1 — boolean node states + NCG output (time series). |
| `bool_input_threshold_case_spec_modified.csv` | Input CSV 2 — PTEN / EGF analog values. |
| `pp_SS_boolean_input_threshold_0.05_validation_v1.csv` | Input CSV 3 — ERK_PP / AKT_PP concentrations. |
| `PCa_NN_shap_summaryplot.png` | Output — SHAP feature ranking for the NN. |
| `PCa_SVM_shap_summaryplot.png` | Output — SHAP feature ranking for the SVM. |
| `PCa_RF_shap_summaryplot.png` | Output — SHAP feature ranking for the RF. |

## Setup (required for both parts)

Python 3.8.

The steps below create a **new** virtual environment named `Shap`. You do not
need any pre-existing environment — `python -m venv Shap` builds it from scratch
in the current directory, and `source Shap/bin/activate` switches into it. (You
only create it once; on later sessions just run the `source` line to reactivate.)

```bash
# 1. Create a fresh virtual environment called "Shap" (run once)
python -m venv Shap

# 2. Activate it (run every new terminal session)
source Shap/bin/activate

# 3. Upgrade pip, then install the pinned dependencies
pip install --upgrade pip
pip install \
  tensorflow==2.13.1 \
  numpy==1.24.3 \
  scikit-learn==1.3.2 \
  shap==0.44.1 \
  imbalanced-learn==0.12.0 \
  pandas==2.0.3 \
  matplotlib==3.7.5 \
  scipy==1.10.1
```

You can name the environment anything — `Shap` is just the name used here.
When you are done, run `deactivate` to leave the environment.

The key packages are pinned to the versions this project was run with —
`tensorflow`, `numpy`, and `scikit-learn` in particular are version-sensitive,
so keeping these pins is the safest way to reproduce the results.

The scripts read the three input CSVs from `~/Shapley_PCa/`, so clone/copy this
folder into your home directory (or adjust the `pd.read_csv(...)` paths at the
top of each file).

---

## Part 1 — Surrogate model training

**File:** `NN_SVM_RF.ipynb`

Trains the three surrogate classifiers (Neural Network, SVM, Random Forest) and
evaluates them across a range of NCG thresholds. For each threshold it reports
the average (cross-validated) and test **balanced accuracy** of each model, and
saves the accuracy bar plots.

> The notebook additionally uses `tabulate` for its in-notebook result tables.
> Install it alongside the Setup dependencies: `pip install tabulate`.

### Run

Open the notebook in Jupyter and run the cells top to bottom:

```bash
source Shap/bin/activate
jupyter notebook NN_SVM_RF.ipynb
```

The notebook loads the three input CSVs, builds the per-cohort feature table,
sweeps the NCG thresholds, and produces the balanced-accuracy plots for the
NN, SVM, and RF.

---

## Part 2 — Feature ranking with Shapley

**File:** `feature_ranking.py`

Trains the NN, SVM, and RF at the fixed NCG threshold of **0.15**, then runs a
SHAP `KernelExplainer` on each model to rank the 22 input features by their
contribution to the prediction. One SHAP summary plot is saved per model.

### Run

```bash
source Shap/bin/activate
python feature_ranking.py
```

This writes the three feature-ranking plots:

- `PCa_NN_shap_summaryplot.png`
- `PCa_SVM_shap_summaryplot.png`
- `PCa_RF_shap_summaryplot.png`

On a cluster, wrap this in your usual job-submission script (activating the
environment and running `python feature_ranking.py`).
