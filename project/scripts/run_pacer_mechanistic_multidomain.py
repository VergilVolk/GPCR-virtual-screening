"""Cross-publication validation of assay-conditional M4 PAM mechanism models."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from lightgbm import LGBMRegressor
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from scipy.stats import spearmanr
from sklearn.cross_decomposition import PLSRegression
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_absolute_error
from sklearn.model_selection import KFold
from sklearn.preprocessing import StandardScaler

from evaluate_multidomain_sar import features as molecular_features
from run_pacer_fs_baselines import features as frozen_features


P = Path(__file__).resolve().parents[1]
DATA = P / "data" / "benchmarks" / "m4_pam_v1" / "mechanistic_chembl" / "mechanistic_multidomain.csv"
HISTORY = P / "data" / "benchmarks" / "m4_pam_v1" / "potency_molecules.csv"
MONASH = P / "data" / "benchmarks" / "m4_pam_v1" / "external_monash_ly2033298" / "external_allostery.csv"
OUT = P / "results" / "pacer_mechanistic_multidomain_v01"
ALPHAS = (0.1, 1.0, 10.0, 100.0)
PATHWAYS = ("GoB", "cAMP", "arrestin")
PARAMETERS = ("log_tauB", "log_alpha_beta")
CONSENSUS = ("pKB", "consensus_log_tauB", "consensus_log_alpha_beta")
SEED = 20260830


def rho(y, prediction):
    value = spearmanr(y, prediction).statistic
    return float(value) if np.isfinite(value) else 0.0


def frozen_prior(history, frames):
    smiles = pd.concat([history.canonical_smiles] + [item.canonical_smiles for item in frames], ignore_index=True)
    X, D, bitvectors = frozen_features(smiles)
    n = len(history)
    selected = np.argsort(X[:n].var(axis=0))[-1024:]
    design = np.c_[X[:, selected], D]
    model = LGBMRegressor(
        n_estimators=350, learning_rate=.03, num_leaves=15, max_depth=5,
        min_child_samples=12, reg_alpha=.1, reg_lambda=1,
        random_state=11, verbosity=-1, n_jobs=6,
    ).fit(design[:n], history.pEC50.to_numpy(float))
    prediction = model.predict(design[n:])
    output, cursor = [], 0
    for frame in frames:
        output.append(prediction[cursor:cursor + len(frame)])
        cursor += len(frame)
    return output


def splits(n):
    return KFold(n_splits=min(4, n), shuffle=True, random_state=SEED).split(np.arange(n))


def fit_ridge(X, y, train, query, alpha):
    scaler = StandardScaler().fit(X[train])
    model = Ridge(alpha=alpha, solver="lsqr").fit(scaler.transform(X[train]), y[train])
    return model.predict(scaler.transform(X[query]))


def select_ridge(X, y):
    best, loss = ALPHAS[0], np.inf
    for alpha in ALPHAS:
        observed, predicted = [], []
        for train, query in splits(len(y)):
            predicted.extend(fit_ridge(X, y, train, query, alpha)); observed.extend(y[query])
        value = mean_absolute_error(observed, predicted)
        if value < loss:
            best, loss = alpha, value
    return best


def fingerprints(smiles):
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    return [generator.GetFingerprint(Chem.MolFromSmiles(item)) for item in smiles]


def nearest_predict(train_smiles, train_y, query_smiles):
    train_fp, query_fp = fingerprints(train_smiles), fingerprints(query_smiles)
    output = []
    for fp in query_fp:
        similarity = DataStructs.BulkTanimotoSimilarity(fp, train_fp)
        output.append(float(train_y[int(np.argmax(similarity))]))
    return np.asarray(output)


def long_parameter(frame, parameter):
    rows = []
    for index, item in frame.iterrows():
        for pathway in PATHWAYS:
            column = f"{pathway}_{parameter}"
            if column in frame.columns and pd.notna(item.get(column)):
                rows.append({"molecule_index": index, "pathway": pathway, "value": float(item[column])})
    return pd.DataFrame(rows)


def conditional_design(molecule_X, long_rows):
    base = molecule_X[long_rows.molecule_index.to_numpy(int)]
    blocks = [base]
    for pathway in PATHWAYS:
        indicator = (long_rows.pathway.to_numpy() == pathway).astype(float)[:, None]
        blocks.append(base * indicator)
    onehot = np.column_stack([(long_rows.pathway.to_numpy() == item).astype(float) for item in PATHWAYS])
    blocks.append(onehot)
    return np.column_stack(blocks)


def molecule_cv_indices(long_rows, n_molecules):
    for molecule_train, molecule_query in splits(n_molecules):
        train = np.where(long_rows.molecule_index.isin(molecule_train))[0]
        query = np.where(long_rows.molecule_index.isin(molecule_query))[0]
        if len(train) and len(query):
            yield train, query


def select_conditional(design, target, long_rows, n_molecules):
    best, loss = ALPHAS[0], np.inf
    for alpha in ALPHAS:
        observed, predicted = [], []
        for train, query in molecule_cv_indices(long_rows, n_molecules):
            predicted.extend(fit_ridge(design, target, train, query, alpha)); observed.extend(target[query])
        value = mean_absolute_error(observed, predicted)
        if value < loss:
            best, loss = alpha, value
    return best


def predict_conditional(train_frame, query_frame, parameter, molecule_X_train, molecule_X_query):
    train_long = long_parameter(train_frame.reset_index(drop=True), parameter)
    design_train = conditional_design(molecule_X_train, train_long)
    alpha = select_conditional(design_train, train_long.value.to_numpy(float), train_long, len(train_frame))
    scaler = StandardScaler().fit(design_train)
    model = Ridge(alpha=alpha, solver="lsqr").fit(scaler.transform(design_train), train_long.value)
    predictions = {}
    for pathway in PATHWAYS:
        query_long = pd.DataFrame({"molecule_index": np.arange(len(query_frame)), "pathway": pathway,
                                   "value": np.zeros(len(query_frame))})
        design_query = conditional_design(molecule_X_query, query_long)
        predictions[pathway] = model.predict(scaler.transform(design_query))
    return predictions, alpha


def fit_pls(X, prior, Y, train, query, components):
    design = np.c_[X, prior]
    scaler = StandardScaler().fit(design[train])
    components = max(1, min(components, len(train) - 1, Y.shape[1]))
    model = PLSRegression(n_components=components, scale=True, max_iter=1000)
    model.fit(scaler.transform(design[train]), Y[train])
    return model.predict(scaler.transform(design[query]))


def select_pls(X, prior, Y):
    scale = Y.std(axis=0); scale[scale < 1e-8] = 1.0
    best, loss = 1, np.inf
    for components in (1, 2, 3):
        values = []
        for train, query in splits(len(Y)):
            prediction = fit_pls(X, prior, Y, train, query, components)
            values.extend(np.mean(np.abs(prediction - Y[query]) / scale, axis=1))
        value = float(np.mean(values))
        if value < loss:
            best, loss = components, value
    return best


def cross_document(frame, history_prior):
    rows, fold_rows = [], []
    for train_document, query_document in (("Jorg2023", "Liu2024"), ("Liu2024", "Jorg2023")):
        train_frame = frame[frame.document == train_document].reset_index(drop=True)
        query_frame = frame[(frame.document == query_document) & (~frame.exact_other_document_overlap)].reset_index(drop=True)
        train_X = molecular_features(train_frame.canonical_smiles)
        query_X = molecular_features(query_frame.canonical_smiles)
        all_X = np.vstack([train_X, query_X])
        train_index = np.arange(len(train_frame)); query_index = np.arange(len(train_frame), len(train_frame) + len(query_frame))
        prior_train = history_prior[frame.document.to_numpy() == train_document]
        prior_query = history_prior[(frame.document.to_numpy() == query_document) & (~frame.exact_other_document_overlap.to_numpy())]

        # Shared consensus tasks and PACER-Triad.
        train_Y = train_frame[list(CONSENSUS)].to_numpy(float)
        components = select_pls(train_X, prior_train, train_Y)
        triad = fit_pls(all_X, np.r_[prior_train, prior_query], train_Y, train_index, query_index, components)
        for endpoint_index, endpoint in enumerate(CONSENSUS):
            y_train = train_frame[endpoint].to_numpy(float)
            alpha = select_ridge(train_X, y_train)
            pooled = fit_ridge(all_X, y_train, train_index, query_index, alpha)
            nearest = nearest_predict(train_frame.canonical_smiles, y_train, query_frame.canonical_smiles)
            methods = {"Historical-Absolute": prior_query, "Molecular-1NN": nearest,
                       "Pooled-Ridge": pooled, "PACER-Triad": triad[:, endpoint_index]}
            for i, item in query_frame.iterrows():
                row = {"train_document": train_document, "query_document": query_document,
                       "molecule_chembl_id": item.molecule_chembl_id, "endpoint": endpoint,
                       "observed": float(item[endpoint])}
                row.update({name: float(value[i]) for name, value in methods.items()})
                rows.append(row)
            fold_rows.append({"train_document": train_document, "query_document": query_document,
                              "endpoint": endpoint, "pooled_alpha": alpha, "triad_components": components})

        # Exact common-pathway transfer with an explicit assay-conditional branch.
        for parameter in PARAMETERS:
            conditional, conditional_alpha = predict_conditional(
                train_frame, query_frame, parameter, train_X, query_X
            )
            for pathway in PATHWAYS:
                endpoint = f"{pathway}_{parameter}"
                if endpoint not in query_frame or query_frame[endpoint].notna().sum() == 0:
                    continue
                y_train_long = long_parameter(train_frame, parameter)
                pooled_X = train_X[y_train_long.molecule_index.to_numpy(int)]
                pooled_y = y_train_long.value.to_numpy(float)
                pooled_alpha = select_ridge(pooled_X, pooled_y)
                design = np.vstack([pooled_X, query_X])
                pooled = fit_ridge(design, pooled_y, np.arange(len(pooled_y)),
                                   np.arange(len(pooled_y), len(pooled_y) + len(query_frame)), pooled_alpha)
                pathway_train = train_frame[endpoint].notna() if endpoint in train_frame else np.zeros(len(train_frame), bool)
                if np.sum(pathway_train):
                    nearest = nearest_predict(train_frame.loc[pathway_train, "canonical_smiles"],
                                              train_frame.loc[pathway_train, endpoint].to_numpy(float),
                                              query_frame.canonical_smiles)
                else:
                    nearest = np.repeat(float(np.mean(pooled_y)), len(query_frame))
                mask = query_frame[endpoint].notna().to_numpy()
                for i in np.where(mask)[0]:
                    rows.append({"train_document": train_document, "query_document": query_document,
                                 "molecule_chembl_id": query_frame.loc[i, "molecule_chembl_id"],
                                 "endpoint": endpoint, "observed": float(query_frame.loc[i, endpoint]),
                                 "Molecular-1NN": float(nearest[i]), "Pooled-Ridge": float(pooled[i]),
                                 "PACER-AssayConditional": float(conditional[pathway][i])})
                fold_rows.append({"train_document": train_document, "query_document": query_document,
                                  "endpoint": endpoint, "pooled_alpha": pooled_alpha,
                                  "conditional_alpha": conditional_alpha})
    return pd.DataFrame(rows), pd.DataFrame(fold_rows)


def metric_table(predictions):
    rows = []
    identity = {"train_document", "query_document", "molecule_chembl_id", "endpoint", "observed"}
    methods = [item for item in predictions.columns if item not in identity]
    for endpoint, block in predictions.groupby("endpoint"):
        for method in methods:
            valid = block[method].notna()
            if valid.sum() < 3:
                continue
            query = block[valid]
            direction = [rho(part.observed, part[method]) for _, part in query.groupby("query_document") if len(part) >= 3]
            rows.append({"endpoint": endpoint, "method": method, "n": len(query),
                         "Spearman": rho(query.observed, query[method]),
                         "MAE": float(mean_absolute_error(query.observed, query[method])),
                         "worst_direction_Spearman": min(direction) if direction else np.nan})
    return pd.DataFrame(rows)


def paired_bootstrap(predictions, endpoint, candidate, reference, n=20000):
    block = predictions[predictions.endpoint == endpoint].dropna(subset=[candidate, reference]).reset_index(drop=True)
    rng = np.random.default_rng(SEED); values = []
    directions = block.query_document.unique()
    for _ in range(n):
        index = []
        for direction in directions:
            available = np.where(block.query_document.to_numpy() == direction)[0]
            index.extend(rng.choice(available, len(available), replace=True))
        sample = block.iloc[index]
        values.append(rho(sample.observed, sample[candidate]) - rho(sample.observed, sample[reference]))
    values = np.asarray(values)
    estimate = rho(block.observed, block[candidate]) - rho(block.observed, block[reference])
    return {"estimate": float(estimate), "ci95": [float(x) for x in np.quantile(values, [.025, .975])],
            "probability_gt_zero": float(np.mean(values > 0)), "n_boot": n}


def monash_exploratory(frame, monash, history_prior_train, history_prior_monash):
    train_X = molecular_features(frame.canonical_smiles); query_X = molecular_features(monash.canonical_smiles)
    all_X = np.vstack([train_X, query_X]); train = np.arange(len(frame)); query = np.arange(len(frame), len(frame) + len(monash))
    Y = frame[list(CONSENSUS)].to_numpy(float)
    components = select_pls(train_X, history_prior_train, Y)
    triad = fit_pls(all_X, np.r_[history_prior_train, history_prior_monash], Y, train, query, components)
    mapping = {"functional_pKB": 0, "log_tauB": 1, "log_alpha_beta": 2}
    rows = []
    eligible = (monash.pam_label == 1) & (~monash.exact_history_overlap)
    for endpoint, index in mapping.items():
        alpha = select_ridge(train_X, Y[:, index])
        pooled = fit_ridge(all_X, Y[:, index], train, query, alpha)
        nearest = nearest_predict(frame.canonical_smiles, Y[:, index], monash.canonical_smiles)
        for method, prediction in {"Historical-Absolute": history_prior_monash,
                                   "Mechanistic-1NN": nearest, "Mechanistic-Pooled-Ridge": pooled,
                                   "PACER-Triad": triad[:, index]}.items():
            rows.append({"endpoint": endpoint, "method": method, "n": int(eligible.sum()),
                         "Spearman": rho(monash.loc[eligible, endpoint], prediction[eligible]),
                         "MAE": float(mean_absolute_error(monash.loc[eligible, endpoint], prediction[eligible]))})
    return pd.DataFrame(rows)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    frame = pd.read_csv(DATA); history = pd.read_csv(HISTORY); monash = pd.read_csv(MONASH)
    prior_frame, prior_monash = frozen_prior(history, [frame, monash])
    predictions, folds = cross_document(frame, prior_frame)
    metrics = metric_table(predictions)
    monash_metrics = monash_exploratory(frame, monash, prior_frame, prior_monash)
    predictions.to_csv(OUT / "cross_document_predictions.csv", index=False)
    metrics.to_csv(OUT / "cross_document_metrics.csv", index=False)
    folds.to_csv(OUT / "folds.csv", index=False)
    monash_metrics.to_csv(OUT / "monash_exploratory_metrics.csv", index=False)

    primary = metrics[metrics.endpoint == "consensus_log_alpha_beta"].sort_values("Spearman", ascending=False)
    best_baseline = primary[primary.method != "PACER-Triad"].iloc[0].method
    delta = paired_bootstrap(predictions, "consensus_log_alpha_beta", "PACER-Triad", best_baseline)
    candidate = primary[primary.method == "PACER-Triad"].iloc[0]
    baseline = primary[primary.method == best_baseline].iloc[0]
    supported = bool(delta["ci95"][0] > 0 and candidate.worst_direction_Spearman >= baseline.worst_direction_Spearman)
    audit = {"status": "cross_publication_mechanistic_development_benchmark",
             "primary_endpoint": "consensus_log_alpha_beta", "primary_method": "PACER-Triad",
             "best_baseline": best_baseline, "paired_bootstrap_delta": delta,
             "preregistered_support": supported, "n_exact_cross_document_queries_excluded": 2,
             "monash_role": "post_failure_exploratory_endpoint_transfer_not_untouched_external",
             "claim_boundary": "Cross-publication computational mechanism transfer; no prospective or wet-lab PAM confirmation."}
    (OUT / "audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    report = ["# PACER-M4 cross-publication mechanistic benchmark", "",
              f"- Primary endpoint: consensus log(alpha*beta).",
              f"- PACER-Triad vs best baseline ({best_baseline}) delta: {delta['estimate']:+.3f} "
              f"(95% CI {delta['ci95'][0]:+.3f} to {delta['ci95'][1]:+.3f}).",
              f"- Preregistered support: **{supported}**.", "",
              "## Cross-publication results", "",
              "| Endpoint | Method | N | Spearman | MAE | Worst direction |", "|---|---|---:|---:|---:|---:|"]
    for row in metrics.sort_values(["endpoint", "Spearman"], ascending=[True, False]).itertuples():
        report.append(f"| {row.endpoint} | {row.method} | {int(row.n)} | {row.Spearman:.3f} | {row.MAE:.3f} | {row.worst_direction_Spearman:.3f} |")
    report += ["", "## Monash exploratory transfer", "",
               "This analysis was motivated by the first Monash failure and is not an untouched external test.", "",
               "| Endpoint | Method | N | Spearman | MAE |", "|---|---|---:|---:|---:|"]
    for row in monash_metrics.sort_values(["endpoint", "Spearman"], ascending=[True, False]).itertuples():
        report.append(f"| {row.endpoint} | {row.method} | {int(row.n)} | {row.Spearman:.3f} | {row.MAE:.3f} |")
    report += ["", "No result here confirms a PAM; all operational parameters remain probe- and pathway-conditional."]
    (OUT / "REPORT.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    print(metrics.to_string(index=False)); print(monash_metrics.to_string(index=False)); print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
