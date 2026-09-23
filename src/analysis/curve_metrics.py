from __future__ import annotations

import numpy as np
from numpy.typing import NDArray


def evaluate_curve(
    scores: NDArray[np.float64],
    y_true: NDArray[np.int64],
) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64], NDArray[np.float64], float, float]:
    """Return deterministic ROC and PR curve points plus AUCs.

    The curve is built from the raw score vector without sklearn.
    PR AUC uses the same stepwise semantics as ``scoring.pr_auc``.
    """
    order = np.argsort(scores, kind="stable")[::-1]
    y_true = y_true[order].astype(np.int64)
    scores = scores[order]

    n_pos = int(y_true.sum())
    n_neg = int((1 - y_true).sum())

    if n_pos == 0 or n_neg == 0:
        fpr = np.array([0.0, 1.0], dtype=np.float64)
        tpr = np.array([0.0, 1.0], dtype=np.float64)
        precision = np.array([1.0, 0.0], dtype=np.float64)
        recall = np.array([0.0, 1.0], dtype=np.float64)
        return fpr, tpr, recall, precision, 0.5, 0.5

    tps = np.cumsum(y_true)
    fps = np.cumsum(1 - y_true)

    tpr = np.r_[0.0, tps / n_pos]
    fpr = np.r_[0.0, fps / n_neg]

    precision = np.r_[1.0, tps / np.arange(1, len(tps) + 1), n_pos / (n_pos + n_neg)]
    recall = np.r_[0.0, tps / n_pos, 1.0]

    roc_auc = float(np.trapz(tpr, fpr))
    pr_area = 0.0
    prev_recall = 0.0
    for index in range(1, len(tps) + 1):
        hits = int(tps[index - 1])
        recall_val = hits / n_pos
        precision_val = hits / index
        pr_area += (recall_val - prev_recall) * precision_val
        prev_recall = recall_val
    return fpr, tpr, recall, precision, roc_auc, pr_area
