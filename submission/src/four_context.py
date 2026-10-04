from __future__ import annotations

import numpy as np


def delta_int(
    apo: np.ndarray,
    probe: np.ndarray,
    candidate: np.ndarray,
    combined: np.ndarray,
) -> np.ndarray:
    """Candidate-by-probe interaction after subtracting both main effects."""
    return combined - probe - candidate + apo


def cosine(left: np.ndarray, right: np.ndarray) -> float:
    denominator = np.linalg.norm(left) * np.linalg.norm(right)
    if denominator == 0:
        return float("nan")
    return float(np.dot(left, right) / denominator)
