"""Safe loading of frozen .npy artifacts.

A subtle hazard: `np.asarray(memmap[...])` returns a view whose buffer is the
memory map itself.  Closing the `.\_mmap` while such a view is still alive is a
use-after-free and crashes the interpreter with an access violation
(0xC0000005).  Every read of a frozen .npy in this package therefore goes
through `load_npy_copy`, which materialises an independent array before the
mapping is released.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np


def load_npy_copy(path: Path | str) -> np.ndarray:
    """Load a .npy fully into owned memory, then release the mapping."""
    array = np.load(Path(path), mmap_mode="r", allow_pickle=False)
    try:
        return np.array(array, dtype=array.dtype, copy=True)
    finally:
        array._mmap.close()


def check_npy(path: Path | str, expected_shape: tuple[int, ...], expected_dtype=np.float32) -> np.ndarray:
    """Load and validate shape/dtype/finiteness; return the owned array."""
    array = load_npy_copy(path)
    if array.shape != expected_shape:
        raise ValueError(f"{path}: shape {array.shape} != {expected_shape}")
    if array.dtype != expected_dtype:
        raise ValueError(f"{path}: dtype {array.dtype} != {expected_dtype}")
    if not np.isfinite(array).all():
        raise ValueError(f"{path}: contains non-finite values")
    return array
