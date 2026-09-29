#!/usr/bin/env python3
"""Apply-only wrapper: run the frozen CGDA LOSO runner on the Science-2026
pickle root. Patches pocket_map to load the pockets dict built by
build_cgda_root_science2026.py instead of the historical LMDB archive layout.
No code inside evaluate_drugclip_cgda_loso.py is modified."""
from __future__ import annotations
import pickle, sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import evaluate_drugclip_cgda_loso as cgda

def _pocket_map(root, archive_path):
    with open(archive_path, "rb") as f:
        return pickle.load(f)

cgda.pocket_map = _pocket_map
cgda.main()
