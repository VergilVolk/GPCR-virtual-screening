"""Manifest-driven, apply-only runners for the frozen PACER-FKG v02 pipeline.

This package adds NO scientific numerics of its own.  Every frozen quantity is
imported from the byte-frozen v02 modules (see `engines.py`); the modules here
only supply:

  * a fail-closed run-spec loader/validator (`spec.py`),
  * an output-path policy that can never touch a protected historical root
    (`paths.py`),
  * a parameterised temporal-block constructor plus bit-identity proofs
    against both frozen implementations (`blocks.py`),
  * the three generic phase drivers (`fkg_phase1.py`, `fkg_phase2a.py`,
    `fkg_phase2b.py`) and a single CLI (`cli.py`),
  * the CM00734 payload-identity regression harness (`regression.py`).
"""

from __future__ import annotations

__all__ = ["spec", "engines", "blocks", "paths", "fkg_phase1", "fkg_phase2a", "fkg_phase2b", "cli", "regression"]

VERSION = "PACER_FKG_V02_FROZEN_RUNNER_v01"
