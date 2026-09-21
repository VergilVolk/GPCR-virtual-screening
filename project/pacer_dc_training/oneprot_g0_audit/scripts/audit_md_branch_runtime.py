#!/usr/bin/env python
"""Assertive runtime audit of the OneProt MD branch load (G0 step 2 of the protocol).

`load_md_branch.py` proved the load was *possible*: it printed
``missing_keys=0 / unexpected_keys=1``.  This script turns that observation into
an assertion, per `START_HERE_ONEPROT_PACER_DC.md`:

    "must check the return value of load_state_dict, require zero missing keys
     and zero shape mismatches under network.md.transformer.*; do not accept
     strict=false merely because the program did not raise."

Provenance of every claim made here:

* ``strict=False`` does not apply a tensor whose shape disagrees with the
  destination parameter, and it reports that key in ``missing_keys``.  PyTorch's
  return object therefore cannot separate "absent" from "shape-mismatched", so
  this script re-derives both from an explicit per-key shape comparison.
* ``network.md.norm.1.log_logit_scale`` comes from ``LearnableLogitScaling`` in
  ``src/models/components/base_encoder.py`` (line 30: a *registered buffer* when
  ``learnable=False``).  The PACER-DC encoder is built with the default
  ``use_logit_scale=False``, so ``_create_normalization`` returns
  ``nn.Sequential(Normalize(dim=-1))`` only and index 1 does not exist locally.
  ``use_logit_scale`` is a contrastive-training device; ``TrajectoryEncoder.forward``
  ends with ``self.norm(projected)``, i.e. the Normalize layer, so this buffer is
  never read during embedding extraction.
* To prove that difference is purely a constructor flag and not a missing weight,
  the audit runs a second configuration with ``use_logit_scale=True`` and asserts
  an *exact* match (zero missing, zero unexpected, zero shape mismatch).

Exit codes: 0 = every load-bearing assertion passed; 2 = at least one failed.
"""
from __future__ import annotations

import argparse
import json
import sys
import types
from pathlib import Path

REPO_DEFAULT = Path("project/tools/oneprot-embeddings")
CKPT_DEFAULT = (
    REPO_DEFAULT / "artifacts" / "Pocket_Text_ST_SG_MD" / "epoch_012_01100-v1.ckpt"
)
OUTER_AUDIT_DEFAULT = Path("runs/oneprot_checkpoint_audit.json")

# The ONLY tolerated difference, with its justification. Anything outside this
# set fails the audit; the set is asserted to be exact, not merely a subset.
KNOWN_UNPROVIDED = {
    "norm.1.log_logit_scale": (
        "LearnableLogitScaling buffer (base_encoder.py:30), created only when "
        "use_logit_scale=True. The formal G0 encoder runs with the default "
        "use_logit_scale=False, so this buffer has no destination. It is a "
        "contrastive-training temperature, never read by TrajectoryEncoder.forward "
        "(which returns self.norm(projected) = Normalize only). It is a buffer, "
        "not a parameter: no learnable weight is lost. The counterfactual "
        "configuration in this same report loads it exactly."
    )
}


def tensor_facts(value: object) -> dict:
    import torch

    if not isinstance(value, torch.Tensor):
        return {"is_tensor": False, "type": type(value).__name__}
    return {
        "is_tensor": True,
        "shape": list(value.shape),
        "dtype": str(value.dtype),
        "numel": int(value.numel()),
        "element_size": int(value.element_size()),
        "bytes": int(value.numel() * value.element_size()),
        "is_floating_point": bool(value.is_floating_point()),
    }


def audit_configuration(encoder, remapped: dict, label: str, use_logit_scale: bool,
                        tolerated: dict) -> dict:
    """Load `remapped` into a fresh encoder and classify every key explicitly."""
    own = encoder.state_dict()
    result = encoder.load_state_dict(remapped, strict=False)

    matched: list[str] = []
    shape_mismatch: list[dict] = []
    absent: list[dict] = []
    for key, value in remapped.items():
        if key not in own:
            absent.append({"key": key, **tensor_facts(value)})
            continue
        target = own[key]
        if tuple(target.shape) != tuple(value.shape):
            shape_mismatch.append({
                "key": key,
                "checkpoint_shape": list(value.shape),
                "encoder_shape": list(target.shape),
            })
        else:
            matched.append(key)

    not_supplied = sorted(set(own) - set(matched))
    absent_keys = sorted(row["key"] for row in absent)
    tolerated_keys = sorted(tolerated)
    unexpected_keys = sorted(result.unexpected_keys)
    missing_keys = sorted(result.missing_keys)
    mtransformer = [k for k in missing_keys if k == "transformer" or k.startswith("transformer.")]
    mmtransformer = [
        row["key"] for row in shape_mismatch
        if row["key"] == "transformer" or row["key"].startswith("transformer.")
    ]
    unclassified = [k for k in absent_keys if k not in tolerated]

    # torch's `missing_keys` means "the encoder has this key and the supplied dict
    # did not" -- it is the mirror image of `not_supplied`, NOT of `absent`. Keep
    # the two directions straight and assert both.
    supplied_keys = set(remapped)
    encoder_side_uncovered = sorted(supplied_keys - set(matched) - set(tolerated))
    encoder_side_unexpected_supplied = sorted(supplied_keys - set(own))

    checks = {
        # a tensor whose shape disagrees is reported by torch as "missing", so both
        # namespaces must be checked separately.
        "zero_missing_keys_under_transformer": not mtransformer,
        "zero_shape_mismatch_under_transformer": not mmtransformer,
        "zero_shape_mismatch_anywhere": not shape_mismatch,
        # direction 1: every key the encoder owns was supplied by the checkpoint
        "every_encoder_key_was_supplied": not not_supplied,
        # direction 2: every checkpoint key landed somewhere, or is explicitly tolerated
        "absent_keys_are_exactly_the_tolerated_set": absent_keys == tolerated_keys,
        "no_unclassified_absent_keys": not unclassified,
        # cross-checks of the torch return value against the explicit classification
        "load_state_dict_missing_keys_match_encoder_side_direction": (
            missing_keys == not_supplied
        ),
        "load_state_dict_unexpected_keys_are_tolerated_set": unexpected_keys == tolerated_keys,
        "every_checkpoint_key_is_matched_or_tolerated": not encoder_side_uncovered,
        "every_extra_supplied_key_is_a_tolerated_key": (
            encoder_side_unexpected_supplied == tolerated_keys
        ),
    }
    return {
        "configuration": label,
        "use_logit_scale": use_logit_scale,
        "encoder_state_dict_key_count": len(own),
        "encoder_parameter_count": sum(1 for _ in encoder.parameters()),
        "checkpoint_keys_considered": len(remapped),
        "matched_count": len(matched),
        "shape_mismatch_count": len(shape_mismatch),
        "shape_mismatch": shape_mismatch,
        "absent_count": len(absent),
        "absent": absent,
        "encoder_keys_not_supplied": not_supplied,
        "encoder_keys_not_covered_by_checkpoint": encoder_side_uncovered,
        "extra_supplied_keys_without_encoder_destination": encoder_side_unexpected_supplied,
        "transformer_checkpoint_keys": len(
            [k for k in remapped if k == "transformer" or k.startswith("transformer.")]
        ),
        "transformer_missing_keys": mtransformer,
        "transformer_shape_mismatches": mmtransformer,
        "load_state_dict_return_value": {
            "strict": False,
            "missing_keys": missing_keys,
            "unexpected_keys": unexpected_keys,
        },
        "tolerated_keys": {
            key: {"reason": tolerated[key], **tensor_facts(remapped[key])}
            for key in tolerated_keys
        },
        "unclassified_absent_keys": unclassified,
        "checks": checks,
        "all_checks_passed": bool(all(checks.values())),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", type=Path, default=REPO_DEFAULT)
    parser.add_argument("--checkpoint", type=Path, default=CKPT_DEFAULT)
    parser.add_argument("--outer-audit", type=Path, default=OUTER_AUDIT_DEFAULT)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--top-level-keys", type=int, default=None,
                        help="Expected len(obj['state_dict']); omit to skip the check.")
    args = parser.parse_args()

    repo = args.repo.resolve()
    if not args.checkpoint.is_file():
        raise FileNotFoundError(args.checkpoint)

    # --- stub packages so src/{models,data}/__init__.py (Lightning) is skipped ---
    for name, rel in [
        ("src", "src"),
        ("src.models", "src/models"),
        ("src.models.components", "src/models/components"),
        ("src.data", "src/data"),
        ("src.data.datasets", "src/data/datasets"),
    ]:
        mod = types.ModuleType(name)
        mod.__path__ = [str(repo / rel)]
        sys.modules[name] = mod

    sys.path.insert(0, str(repo))
    sys.path.insert(0, str(repo / "external" / "mdgen"))

    import torch  # noqa: E402
    from src.models.components.md_encoder import TrajectoryEncoder  # noqa: E402

    report: dict = {
        "method": "assertive runtime audit of the OneProt MD branch load",
        "checkpoint": str(args.checkpoint),
        "file_size_bytes": args.checkpoint.stat().st_size,
        "torch_version": torch.__version__,
        "protocol_requirement": (
            "pretrained=false / model_path=null instantiation, load the outer "
            "checkpoint, require zero missing keys and zero shape mismatches "
            "under network.md.transformer.*"
        ),
    }

    if args.outer_audit.is_file():
        outer = json.loads(args.outer_audit.read_text(encoding="utf-8"))
        report["outer_audit"] = {
            "path": str(args.outer_audit),
            "sha256": outer.get("sha256"),
            "embedded_md_transformer_candidate": outer.get("embedded_md_transformer_candidate"),
            "decision": outer.get("decision"),
            "md_transformer_key_count": outer.get("md_transformer_key_count"),
            "md_transformer_tensor_bytes": outer.get("md_transformer_tensor_bytes"),
        }

    payload = torch.load(args.checkpoint, map_location="cpu", weights_only=False)
    top = payload["state_dict"] if isinstance(payload, dict) and "state_dict" in payload else payload
    if not isinstance(top, dict):
        raise TypeError(f"Unsupported checkpoint payload: {type(top).__name__}")
    md = {str(k): v for k, v in top.items() if str(k).startswith("network.md.")}
    remapped = {k[len("network.md."):]: v for k, v in md.items()}
    report["checkpoint_top_level_key_count"] = len(top)
    report["network_md_key_count"] = len(md)
    report["network_md_prefix_counts"] = {
        prefix: sum(1 for k in md if k.startswith(prefix))
        for prefix in ("network.md.transformer", "network.md.proj", "network.md.norm")
    }

    def build(**overrides):
        kwargs = dict(
            output_dim=1024, hidden_size=21, num_layers=4, num_heads=8,
            pretrained=False, frozen=True, proj_type="mlp", num_frames=100,
            suffix="_i1",
        )
        kwargs.update(overrides)
        return TrajectoryEncoder(**kwargs)

    # ---- configuration A: exactly what the formal G0 smoke used -------------
    report["configurations"] = [
        audit_configuration(
            build(),
            remapped,
            label="A_formal_g0_encoder_use_logit_scale_false",
            use_logit_scale=False,
            tolerated=KNOWN_UNPROVIDED,
        ),
        audit_configuration(
            build(use_logit_scale=True),
            remapped,
            label="B_counterfactual_use_logit_scale_true",
            use_logit_scale=True,
            tolerated={},
        ),
    ]
    # `model_path` is never passed: the protocol demands pretrained=false/model_path=null.
    report["encoder_constructed_without_forward_sim_ckpt"] = True
    report["model_path_argument_passed"] = None

    a, b = report["configurations"]
    report["configuration_equivalence"] = {
        "A_absent_keys": sorted(row["key"] for row in a["absent"]),
        "B_absent_keys": sorted(row["key"] for row in b["absent"]),
        "difference_is_exactly_the_logit_scale_buffer": (
            sorted(row["key"] for row in a["absent"]) == sorted(KNOWN_UNPROVIDED)
            and not b["absent"]
            and not b["shape_mismatch"]
            and b["matched_count"] == b["checkpoint_keys_considered"]
        ),
    }

    checks = {
        "A_all_checks_passed": a["all_checks_passed"],
        "B_exact_match_zero_difference": b["all_checks_passed"],
        "difference_is_exactly_the_logit_scale_buffer": (
            report["configuration_equivalence"]["difference_is_exactly_the_logit_scale_buffer"]
        ),
        "both_configurations_supply_every_transformer_key": (
            not a["transformer_missing_keys"] and not b["transformer_missing_keys"]
        ),
        "no_shape_mismatch_in_either_configuration": (
            not a["shape_mismatch"] and not b["shape_mismatch"]
        ),
    }
    if args.top_level_keys is not None:
        checks["checkpoint_top_level_key_count_as_expected"] = (
            len(top) == args.top_level_keys
        )
    report["checks"] = checks
    report["all_checks_passed"] = bool(all(checks.values()))
    report["claim_boundary"] = (
        "G0 technical load audit only: proves the outer checkpoint carries a complete, "
        "shape-compatible MD transformer for the pretrained=false encoder. It is not "
        "PAM efficacy evidence and says nothing about embedding quality. The 4-byte "
        "norm.1.log_logit_scale buffer is contrastive-only and is not read by forward()."
    )

    text = json.dumps(report, ensure_ascii=False, indent=2)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(text, encoding="utf-8")
    print(text)

    if not report["all_checks_passed"]:
        print("\nRUNTIME_AUDIT_FAILED")
        raise SystemExit(2)
    print("\nRUNTIME_AUDIT_PASSED")


if __name__ == "__main__":
    main()
