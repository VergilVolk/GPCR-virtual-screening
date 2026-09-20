"""Export official DeepRLI scores and the pretrained 64D complex representation z."""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch

from deeprli.datasets import HeavyDatasetForInfer
from deeprli.model import DeepRLI


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-root", required=True)
    ap.add_argument("--index", default="index/complexes.processed.csv")
    ap.add_argument("--compiled", default="compiled/complexes.processed.pkl")
    ap.add_argument("--weights", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--batch", type=int, default=8)
    ap.add_argument("--chunk-size", type=int, default=64)
    args = ap.parse_args()
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    dataset = HeavyDatasetForInfer(root=args.data_root, data_index=args.index, data_file=args.compiled)
    model = DeepRLI(f_dropout_rate=0, g_dropout_rate=0, hidden_dim=64,
                    num_attention_heads=8, use_layer_norm=False,
                    use_batch_norm=True, use_residual=True,
                    use_envelope=True, use_multi_obj=True)
    state = torch.load(args.weights, map_location="cpu")
    model.load_state_dict(state, strict=True); model.eval()
    chunk_dir = out / "chunks"; chunk_dir.mkdir(exist_ok=True)
    t0 = time.time()
    for start in range(0, len(dataset), args.chunk_size):
        end = min(start + args.chunk_size, len(dataset))
        chunk_path = chunk_dir / f"chunk_{start:04d}_{end:04d}.npz"
        if chunk_path.exists():
            saved = np.load(chunk_path)
            if saved["z"].shape == (end - start, 64):
                print(f"resume {end}/{len(dataset)}", flush=True)
                continue
        subset = torch.utils.data.Subset(dataset, range(start, end))
        loader = torch.utils.data.DataLoader(subset, batch_size=args.batch,
                                             collate_fn=HeavyDatasetForInfer.collate_fn, num_workers=0)
        captured, local_scores = [], [[], [], []]
        handle = model.readout1.register_forward_pre_hook(
            lambda module, inputs: captured.append(inputs[0].detach().cpu().numpy().copy()))
        with torch.no_grad():
            for batch in loader:
                result = model(batch)
                for i in range(3):
                    local_scores[i].extend(result[i].cpu().numpy().tolist())
        handle.remove()
        local_z = np.vstack(captured)
        local_scores = np.asarray(local_scores, dtype=np.float32).T
        if local_z.shape != (end - start, 64) or local_scores.shape != (end - start, 3):
            raise RuntimeError(f"bad chunk {start}:{end} z={local_z.shape} scores={local_scores.shape}")
        np.savez_compressed(chunk_path, z=local_z, scores=local_scores)
        print(f"completed {end}/{len(dataset)} elapsed_s={time.time()-t0:.1f}", flush=True)
    all_z, all_scores = [], []
    for start in range(0, len(dataset), args.chunk_size):
        end = min(start + args.chunk_size, len(dataset))
        saved = np.load(chunk_dir / f"chunk_{start:04d}_{end:04d}.npz")
        all_z.append(saved["z"]); all_scores.append(saved["scores"])
    z = np.vstack(all_z); score_matrix = np.vstack(all_scores)
    if z.shape != (len(dataset), 64) or not np.isfinite(z).all():
        raise RuntimeError(f"invalid embedding matrix: {z.shape}")
    index = pd.read_csv(Path(args.data_root) / args.index)
    scores = index.copy()
    scores["scoring_score"] = score_matrix[:, 0]
    scores["docking_score"] = score_matrix[:, 1]
    scores["screening_score"] = score_matrix[:, 2]
    scores.to_csv(out / "deeprli_scores.csv", index=False)
    np.save(out / "deeprli_z64.npy", z)
    columns = [f"z{i:02d}" for i in range(64)]
    pd.concat([scores[["complex_path"]], pd.DataFrame(z, columns=columns)], axis=1).to_csv(
        out / "deeprli_z64.csv", index=False)
    audit = {"complexes": len(dataset), "embedding_shape": list(z.shape),
             "finite_scores": bool(np.isfinite(score_matrix).all()),
             "weights_loaded_strict": True,
             "resumable_chunks": len(all_z),
             "embedding_source": "input z to official DeepRLI readout1 after 10 graph-transformer layers"}
    (out / "audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
