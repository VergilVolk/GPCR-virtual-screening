#!/usr/bin/env python3
"""READ-ONLY Phase-1 audit: what four-context material actually exists."""
import json
import struct
from pathlib import Path

W = Path("/root/work/GPCR-virtual-screening")
PROD = W / "project/results/pacer_dc_production_v01"
REF = W / "project/results/pacer_dc_membrane_reference_v01"

SYS = ["apo", "probe_only", "LY2119620__candidate_probe", "LY2119620__candidate_no_probe"]


def dcd_frames(p):
    if not p.is_file():
        return None
    with p.open("rb") as fh:
        head = fh.read(16)
    if head[4:8] == b"CORD":
        return struct.unpack("<i", head[8:12])[0]
    if head[:4] == b"CORD":
        return struct.unpack("<i", head[4:8])[0]
    return None


def state_rows(p):
    if not p.is_file():
        return None, None, None
    rows = []
    for line in p.open():
        parts = line.split(",")
        if parts[0].strip().isdigit():
            rows.append((int(float(parts[0])), float(parts[1])))
    if not rows:
        return 0, None, None
    dt = (rows[-1][1] - rows[0][1]) / (len(rows) - 1) if len(rows) > 1 else None
    return len(rows), rows[0][1], dt


print("=" * 100)
print("Phase 1 audit — PACER-DC production material (read-only)")
print("=" * 100)
print(f"{'system':<34}{'rep':<5}{'status':<10}{'ns':<7}{'dcd_frames':<12}{'csv_rows':<10}{'ps/frame':<10}{'topology'}")
print("-" * 100)
for s in SYS:
    for r in (1, 2, 3):
        d = PROD / s / f"replica_{r:02d}"
        ref = REF / s
        if not d.is_dir():
            print(f"{s:<34}{r:<5}{'MISSING':<10}{'-':<7}{'-':<12}{'-':<10}{'-':<10}"
                  f"{'ok' if (ref/'minimized.pdb').is_file() else 'NO'}")
            continue
        led = d / "progress.json"
        L = json.loads(led.read_text()) if led.is_file() else {}
        fr = dcd_frames(d / "trajectory.dcd")
        rows, t0, dt = state_rows(d / "state.csv")
        print(f"{s:<34}{r:<5}{str(L.get('status')):<10}{str(L.get('completed_ns')):<7}{str(fr):<12}"
              f"{str(rows):<10}{str(dt):<10}{'ok' if (ref/'minimized.pdb').is_file() else 'NO'}")

print()
print("=== frame spacing sanity (apo r1) ===")
rows, t0, dt = state_rows(PROD / "apo/replica_01/state.csv")
print(f"  rows={rows} first_time_ps={t0} spacing_ps={dt} -> span={(rows-1)*dt if rows else None} ps")

print()
print("=== reference topology present for each system ===")
for s in SYS:
    ref = REF / s
    files = sorted(p.name for p in ref.iterdir()) if ref.is_dir() else []
    print(f"  {s:<34} {files}")
