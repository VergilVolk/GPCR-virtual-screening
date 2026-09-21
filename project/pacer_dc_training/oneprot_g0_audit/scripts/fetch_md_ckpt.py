#!/usr/bin/env python
"""Resumably extract a single member from a remote ZIP via HTTP range requests.

Why: checkpoints.zip on Zenodo is 14.73 GB and holds seven ~3.5 GB models, but
only the MD-extended one is needed.  Zenodo honours Range requests even though
it does not advertise Accept-Ranges, so the deflate stream for one member can be
fetched on its own.  The member payload starts at DATA_OFFSET and is
COMPRESS_SIZE bytes long; the enclosing ZIP local header is skipped.

Stages are idempotent: the compressed payload is kept as *.deflate.part and the
download resumes from whatever is already on disk, so an interrupted run can be
restarted without refetching.
"""
from __future__ import annotations

import os
import time
import urllib.request
import zlib

URL = "https://zenodo.org/api/records/20997998/files/checkpoints.zip/content"

# From the ZIP central directory (read separately, no full download needed).
DATA_OFFSET = 11_350_667_679
COMPRESS_SIZE = 2_469_365_696
FILE_SIZE = 3_614_862_456
CRC32 = 0x2F07819F

OUT = (
    "project/tools/oneprot-embeddings/artifacts/"
    "Pocket_Text_ST_SG_MD/epoch_012_01100-v1.ckpt"
)
PART = OUT + ".deflate.part"
CHUNK = 8 * 1024 * 1024


def fetch_range(start: int, end: int, retries: int = 6) -> bytes:
    req = urllib.request.Request(URL, headers={"Range": f"bytes={start}-{end}"})
    last = None
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(req, timeout=300) as resp:
                return resp.read()
        except Exception as exc:  # network flake: back off and retry
            last = exc
            print(f"    retry {attempt + 1}/{retries}: {exc}", flush=True)
            time.sleep(5 * (attempt + 1))
    raise RuntimeError(f"range fetch failed: {last}")


def main() -> None:
    os.makedirs(os.path.dirname(OUT), exist_ok=True)

    print("=== stage 1: fetch compressed member payload ===", flush=True)
    have = os.path.getsize(PART) if os.path.exists(PART) else 0
    if have > COMPRESS_SIZE:
        print(f"  partial file larger than expected ({have}); removing", flush=True)
        os.remove(PART)
        have = 0
    print(f"  resuming at {have:,} / {COMPRESS_SIZE:,}", flush=True)

    started = time.time()
    baseline = have
    with open(PART, "ab") as fh:
        while have < COMPRESS_SIZE:
            n = min(CHUNK, COMPRESS_SIZE - have)
            data = fetch_range(DATA_OFFSET + have, DATA_OFFSET + have + n - 1)
            if not data:
                raise RuntimeError("empty range response")
            fh.write(data)
            fh.flush()
            have += len(data)
            elapsed = time.time() - started
            got = have - baseline
            rate = got / elapsed / 1024 if elapsed > 0 else 0
            eta = (COMPRESS_SIZE - have) / (got / elapsed) if got > 0 else 0
            print(
                f"  {have:>13,}/{COMPRESS_SIZE:,}  {100*have/COMPRESS_SIZE:6.2f}%"
                f"  {rate:7.1f} KB/s  ETA {eta/60:6.1f} min",
                flush=True,
            )

    print("=== stage 2: inflate and verify ===", flush=True)
    decomp = zlib.decompressobj(-15)
    crc = 0
    total = 0
    with open(PART, "rb") as src, open(OUT, "wb") as dst:
        while True:
            block = src.read(CHUNK)
            if not block:
                break
            out = decomp.decompress(block)
            dst.write(out)
            crc = zlib.crc32(out, crc)
            total += len(out)
        tail = decomp.flush()
        dst.write(tail)
        crc = zlib.crc32(tail, crc)
        total += len(tail)

    print(f"  uncompressed : {total:,}  (expected {FILE_SIZE:,})  match={total == FILE_SIZE}", flush=True)
    print(f"  crc32        : {crc & 0xFFFFFFFF:08x}  (expected {CRC32:08x})"
          f"  match={(crc & 0xFFFFFFFF) == CRC32}", flush=True)
    if total != FILE_SIZE or (crc & 0xFFFFFFFF) != CRC32:
        raise SystemExit("verification FAILED; keeping .part for inspection")
    print(f"  OK -> {OUT}", flush=True)
    print("MD_CHECKPOINT_READY", flush=True)


if __name__ == "__main__":
    main()
