# -*- coding: utf-8 -*-
"""Resumable download and size verification for Figshare article 33283491."""
from __future__ import annotations

import json
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import requests


ARTICLE = 33283491
P = Path(__file__).resolve().parents[1]
OUT = P / "data" / "m4_gamd_figshare_33283491"
OUT.mkdir(parents=True, exist_ok=True)


def download(meta):
    path = OUT / meta["name"]
    expected = int(meta["size"])
    if path.exists() and path.stat().st_size == expected:
        return meta["name"], "cached", expected
    tmp = path.with_suffix(path.suffix + ".part")
    have = tmp.stat().st_size if tmp.exists() else 0
    headers = {"Range": f"bytes={have}-"} if have else {}
    with requests.get(meta["download_url"], headers=headers, stream=True, timeout=120) as r:
        if have and r.status_code != 206:
            have = 0
            tmp.unlink(missing_ok=True)
        r.raise_for_status()
        with tmp.open("ab" if have else "wb") as fh:
            for chunk in r.iter_content(1024 * 1024):
                if chunk:
                    fh.write(chunk)
    if tmp.stat().st_size != expected:
        raise IOError(f"size mismatch {meta['name']}: {tmp.stat().st_size} != {expected}")
    tmp.replace(path)
    return meta["name"], "downloaded", expected


def main():
    article = requests.get(f"https://api.figshare.com/v2/articles/{ARTICLE}", timeout=60).json()
    files = article["files"]
    (OUT / "figshare_manifest.json").write_text(json.dumps({
        "article_id": ARTICLE, "title": article["title"], "files": files,
        "total_bytes": sum(int(x["size"]) for x in files),
    }, indent=2), encoding="utf-8")
    rows = []
    with ThreadPoolExecutor(max_workers=4) as ex:
        jobs = [ex.submit(download, x) for x in files]
        for job in as_completed(jobs):
            row = job.result()
            rows.append(row)
            print(*row, flush=True)
    print("verified_files", len(rows), "verified_bytes", sum(x[2] for x in rows), flush=True)


if __name__ == "__main__":
    main()
