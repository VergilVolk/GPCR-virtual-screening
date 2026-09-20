"""Recover the compound-table images and OCR text from WO2025099660A1.

This script does not create a benchmark or run a model.  It only builds an
auditable source ledger from the Google Patents HTML.  The functional labels
and structures are frozen in a separate step after manual/chemical QC.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import subprocess
import urllib.request
from pathlib import Path

import numpy as np
from PIL import Image


ROOT = Path(__file__).resolve().parents[2]
PAPER_DIR = ROOT / "project" / "data" / "papers" / "WO2025099660A1"
HTML_PATH = PAPER_DIR / "WO2025099660A1.html"
OUT_DIR = PAPER_DIR / "compound_tables"
TESSERACT = Path(r"C:\Program Files\Tesseract-OCR\tesseract.exe")


IMG_RE = re.compile(
    r'<img id="(?P<id>imgf(?P<page>\d{6})_(?P<seq>\d{4}))"[^>]*'
    r'img-content="(?P<content>[^"]+)"[^>]*src="(?P<url>[^"]+)"',
    re.IGNORECASE,
)


def extract_manifest(html: str) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for match in IMG_RE.finditer(html):
        page = int(match.group("page"))
        # Compound examples 3--126 occupy the source table images 36--84.
        if 36 <= page <= 84 and match.group("content").lower() == "table":
            rows.append(
                {
                    "image_id": match.group("id"),
                    "source_page": page,
                    "sequence": int(match.group("seq")),
                    "url": match.group("url"),
                }
            )
    return rows


def download(url: str, destination: Path) -> None:
    if destination.exists() and destination.stat().st_size > 1_000:
        return
    request = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(request, timeout=60) as response:
        destination.write_bytes(response.read())


def ocr(image_path: Path, output_path: Path) -> str:
    if output_path.exists() and output_path.stat().st_size > 20:
        return output_path.read_text(encoding="utf-8", errors="replace")
    result = subprocess.run(
        [str(TESSERACT), str(image_path), "stdout", "--psm", "6"],
        check=True,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
    )
    output_path.write_text(result.stdout, encoding="utf-8")
    return result.stdout


def _line_centers(mask_fraction: np.ndarray, cutoff: float = 0.85) -> list[int]:
    indices = np.flatnonzero(mask_fraction >= cutoff).tolist()
    groups: list[list[int]] = []
    for index in indices:
        if not groups or index > groups[-1][-1] + 1:
            groups.append([index])
        else:
            groups[-1].append(index)
    return [int(round(sum(group) / len(group))) for group in groups]


def _ocr_cell(image: Image.Image, psm: int, digits_only: bool = False) -> str:
    command = [str(TESSERACT), "stdin", "stdout", "--psm", str(psm)]
    if digits_only:
        command += ["-c", "tessedit_char_whitelist=0123456789"]
    import io

    buffer = io.BytesIO()
    image.save(buffer, format="PNG")
    result = subprocess.run(
        command,
        input=buffer.getvalue(),
        check=True,
        capture_output=True,
    )
    return result.stdout.decode("utf-8", errors="replace").strip()


def segment_compound_rows(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    """OCR ID/name cells using the patent's explicit table grid.

    The chemical name remains provisional OCR evidence.  Row segmentation is
    deterministic and preserves the source image and pixel interval for audit.
    """
    cells: list[dict[str, object]] = []
    for source in rows:
        image_path = ROOT / str(source["local_image"])
        image = Image.open(image_path).convert("L")
        array = np.asarray(image)
        horizontal = _line_centers((array < 80).mean(axis=1))
        vertical = _line_centers((array < 80).mean(axis=0))
        if len(vertical) < 4:
            continue
        id_left, id_right = vertical[0] + 4, vertical[1] - 4
        name_left, name_right = vertical[2] + 4, vertical[3] - 4
        data_left, data_right = vertical[3] + 4, vertical[4] - 4
        for top, bottom in zip(horizontal, horizontal[1:]):
            if bottom - top < 90:
                continue
            id_crop = image.crop((id_left, top + 4, id_right, bottom - 4))
            compound_text = _ocr_cell(id_crop, psm=10, digits_only=True)
            matches = re.findall(r"\d{1,3}", compound_text)
            if not matches:
                continue
            compound_id = int(matches[0])
            if not 3 <= compound_id <= 126:
                continue
            name_crop = image.crop((name_left, top + 4, name_right, bottom - 4))
            name_ocr = _ocr_cell(name_crop, psm=6)
            data_crop = image.crop((data_left, top + 4, data_right, bottom - 4))
            characterization_ocr = _ocr_cell(data_crop, psm=6)
            mass_matches = re.findall(
                r"Mass\s*\(m/z\)\s*:\s*([0-9]+(?:\.[0-9]+)?)",
                characterization_ocr,
                flags=re.IGNORECASE,
            )
            cells.append(
                {
                    "compound_id": compound_id,
                    "image_id": source["image_id"],
                    "source_page": source["source_page"],
                    "row_top_px": top,
                    "row_bottom_px": bottom,
                    "name_ocr": " ".join(name_ocr.split()),
                    "characterization_ocr": " ".join(characterization_ocr.split()),
                    "reported_mh_plus": float(mass_matches[-1]) if mass_matches else "",
                    "structure_status": "not_reconstructed",
                    "qc_status": "raw_ocr",
                }
            )
    return cells


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--skip-ocr", action="store_true")
    args = parser.parse_args()

    if not HTML_PATH.exists():
        raise FileNotFoundError(HTML_PATH)
    if not args.skip_ocr and not TESSERACT.exists():
        raise FileNotFoundError(TESSERACT)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rows = extract_manifest(HTML_PATH.read_text(encoding="utf-8", errors="replace"))
    for row in rows:
        image_path = OUT_DIR / f"{row['image_id']}.png"
        text_path = OUT_DIR / f"{row['image_id']}.txt"
        download(str(row["url"]), image_path)
        row["local_image"] = image_path.relative_to(ROOT).as_posix()
        if not args.skip_ocr:
            text = ocr(image_path, text_path)
            row["local_ocr"] = text_path.relative_to(ROOT).as_posix()
            row["ocr_characters"] = len(text)

    with (OUT_DIR / "image_manifest.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    cells = [] if args.skip_ocr else segment_compound_rows(rows)
    if cells:
        with (OUT_DIR / "compound_row_ocr.csv").open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(cells[0]))
            writer.writeheader()
            writer.writerows(cells)
    (OUT_DIR / "audit.json").write_text(
        json.dumps(
            {
                "patent": "WO2025099660A1",
                "status": "raw_source_recovery_only",
                "n_table_images": len(rows),
                "n_segmented_compound_rows": len(cells),
                "n_unique_segmented_ids": len({row["compound_id"] for row in cells}),
                "image_page_range": [36, 84],
                "benchmark_frozen": False,
                "model_predictions_run": False,
                "warning": "OCR is evidence assistance, not ground truth; every structure requires chemical QC.",
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    print(
        json.dumps(
            {
                "images": len(rows),
                "segmented_rows": len(cells),
                "unique_ids": len({row["compound_id"] for row in cells}),
                "output": str(OUT_DIR),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
