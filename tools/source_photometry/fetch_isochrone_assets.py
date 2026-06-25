#!/usr/bin/env python3
"""Fetch MIST isochrone assets used for source-property systematics tests."""

from __future__ import annotations

import argparse
import hashlib
import json
import urllib.request
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
SOURCE_PHOTOMETRY = REPO_ROOT / "input_files" / "source_photometry"
MIST_BASE = "https://mist.science"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def download(url: str, output: Path, *, force: bool = False) -> dict:
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.exists() and not force:
        return {
            "url": url,
            "path": str(output.relative_to(REPO_ROOT)),
            "bytes": output.stat().st_size,
            "sha256": sha256(output),
            "status": "reused",
        }
    print(f"download {url} -> {output}", flush=True)
    with urllib.request.urlopen(url) as response, output.open("wb") as handle:
        while True:
            chunk = response.read(1024 * 1024)
            if not chunk:
                break
            handle.write(chunk)
    return {
        "url": url,
        "path": str(output.relative_to(REPO_ROOT)),
        "bytes": output.stat().st_size,
        "sha256": sha256(output),
        "status": "downloaded",
    }


def planned_download_record(url: str, output: Path) -> dict:
    status = "would_reuse" if output.exists() else "would_download"
    return {
        "url": url,
        "path": str(output.relative_to(REPO_ROOT)),
        "bytes": output.stat().st_size if output.exists() else 0,
        "sha256": sha256(output) if output.exists() else None,
        "status": status,
    }


def fetch_mist(args) -> list[dict]:
    records = []
    for product in args.mist_product:
        url = f"{MIST_BASE}/data/tarballs_v2.5/isos/{product}.txz"
        output = SOURCE_PHOTOMETRY / "mist" / "v2.5" / "raw" / "isos" / f"{product}.txz"
        records.append(
            planned_download_record(url, output)
            if args.dry_run
            else download(url, output, force=args.force)
        )
    return records


def write_manifest(records: list[dict], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps({
        "schema": "genulens.source_photometry.download_manifest.v1",
        "records": records,
    }, indent=2) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--mist-product", action="append")
    parser.add_argument("--manifest", type=Path, default=SOURCE_PHOTOMETRY / "download_manifest.mist.json")
    args = parser.parse_args()
    args.mist_product = args.mist_product or ["Roman", "UBVRIplus"]

    records = fetch_mist(args)
    write_manifest(records, args.manifest)
    if args.dry_run:
        counts = {}
        for record in records:
            counts[record["status"]] = counts.get(record["status"], 0) + 1
        print(f"dry-run records: {counts}")
    print(f"wrote {args.manifest}")


if __name__ == "__main__":
    main()
