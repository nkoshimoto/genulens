#!/usr/bin/env python3
"""Compare PARSEC/CMD normalized tables with copied genstars photometry tables."""

from __future__ import annotations

import argparse
from pathlib import Path


def read_table(path: Path) -> list[dict[str, float]]:
    columns: list[str] | None = None
    rows: list[dict[str, float]] = []
    for line in path.read_text().splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith("#"):
            tokens = stripped[1:].split()
            if tokens and tokens[0] in {"Mini", "MPD", "Mass"}:
                columns = tokens
            continue
        if columns is None:
            raise ValueError(f"missing header in {path}")
        values = [float(token) for token in stripped.split()]
        rows.append(dict(zip(columns, values)))
    return rows


def parsec_table(path: Path) -> list[dict[str, float]]:
    rows = []
    for line in path.read_text().splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith("#"):
            if stripped.startswith("# Mini"):
                columns = stripped[1:].split()
            continue
        values = [float(token) for token in stripped.split()]
        rows.append(dict(zip(columns, values)))
    return rows


def nearest(rows: list[dict[str, float]], mini: float, mass: float) -> dict[str, float] | None:
    best = None
    best_diff = 1e99
    for row in rows:
        diff = abs(row["Mini"] - mini) + abs(row["Mass"] - mass)
        if diff < best_diff:
            best = row
            best_diff = diff
    if best_diff > 5e-3:
        return None
    return best


def compare(label: str, legacy_path: Path, parsec_path: Path, columns: list[tuple[str, str]], min_mass: float, max_mass: float) -> list[str]:
    legacy = [
        row
        for row in read_table(legacy_path)
        if row["Mini"] >= min_mass
        and row["Mini"] <= max_mass
        and abs(row["MPD"] - row["Mini"]) <= 1e-2
        and all(row[c[0]] < 90 for c in columns)
    ]
    parsec = parsec_table(parsec_path)
    lines = [
        f"### {label}",
        "",
        f"Mass range: `{min_mass} <= Mini <= {max_mass}` and `abs(MPD - Mini) <= 0.01`.",
        "",
        "| column | n | median_abs_diff | max_abs_diff |",
        "| --- | ---: | ---: | ---: |",
    ]
    for old_col, new_col in columns:
        diffs = []
        for old in legacy:
            new = nearest(parsec, old["Mini"], old["MPD"])
            if new is not None:
                diffs.append(abs(old[old_col] - new[new_col]))
        diffs.sort()
        if diffs:
            median = diffs[len(diffs) // 2]
            max_diff = max(diffs)
            lines.append(f"| {old_col} vs {new_col} | {len(diffs)} | {median:.6g} | {max_diff:.6g} |")
    lines.append("")
    return lines


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-dir", type=Path, default=Path("input_files/source_photometry"))
    parser.add_argument("--output", type=Path, default=Path("input_files/source_photometry/parsec_cmd/genstars_compatibility.md"))
    args = parser.parse_args()

    out = [
        "# genstars Compatibility Check",
        "",
        "This compares the regenerated PARSEC/CMD normalized tables against the copied",
        "`genstars` processed source-photometry tables. Differences are expected for the",
        "low-mass empirical/Baraffe parts of the legacy tables, so the primary checks use",
        "the PARSEC-dominated ranges.",
        "",
    ]

    components = ["thin1", "thin2", "thin3", "thin4", "thin5", "thin6", "thin7", "bar", "NSD"]
    for component in components:
        out.extend(
            compare(
                f"prime {component}",
                args.source_dir / "prime_hybrid" / f"isoemp_{component}.dat",
                args.source_dir / "parsec_cmd" / "normalized" / f"{component}_prime_parsec.dat",
                [("MV_j", "Vmag"), ("MI_c", "Imag"), ("MJ_2M", "Jmag_2mass"), ("MH_2M", "Hmag_2mass"), ("MK_2M", "Ksmag_2mass")],
                0.55,
                1.0,
            )
        )
        out.extend(
            compare(
                f"roman {component}",
                args.source_dir / "roman_isochrone" / f"isochrone_{component}.dat",
                args.source_dir / "parsec_cmd" / "normalized" / f"{component}_roman_parsec.dat",
                [("MJ_2M", "Jmag_2mass"), ("MH_2M", "Hmag_2mass"), ("MK_2M", "Ksmag_2mass"), ("MZ087", "F087mag"), ("MW146", "F146mag"), ("MF213", "F213mag")],
                0.09,
                1.0,
            )
        )

    out.extend(
        compare(
            "prime thick",
            args.source_dir / "prime_hybrid" / "isoemp_thick2.dat",
            args.source_dir / "parsec_cmd" / "normalized" / "thick_prime_parsec.dat",
            [("MV_j", "Vmag"), ("MI_c", "Imag"), ("MJ_2M", "Jmag_2mass"), ("MH_2M", "Hmag_2mass"), ("MK_2M", "Ksmag_2mass")],
            0.55,
            1.0,
        )
    )
    out.extend(
        compare(
            "roman thick",
            args.source_dir / "roman_isochrone" / "isochrone_thick.dat",
            args.source_dir / "parsec_cmd" / "normalized" / "thick_roman_parsec.dat",
            [("MJ_2M", "Jmag_2mass"), ("MH_2M", "Hmag_2mass"), ("MK_2M", "Ksmag_2mass"), ("MZ087", "F087mag"), ("MW146", "F146mag"), ("MF213", "F213mag")],
            0.09,
            1.0,
        )
    )

    args.output.write_text("\n".join(out) + "\n")


if __name__ == "__main__":
    main()
