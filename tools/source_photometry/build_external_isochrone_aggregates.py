#!/usr/bin/env python3
"""Build normalized MIST aggregate tables for source-forward tests."""

from __future__ import annotations

import argparse
import hashlib
import io
import math
import tarfile
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
SOURCE_PHOTOMETRY = REPO_ROOT / "input_files" / "source_photometry"

ROMAN_COLUMNS = [
    "component", "component_index", "family", "logAge", "MH", "Zini",
    "Mini", "Mass", "radius_rsun", "Teff_K", "logg",
    "Jmag_2mass", "Hmag_2mass", "Ksmag_2mass",
    "F087mag", "F146mag", "F213mag",
]
PRIME_COLUMNS = [
    "component", "component_index", "family", "logAge", "MH", "Zini",
    "Mini", "Mass", "radius_rsun", "Teff_K", "logg",
    "Vmag", "Imag", "Jmag_2mass", "Hmag_2mass", "Ksmag_2mass",
]

COMPONENTS = [
    {"name": "thin1", "component_index": 0, "family": "thin_disk", "ages_gyr": [0.05, 0.10, 0.15]},
    {"name": "thin2", "component_index": 1, "family": "thin_disk", "ages_gyr": [0.15, 0.586449, 1.0]},
    {"name": "thin3", "component_index": 2, "family": "thin_disk", "ages_gyr": [1.0, 1.516357, 2.0]},
    {"name": "thin4", "component_index": 3, "family": "thin_disk", "ages_gyr": [2.0, 2.516884, 3.0]},
    {"name": "thin5", "component_index": 4, "family": "thin_disk", "ages_gyr": [3.0, 4.068387, 5.0]},
    {"name": "thin6", "component_index": 5, "family": "thin_disk", "ages_gyr": [5.0, 6.069263, 7.0]},
    {"name": "thin7", "component_index": 6, "family": "thin_disk", "ages_gyr": [7.0, 8.656024, 10.0]},
    {"name": "thick", "component_index": 7, "family": "thick_disk", "ages_gyr": [12.0]},
    {"name": "bar", "component_index": 8, "family": "bar", "ages_gyr": [8.0, 9.0, 10.0]},
    {"name": "NSD", "component_index": 9, "family": "NSD", "ages_gyr": [6.0, 7.0, 8.0]},
    {"name": "halo", "component_index": 10, "family": "stellar_halo", "ages_gyr": [12.0]},
]

MIST_FEH = [-1.25, -1.00, -0.75, -0.50, -0.25, 0.00, 0.25]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def write_table(path: Path, columns: list[str], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        handle.write("# " + " ".join(columns) + "\n")
        for row in rows:
            handle.write(" ".join(format_value(row[col]) for col in columns) + "\n")


def format_value(value) -> str:
    if isinstance(value, str):
        return value
    if isinstance(value, int):
        return str(value)
    return f"{value:.10g}"


def merge_by_initial_mass(*tables: list[dict[str, float]], key: str) -> list[list[dict[str, float]]]:
    indexes = [0 for _ in tables]
    merged = []
    for row in tables[0]:
        mass = row[key]
        group = [row]
        ok = True
        for table_index, table in enumerate(tables[1:], start=1):
            idx = indexes[table_index]
            while idx < len(table) and table[idx][key] < mass - 1e-7:
                idx += 1
            indexes[table_index] = idx
            if idx >= len(table) or abs(table[idx][key] - mass) > 1e-7:
                ok = False
                break
            group.append(table[idx])
        if ok:
            merged.append(group)
    return merged


def mist_feh_token(feh: float) -> str:
    prefix = "p" if feh >= 0.0 else "m"
    return f"feh_{prefix}{int(round(abs(feh) * 100)):03d}"


def mist_afe_token(alpha_fe: float) -> str:
    prefix = "p" if alpha_fe >= 0.0 else "m"
    return f"afe_{prefix}{int(round(abs(alpha_fe) * 10))}"


def build_mist(args) -> None:
    roman_archive = SOURCE_PHOTOMETRY / "mist" / "v2.5" / "raw" / "isos" / "Roman.txz"
    ubvri_archive = SOURCE_PHOTOMETRY / "mist" / "v2.5" / "raw" / "isos" / "UBVRIplus.txz"
    for alpha_fe, label, subdir in [
        (0.0, "solar", "solar_scaled"),
        (0.4, "alpha_afe_p0p40", "alpha_enhanced/afe_p0p40"),
    ]:
        roman_tables = load_mist_tables(roman_archive, "Roman", alpha_fe)
        ubvri_tables = load_mist_tables(ubvri_archive, "UBVRIplus", alpha_fe)
        roman_rows = []
        prime_rows = []
        for component in COMPONENTS:
            for target_age_gyr in component["ages_gyr"]:
                target_log_age = math.log10(target_age_gyr * 1.0e9)
                for feh in MIST_FEH:
                    roman = roman_tables.get(feh)
                    ubvri = ubvri_tables.get(feh)
                    if roman is None or ubvri is None:
                        continue
                    roman_age = nearest_age_rows(roman, target_log_age)
                    ubvri_age = nearest_age_rows(ubvri, target_log_age)
                    for rrow, urow in merge_by_initial_mass(roman_age, ubvri_age, key="initial_mass"):
                        base = {
                            "component": component["name"],
                            "component_index": component["component_index"],
                            "family": component["family"],
                            "logAge": rrow["log10_isochrone_age_yr"],
                            "MH": rrow.get("[Fe/H]", rrow.get("[Fe/H]_init", feh)),
                            "Zini": rrow["Zinit"],
                            "Mini": rrow["initial_mass"],
                            "Mass": rrow["star_mass"],
                            "radius_rsun": 10.0 ** rrow["log_R"],
                            "Teff_K": 10.0 ** rrow["log_Teff"],
                            "logg": rrow["log_g"],
                            "Jmag_2mass": urow["2MASS_J"],
                            "Hmag_2mass": urow["2MASS_H"],
                            "Ksmag_2mass": urow["2MASS_Ks"],
                            "F087mag": rrow["Roman_F087"],
                            "F146mag": rrow["Roman_F146"],
                            "F213mag": rrow["Roman_F213"],
                        }
                        roman_rows.append(base)
                        prow = dict(base)
                        prow["Vmag"] = urow["Bessell_V"]
                        prow["Imag"] = urow["Bessell_I"]
                        prime_rows.append(prow)

        outdir = SOURCE_PHOTOMETRY / "mist" / "v2.5" / subdir / "normalized"
        write_table(outdir / f"all_roman_mist_{label}.dat", ROMAN_COLUMNS, roman_rows)
        write_table(outdir / f"all_prime_mist_{label}.dat", PRIME_COLUMNS, prime_rows)


def load_mist_tables(archive: Path, suffix: str, alpha_fe: float) -> dict[float, list[dict[str, float]]]:
    targets = {
        f"{mist_feh_token(feh)}_{mist_afe_token(alpha_fe)}_vvcrit0.0_full.iso.{suffix}": feh
        for feh in MIST_FEH
    }
    tables: dict[float, list[dict[str, float]]] = {}
    with tarfile.open(archive, "r:*") as tar:
        for member in tar:
            if not member.isfile() or member.name not in targets:
                continue
            fileobj = tar.extractfile(member)
            if fileobj is None:
                continue
            text = io.TextIOWrapper(fileobj, encoding="utf-8", errors="replace").read()
            _, rows = parse_mist_table(text)
            tables[targets[member.name]] = rows
            if len(tables) == len(targets):
                break
    return tables


def parse_mist_table(text: str) -> tuple[dict[str, float], list[dict[str, float]]]:
    metadata: dict[str, float] = {}
    columns: list[str] | None = None
    rows: list[dict[str, float]] = []
    lines = text.splitlines()
    for i, line in enumerate(lines):
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith("#"):
            tokens = stripped[1:].split()
            if tokens[:5] == ["Yinit", "Zinit", "[Fe/H]", "[a/Fe]", "v/vcrit"] and i + 1 < len(lines):
                values = lines[i + 1].strip().lstrip("#").split()
                metadata["Yinit"] = float(values[0])
                metadata["Zinit"] = float(values[1])
            if tokens and tokens[0] == "EEP":
                columns = tokens
            continue
        if columns is None:
            continue
        values = stripped.split()
        if len(values) != len(columns):
            continue
        row = {key: float(value) for key, value in zip(columns, values)}
        row["Zinit"] = metadata.get("Zinit", math.nan)
        rows.append(row)
    return metadata, rows


def nearest_age_rows(rows: list[dict[str, float]], target_log_age: float) -> list[dict[str, float]]:
    ages = sorted({row["log10_isochrone_age_yr"] for row in rows})
    nearest = min(ages, key=lambda age: abs(age - target_log_age))
    return [row for row in rows if row["log10_isochrone_age_yr"] == nearest]


def write_checksums() -> None:
    files = sorted(
        path for root in [
            SOURCE_PHOTOMETRY / "mist" / "v2.5" / "solar_scaled",
            SOURCE_PHOTOMETRY / "mist" / "v2.5" / "alpha_enhanced",
        ]
        if root.exists()
        for path in root.rglob("*.dat")
    )
    out = SOURCE_PHOTOMETRY / "external_isochrone_checksums.sha256"
    with out.open("w") as handle:
        for path in files:
            handle.write(f"{sha256(path)}  {path.relative_to(SOURCE_PHOTOMETRY)}\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.parse_args()
    build_mist(None)
    write_checksums()


if __name__ == "__main__":
    main()
