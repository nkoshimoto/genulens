#!/usr/bin/env python3
"""Normalize external isochrone tables into the genulens source-prior schema.

This importer is intentionally model-agnostic.  PARSEC/CMD and MIST
tables use different column names, so the caller can either rely on common
aliases or pass explicit column mappings.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path


TSUN_K = 5772.0
LOGG_SUN_CGS = 4.438
CORE_COLUMNS = [
    "component",
    "component_index",
    "family",
    "logAge",
    "MH",
    "Zini",
    "Mini",
    "Mass",
    "radius_rsun",
    "Teff_K",
    "logg",
]

ALIASES = {
    "Mini": ["Mini", "initial_mass", "M_ini", "M/Mo(ini)", "M_ini/Mo"],
    "Mass": ["Mass", "star_mass", "current_mass", "M_act", "M/Mo(fin)", "M/Mo", "M"],
    "logL": ["logL", "log_L", "log_L/Lo", "log(L/Lo)"],
    "logTe": ["logTe", "log_Teff", "logTeff", "log_Te"],
    "Teff_K": ["Teff_K", "Teff", "T_eff", "effective_temperature"],
    "logg": ["logg", "log_g", "logg_surface"],
    "Zini": ["Zini", "Z", "Zini"],
    "MH": ["MH", "M_H", "FeH", "[Fe/H]", "feh"],
    "logAge": ["logAge", "log_age", "log10_age", "log10_isochrone_age_yr"],
}

HEADER_MARKERS = {
    "EEP",
    "Mini",
    "M/Mo(ini)",
    "initial_mass",
    "log10_isochrone_age_yr",
    "logTe",
    "log_Teff",
    "log(L/Lo)",
}


def split_mapping(values: list[str]) -> dict[str, str]:
    mapping: dict[str, str] = {}
    for value in values:
        if "=" not in value:
            raise SystemExit(f"mapping must be OUT=IN, got {value!r}")
        out, raw = value.split("=", 1)
        mapping[out.strip()] = raw.strip()
    return mapping


def read_table(path: Path) -> tuple[list[str], list[dict[str, float]]]:
    columns: list[str] | None = None
    rows: list[dict[str, float]] = []
    for line in path.read_text().splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith("#"):
            candidate = stripped[1:].split()
            if candidate and any(token in HEADER_MARKERS for token in candidate):
                columns = candidate
            continue
        if columns is None:
            columns = stripped.split()
            continue
        values = stripped.split()
        if len(values) != len(columns):
            raise SystemExit(
                f"{path}: row has {len(values)} values but header has {len(columns)} columns"
            )
        row: dict[str, float] = {}
        for key, value in zip(columns, values):
            try:
                row[key] = float(value)
            except ValueError:
                # Evolutionary phase labels or comments are irrelevant here.
                row[key] = math.nan
        rows.append(row)
    if columns is None or not rows:
        raise SystemExit(f"{path}: no table rows found")
    return columns, rows


def pick(row: dict[str, float], logical: str, explicit: dict[str, str], default=None):
    if logical in explicit:
        name = explicit[logical]
        if name not in row:
            raise SystemExit(f"mapped column {name!r} for {logical!r} is missing")
        return row[name]
    for name in ALIASES.get(logical, []):
        if name in row:
            return row[name]
    if default is not None:
        return default
    raise SystemExit(f"could not infer required column {logical!r}; pass --column {logical}=RAW")


def radius_from_row(row: dict[str, float], explicit: dict[str, str]) -> float:
    if "radius_rsun" in explicit:
        return pick(row, "radius_rsun", explicit)
    for name in ["radius_rsun", "R", "R_Rsun", "radius", "R/RSun"]:
        if name in row:
            return row[name]
    for name in ["log_R", "logR"]:
        if name in row:
            return 10.0 ** row[name]
    log_l = pick(row, "logL", explicit)
    log_te = pick(row, "logTe", explicit)
    teff = 10.0**log_te
    return math.sqrt(10.0**log_l) * (TSUN_K / teff) ** 2


def teff_from_row(row: dict[str, float], explicit: dict[str, str]) -> float:
    try:
        return pick(row, "Teff_K", explicit)
    except SystemExit:
        return 10.0 ** pick(row, "logTe", explicit)


def logg_from_row(row: dict[str, float], explicit: dict[str, str], mass_msun: float, radius_rsun: float) -> float:
    try:
        return pick(row, "logg", explicit)
    except SystemExit:
        return LOGG_SUN_CGS + math.log10(mass_msun) - 2.0 * math.log10(radius_rsun)


def fmt(value) -> str:
    if isinstance(value, str):
        return value
    if isinstance(value, int):
        return str(value)
    return f"{value:.10g}"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--component", required=True)
    parser.add_argument("--component-index", type=int, required=True)
    parser.add_argument("--family", required=True)
    parser.add_argument("--log-age", type=float)
    parser.add_argument("--mh", type=float)
    parser.add_argument("--zini", type=float)
    parser.add_argument("--column", action="append", default=[], help="Map normalized logical column to raw column, e.g. Mini=initial_mass")
    parser.add_argument("--band", action="append", default=[], help="Map output band to raw column, e.g. F146mag=F146W")
    args = parser.parse_args()

    explicit = split_mapping(args.column)
    bands = split_mapping(args.band)
    if not bands:
        raise SystemExit("at least one --band OUT=RAW mapping is required")

    _, raw_rows = read_table(args.input)
    out_rows = []
    for row in raw_rows:
        log_age = args.log_age if args.log_age is not None else pick(row, "logAge", explicit)
        mh = args.mh if args.mh is not None else pick(row, "MH", explicit)
        zini = args.zini if args.zini is not None else pick(row, "Zini", explicit)
        mass = pick(row, "Mass", explicit)
        radius = radius_from_row(row, explicit)
        out = {
            "component": args.component,
            "component_index": args.component_index,
            "family": args.family,
            "logAge": log_age,
            "MH": mh,
            "Zini": zini,
            "Mini": pick(row, "Mini", explicit),
            "Mass": mass,
            "radius_rsun": radius,
            "Teff_K": teff_from_row(row, explicit),
            "logg": logg_from_row(row, explicit, mass, radius),
        }
        for out_band, raw_band in bands.items():
            if raw_band not in row:
                raise SystemExit(f"raw band column {raw_band!r} is missing")
            out[out_band] = row[raw_band]
        out_rows.append(out)

    columns = CORE_COLUMNS + list(bands.keys())
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w") as handle:
        handle.write("# " + " ".join(columns) + "\n")
        for row in out_rows:
            handle.write(" ".join(fmt(row[col]) for col in columns) + "\n")


if __name__ == "__main__":
    main()
