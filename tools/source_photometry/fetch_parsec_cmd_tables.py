#!/usr/bin/env python3
"""Fetch PARSEC/CMD source-photometry tables and build normalized tables."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path


TSUN_K = 5772.0


def read_config(path: Path) -> dict:
    return json.loads(path.read_text())


def radius_from_logl_logte(log_l: float, log_te: float) -> float:
    teff = 10.0**log_te
    return math.sqrt(10.0**log_l) * (TSUN_K / teff) ** 2


def write_table(path: Path, columns: list[str], rows: list[dict[str, float]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as out:
        out.write("# " + " ".join(columns) + "\n")
        for row in rows:
            out.write(" ".join(format_value(row[col]) for col in columns) + "\n")


def format_value(value: float) -> str:
    if isinstance(value, int):
        return str(value)
    return f"{value:.10g}"


def save_raw_dataframe(df, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as out:
        out.write("# downloaded from PARSEC/CMD via ezpadova\n")
        out.write("# " + " ".join(str(col) for col in df.columns) + "\n")
        df.to_csv(out, sep=" ", index=False, header=False, float_format="%.10g")


def table_rows(df) -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    for rec in df.to_dict(orient="records"):
        rows.append({str(k): float(v) for k, v in rec.items()})
    return rows


def read_raw_table(path: Path) -> list[dict[str, float]]:
    columns: list[str] | None = None
    rows: list[dict[str, float]] = []
    for line in path.read_text().splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith("#"):
            tokens = stripped[1:].split()
            if tokens and tokens[0] == "Zini":
                columns = tokens
            continue
        if columns is None:
            raise ValueError(f"missing raw header in {path}")
        rows.append(dict(zip(columns, (float(token) for token in stripped.split()))))
    return rows


def load_existing_raw(config: dict, output_dir: Path, components_to_update: set[str]) -> dict[tuple[str, str], list[dict[str, float]]]:
    raw: dict[tuple[str, str], list[dict[str, float]]] = {}
    for component in config["components"]:
        name = component["name"]
        if name not in components_to_update:
            continue
        for phot_name in config["photometric_systems"]:
            raw[(name, phot_name)] = read_raw_table(output_dir / "raw" / f"{name}_{phot_name}.dat")
    return raw


def merge_by_mini(left: list[dict[str, float]], right: list[dict[str, float]], *, tol: float = 1e-7) -> list[tuple[dict[str, float], dict[str, float]]]:
    pairs = []
    j = 0
    for lrow in left:
        mini = lrow["Mini"]
        while j < len(right) and right[j]["Mini"] < mini - tol:
            j += 1
        if j < len(right) and abs(right[j]["Mini"] - mini) <= tol:
            pairs.append((lrow, right[j]))
    return pairs


def build_normalized(config: dict, raw: dict[tuple[str, str], list[dict[str, float]]], output_dir: Path, allowed_components: set[str] | None = None) -> None:
    for component in config["components"]:
        name = component["name"]
        if allowed_components and name not in allowed_components:
            continue
        ubv = raw[(name, "ubvrijhk")]
        twomass = raw[(name, "twomass")]
        roman = raw[(name, "roman2021")]

        prime_rows = []
        for urow, trow in merge_by_mini(ubv, twomass):
            radius = radius_from_logl_logte(trow["logL"], trow["logTe"])
            prime_rows.append(
                {
                    "Mini": trow["Mini"],
                    "Mass": trow["Mass"],
                    "radius_rsun": radius,
                    "Teff_K": 10.0 ** trow["logTe"],
                    "logg": trow["logg"],
                    "Vmag": urow["Vmag"],
                    "Imag": urow["Imag"],
                    "Jmag_2mass": trow["Jmag"],
                    "Hmag_2mass": trow["Hmag"],
                    "Ksmag_2mass": trow["Ksmag"],
                }
            )
        write_table(output_dir / "normalized" / f"{name}_prime_parsec.dat", config["normalized_tables"]["prime_parsec"]["columns"], prime_rows)

        roman_rows = []
        for trow, rrow in merge_by_mini(twomass, roman):
            radius = radius_from_logl_logte(trow["logL"], trow["logTe"])
            roman_rows.append(
                {
                    "Mini": trow["Mini"],
                    "Mass": trow["Mass"],
                    "radius_rsun": radius,
                    "Teff_K": 10.0 ** trow["logTe"],
                    "logg": trow["logg"],
                    "Jmag_2mass": trow["Jmag"],
                    "Hmag_2mass": trow["Hmag"],
                    "Ksmag_2mass": trow["Ksmag"],
                    "F087mag": rrow["F087mag"],
                    "F146mag": rrow["F146mag"],
                    "F213mag": rrow["F213mag"],
                }
            )
        write_table(output_dir / "normalized" / f"{name}_roman_parsec.dat", config["normalized_tables"]["roman_parsec"]["columns"], roman_rows)


def write_checksums(output_dir: Path) -> None:
    files = sorted(p for p in output_dir.rglob("*.dat") if p.is_file())
    with (output_dir / "checksums.sha256").open("w") as out:
        for path in files:
            digest = hashlib.sha256(path.read_bytes()).hexdigest()
            out.write(f"{digest}  {path.relative_to(output_dir)}\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=Path("input_files/source_photometry/parsec_cmd/requests.json"))
    parser.add_argument("--output-dir", type=Path, default=Path("input_files/source_photometry/parsec_cmd"))
    parser.add_argument("--component", action="append", help="Fetch only the named component. Can be repeated.")
    args = parser.parse_args()

    try:
        import ezpadova
    except ImportError as exc:
        raise SystemExit("Install ezpadova first: python -m pip install git+https://github.com/mfouesneau/ezpadova") from exc

    config = read_config(args.config)
    raw_tables: dict[tuple[str, str], list[dict[str, float]]] = {}

    allowed_components = set(args.component or [])
    components = [component for component in config["components"] if not allowed_components or component["name"] in allowed_components]

    for component in components:
        name = component["name"]
        log_age = component["logAge"]
        metallicity = component.get("MH", config["defaults"]["MH"])
        for phot_name, phot_cfg in config["photometric_systems"].items():
            print(f"fetch {name} {phot_name} logAge={log_age} MH={metallicity}")
            df = ezpadova.get_isochrones(
                logage=(log_age, log_age, 0.0),
                MH=(metallicity, metallicity, 0.0),
                photsys_file=phot_cfg["photsys_file"],
                track_parsec=config["defaults"]["track_parsec"],
                track_colibri=config["defaults"]["track_colibri"],
                photsys_version=config["defaults"]["photsys_version"],
                kind_interp=config["defaults"]["kind_interp"],
                kind_mag=config["defaults"]["kind_mag"],
                output_kind=config["defaults"]["output_kind"],
                output_evstage=config["defaults"]["output_evstage"],
            )
            save_raw_dataframe(df, args.output_dir / "raw" / f"{name}_{phot_name}.dat")
            raw_tables[(name, phot_name)] = table_rows(df)

    if allowed_components:
        existing_raw = load_existing_raw(config, args.output_dir, allowed_components)
        existing_raw.update(raw_tables)
        raw_tables = existing_raw

    build_normalized(config, raw_tables, args.output_dir, allowed_components)
    write_checksums(args.output_dir)


if __name__ == "__main__":
    main()
