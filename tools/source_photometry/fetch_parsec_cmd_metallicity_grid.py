#!/usr/bin/env python3
"""Fetch a PARSEC/CMD metallicity grid for source-prior sampling."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path


TSUN_K = 5772.0


def read_config(path: Path) -> dict:
    return json.loads(path.read_text())


def mh_token(mh: float) -> str:
    prefix = "p" if mh >= 0 else "m"
    value = abs(mh)
    return prefix + f"{value:.2f}".replace(".", "p")


def age_token(log_age: float) -> str:
    return "age" + f"{log_age:.5f}".replace(".", "p")


def component_log_ages(component: dict) -> list[float]:
    if "logAge_grid" in component:
        return [float(value) for value in component["logAge_grid"]]
    return [float(component["logAge"])]


def radius_from_logl_logte(log_l: float, log_te: float) -> float:
    teff = 10.0**log_te
    return math.sqrt(10.0**log_l) * (TSUN_K / teff) ** 2


def format_value(value) -> str:
    if isinstance(value, str):
        return value
    if isinstance(value, int):
        return str(value)
    return f"{value:.10g}"


def write_table(path: Path, columns: list[str], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as out:
        out.write("# " + " ".join(columns) + "\n")
        for row in rows:
            out.write(" ".join(format_value(row[col]) for col in columns) + "\n")


def save_raw_dataframe(df, path: Path, *, component: dict, mh: float, phot_name: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as out:
        out.write("# downloaded from PARSEC/CMD via ezpadova\n")
        out.write(f"# component={component['name']} family={component['family']} component_index={component['component_index']} logAge={component['logAge']} MH={mh} photometric_system={phot_name}\n")
        out.write("# " + " ".join(str(col) for col in df.columns) + "\n")
        df.to_csv(out, sep=" ", index=False, header=False, float_format="%.10g")


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


def table_rows(df) -> list[dict[str, float]]:
    return [{str(k): float(v) for k, v in rec.items()} for rec in df.to_dict(orient="records")]


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


def metadata(component: dict, row: dict[str, float], mh: float) -> dict:
    return {
        "component": component["name"],
        "component_index": component["component_index"],
        "family": component["family"],
        "logAge": component["logAge"],
        "MH": mh,
        "Zini": row["Zini"],
    }


def grid_paths(output_dir: Path, component: dict, mh: float, log_age: float, multi_age: bool, phot_name: str | None = None) -> dict[str, Path]:
    metallicity_token = mh_token(mh)
    if multi_age:
        age_part = age_token(log_age)
        paths = {
            "normalized_prime": output_dir / "normalized" / component["name"] / age_part / f"{component['name']}_{age_part}_{metallicity_token}_prime_parsec.dat",
            "normalized_roman": output_dir / "normalized" / component["name"] / age_part / f"{component['name']}_{age_part}_{metallicity_token}_roman_parsec.dat",
        }
        if phot_name is not None:
            paths["raw"] = output_dir / "raw" / "grid" / age_part / f"{metallicity_token}_{phot_name}.dat"
        return paths
    paths = {
        "normalized_prime": output_dir / "normalized" / component["name"] / f"{component['name']}_{metallicity_token}_prime_parsec.dat",
        "normalized_roman": output_dir / "normalized" / component["name"] / f"{component['name']}_{metallicity_token}_roman_parsec.dat",
    }
    if phot_name is not None:
        paths["raw"] = output_dir / "raw" / component["name"] / f"{metallicity_token}_{phot_name}.dat"
    return paths


def build_rows(component: dict, mh: float, raw: dict[str, list[dict[str, float]]]) -> tuple[list[dict], list[dict]]:
    ubv = raw["ubvrijhk"]
    twomass = raw["twomass"]
    roman = raw["roman2021"]

    prime_rows: list[dict] = []
    for urow, trow in merge_by_mini(ubv, twomass):
        base = metadata(component, trow, mh)
        base.update(
            {
                "Mini": trow["Mini"],
                "Mass": trow["Mass"],
                "radius_rsun": radius_from_logl_logte(trow["logL"], trow["logTe"]),
                "Teff_K": 10.0 ** trow["logTe"],
                "logg": trow["logg"],
                "Vmag": urow["Vmag"],
                "Imag": urow["Imag"],
                "Jmag_2mass": trow["Jmag"],
                "Hmag_2mass": trow["Hmag"],
                "Ksmag_2mass": trow["Ksmag"],
            }
        )
        prime_rows.append(base)

    roman_rows: list[dict] = []
    for trow, rrow in merge_by_mini(twomass, roman):
        base = metadata(component, trow, mh)
        base.update(
            {
                "Mini": trow["Mini"],
                "Mass": trow["Mass"],
                "radius_rsun": radius_from_logl_logte(trow["logL"], trow["logTe"]),
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
        roman_rows.append(base)

    return prime_rows, roman_rows


def write_checksums(output_dir: Path) -> None:
    files = sorted(p for p in output_dir.rglob("*.dat") if p.is_file())
    with (output_dir / "checksums.sha256").open("w") as out:
        for path in files:
            digest = hashlib.sha256(path.read_bytes()).hexdigest()
            out.write(f"{digest}  {path.relative_to(output_dir)}\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=Path("input_files/source_photometry/parsec_cmd/metallicity_grid/requests.json"))
    parser.add_argument("--output-dir", type=Path, default=Path("input_files/source_photometry/parsec_cmd/metallicity_grid"))
    parser.add_argument("--component", action="append", help="Fetch only the named component. Can be repeated.")
    parser.add_argument("--mh", type=float, action="append", help="Fetch only this [M/H]. Can be repeated.")
    parser.add_argument("--force", action="store_true", help="Refetch raw tables even when they already exist.")
    args = parser.parse_args()

    try:
        import ezpadova
    except ImportError as exc:
        raise SystemExit("Install ezpadova first: python -m pip install git+https://github.com/mfouesneau/ezpadova") from exc

    config = read_config(args.config)
    allowed_components = set(args.component or [])
    allowed_mh = set(args.mh or [])
    all_prime_rows: list[dict] = []
    all_roman_rows: list[dict] = []
    fetched_raw_paths: set[Path] = set()

    for component in config["components"]:
        if allowed_components and component["name"] not in allowed_components:
            continue
        log_ages = component_log_ages(component)
        multi_age = "logAge_grid" in component
        if multi_age and allowed_components:
            print(f"{component['name']} has {len(log_ages)} logAge grid points")
        for mh in component["mh_grid"]:
            if allowed_mh and mh not in allowed_mh:
                continue
            for log_age in log_ages:
                component_at_age = {**component, "logAge": log_age}
                raw_tables: dict[str, list[dict[str, float]]] = {}
                for phot_name, phot_cfg in config["photometric_systems"].items():
                    raw_path = grid_paths(args.output_dir, component, mh, log_age, multi_age, phot_name)["raw"]
                    if raw_path.exists() and (not args.force or raw_path in fetched_raw_paths):
                        print(f"reuse {component['name']} {phot_name} logAge={log_age} MH={mh}")
                        raw_tables[phot_name] = read_raw_table(raw_path)
                        continue
                    print(f"fetch {component['name']} {phot_name} logAge={log_age} MH={mh}", flush=True)
                    df = ezpadova.get_isochrones(
                        logage=(log_age, log_age, 0.0),
                        MH=(mh, mh, 0.0),
                        photsys_file=phot_cfg["photsys_file"],
                        track_parsec=config["defaults"]["track_parsec"],
                        track_colibri=config["defaults"]["track_colibri"],
                        photsys_version=config["defaults"]["photsys_version"],
                        kind_interp=config["defaults"]["kind_interp"],
                        kind_mag=config["defaults"]["kind_mag"],
                        output_kind=config["defaults"]["output_kind"],
                        output_evstage=config["defaults"]["output_evstage"],
                    )
                    save_raw_dataframe(df, raw_path, component=component_at_age, mh=mh, phot_name=phot_name)
                    fetched_raw_paths.add(raw_path)
                    raw_tables[phot_name] = table_rows(df)

                prime_rows, roman_rows = build_rows(component_at_age, mh, raw_tables)
                prime_columns = config["normalized_tables"]["prime_parsec"]["columns"]
                roman_columns = config["normalized_tables"]["roman_parsec"]["columns"]
                paths = grid_paths(args.output_dir, component, mh, log_age, multi_age)
                write_table(paths["normalized_prime"], prime_columns, prime_rows)
                write_table(paths["normalized_roman"], roman_columns, roman_rows)
                all_prime_rows.extend(prime_rows)
                all_roman_rows.extend(roman_rows)

    if all_prime_rows:
        write_table(args.output_dir / "normalized" / "all_prime_parsec.dat", config["normalized_tables"]["prime_parsec"]["columns"], all_prime_rows)
    if all_roman_rows:
        write_table(args.output_dir / "normalized" / "all_roman_parsec.dat", config["normalized_tables"]["roman_parsec"]["columns"], all_roman_rows)
    write_checksums(args.output_dir)


if __name__ == "__main__":
    main()
