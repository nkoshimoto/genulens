from pathlib import Path
import sys

import numpy as np
import pandas as pd


def import_genulens():
    repo_root = Path.cwd().resolve()
    for path in [repo_root, *repo_root.parents]:
        if (path / "CMakeLists.txt").exists() and (path / "build").exists():
            repo_root = path
            break
    build_dir = repo_root / "build"
    if build_dir.exists():
        sys.path.insert(0, str(build_dir))
    import genulens

    return genulens


def as_dataframe(result):
    return pd.DataFrame(result.to_numpy(), columns=result.columns)


def base_config(genulens):
    cfg = genulens.Config(l=1.0, b=-3.9, n_simu=5, seed=2026)
    cfg.forward_source.enabled = 1
    cfg.forward_source.photometry = "roman"
    cfg.forward_source.min_initial_mass_msun = 0.1
    cfg.forward_source.max_initial_mass_msun = 2.0
    return cfg


def parsec_solar_config(genulens):
    cfg = base_config(genulens)
    cfg.forward_source.isochrone_family = "parsec"
    cfg.forward_source.isochrone_abundance = "solar_scaled"
    return cfg


def mist_alpha_config(genulens):
    cfg = base_config(genulens)
    cfg.forward_source.isochrone_family = "mist"
    cfg.forward_source.isochrone_abundance = "alpha_enhanced"
    cfg.forward_source.isochrone_alpha_fe = 0.4
    return cfg


def main():
    genulens = import_genulens()
    source_columns = [
        "D_S",
        "MH_S",
        "M_S_ini",
        "M_S",
        "R_S",
        "teff_S",
        "logg_S",
        "theta_S",
        "M_F146mag_S",
    ]
    configs = {
        "parsec_solar": parsec_solar_config(genulens),
        "mist_alpha": mist_alpha_config(genulens),
    }

    rows = []
    for name, cfg in configs.items():
        df = as_dataframe(genulens.simulate(cfg))
        assert len(df) == cfg.n_simu, name
        assert np.isfinite(df[source_columns].to_numpy()).all(), name
        rows.append({
            "run": name,
            "rows": len(df),
            "median_D_S": df["D_S"].median(),
            "median_teff_S": df["teff_S"].median(),
            "median_R_S": df["R_S"].median(),
            "median_M_F146mag_S": df["M_F146mag_S"].median(),
        })

    spec = genulens.IsochroneLibrarySpec()
    spec.family = "parsec"
    spec.photometry = "roman"
    spec.abundance = "alpha_enhanced"
    try:
        genulens.default_isochrone_table_path(spec)
    except RuntimeError as exc:
        assert "PARSEC/CMD alpha-enhanced" in str(exc)
    else:
        raise AssertionError("PARSEC alpha unexpectedly resolved to a table path")

    print(pd.DataFrame(rows))


if __name__ == "__main__":
    main()
