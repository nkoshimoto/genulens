"""Compare source-forward predictions for alternative isochrone tables.

This script is intended for alpha-enhancement systematics tests. It runs the
same source prior with different fractions of a secondary isochrone table and
summarizes the source-property shifts. By default it uses PRIME photometry and
an apparent `Imag` source-selection range when `--use-selection` is supplied.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import numpy as np
import pandas as pd


def import_genulens():
    repo_root = Path(__file__).resolve().parents[1]
    build_dir = repo_root / "build"
    if build_dir.exists():
        sys.path.insert(0, str(build_dir))
    import genulens

    return genulens


def as_dataframe(result):
    return pd.DataFrame(result.to_numpy(), columns=result.columns)


def make_config(genulens, args, alpha_fraction):
    cfg = genulens.Config(l=args.l, b=args.b, n_simu=args.n_simu, seed=args.seed)
    # Keep this example focused on source selection and source-property
    # predictions; auto I-band extinction would otherwise also activate the
    # default lens-flux IL constraint.
    cfg.observation.IL_err = 0.0
    cfg.source.mode = "isochrone"
    cfg.source.photometry = args.photometry
    cfg.source.isochrone_model = "alpha_mixture" if alpha_fraction > 0.0 else "parsec_solar_scaled"
    cfg.source.isochrone_family = args.primary_family
    cfg.source.isochrone_abundance = args.primary_abundance
    cfg.source.isochrone_alpha_fe = args.primary_alpha_fe
    cfg.source.secondary_isochrone_family = args.secondary_family
    cfg.source.secondary_isochrone_abundance = args.secondary_abundance
    cfg.source.secondary_isochrone_alpha_fe = args.secondary_alpha_fe
    cfg.source.min_initial_mass_msun = args.min_mass
    cfg.source.max_initial_mass_msun = args.max_mass
    cfg.source.use_magnitude_selection = 1 if args.use_selection else 0
    if args.use_selection:
        cfg.source.band = args.selection_band
        cfg.source.min_magnitude = args.selection_min
        cfg.source.max_magnitude = args.selection_max
        cfg.source.apparent_magnitude = 1 if args.apparent else 0
    if args.primary_table:
        cfg.source.isochrone_table_path = args.primary_table
    if alpha_fraction > 0.0:
        cfg.source.alpha_enhanced_fraction = alpha_fraction
        if args.secondary_table:
            cfg.source.secondary_isochrone_table_path = args.secondary_table
    if args.alpha_components:
        cfg.source.alpha_enhanced_fraction = 0.0
        cfg.source.alpha_enhanced_components = [
            int(item) for item in args.alpha_components.split(",") if item
        ]
        cfg.source.alpha_enhanced_component_fractions = [
            alpha_fraction for _ in cfg.source.alpha_enhanced_components
        ]
        if alpha_fraction > 0.0:
            cfg.source.isochrone_model = "alpha_mixture"

    if args.use_selection and args.apparent:
        cfg.source.extinction_mode = "genstars"
        cfg.source.extinction_law = args.extinction_law
        cfg.source.ejk_rc = args.ejk_rc
        cfg.source.dm_rc = args.dm_rc
    return cfg


def weighted_quantile(values, weights, quantiles):
    values = np.asarray(values)
    weights = np.asarray(weights)
    mask = np.isfinite(values) & np.isfinite(weights) & (weights > 0)
    values = values[mask]
    weights = weights[mask]
    if len(values) == 0:
        return [np.nan for _ in quantiles]
    order = np.argsort(values)
    values = values[order]
    weights = weights[order]
    cdf = np.cumsum(weights)
    cdf /= cdf[-1]
    return np.interp(quantiles, cdf, values)


def summarize(df, columns):
    rows = []
    weights = df["wtj"].to_numpy()
    for column in columns:
        if column not in df:
            continue
        q16, q50, q84 = weighted_quantile(df[column].to_numpy(), weights, [0.16, 0.50, 0.84])
        rows.append({"column": column, "q16": q16, "q50": q50, "q84": q84})
    return pd.DataFrame(rows)


def write_hr_plot(sample_df, output_path):
    import matplotlib.pyplot as plt

    if not {"teff_S", "M_Imag_S", "wtj", "alpha_enhanced_fraction"}.issubset(sample_df.columns):
        return False
    plot_mask = (
        np.isfinite(sample_df["teff_S"]) &
        np.isfinite(sample_df["M_Imag_S"]) &
        np.isfinite(sample_df["wtj"]) &
        (sample_df["wtj"] > 0)
    )
    plot_source = sample_df.loc[plot_mask].copy()
    if plot_source.empty:
        return False
    fractions = list(dict.fromkeys(plot_source["alpha_enhanced_fraction"].tolist()))
    if not fractions:
        return False
    fig, axes = plt.subplots(
        1,
        len(fractions),
        figsize=(5 * len(fractions), 4),
        sharex=True,
        sharey=True,
        constrained_layout=True,
    )
    if len(fractions) == 1:
        axes = [axes]
    scatter = None
    for ax, fraction in zip(axes, fractions):
        plot_df = plot_source.loc[
            plot_source["alpha_enhanced_fraction"] == fraction,
            ["teff_S", "M_Imag_S", "wtj"],
        ].copy()
        if len(plot_df) > 4000:
            plot_df = plot_df.sample(4000, random_state=0)
        color = np.log10(plot_df["wtj"].to_numpy())
        scatter = ax.scatter(
            plot_df["teff_S"],
            plot_df["M_Imag_S"],
            c=color,
            s=8,
            alpha=0.45,
            cmap="viridis",
            linewidths=0,
        )
        ax.invert_xaxis()
        ax.invert_yaxis()
        ax.set_title(f"alpha fraction = {fraction:g}")
        ax.set_xlabel("teff_S [K]")
        ax.set_ylabel("M_Imag_S")
    if scatter is not None:
        fig.colorbar(scatter, ax=axes, label="log10(wtj)")
    fig.savefig(output_path, dpi=160)
    plt.close(fig)
    return True


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--secondary-table", default="", help="Optional normalized secondary table path.")
    parser.add_argument("--alpha-table", default="", help="Legacy alias for --secondary-table.")
    parser.add_argument("--primary-table", default="", help="Optional normalized primary table path.")
    parser.add_argument("--primary-family", default="parsec", choices=["parsec", "mist"])
    parser.add_argument("--primary-abundance", default="solar_scaled", choices=["solar_scaled", "alpha_enhanced"])
    parser.add_argument("--primary-alpha-fe", type=float, default=0.0)
    parser.add_argument("--secondary-family", default="mist", choices=["parsec", "mist"])
    parser.add_argument("--secondary-abundance", default="alpha_enhanced", choices=["solar_scaled", "alpha_enhanced"])
    parser.add_argument("--secondary-alpha-fe", type=float, default=0.4)
    parser.add_argument("--photometry", default="prime", choices=["roman", "prime"])
    parser.add_argument("--fractions", default="0,0.5,1", help="Comma-separated secondary fractions.")
    parser.add_argument(
        "--alpha-components",
        default="8,9",
        help="Optional comma-separated component indices that receive the tested alpha fraction, e.g. 8,9 for bar/NSD.",
    )
    parser.add_argument("--n-simu", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--l", type=float, default=1.0)
    parser.add_argument("--b", type=float, default=-3.9)
    parser.add_argument("--min-mass", type=float, default=0.1)
    parser.add_argument("--max-mass", type=float, default=2.0)
    parser.add_argument("--selection-band", default="Imag")
    parser.add_argument("--selection-min", type=float, default=12.0)
    parser.add_argument("--selection-max", type=float, default=21.0)
    parser.add_argument("--use-selection", action="store_true", help="Fold the band cut into the source-distance prior. This can be slow for MIST.")
    parser.add_argument("--absolute", action="store_true", help="Interpret the selection as absolute magnitudes.")
    parser.add_argument("--dm-rc", type=float, default=14.5)
    parser.add_argument("--ejk-rc", type=float, default=1.0)
    parser.add_argument("--extinction-law", type=int, default=1)
    parser.add_argument("--output", type=Path, default=Path("source_isochrone_systematics_summary.csv"))
    parser.add_argument(
        "--samples-output",
        type=Path,
        default=Path("source_isochrone_systematics_samples.csv"),
        help="Write per-event source-property samples for direct inspection.",
    )
    parser.add_argument(
        "--hr-output",
        type=Path,
        default=Path("source_isochrone_systematics_hr.png"),
        help="Write an HR-diagram scatter plot. Use an empty string to disable.",
    )
    args = parser.parse_args()
    args.apparent = not args.absolute
    if args.alpha_table and not args.secondary_table:
        args.secondary_table = args.alpha_table
    genulens = import_genulens()
    fractions = [float(item) for item in args.fractions.split(",") if item]
    magnitude_column = f"M_{args.selection_band}_S"
    photometry_magnitude_columns = (
        ["M_Vmag_S", "M_Imag_S", "M_Jmag_2mass_S", "M_Hmag_2mass_S", "M_Ksmag_2mass_S"]
        if args.photometry == "prime"
        else ["M_F087mag_S", "M_F146mag_S", "M_F213mag_S"]
    )
    columns = [
        "D_S",
        "logage_S",
        "MH_S",
        "M_S_ini",
        "M_S",
        "R_S",
        "teff_S",
        "logg_S",
        "theta_S",
        magnitude_column,
        *photometry_magnitude_columns,
        "VI_S",
    ]
    columns = list(dict.fromkeys(columns))

    summaries = []
    sample_tables = []
    for fraction in fractions:
        cfg = make_config(genulens, args, fraction)
        df = as_dataframe(genulens.simulate(cfg))
        if {"M_Vmag_S", "M_Imag_S"}.issubset(df.columns):
            df["VI_S"] = df["M_Vmag_S"] - df["M_Imag_S"]
        summary = summarize(df, columns)
        summary.insert(0, "alpha_enhanced_fraction", fraction)
        summaries.append(summary)
        sample_columns = ["iS", "D_S", "wtj", *[column for column in columns if column != "D_S"]]
        sample_columns = [column for column in sample_columns if column in df]
        samples = df[sample_columns].copy()
        samples.insert(0, "alpha_enhanced_fraction", fraction)
        sample_tables.append(samples)

    out = pd.concat(summaries, ignore_index=True)
    out.to_csv(args.output, index=False)
    sample_out = pd.concat(sample_tables, ignore_index=True)
    sample_out.to_csv(args.samples_output, index=False)
    wrote_hr = False
    if args.hr_output:
        wrote_hr = write_hr_plot(sample_out, args.hr_output)
    print(out)
    print(f"wrote {args.output}")
    print(f"wrote {args.samples_output}")
    if wrote_hr:
        print(f"wrote {args.hr_output}")


if __name__ == "__main__":
    main()
