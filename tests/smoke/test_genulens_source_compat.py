import os
import re
import subprocess
from pathlib import Path

import genulens
import pytest


LEGACY_DIR_ENV = os.environ.get("GENULENS_SOURCE_COMPAT_DIR")
LEGACY_DIR = Path(LEGACY_DIR_ENV) if LEGACY_DIR_ENV else None
LEGACY_EXE = LEGACY_DIR / "genulens_source" if LEGACY_DIR is not None else None


def relative_difference(new, reference):
    return abs(new / reference - 1.0)


def parse_legacy_summary(stdout):
    density_match = re.search(r"Source number density=\s*([0-9.eE+-]+)", stdout)
    rate_match = re.search(
        r"avetE=\s*([0-9.eE+-]+)\s+days,"
        r"\s*medtE=\s*([0-9.eE+-]+)\s+days,"
        r"\s*tau=\s*([0-9.eE+-]+)\s*,"
        r"\s*event_rate=\s*([0-9.eE+-]+)\s*/star/yr"
        r"\s*or\s*([0-9.eE+-]+)\s*/deg\^2/yr",
        stdout,
    )
    if density_match is None or rate_match is None:
        raise AssertionError("could not parse genulens_source rate summary")
    return {
        "source_density_arcmin2": float(density_match.group(1)),
        "mean_tE_days": float(rate_match.group(1)),
        "median_tE_days": float(rate_match.group(2)),
        "tau": float(rate_match.group(3)),
        "event_rate_per_star_per_year": float(rate_match.group(4)),
        "event_rate_per_deg2_per_year": float(rate_match.group(5)),
    }


def run_legacy_genulens_source(l_deg, b_deg, n_simu=1000):
    if LEGACY_DIR is None or LEGACY_EXE is None:
        raise RuntimeError("set GENULENS_SOURCE_COMPAT_DIR to run this compatibility test")
    command = [
        str(LEGACY_EXE),
        "l",
        str(l_deg),
        "b",
        str(b_deg),
        "Nsimu",
        str(n_simu),
        "Magrange",
        "14",
        "21",
        "AIrc",
        "1.0",
        "iMag",
        "1",
        "EXTLAW",
        "1",
        "VERBOSITY",
        "0",
        "seed",
        "2026",
    ]
    completed = subprocess.run(
        command,
        cwd=LEGACY_DIR,
        text=True,
        capture_output=True,
        check=True,
    )
    return parse_legacy_summary(completed.stdout)


def run_current_genulens(l_deg, b_deg, n_simu=1000):
    cfg = genulens.Config(l=l_deg, b=b_deg, n_simu=n_simu, seed=2026)
    cfg.use_classic_source(i_min=14.0, i_max=21.0)
    cfg.use_manual_extinction(ai_rc=1.0, dm_rc=0.0)
    cfg.observation.IL_err = 0.0
    summary = genulens.compute_rate_summary(cfg)
    return {
        "source_density_arcmin2": summary.source_density_arcmin2,
        "mean_tE_days": summary.mean_tE_days,
        "median_tE_days": summary.median_tE_days,
        "tau": summary.tau,
        "event_rate_per_star_per_year": summary.event_rate_per_star_per_year,
        "event_rate_per_deg2_per_year": summary.event_rate_per_deg2_per_year,
    }


@pytest.mark.skipif(
    LEGACY_EXE is None or not LEGACY_EXE.exists(),
    reason=(
        "legacy genulens_source executable is not available; "
        "set GENULENS_SOURCE_COMPAT_DIR to enable this optional compatibility test"
    ),
)
def test_genulens_source_rate_compatibility_over_wide_field():
    grid = [
        (l_deg, b_deg)
        for b_deg in (-6.0, -4.0, -2.0)
        for l_deg in (-4.0, 0.0, 4.0)
    ]

    for l_deg, b_deg in grid:
        legacy = run_legacy_genulens_source(l_deg, b_deg)
        current = run_current_genulens(l_deg, b_deg)

        assert relative_difference(current["source_density_arcmin2"],
                                   legacy["source_density_arcmin2"]) < 0.02
        assert relative_difference(current["tau"], legacy["tau"]) < 0.04

        # The event rate uses the Monte Carlo mean tE. Keep this looser than
        # tau/source density, but still tight enough to catch rate-regression
        # bugs in the refactored sampler.
        assert relative_difference(current["event_rate_per_star_per_year"],
                                   legacy["event_rate_per_star_per_year"]) < 0.08
        assert relative_difference(current["event_rate_per_deg2_per_year"],
                                   legacy["event_rate_per_deg2_per_year"]) < 0.08
