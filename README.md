[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-yellow.svg)](http://www.gnu.org/licenses/gpl-3.0)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4784948.svg)](https://doi.org/10.5281/zenodo.4784948)
[![arXiv](http://img.shields.io/badge/arXiv-2104.03306-orange.svg?style=flat)](https://arxiv.org/abs/2104.03306)

# genulens v2

`genulens` ("generate microlensing") simulates microlensing events with the
Galactic model of [Koshimoto, Baba & Bennett (2021), ApJ, 917, 78](https://ui.adsabs.harvard.edu/abs/2021ApJ...917...78K/abstract).
The model is optimized for bulge sightlines and is most appropriate around
`|l| < 10 deg` and `|b| < 7 deg`.

Version 2 keeps the command-line simulator but refactors the C++ core into
shared model and simulation objects, adds a direct Python API, and exposes new
source-population and rate-map workflows.

Please cite Koshimoto, Baba & Bennett (2021) and
[Koshimoto & Ranc (2021), Zenodo.4784948](http://doi.org/10.5281/zenodo.4784948)
if you use this code in your research. A separate star simulator,
[`genstars`](https://github.com/nkoshimoto/genstars), is also available.

The copyright of the included supplementary code `option.cpp` belongs to
Ian A. Bond and Takahiro Sumi.

## What Is New

- Refactored C++ simulation core with reusable objects under `src/genulens/`.
- Direct Python bindings that call the C++ core without running `./genulens` as
  a subprocess or parsing CLI stdout.
- Typed Python configuration through `genulens.Config`.
- NumPy-compatible simulation results with CLI-style event labels.
- Python custom likelihood functions evaluated inside the Monte Carlo event
  loop.
- Source-forward mode that attaches isochrone-based source properties such as
  source mass, radius, effective temperature, surface gravity, angular radius,
  metallicity metadata, and absolute magnitudes.
- PARSEC and MIST source-property tables, including MIST alpha-enhanced
  systematics experiments.
- Bundled genstars-style `E(J-Ks)` extinction map support with adjustable
  extinction scale and offset.
- Direct optical-depth and event-rate summary APIs for single sightlines and
  maps.
- Refactored `pre_gapmoe` helper tools linked against the same Galactic model
  core.

## Installation

`genulens` requires a C++ compiler, CMake, and GSL.

Check that GSL is visible:

```bash
gsl-config --libs
```

Clone and build the command-line program:

```bash
git clone https://github.com/nkoshimoto/genulens.git
cd genulens
make
./genulens
```

If GSL is installed in a non-standard prefix:

```bash
GSL_ROOT=/path/to/gsl make
```

Additional build targets:

```bash
make python
make pre_gapmoe
make test
make clean
```

For Python use from a source checkout:

```bash
make python
PYTHONPATH=build python -c "import genulens; print(genulens.__file__)"
```

The Python extension can also be built with:

```bash
pip install .
```

Editable installs are supported in environments with `scikit-build-core` and
`pybind11`:

```bash
pip install -e .
```

## Quick Start

Run a plain simulation from Python:

```python
import pandas as pd
import genulens

cfg = genulens.Config(l=1.0, b=-3.9, n_simu=100_000, seed=42)
result = genulens.simulate(cfg)
df = pd.DataFrame(result.to_numpy(), columns=result.columns)
```

The default Python result uses a `VERBOSITY=3`-style event layout:

```text
wtj, M_L, D_L, D_S, t_E, theta_E, pi_E, pi_EN, pi_EE,
mu_rel, mu_rel_N, mu_rel_E, mu_Sl, mu_Sb, I_L, K_L, iS, iL, fREM
```

Unlike the historical CLI output, the Python table includes `mu_rel_N` and
`mu_rel_E` directly.

Configure observation constraints and model options through typed objects:

```python
cfg.observation.tE_obs = 54.5
cfg.observation.tE_err = 5.0
cfg.sampling.small_gamma = 1
cfg.model.imf.alpha2 = -1.35
```

Pass a Python likelihood function:

```python
def likelihood(event):
    return 1.0 if 10.0 < event.t_E < 100.0 else 0.0

result = genulens.simulate(cfg, likelihood=likelihood)
```

## Source-Forward Mode

Classic source mode keeps the historical luminosity-function source selection:

```python
cfg = genulens.Config(l=1.0, b=-3.9, n_simu=20_000, seed=42)
cfg.use_classic_source(i_min=14.0, i_max=21.0)
```

Isochrone source mode samples a concrete source star and appends source
properties to each simulated event:

```python
cfg = genulens.Config(l=1.0, b=-3.9, n_simu=20_000, seed=42)
cfg.use_isochrone_source(
    i_min=14.0,
    i_max=21.0,
    band="Imag",
    photometry="prime",
    apparent=True,
    min_mass=0.09,
    max_mass=2.0,
)
cfg.use_genstars_extinction_map(extinction_law=1, extinction_map=1)

result = genulens.simulate(cfg)
df = pd.DataFrame(result.to_numpy(), columns=result.columns)
```

Isochrone mode appends columns such as:

```text
logage_S, MH_S, M_S_ini, M_S, R_S, teff_S, logg_S,
theta_S, M_Vmag_S, M_Imag_S, M_Hmag_2mass_S, ...
```

Apparent-magnitude source cuts require extinction. You can provide manual band
extinctions or use the genstars-style extinction law/map:

```python
cfg.use_genstars_extinction(dm_rc=14.5, ejk_rc=1.0)

cfg.use_genstars_extinction_map(
    extinction_law=1,
    extinction_map=1,
    ejk_scale=1.0,
    ejk_offset=0.0,
)
```

## Isochrone Systematics

Source-property predictions can use PARSEC or MIST normalized isochrone tables:

```python
cfg.source.mode = "isochrone"
cfg.source.photometry = "roman"
cfg.source.isochrone_family = "mist"
cfg.source.isochrone_abundance = "solar_scaled"
```

For alpha-enhancement systematics, use MIST alpha-enhanced tables or a mixture:

```python
cfg.source.isochrone_model = "alpha_mixture"
cfg.source.secondary_isochrone_family = "mist"
cfg.source.secondary_isochrone_abundance = "alpha_enhanced"
cfg.source.secondary_isochrone_alpha_fe = 0.4
cfg.source.alpha_enhanced_fraction = 0.5
```

PARSEC alpha-enhanced tables are not available in the bundled table set; trying
to select PARSEC alpha-enhanced models raises an explicit error. Treat the alpha
mixture as a systematics knob, not as a full chemical-evolution model.

## Optical Depth and Event-Rate Summaries

The Python API can return optical-depth and event-rate summaries directly from
the C++ simulation state:

```python
cfg = genulens.Config(l=1.0, b=-3.9, n_simu=10_000, seed=42)
cfg.use_classic_source(i_min=14.0, i_max=21.0)
cfg.use_genstars_extinction_map(extinction_law=1, extinction_map=1)
cfg.observation.IL_err = 0.0

summary = genulens.compute_rate_summary(cfg)
print(summary.tau)
print(summary.event_rate_per_star_per_year)
print(summary.event_rate_per_deg2_per_year)
```

For maps, evaluate one config per sightline:

```python
configs = []
for b in b_values:
    for l in l_values:
        cfg_lb = genulens.Config(l=float(l), b=float(b), n_simu=10_000, seed=42)
        cfg_lb.use_isochrone_source(
            i_min=14.0,
            i_max=21.0,
            band="Imag",
            photometry="prime",
            apparent=True,
        )
        cfg_lb.use_genstars_extinction_map(extinction_law=1, extinction_map=1)
        cfg_lb.observation.IL_err = 0.0
        configs.append(cfg_lb)

rate_map = genulens.compute_rate_summaries(configs)
df = pd.DataFrame(rate_map.to_numpy(), columns=rate_map.columns)
```

## Documentation and Examples

Project documentation is under [`docs/`](docs/):

- [Python API](docs/python_api.md): Python configuration, result tables, custom
  likelihoods, source modes, extinction, isochrone systematics, and rate
  summaries.
- [Architecture](docs/architecture.md): current object boundaries and
  collaborator-facing design notes.
- [pre_gapmoe tools](docs/pre_gapmoe.md): helper-tool build and usage notes.

Notebook examples:

- [`examples/python_binding.ipynb`](examples/python_binding.ipynb): basic Python
  simulation, result tables, plotting, and custom likelihoods.
- [`examples/source_isochrone_systematics.ipynb`](examples/source_isochrone_systematics.ipynb):
  source-property and HR-diagram comparisons across isochrone assumptions.
- [`examples/rate_summary_map.ipynb`](examples/rate_summary_map.ipynb): optical
  depth, event-rate, and extinction maps for I- and H-band source selections.

Script examples:

- [`examples/source_isochrone_systematics.py`](examples/source_isochrone_systematics.py):
  fast source-property systematics check; pass `--use-selection` to fold the
  band cut into the source-distance prior.

## Scientific Notes

- The Galactic model was developed for bulge microlensing applications. Be
  careful outside the recommended `|l| < 10 deg`, `|b| < 7 deg` range.
- Source-forward mode is intended for physically motivated prior construction
  and systematics exploration. Its output depends on the isochrone library,
  abundance assumptions, IMF, source-population prior, extinction law, and
  extinction map.
- MIST/PARSEC differences include stellar evolution, atmosphere/bolometric
  correction choices, passbands, and zero-points. Do not interpret differences
  between libraries as only alpha-enhancement effects.
- Apparent source selections require an extinction model. For maps, use
  line-of-sight dependent extinction such as `use_genstars_extinction_map`
  rather than reusing one fixed red-clump extinction value everywhere.
- Building source-forward selected-density grids can take time. The simulator
  prints progress to `stderr`, and repeated runs in the same Python process can
  reuse cached selection probabilities.

## Legacy CLI and pre_gapmoe

The historical CLI remains available:

```bash
./genulens
```

The original usage guide is still useful for legacy command-line options:

- [Usage.pdf](Usage.pdf)

The refactored `pre_gapmoe` helpers are built with:

```bash
make pre_gapmoe
```

They are written to `build/pre_gapmoe/` when using the CMake build. See
[docs/pre_gapmoe.md](docs/pre_gapmoe.md) for details.

## Release History

- v2: refactored C++ core, direct Python API, source-forward isochrone support,
  extinction-map support, custom Python likelihoods, and rate-summary APIs.
- v1.2, June-July 2022: importance sampling, NSD component, updated Galactic
  Center position, revised usage documentation, and related `genstars` release.
- v1.1, June 2021: switched to the GSL random number generator.
- v1.0, May 2021: initial public release.
