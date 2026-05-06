# Python API

The Python API calls the shared C++ simulation core directly. It does not run
`./genulens` as a subprocess and it does not parse command-line stdout.

The main entry point is:

```python
import genulens

cfg = genulens.Config(l=1.0, b=-3.9, n_simu=100_000, seed=42)
result = genulens.simulate(cfg)
arr = result.to_numpy()
```

`result` is a `SimulationResult`. `result.columns` gives the column order used
by `to_numpy()`.

## Build and import

From a source checkout:

```bash
cmake --build build --target genulens_python
PYTHONPATH=build python -c "import genulens; print(genulens.__file__)"
```

For editable installs:

```bash
pip install -e .
```

The package currently exposes the Python extension built from the C++ core. The
Python API should be treated as the first public milestone of the refactor, not
as a fully frozen long-term interface.

## Basic configuration

```python
cfg = genulens.Config()
cfg.l = 1.0
cfg.b = -3.9
cfg.n_simu = 20_000
cfg.seed = 42

result = genulens.simulate(cfg)
```

`n_simu` is the target number of accepted events returned by the simulation.

`n_like_min` is separate:

```python
cfg.sampling.n_like_min = 5_000
```

This is equivalent to CLI `NlikeMIN`. It allows the internal run length to grow
when too few events pass the likelihood. Therefore `n_simu` and `n_like_min` can
be used together, but the final internal number of generated events is no longer
fixed when `n_like_min > 0`.

## Result columns

The default event result uses the `VERBOSITY=3` style layout:

- `wtj`
- `M_L`
- `D_L`
- `D_S`
- `t_E`
- `theta_E`
- `pi_E`
- `pi_EN`
- `pi_EE`
- `mu_rel`
- `mu_rel_N`
- `mu_rel_E`
- `mu_Sl`
- `mu_Sb`
- `I_L`
- `K_L`
- `iS`
- `iL`
- `fREM`

Example conversion to pandas:

```python
import pandas as pd

df = pd.DataFrame(result.to_numpy(), columns=result.columns)
```

Use the `wtj` column for histograms, quantiles, and population fractions.

`cfg.sampling.verbosity` selects the event layout:

```python
cfg.sampling.verbosity = 3
result = genulens.simulate(cfg)
print(result.columns)
```

When left at the Python default, the result uses the same layout as
`verbosity = 3`. The Python result includes `mu_rel_N` and `mu_rel_E` directly,
so callers do not need to reconstruct the relative proper-motion components from
parallax components. The older descriptive result labels such as
`lens_mass_msun`, `lens_distance_pc`, and `source_distance_pc` are not part of
the public result table.

## Observation constraints

Observation likelihoods are configured through `cfg.observation`:

```python
cfg = genulens.Config(l=1.0, b=-3.9, n_simu=20_000, seed=42)
cfg.observation.tE_obs = 54.5
cfg.observation.tE_err = 5.0
cfg.observation.thetaE_obs = 0.55
cfg.observation.thetaE_err = 0.15

result = genulens.simulate(cfg)
```

These values are passed as typed C++ fields. They are not converted into CLI
option strings.

## Source, model, and sampling options

Common configuration groups are exposed as nested objects:

```python
cfg.source.i_min = 14.0
cfg.source.i_max = 21.0
cfg.source.ai_rc = 1.5
cfg.source.evi_rc = 1.2

cfg.model.imf.alpha2 = -1.35
cfg.model.density.stellar_halo = 0
cfg.model.nsd.enabled = 0
cfg.model.kinematics.omega_p = 45.0

cfg.sampling.binary = 1
cfg.sampling.remnant = 1
cfg.sampling.small_gamma = 1
cfg.sampling.verbosity = 3
```

Important sampling option mappings:

| Python field | CLI option |
|---|---|
| `cfg.sampling.n_like_min` | `NlikeMIN` |
| `cfg.sampling.small_gamma` | `SMALLGAMMA` |
| `cfg.sampling.no_gamma_importance_sampling` | `NoGAMMAIS` |
| `cfg.sampling.weight_lens_mass` | `wtM_L` |
| `cfg.sampling.weight_lens_distance` | `wtD_L` |
| `cfg.sampling.gamma_ds` | `gammaDs` |
| `cfg.sampling.binary` | `BINARY` |
| `cfg.sampling.remnant` | `REMNANT` |
| `cfg.sampling.only_white_dwarf` | `onlyWD` |
| `cfg.sampling.uniform_likelihood` | `UNIFORM` |
| `cfg.sampling.verbosity` | `VERBOSITY` |

`small_gamma = 1` keeps low-Gamma events and carries Gamma through the event
weight. With `small_gamma = 0`, low-Gamma events can be rejected by the usual
Gamma rejection step.

## Forward source annotations

The default simulator keeps the historical LF/CMF source-selection and
event-rate weighting. To append stellar properties drawn from the shared
isochrone lookup as annotations, enable forward source mode:

```python
cfg = genulens.Config(l=1.0, b=-3.9, n_simu=10_000, seed=42)
cfg.forward_source.enabled = 1
cfg.forward_source.photometry = "roman"
cfg.forward_source.min_initial_mass_msun = 0.1
cfg.forward_source.max_initial_mass_msun = 1.0

result = genulens.simulate(cfg)
df = pd.DataFrame(result.to_numpy(), columns=result.columns)
```

This mode appends source-property columns after the event columns:

- `source_log_age`
- `source_metallicity_mh`
- `source_zini`
- `source_initial_mass_msun`
- `source_current_mass_msun`
- `source_radius_rsun`
- `source_teff_k`
- `source_logg`
- `source_angular_radius_microarcsec`
- `source_abs_<band>mag`

These columns are annotations for forward-prior studies. They do not yet replace
the legacy LF/CMF source-selection machinery used by the event sampler. Rows
without a matched stellar source entry can contain `NaN` source annotations.

For standalone source-population tests, use `ForwardSourceGenerator` directly:

```python
generator = genulens.ForwardSourceGenerator()
samples = generator.sample_many(n=1000, component=8, distance_pc=8000.0, seed=3)
source_df = pd.DataFrame(samples.to_numpy(), columns=samples.columns)
```

## Custom likelihood

Pass a Python callable as `likelihood=`:

```python
import math

def likelihood(event):
    if event.D_L >= event.D_S:
        return 0.0
    z = (event.t_E - 54.5) / 8.0
    return math.exp(-0.5 * z * z)

result = genulens.simulate(cfg, likelihood=likelihood)
```

The callable receives an `Event` object and returns a multiplicative likelihood.
A return value less than or equal to zero rejects the event.

The `Event` object uses the same event labels as the default result table, such
as `event.t_E`, `event.M_L`, `event.D_L`, `event.D_S`, `event.pi_EN`, and
`event.mu_rel_N`.

Avoid extremely narrow hard cuts unless you also understand the runtime
implications. The sampler will keep drawing until it has enough accepted events,
so a likelihood that rejects nearly everything can run for a long time.

## Shortcut API

For quick experiments:

```python
result = genulens.ruc(l=0.0, b=-2.0, n_simu=10_000, seed=1)
```

Use `Config` for any analysis that needs observation constraints, model
parameters, source selection, or sampling options.
