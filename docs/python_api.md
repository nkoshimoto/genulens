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

The default event result contains:

- `weight`
- `tE`
- `thetaE`
- `piE`
- `lens_distance_pc`
- `source_distance_pc`
- `lens_mass_msun`
- `mu_rel_masyr`
- `lens_component`
- `source_component`

Example conversion to pandas:

```python
import pandas as pd

df = pd.DataFrame(result.to_numpy(), columns=result.columns)
```

Use the `weight` column for histograms, quantiles, and population fractions.

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

`small_gamma = 1` keeps low-Gamma events and carries Gamma through the event
weight. With `small_gamma = 0`, low-Gamma events can be rejected by the usual
Gamma rejection step.

## Custom likelihood

Pass a Python callable as `likelihood=`:

```python
import math

def likelihood(event):
    if event.lens_distance_pc >= event.source_distance_pc:
        return 0.0
    z = (event.tE - 54.5) / 8.0
    return math.exp(-0.5 * z * z)

result = genulens.simulate(cfg, likelihood=likelihood)
```

The callable receives an `Event` object and returns a multiplicative likelihood.
A return value less than or equal to zero rejects the event.

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
