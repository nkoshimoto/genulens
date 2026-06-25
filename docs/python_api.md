# Python API

The Python API calls the shared C++ simulation core directly. It does not run
`./genulens` as a subprocess and does not parse command-line stdout.

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
make python
PYTHONPATH=build python -c "import genulens; print(genulens.__file__)"
```

For editable installs:

```bash
pip install -e .
```

The package exposes a Python extension built from the C++ core. The v2 Python
API is intended to be the public user surface for scripted analyses, notebooks,
custom likelihoods, source-forward simulations, and rate-map workflows.
Collaborators should still treat lower-level implementation details as
evolving.

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
when too few events pass the likelihood. Therefore `n_simu` and `n_like_min`
can be used together, but the final internal number of generated events is no
longer fixed when `n_like_min > 0`.

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
so callers do not need to reconstruct the relative proper-motion components
from parallax components. The older descriptive result labels such as
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

## Configuration groups

Common configuration groups are exposed as nested objects:

```python
cfg.source.mode = "classic"
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

Avoid extremely narrow hard cuts unless you understand the runtime implications.
The sampler will keep drawing until it has enough accepted events, so a
likelihood that rejects nearly everything can run for a long time.

## Shortcut API

For quick experiments:

```python
result = genulens.ruc(l=0.0, b=-2.0, n_simu=10_000, seed=1)
```

Use `Config` for any analysis that needs observation constraints, model
parameters, source selection, or sampling options.

## Related pages

- [Source-forward mode](source_forward.md)
- [Extinction](extinction.md)
- [Optical-depth and event-rate maps](rate_maps.md)
- [Isochrone systematics](isochrone_systematics.md)
