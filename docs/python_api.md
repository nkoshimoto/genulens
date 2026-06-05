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

## Rate Summaries

The Python API can return line-of-sight optical-depth and event-rate summaries
directly from the C++ simulation state. These values are not parsed from CLI
stdout.

```python
cfg = genulens.Config(l=1.0, b=-3.9, n_simu=10_000, seed=42)
cfg.source.mode = "classic"
cfg.source.i_min = 14.0
cfg.source.i_max = 21.0
cfg.source.ai_rc = 1.0
cfg.source.dm_rc = 14.5
cfg.observation.IL_err = 0.0

summary = genulens.compute_rate_summary(cfg)
print(summary.tau)
print(summary.event_rate_per_star_per_year)
print(summary.event_rate_per_deg2_per_year)
```

`compute_rate_summary` copies the config, enables the optical-depth calculation
internally, runs the same event sampler used by `simulate`, and returns:

- `source_density_arcmin2`
- `source_density_raw_arcmin2`
- `tau`
- `mean_tE_days`
- `median_tE_days`
- `event_rate_per_star_per_year`
- `event_rate_per_deg2_per_year`
- `sum_gamma`
- `sum_tE_gamma`

The event rate follows the existing CLI summary convention:

```text
event_rate = 2 * tau / pi / mean_tE_days * 365.25
```

For quick maps, evaluate the same summary on an `l,b` grid:

```python
rate_map = genulens.compute_rate_map(
    l_values=np.linspace(-5.0, 5.0, 21),
    b_values=np.linspace(-6.0, -2.0, 9),
    base_config=cfg,
)
df = pd.DataFrame(rate_map.to_numpy(), columns=rate_map.columns)
```

`compute_rate_map` only changes `cfg.l` and `cfg.b`. It does not infer a
line-of-sight extinction map. If `cfg.source.ai_rc`, `cfg.source.dm_rc`, or the
genstars-style reference extinction values are fixed in `base_config`, the
visible source cut uses those same values at every map point.

For realistic visible-band maps, build one config per sightline and set the
red-clump extinction and distance modulus for that sightline explicitly:

```python
configs = []
for b in b_values:
    for l in l_values:
        cfg_lb = genulens.Config(l=l, b=b, n_simu=10_000, seed=42)
        cfg_lb.source.mode = "classic"
        cfg_lb.source.i_min = 14.0
        cfg_lb.source.i_max = 21.0
        cfg_lb.source.ai_rc = extinction_map_ai_rc(l, b)
        cfg_lb.source.dm_rc = red_clump_dm(l, b)
        cfg_lb.observation.IL_err = 0.0
        configs.append(cfg_lb)

rate_map = genulens.compute_rate_summaries(configs)
df = pd.DataFrame(rate_map.to_numpy(), columns=rate_map.columns)
```

The summary depends on the active source selection, extinction, source-density
model, and Monte Carlo settings. In particular, apparent source cuts require
valid extinction parameters, and event rates use the Monte Carlo estimate of
`mean_tE_days`.

## Source Modes

The Python API exposes two source-selection modes through `cfg.source`.

Classic mode keeps the historical luminosity-function/color-magnitude-function
source selection:

```python
cfg = genulens.Config(l=1.0, b=-3.9, n_simu=10_000, seed=42)
cfg.source.mode = "classic"
cfg.source.i_min = 14.0
cfg.source.i_max = 21.0
```

Isochrone mode uses the shared isochrone/IMF source prior and appends source
stellar properties to each simulated event:

```python
cfg = genulens.Config(l=1.0, b=-3.9, n_simu=10_000, seed=42)
cfg.source.mode = "isochrone"
cfg.source.photometry = "prime"
cfg.source.band = "Imag"
cfg.source.min_magnitude = 12.0
cfg.source.max_magnitude = 21.0
cfg.source.min_initial_mass_msun = 0.1
cfg.source.max_initial_mass_msun = 2.0
cfg.source.extinction_mode = "genstars"
cfg.source.extinction_law = 1
cfg.source.ejk_rc = 1.0
cfg.source.dm_rc = 14.5

result = genulens.simulate(cfg)
df = pd.DataFrame(result.to_numpy(), columns=result.columns)
```

Convenience methods are available for common setup:

```python
cfg.use_classic_source(i_min=14.0, i_max=21.0)
cfg.use_isochrone_source(i_min=12.0, i_max=21.0, band="Imag", photometry="prime")
cfg.use_genstars_extinction(dm_rc=14.5, ejk_rc=1.0)
```

Isochrone mode appends source-property columns after the event columns:

- `iS`
- `D_S`
- `logage_S`
- `MH_S`
- `M_S_ini`
- `M_S`
- `R_S`
- `teff_S`
- `logg_S`
- `theta_S`
- `M_<band>_S`

With `cfg.source.use_magnitude_selection = 1`, the simulator folds the
isochrone/IMF selection probability into the source-distance density grid
before sampling events. It then samples a concrete source star from the
selected initial-mass intervals for the sampled source distance and component.
For speed, the event-level source-star proposal reuses conservative
source-distance bins when finding allowed initial-mass intervals, then applies
the exact apparent or absolute magnitude cut at the sampled `D_S` before
accepting the source. This keeps large `n_simu` runs from repeating the same
isochrone interval search for every event while preserving the requested
selection support.
Set `cfg.source.use_magnitude_selection = 0` to append source properties
without conditioning the source-distance prior on a magnitude cut.

`cfg.source.min_initial_mass_msun` and `cfg.source.max_initial_mass_msun`
define the source-star IMF range used both for the selected source-density
normalization and for drawing the concrete source star. Treat them as part of
the physical source prior, not only as proposal bounds.

Cuts are interpreted as apparent magnitudes by default. Set
`cfg.source.apparent_magnitude = 0` for absolute-magnitude cuts. Apparent cuts
require extinction for the selected band. Manual mode uses direct band
extinctions at the reference red-clump distance:

```python
cfg.source.extinction_mode = "manual"
cfg.source.dm_rc = 14.5
cfg.source.av_rc = 2.0   # Vmag
cfg.source.ai_rc = 1.4   # Imag
cfg.source.aj_rc = 0.5   # Jmag / Jmag_2mass
cfg.source.ah_rc = 0.2   # Hmag / Hmag_2mass
cfg.source.ak_rc = 0.1   # Kmag / Ksmag_2mass / K_L
```

Automatic extinction follows the `genstars` law from `E(J-K)`:

```python
cfg.source.extinction_mode = "genstars"
cfg.source.extinction_law = 1
cfg.source.ejk_rc = 1.0
cfg.source.dm_rc = 14.5
```

This uses a single `E(J-Ks)` value for the current sightline. To read the
genstars extinction map directly, use:

```python
cfg.use_genstars_extinction_map(extinction_law=1, extinction_map=1)
```

You can adjust the map without editing the input file. The effective value is
`E(J-Ks)_eff = E(J-Ks)_map * ejk_scale + ejk_offset`:

```python
cfg.use_genstars_extinction_map(
    extinction_law=1,
    extinction_map=1,
    ejk_scale=1.15,
    ejk_offset=0.05,
)
```

or equivalently:

```python
cfg.source.extinction_mode = "genstars_map"
cfg.source.extinction_law = 1
cfg.source.extinction_map = 1
cfg.source.ejk_scale = 1.15
cfg.source.ejk_offset = 0.05
cfg.source.dm_rc = 0.0
```

`extinction_map = 1` uses the checked-in low-resolution
`EJK_G12_S20_LR.dat` map with subgrid values where available. `extinction_map =
2` uses only the 0.025 deg cell mean. The public genstars high-resolution map
is not bundled here; if a local `EJK_G12_S20.dat` is available,
`extinction_map = 0` can resolve it through the normal input-file search path.
When `dm_rc = 0`, genulens uses the same Nataf-style red-clump distance-modulus
formula used by the existing simulator path.

You can inspect the map and the law conversion directly:

```python
emap = genulens.GenstarsExtinctionMap.load_default(extinction_map=1)
sample = emap.lookup(l=1.0, b=-3.9)
ref = genulens.genstars_reference_extinction(1.0, -3.9, sample.ejk, 1)
print(sample.ejk, ref.Imag)
```

This mode supports `Vmag`, `Imag`, `Jmag_2mass`, `Hmag_2mass`,
`Ksmag_2mass`, and the Roman bands `F087mag`, `F146mag`, and `F213mag`.
`evi_rc` is retained for the legacy `V-I` source-selection path. For
isochrone source cuts, prefer direct band extinctions or
`extinction_mode = "genstars"` / `"genstars_map"`.

Age and metallicity are represented by a conservative discrete source-population
prior. The simulator marginalizes source-selection probability over the
component's `(logAge, [M/H])` prior points, then draws the concrete source star
from the same selected mixture. `IsochroneGrid` currently uses nearest-grid
sequence lookup; it does not interpolate in age or metallicity.

The older `cfg.forward_source` object still exists as an internal compatibility
layer, but new Python examples should use `cfg.source`.

### Isochrone Systematics

The default source-property tables are solar-scaled PARSEC/CMD tables. The
active library can be selected by family and abundance:

```python
cfg.source.mode = "isochrone"
cfg.source.photometry = "roman"  # "roman" or "prime"
cfg.source.isochrone_family = "parsec"  # "parsec" or "mist"
cfg.source.isochrone_abundance = "solar_scaled"
```

The built-in resolver returns the normalized table path expected by genulens:

```python
spec = genulens.IsochroneLibrarySpec()
spec.family = "mist"
spec.photometry = "roman"
spec.abundance = "alpha_enhanced"
spec.alpha_fe = 0.4
path = genulens.default_isochrone_table_path(spec)
```

The checked-in normalized external tables include the PARSEC solar-scaled
baseline plus MIST solar-scaled and alpha-enhanced systematics tables. MIST
currently covers all source components (`0-10`). PARSEC alpha-enhanced tables
are not available; selecting PARSEC with `alpha_enhanced` raises an explicit
error. Use MIST for alpha systematics. Raw MIST archives are intentionally not
required in the repository; recreate them with:

```bash
python tools/source_photometry/fetch_isochrone_assets.py
python tools/source_photometry/build_external_isochrone_aggregates.py
```

`tools/source_photometry/normalize_isochrone_table.py` remains available for
one-off conversion of additional upstream tables.

For bulge-systematics tests, isochrone mode can mix a second normalized
isochrone table, for example an alpha-enhanced table converted to the same
schema and photometric bands:

```python
cfg.source.mode = "isochrone"
cfg.source.photometry = "roman"
cfg.source.isochrone_model = "alpha_mixture"
cfg.source.secondary_isochrone_family = "mist"
cfg.source.secondary_isochrone_abundance = "alpha_enhanced"
cfg.source.secondary_isochrone_alpha_fe = 0.4
cfg.source.alpha_enhanced_fraction = 0.5
```

The primary table remains the default solar-scaled table unless
`isochrone_table_path` is set. `alpha_enhanced_fraction` is interpreted as the
global prior fraction of the secondary table in the source-population model.
If source-selection cuts are active, the mixture is folded into the selected
source-density grid, so both `D_S` sampling and emitted source properties use
the same isochrone mixture. Building this selected source-density grid can take
tens of seconds for MIST tables; genulens prints progress messages to `stderr`
while the grid is being built. Selection probabilities are cached in memory
within the current process, so repeated simulations with the same isochrone,
IMF, extinction, mass range, and selection cuts reuse the expensive part of the
grid build. The secondary table must have the same normalized columns and
magnitude bands as the active primary table. Use
`secondary_isochrone_table_path` or the legacy alias `alpha_enhanced_table_path`
to point at a specific converted table.

For component-specific experiments, override the mixture fraction for selected
Galactic components:

```python
# Mix alpha-enhanced tables only for bar and NSD sources.
cfg.source.isochrone_model = "alpha_mixture"
cfg.source.alpha_enhanced_fraction = 0.0
cfg.source.alpha_enhanced_components = [8, 9]
cfg.source.alpha_enhanced_component_fractions = [1.0, 1.0]
```

Component indices follow the event output convention: `0-6` thin disk, `7`
thick disk, `8` bar, `9` NSD, and `10` stellar halo. Halo source properties
currently use the thick-disk isochrone proxy internally, so component `10` maps
to the thick sequence for source-property lookup.

Scientific caveats:

- `alpha_enhanced_fraction` is a systematics knob, not a chemical-evolution
  model.
- MIST and PARSEC differ in stellar physics, bolometric corrections,
  passbands, and zero-points; differences should not be interpreted as only
  `[alpha/Fe]` effects.

The source-isochrone systematics example defaults to a fast annotation-only
comparison. Add `--use-selection` to fold a band cut into the source-distance
prior:

```bash
PYTHONPATH=build python examples/source_isochrone_systematics.py --use-selection
```

For standalone source-population tests, use `ForwardSourceGenerator` directly:

```python
generator = genulens.ForwardSourceGenerator()
samples = generator.sample_many(n=1000, component=8, distance_pc=8000.0, seed=3)
source_df = pd.DataFrame(samples.to_numpy(), columns=samples.columns)
```

`ForwardSourceResult` uses the same source labels as the appended simulation
columns, with `iS` and `D_S` prepended for the sampled source component and
distance.

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
