# Optical-depth and event-rate maps

The Python API can return line-of-sight optical-depth and event-rate summaries
directly from the C++ simulation state. These values are not parsed from CLI
stdout.

## Single sightline

```python
import genulens

cfg = genulens.Config(l=1.0, b=-3.9, n_simu=10_000, seed=42)
cfg.use_classic_source(i_min=14.0, i_max=21.0)
cfg.use_genstars_extinction_map(extinction_law=1, extinction_map=1)
cfg.observation.IL_err = 0.0

summary = genulens.compute_rate_summary(cfg)
print(summary.tau)
print(summary.event_rate_per_star_per_year)
print(summary.event_rate_per_deg2_per_year)
```

`compute_rate_summary` copies the config, enables optical-depth calculation
internally, runs the same event sampler used by `simulate`, and returns a
`RateSummary`.

Fields include:

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

## Multiple sightlines

For maps, build one config per sightline and call
`compute_rate_summaries`:

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

Use one config per sightline when the extinction, source selection, or model
options depend on `l,b`.

## Convenience map API

`compute_rate_map` changes only `cfg.l` and `cfg.b` on a base config:

```python
rate_map = genulens.compute_rate_map(
    l_values=np.linspace(-5.0, 5.0, 21),
    b_values=np.linspace(-6.0, -2.0, 9),
    base_config=cfg,
)
```

This is useful for quick tests. It does not infer a line-of-sight extinction
map. If the base config contains fixed manual extinction values, those values
are reused at every map point.

For realistic visible-band maps, prefer `compute_rate_summaries` with
`use_genstars_extinction_map` or another per-sightline extinction model.

## Classic versus source-forward maps

Classic source selection uses the historical luminosity-function path.

Source-forward isochrone selection folds a band and magnitude cut into the
source-distance density grid. This is the mode to use when the event-rate map
must be tied to source properties or to bands beyond the historical `I`/`V-I`
selection.

The current I-band classic rate path has been checked against the legacy
`genulens_source` summary calculation. H-band maps in the public examples use
the new isochrone source-selection model, so they should be treated as
model-dependent forward predictions rather than legacy-compatible rates.

Do not compare event rates from different source-selection models as if they
were only band differences. For example, comparing classic `I` selection to
isochrone `H` selection mixes source-density normalizations and selection
models. For a cleaner I/H comparison, use isochrone mode for both bands.

## Example notebook

See [`examples/rate_summary_map.ipynb`](../examples/rate_summary_map.ipynb) for
I- and H-band maps, extinction maps, per-star rates, per-square-degree rates,
latitude profiles integrated over a longitude range, and an I-band comparison
with the central-region MOA-II 9 yr empirical fits from Nunota et al. (2025).

The Nunota et al. comparison uses the published `|l| < 5 deg` fits:

```text
tau = 1.75e-6 * exp(0.34 * (3 - |b|))
Gamma = 16.08e-6 * exp(0.44 * (3 - |b|)) star^-1 yr^-1
```

These fits are for `I_S < 21.4 mag`, `t_E < 760 d`, and the measured southern
bulge range. The notebook overlays them only as a diagnostic reference for
`tau` and per-star event rate. It does not compare per-square-degree rates
because that also requires matching the observed source-density normalization.
