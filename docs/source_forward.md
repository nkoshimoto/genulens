# Source-forward mode

Source-forward mode is the v2 path for attaching physical source-star
properties to simulated microlensing events.

The core idea is:

1. Sample the source component and source distance from the Galactic model.
2. Use the component-dependent source-population prior for age and metallicity.
3. Draw an initial source mass from the IMF.
4. Query an isochrone table for radius, temperature, surface gravity, current
   mass, and absolute magnitudes.
5. Optionally fold a magnitude selection into the source-distance prior.

## Classic versus isochrone source modes

Classic source mode keeps the historical luminosity-function and
color-magnitude-function source selection:

```python
cfg = genulens.Config(l=1.0, b=-3.9, n_simu=20_000, seed=42)
cfg.use_classic_source(i_min=14.0, i_max=21.0)
```

Isochrone mode uses the shared isochrone/IMF source prior:

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
```

The convenience method sets:

- `cfg.source.mode = "isochrone"`
- `cfg.source.band`
- `cfg.source.photometry`
- `cfg.source.min_magnitude`
- `cfg.source.max_magnitude`
- `cfg.source.apparent_magnitude`
- `cfg.source.min_initial_mass_msun`
- `cfg.source.max_initial_mass_msun`
- `cfg.source.use_magnitude_selection = 1`

Isochrone mode is a new source-population model, not a byte-for-byte
replacement for the historical luminosity-function selection. It is the right
mode for source physical properties and multi-band forward-prior experiments,
but absolute source densities and event rates depend on the chosen isochrone
family, age/metallicity prior, abundance assumptions, IMF, passbands, and
extinction model.

## Output columns

Isochrone mode appends source-property columns after the event columns:

- `logage_S`
- `MH_S`
- `M_S_ini`
- `M_S`
- `R_S`
- `teff_S`
- `logg_S`
- `theta_S`
- `M_<band>_S`

The exact magnitude columns depend on the active photometry table. For PRIME
tables, examples include `M_Vmag_S`, `M_Imag_S`, `M_Jmag_2mass_S`,
`M_Hmag_2mass_S`, and `M_Ksmag_2mass_S`. Roman tables include Roman band
magnitudes such as `M_F146mag_S`.

`theta_S` is the source angular radius in microarcseconds. Remnants and
unsupported source states can emit `NaN` for stellar-property fields.

## Magnitude selection

With `cfg.source.use_magnitude_selection = 1`, genulens folds the isochrone/IMF
selection probability into the source-distance density grid before event
sampling. This means the sampled `D_S` distribution changes when the selected
band and magnitude range change.

Set:

```python
cfg.source.use_magnitude_selection = 0
```

to append source properties without conditioning the source-distance prior on a
magnitude cut.

Cuts are apparent magnitudes by default. Set:

```python
cfg.source.apparent_magnitude = 0
```

for absolute-magnitude cuts.

Apparent cuts require extinction for the selected band. See
[Extinction](extinction.md).

## Initial mass range

`cfg.source.min_initial_mass_msun` and `cfg.source.max_initial_mass_msun` define
the source-star IMF range used for both:

- selected source-density normalization
- concrete source-star sampling

Treat this range as part of the physical source prior, not only as a proposal
optimization.

## Age and metallicity

Age and metallicity are represented by a conservative discrete
source-population prior. The simulator marginalizes source-selection
probability over the component's `(logAge, [M/H])` prior points, then draws the
concrete source star from the same selected mixture.

`IsochroneGrid` currently uses nearest-grid sequence lookup. It does not
interpolate in age or metallicity.

Component indices follow the event output convention:

- `0-6`: thin disk
- `7`: thick disk
- `8`: bar
- `9`: nuclear stellar disk
- `10`: stellar halo

Halo source properties currently use the thick-disk isochrone proxy internally,
so component `10` maps to the thick sequence for source-property lookup.

## Performance

Building a source-forward selected-density grid can take tens of seconds for
large isochrone tables and map grids. Genulens prints progress messages to
`stderr`, and repeated runs in the same Python process can reuse cached
selection probabilities when the isochrone, IMF, extinction, mass range, and
selection cuts are unchanged.

For notebooks, see
[`examples/rate_summary_map.ipynb`](../examples/rate_summary_map.ipynb) for a
`tqdm`-based progress display that hides detailed native progress logs.

## Related pages

- [Extinction](extinction.md)
- [Isochrone systematics](isochrone_systematics.md)
- [Optical-depth and event-rate maps](rate_maps.md)
