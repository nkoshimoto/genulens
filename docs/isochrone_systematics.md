# Isochrone systematics

Source-forward mode can test how source-property predictions change under
different isochrone assumptions.

The current bundled table set supports:

- PARSEC solar-scaled baseline tables.
- MIST solar-scaled tables.
- MIST alpha-enhanced systematics tables.

PARSEC alpha-enhanced tables are not available in the bundled table set.
Selecting PARSEC with `alpha_enhanced` raises an explicit error. Use MIST for
alpha-enhancement systematics.

## Select a library

```python
cfg.source.mode = "isochrone"
cfg.source.photometry = "roman"  # "roman" or "prime"
cfg.source.isochrone_family = "mist"  # "parsec" or "mist"
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

## Alpha-enhanced mixture

For bulge-systematics tests, isochrone mode can mix a second normalized
isochrone table:

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
the same isochrone mixture.

The secondary table must have the same normalized columns and magnitude bands
as the active primary table. Use `secondary_isochrone_table_path` or the legacy
alias `alpha_enhanced_table_path` to point at a specific converted table.

## Component-specific mixtures

Override the mixture fraction for selected Galactic components:

```python
# Mix alpha-enhanced tables only for bar and NSD sources.
cfg.source.isochrone_model = "alpha_mixture"
cfg.source.alpha_enhanced_fraction = 0.0
cfg.source.alpha_enhanced_components = [8, 9]
cfg.source.alpha_enhanced_component_fractions = [1.0, 1.0]
```

Component indices follow the event output convention:

- `0-6`: thin disk
- `7`: thick disk
- `8`: bar
- `9`: nuclear stellar disk
- `10`: stellar halo

Halo source properties currently use the thick-disk isochrone proxy internally.

## Rebuilding normalized external tables

Raw MIST archives are intentionally not required in the repository. Recreate
the normalized external tables with:

```bash
python tools/source_photometry/fetch_isochrone_assets.py
python tools/source_photometry/build_external_isochrone_aggregates.py
```

`tools/source_photometry/normalize_isochrone_table.py` remains available for
one-off conversion of additional upstream tables.

## Scientific caveats

- `alpha_enhanced_fraction` is a systematics knob, not a chemical-evolution
  model.
- MIST and PARSEC differ in stellar physics, bolometric corrections, passbands,
  and zero-points. Differences should not be interpreted as only `[alpha/Fe]`
  effects.
- Cool-star spectra and molecular absorption features can make inferred
  source properties sensitive to the chosen isochrone/atmosphere model.
- Source-forward outputs depend on the isochrone library, abundance
  assumptions, IMF, source-population prior, extinction law, and extinction map.

## Examples

The source-isochrone systematics example defaults to a fast annotation-only
comparison:

```bash
PYTHONPATH=build python examples/source_isochrone_systematics.py
```

Add `--use-selection` to fold a band cut into the source-distance prior:

```bash
PYTHONPATH=build python examples/source_isochrone_systematics.py --use-selection
```

The notebook version is
[`examples/source_isochrone_systematics.ipynb`](../examples/source_isochrone_systematics.ipynb).
