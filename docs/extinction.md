# Extinction

Genulens supports manual band extinction values and genstars-style extinction
laws. Apparent source selections require an extinction model for the selected
band.

## Manual extinction

Manual mode uses direct band extinctions at the reference red-clump distance:

```python
cfg.source.extinction_mode = "manual"
cfg.source.dm_rc = 14.5
cfg.source.av_rc = 2.0   # Vmag
cfg.source.ai_rc = 1.4   # Imag
cfg.source.aj_rc = 0.5   # Jmag / Jmag_2mass
cfg.source.ah_rc = 0.2   # Hmag / Hmag_2mass
cfg.source.ak_rc = 0.1   # Kmag / Ksmag_2mass / K_L
```

Use manual mode when you want full control over the band extinctions used by a
single sightline.

## genstars law from one E(J-Ks)

Automatic extinction follows the genstars law from a single `E(J-Ks)` value:

```python
cfg.use_genstars_extinction(
    dm_rc=14.5,
    ejk_rc=1.0,
    extinction_law=1,
)
```

This is equivalent to setting:

```python
cfg.source.extinction_mode = "genstars"
cfg.source.extinction_law = 1
cfg.source.ejk_rc = 1.0
cfg.source.dm_rc = 14.5
```

Use this mode for a single sightline when you already know the appropriate
red-clump `E(J-Ks)`.

## genstars extinction map

For maps, use the bundled genstars-style extinction map:

```python
cfg.use_genstars_extinction_map(
    extinction_law=1,
    extinction_map=1,
)
```

`extinction_map = 1` uses the checked-in low-resolution
`input_files/EJK_G12_S20_LR.dat` map with subgrid values where available.
`extinction_map = 2` uses only the 0.025 deg cell mean. The public genstars
high-resolution map is not bundled here; if a local `EJK_G12_S20.dat` is
available, `extinction_map = 0` can resolve it through the normal input-file
search path.

When `dm_rc = 0`, genulens uses the same Nataf-style red-clump
distance-modulus formula used by the existing simulator path.

## Map adjustment knobs

You can adjust the map without editing the input file:

```python
cfg.use_genstars_extinction_map(
    extinction_law=1,
    extinction_map=1,
    ejk_scale=1.15,
    ejk_offset=0.05,
)
```

The effective value is:

```text
E(J-Ks)_eff = E(J-Ks)_map * ejk_scale + ejk_offset
```

Use this for sensitivity checks or field-level calibration experiments.

## Inspecting the map

```python
emap = genulens.GenstarsExtinctionMap.load_default(extinction_map=1)
sample = emap.lookup(l=1.0, b=-3.9)
ref = genulens.genstars_reference_extinction(1.0, -3.9, sample.ejk, 1)

print(sample.ejk)
print(ref.Imag, ref.Hmag, ref.Kmag)
```

Supported source-selection bands include `Vmag`, `Imag`, `Jmag_2mass`,
`Hmag_2mass`, `Ksmag_2mass`, and Roman bands `F087mag`, `F146mag`, and
`F213mag`.

`evi_rc` is retained for the legacy `V-I` source-selection path. For isochrone
source cuts, prefer direct band extinctions or `extinction_mode = "genstars"` /
`"genstars_map"`.

## Map workflows

For realistic event-rate maps, use one config per sightline and call
`use_genstars_extinction_map` for each config. Do not reuse one fixed
red-clump extinction value over a large `l,b` grid.
