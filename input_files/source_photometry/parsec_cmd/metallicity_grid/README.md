# PARSEC/CMD Metallicity Grid

This directory contains the new-mode PARSEC/CMD source-photometry grid. Unlike
the compatibility anchor in `../normalized`, this table set keeps metallicity as
an explicit grid dimension.

## Layout

```text
requests.json
raw/<component>/<mh-token>_<photometric-system>.dat
raw/grid/<age-token>/<mh-token>_<photometric-system>.dat
normalized/<component>/<component>_<mh-token>_prime_parsec.dat
normalized/<component>/<component>_<mh-token>_roman_parsec.dat
normalized/<component>/<age-token>/<component>_<age-token>_<mh-token>_prime_parsec.dat
normalized/<component>/<age-token>/<component>_<age-token>_<mh-token>_roman_parsec.dat
normalized/all_prime_parsec.dat
normalized/all_roman_parsec.dat
checksums.sha256
```

The checked-in import currently contains the multi-age, multi-metallicity grid:

```text
113 component-age-metallicity grid points
77 unique age-metallicity raw grid points
231 shared raw CMD files
226 per-grid-point normalized tables
2 aggregate normalized tables
```

`mh-token` encodes `[M/H]`, for example:

```text
m0p25  -> [M/H] = -0.25
p0p25  -> [M/H] = +0.25
```

`age-token` encodes `logAge`, for example:

```text
age9p95424 -> logAge = 9.95424
```

Age-token raw subdirectories are used only when `requests.json` supplies
`logAge_grid`. The raw CMD files are shared across components because PARSEC
isochrones are functions of `logAge`, `[M/H]`, and photometric system; component
metadata is added only when normalized component tables are written. Single-age
compatibility configs keep the old component-local raw paths.

## Grid

The metallicity grid is intentionally modest:

```text
thin disk:  [-0.5, -0.25, 0.0, +0.25]
thick disk: [-1.2, -1.0, -0.8, -0.6, -0.4]
bar/NSD:    [-0.5, -0.25, 0.0, +0.25]
```

It is meant to represent metallicity uncertainty while keeping the first
implementation auditable. A sampler can start with nearest-neighbor metallicity
selection and later add interpolation where it is well behaved.

The age grid is also finite and intentionally hard-discrete:

```text
thin disk:  three representative points per thin age bin
thick disk: old single population at 12 Gyr
bar:        8, 9, 10 Gyr around the legacy Gaussian SFR center
NSD:        6, 7, 8 Gyr around the legacy Gaussian SFR center
```

This mirrors the age handling in the Koshimoto/genstars luminosity-function
construction without introducing age interpolation. The sampler marginalizes
over these grid points; `IsochroneGrid` then uses nearest-grid sequence lookup.

## Normalized Columns

The normalized tables retain metadata needed by the source-prior sampler:

```text
component  component_index  family  logAge  MH  Zini  Mini  Mass  radius_rsun  Teff_K  logg  ...
```

`radius_rsun` is derived from CMD `logL` and `logTe` using `T_sun = 5772 K`.

## Difference from genstars

`genstars` uses one fixed source-photometry table per Galactic component, with
the age distribution folded into the precomputed luminosity functions. This grid
keeps the same component organization but exposes both `logAge` and `[M/H]` as
explicit dimensions. This makes it a direct extension of the `genstars` table
model, not a separate stellar-population prescription.
