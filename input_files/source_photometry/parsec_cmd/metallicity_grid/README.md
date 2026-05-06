# PARSEC/CMD Metallicity Grid

This directory contains the new-mode PARSEC/CMD source-photometry grid. Unlike
the compatibility anchor in `../normalized`, this table set keeps metallicity as
an explicit grid dimension.

## Layout

```text
requests.json
raw/<component>/<mh-token>_<photometric-system>.dat
normalized/<component>/<component>_<mh-token>_prime_parsec.dat
normalized/<component>/<component>_<mh-token>_roman_parsec.dat
normalized/all_prime_parsec.dat
normalized/all_roman_parsec.dat
checksums.sha256
```

Current import contents:

```text
41 component-metallicity grid points
123 raw CMD files
82 per-grid-point normalized tables
2 aggregate normalized tables
```

`mh-token` encodes `[M/H]`, for example:

```text
m0p25  -> [M/H] = -0.25
p0p25  -> [M/H] = +0.25
```

## Grid

The initial grid is intentionally modest:

```text
thin disk:  [-0.5, -0.25, 0.0, +0.25]
thick disk: [-1.2, -1.0, -0.8, -0.6, -0.4]
bar/NSD:    [-0.5, -0.25, 0.0, +0.25]
```

It is meant to represent metallicity uncertainty while keeping the first
implementation auditable. A sampler can start with nearest-neighbor metallicity
selection and later add interpolation where it is well behaved.

## Normalized Columns

The normalized tables retain metadata needed by the source-prior sampler:

```text
component  component_index  family  logAge  MH  Zini  Mini  Mass  radius_rsun  Teff_K  logg  ...
```

`radius_rsun` is derived from CMD `logL` and `logTe` using `T_sun = 5772 K`.

## Difference from genstars

`genstars` uses one fixed source-photometry table per Galactic component. This
grid keeps that component/age organization but adds `[M/H]` as an explicit
dimension. This makes it a direct extension of the `genstars` table model, not a
separate stellar-population prescription.
