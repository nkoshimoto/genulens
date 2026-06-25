# PARSEC/CMD Source Photometry Tables

This directory contains regenerated PARSEC/CMD source-photometry tables with
stellar parameters needed by the source-prior sampler.

## Provenance

Tables were downloaded from the Padova PARSEC/CMD web service:

```text
https://stev.oapd.inaf.it/cgi-bin/cmd
```

The observed input form was CMD 3.9. The download script used `ezpadova` only as
an acquisition helper; `genulens` does not depend on `ezpadova` at runtime.

The exact request settings are recorded in `requests.json`. The important
choices are:

- PARSEC track set: `parsec_CAF09_v1.2S`
- COLIBRI set: `parsec_CAF09_v1.2S_S_LMC_08_web`
- metallicity: `[M/H] = 0` except `thick`, which uses `[M/H] = -0.8`
- atmosphere/magnitude system: `YBCnewVega`
- component ages copied from the existing `genstars` table comments
- thick-disk metallicity selected from compatibility checks against the existing
  `genstars` thick tables

This root import is the `genstars`-compatibility anchor: it uses one
representative metallicity per component to reproduce the current fixed-table
behavior as closely as possible.

The new source-prior metallicity grid is stored in `metallicity_grid/`.

## Photometric Systems

Three raw CMD outputs are stored per component:

- `*_ubvrijhk.dat`: used for Johnson-Cousins `V` and `I`
- `*_twomass.dat`: used for 2MASS `J`, `H`, and `Ks`
- `*_roman2021.dat`: used for Roman `F087`, `F146`, and `F213`

`Roman2021` was selected because it matches the existing `genstars`
`Z087/W146/F213` processed columns more closely than the newer Roman2024 filter
set.

## Directory Layout

```text
raw/
  <component>_<photometric-system>.dat
normalized/
  <component>_prime_parsec.dat
  <component>_roman_parsec.dat
```

The normalized PRIME-style tables combine `V/I` from `ubvrijhk` and `J/H/Ks`
from `twomass`.

The normalized Roman-style tables combine `J/H/Ks` from `twomass` and
`F087/F146/F213` from `roman2021`.

Normalized columns include:

```text
Mini  Mass  radius_rsun  Teff_K  logg  ...
```

The `metallicity_grid/` tables retain explicit component/age/metallicity
metadata in normalized outputs:

```text
component  logAge  MH  Zini  Mini  Mass  radius_rsun  Teff_K  logg  ...
```

`radius_rsun` is derived from CMD `logL` and `logTe` using
`T_sun = 5772 K`:

```text
R/Rsun = sqrt(10^logL) * (T_sun / 10^logTe)^2
```

## Compatibility Check

`genstars_compatibility.md` compares the normalized tables against the copied
`genstars` processed tables in `../prime_hybrid` and `../roman_isochrone`.

The comparison intentionally focuses on unevolved, PARSEC-dominated mass ranges:

- PRIME/hybrid: `0.55 <= Mini <= 1.0`
- Roman/isochrone: `0.09 <= Mini <= 1.0`

The lower-mass legacy rows are partly empirical/Baraffe, not pure PARSEC. The
evolved high-mass rows are more sensitive to isochrone sampling and historical
table-processing choices.

## Difference from genstars

`genstars` treats age and metallicity as properties baked into one preprocessed
mass-luminosity table per component. It does not sample metallicity at runtime.

The `genulens` source-prior layer keeps that behavior available as a
compatibility mode, and adds a new mode where metallicity is sampled or
conditioned explicitly and the stellar properties are read from
`metallicity_grid/`.

## Citations

When these tables are used in published work, cite the PARSEC/CMD service and
the relevant model/filter papers listed by the CMD output and service page,
including at least:

- Bressan et al. 2012, MNRAS, 427, 127
- Chen et al. 2014, MNRAS, 444, 2525
- Chen et al. 2015, MNRAS, 452, 1068
- Tang et al. 2014, MNRAS, 445, 4287
- Marigo et al. 2017, ApJ, 835, 77
- Pastorelli et al. 2019/2020 for COLIBRI TP-AGB ingredients
- Chen et al. 2019 for the YBC bolometric-correction database

See `CITATION.md` for the reference list recorded with this import.
