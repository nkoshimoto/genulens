# Source photometry tables

These tables are copied from `genstars/input_files` so that `genulens` can build
source-star priors without depending on a sibling checkout of `genstars`.

The tables are grouped by photometry relation, not only by output band set.

The production source-prior work should use `parsec_cmd/` as the canonical
external table family because it includes `Teff_K` and `logg`. The copied
`prime_hybrid/` and `roman_isochrone/` tables are retained as compatibility
references for the current `genstars` behavior.

The intended extension path is:

```text
genstars-compatible:
  component -> one fixed source-photometry table

new source-prior:
  component -> metallicity prior -> PARSEC/CMD table grid
```

This keeps the current `genstars` behavior reproducible while making
metallicity an explicit source-prior dimension.

## `parsec_cmd/`

Regenerated PARSEC/CMD tables with raw downloads, normalized source-prior
tables, checksums, request settings, and compatibility checks.

The root `parsec_cmd/normalized` import is a compatibility anchor with one
representative metallicity per component. `parsec_cmd/metallicity_grid/` stores
the new-mode finite `[M/H]` grid with `[M/H]`/`Zini` explicit in the normalized
tables.

Normalized table families:

```text
<component>_prime_parsec.dat
<component>_roman_parsec.dat
```

See `parsec_cmd/README.md` for the acquisition settings and citation notes.

## Isochrone Library Registry

`isochrone_libraries.json` records the normalized table paths understood by
the source-forward resolver. The code can select `parsec` or `mist`
families. PARSEC is solar-scaled only in this repository; PARSEC
`alpha_enhanced` requests intentionally raise an error so users do not mistake a
missing catalog for a supported model. Use MIST for alpha-enhanced systematics.
The checked-in normalized external tables include the PARSEC solar-scaled
baseline plus MIST comparison tables.

Use `tools/source_photometry/normalize_isochrone_table.py` for one-off
conversion of external model tables. It accepts explicit column mappings, so it
can ingest PARSEC/CMD or MIST tables without hard-coding every upstream header
variant.

Raw external downloads for MIST are tracked by manifest files in this
directory. See `.note/isochrone_download_status.md` for the current downloaded
asset list and normalization notes.

Large raw archives are not intended to be committed. Recreate them and rebuild
the normalized aggregates with:

```bash
python tools/source_photometry/fetch_isochrone_assets.py
python tools/source_photometry/build_external_isochrone_aggregates.py
```

The generated MIST aggregates cover component indices `0-10`.

## `prime_hybrid/`

Files:

```text
isoemp_thin1.dat ... isoemp_thin7.dat
isoemp_thick2.dat
isoemp_bar.dat
isoemp_NSD.dat
```

Columns:

```text
Mini  MPD  Rad  MV_j  MI_c  MJ_2M  MH_2M  MK_2M
```

This is the `genstars` non-Roman / PRIME-style relation. It uses a hybrid
mass-luminosity relation: empirical at low mass and isochrone-based at higher
mass.

## `roman_isochrone/`

Files:

```text
isochrone_thin1.dat ... isochrone_thin7.dat
isochrone_thick.dat
isochrone_bar.dat
isochrone_NSD.dat
```

Columns:

```text
Mini  MPD  Rad  MJ_2M  MH_2M  MK_2M  MZ087  MW146  MF213
```

This is the `genstars` Roman-style relation. It uses isochrone-based values for
the full mass range.

## Notes

- Common `J/H/Ks` bands are intentionally kept in separate relation directories.
  They can differ because the underlying mass-luminosity relation differs.
- The current copied `genstars` tables do not include a stellar halo component.
  A source-prior implementation must either add halo tables, fall back to an
  existing relation, or mark halo photometry unavailable.
- The current copied tables include radius and magnitudes, but not `Teff`,
  `logg`, or microturbulent velocity.
