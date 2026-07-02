# genulens v2 documentation

This directory is the user and collaborator guide for the refactored
`genulens` v2 codebase.

Start with the repository README for a high-level overview. Use these pages for
the working details.

## User guides

- [Quick start](quickstart.md): build the project, import the Python extension,
  run a simulation, and convert results to pandas.
- [Python API](python_api.md): `Config`, `SimulationResult`, event labels,
  observation constraints, sampling options, and custom likelihoods.
- [Source-forward mode](source_forward.md): classic source selection versus
  isochrone source sampling, source-property columns, and magnitude cuts.
- [Extinction](extinction.md): manual extinction, genstars-style extinction
  laws, bundled `E(J-Ks)` map support, and map adjustment knobs.
- [Optical-depth and event-rate maps](rate_maps.md): direct rate summary APIs
  and map workflows.
- [Isochrone systematics](isochrone_systematics.md): PARSEC/MIST tables,
  alpha-enhanced MIST experiments, and scientific caveats.
- [`pre_gapmoe` tools](pre_gapmoe.md): build and usage notes for the helper
  command-line tools.

## Developer notes

- [Architecture](architecture.md): current object boundaries, simulation flow,
  transitional pieces, and suggested next refactor steps.
- [Releasing Python wheels](releasing.md): wheel build, TestPyPI, and PyPI
  publication workflow.

## Examples

- [`examples/python_binding.ipynb`](../examples/python_binding.ipynb): basic
  Python simulation, result tables, plotting, and custom likelihoods.
- [`examples/source_isochrone_systematics.ipynb`](../examples/source_isochrone_systematics.ipynb):
  source-property and HR-diagram comparisons across isochrone assumptions.
- [`examples/source_isochrone_systematics.py`](../examples/source_isochrone_systematics.py):
  script version of the isochrone-systematics comparison.
- [`examples/rate_summary_map.ipynb`](../examples/rate_summary_map.ipynb):
  optical-depth, event-rate, and extinction maps for I- and H-band source
  selections.
- [`examples/cli_examples.ipynb`](../examples/cli_examples.ipynb): command-line
  examples for legacy workflows.

## Scientific scope

The Galactic model was developed for bulge microlensing applications and is
most appropriate around `|l| < 10 deg` and `|b| < 7 deg`.

The source-forward and isochrone-systematics workflows are intended for
physically motivated prior construction and systematic-error exploration. They
depend on the selected isochrone library, abundance assumptions, IMF,
source-population prior, extinction law, and extinction map.
