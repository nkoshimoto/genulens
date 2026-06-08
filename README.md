[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-yellow.svg)](http://www.gnu.org/licenses/gpl-3.0)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4784948.svg)](https://doi.org/10.5281/zenodo.4784948)
[![arXiv](http://img.shields.io/badge/arXiv-2104.03306-orange.svg?style=flat)](https://arxiv.org/abs/2104.03306)

# genulens

`genulens` ("generate microlensing") simulates microlensing events with the
Galactic model of [Koshimoto, Baba & Bennett (2021), ApJ, 917, 78](https://ui.adsabs.harvard.edu/abs/2021ApJ...917...78K/abstract).
The model is optimized for bulge sightlines and is most appropriate around
`|l| < 10 deg` and `|b| < 7 deg`.

The v2 branch keeps the historical command-line simulator and adds a refactored
C++ core, direct Python bindings, source-forward isochrone workflows,
genstars-style extinction-map support, custom Python likelihoods, and
optical-depth/event-rate summary APIs.

For the full guide, start from [`docs/`](docs/). GitHub displays
[`docs/README.md`](docs/README.md) automatically when opening that directory.

## Installation

`genulens` requires a C++ compiler, CMake, and GSL.

Check that GSL is visible:

```bash
gsl-config --libs
```

Clone and build:

```bash
git clone https://github.com/nkoshimoto/genulens.git
cd genulens
make
```

If GSL is installed in a non-standard prefix:

```bash
GSL_ROOT=/path/to/gsl make
```

The same GSL shared libraries must be visible at runtime. On Linux with a
non-standard GSL prefix, set for example:

```bash
export LD_LIBRARY_PATH=/path/to/gsl/lib:$LD_LIBRARY_PATH
```

Build the Python extension:

```bash
make python
PYTHONPATH=build python -c "import genulens; print(genulens.__file__)"
```

The Python extension can also be built with:

```bash
pip install .
```

Editable installs are supported in environments with `scikit-build-core` and
`pybind11`:

```bash
pip install -e .
```

Other build targets:

```bash
make pre_gapmoe
make test
make clean
```

## CLI and Python

The historical CLI remains available:

```bash
./genulens
```

The Python API calls the same C++ simulation core directly. It does not run
`./genulens` as a subprocess and does not parse CLI stdout.

```python
import pandas as pd
import genulens

cfg = genulens.Config(l=1.0, b=-3.9, n_simu=20_000, seed=42)
result = genulens.simulate(cfg)
df = pd.DataFrame(result.to_numpy(), columns=result.columns)
```

See [docs/quickstart.md](docs/quickstart.md) and
[docs/python_api.md](docs/python_api.md) for details.

## Documentation

- [Documentation index](docs/README.md)
- [Quick start](docs/quickstart.md)
- [Python API](docs/python_api.md)
- [Source-forward mode](docs/source_forward.md)
- [Extinction](docs/extinction.md)
- [Optical-depth and event-rate maps](docs/rate_maps.md)
- [Isochrone systematics](docs/isochrone_systematics.md)
- [Architecture notes](docs/architecture.md)
- [pre_gapmoe helper tools](docs/pre_gapmoe.md)

Notebook examples:

- [`examples/python_binding.ipynb`](examples/python_binding.ipynb)
- [`examples/source_isochrone_systematics.ipynb`](examples/source_isochrone_systematics.ipynb)
- [`examples/rate_summary_map.ipynb`](examples/rate_summary_map.ipynb)

The original command-line usage guide remains available as [Usage.pdf](Usage.pdf).

## Citation

Please cite Koshimoto, Baba & Bennett (2021) and
[Koshimoto & Ranc (2021), Zenodo.4784948](http://doi.org/10.5281/zenodo.4784948)
if you use this code in your research.

A separate star simulator, [`genstars`](https://github.com/nkoshimoto/genstars),
is also available.

The copyright of the included supplementary code `option.cpp` belongs to
Ian A. Bond and Takahiro Sumi.

## Release History

- v2: refactored C++ core, direct Python API, source-forward isochrone support,
  extinction-map support, custom Python likelihoods, and rate-summary APIs.
- v1.2, June-July 2022: importance sampling, NSD component, updated Galactic
  Center position, revised usage documentation, and related `genstars` release.
- v1.1, June 2021: switched to the GSL random number generator.
- v1.0, May 2021: initial public release.
