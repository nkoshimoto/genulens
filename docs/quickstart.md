# Quick start

This page shows the shortest path from a source checkout to a Python
simulation.

## Install

Released Linux x86_64 and macOS arm64 wheels can be installed with:

```bash
pip install genulens
```

The wheels include the compiled extension, command-line executable, bundled
input tables, and the GSL shared libraries needed by the extension. macOS
x86_64 may fall back to a source build and therefore requires a system GSL
installation.

Source builds, including `pip install --no-binary genulens genulens` and
`pip install git+https://github.com/nkoshimoto/genulens.git`, require a system
GSL installation.

## Build From Source

To build from a source checkout, `genulens` requires a C++ compiler, CMake, and
GSL.

Check that GSL is visible:

```bash
gsl-config --libs
```

Build the command-line program and Python extension:

```bash
make
make python
```

If GSL is installed in a non-standard prefix:

```bash
GSL_ROOT=/path/to/gsl make
GSL_ROOT=/path/to/gsl make python
```

The same GSL shared libraries must be visible at runtime. On Linux with a
non-standard GSL prefix, set for example:

```bash
export LD_LIBRARY_PATH=/path/to/gsl/lib:$LD_LIBRARY_PATH
```

From a source checkout, import the extension with `PYTHONPATH=build`:

```bash
PYTHONPATH=build python -c "import genulens; print(genulens.__file__)"
```

You can also build/install the Python extension from source with:

```bash
pip install .
```

For a non-standard GSL prefix, pass `GSL_ROOT` to pip:

```bash
GSL_ROOT=/path/to/gsl pip install .
```

For source installs from a non-standard GSL prefix, the built extension records
the linked GSL path in its install RPATH. If your platform strips or ignores
that RPATH, set `LD_LIBRARY_PATH` or the platform equivalent at runtime.

Editable installs are supported in environments with `scikit-build-core` and
`pybind11`:

```bash
pip install -e .
```

## Run a simulation

```python
import pandas as pd
import genulens

cfg = genulens.Config(l=1.0, b=-3.9, n_simu=20_000, seed=42)
result = genulens.simulate(cfg)
df = pd.DataFrame(result.to_numpy(), columns=result.columns)
```

`result.columns` gives the column order. Use the `wtj` column for histograms,
weighted quantiles, and population fractions.

## Add observation constraints

```python
cfg.observation.tE_obs = 54.5
cfg.observation.tE_err = 5.0
cfg.observation.thetaE_obs = 0.55
cfg.observation.thetaE_err = 0.15

result = genulens.simulate(cfg)
```

## Add a custom likelihood

```python
def likelihood(event):
    if event.D_L >= event.D_S:
        return 0.0
    return 1.0 if 10.0 < event.t_E < 100.0 else 0.0

result = genulens.simulate(cfg, likelihood=likelihood)
```

The callable is evaluated inside the Monte Carlo event loop. It receives an
`Event` with labels matching the Python result table.

## Next steps

- [Python API](python_api.md) for result labels and configuration fields.
- [Source-forward mode](source_forward.md) for source stellar properties.
- [Optical-depth and event-rate maps](rate_maps.md) for `tau` and event-rate
  summaries.
