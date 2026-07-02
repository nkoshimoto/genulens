# Releasing Python Wheels

This project uses `cibuildwheel` and PyPI trusted publishing to build and
publish binary wheels.

## Supported Wheels

The release workflow builds wheels for:

- Linux x86_64
- macOS arm64, deployment target 14.0
- CPython 3.9 through 3.13

macOS x86_64 and Windows wheels are not built in the first PyPI release
series. Those platforms may fall back to source builds requiring a system GSL
installation.

The wheels bundle:

- the `genulens` Python extension
- the `genulens` command-line executable
- required `input_files`
- GSL shared libraries used by the extension

Source builds do not bundle GSL. `pip install --no-binary genulens genulens`,
GitHub installs, and editable installs require a system GSL installation. Use
`GSL_ROOT=/path/to/gsl` for non-standard prefixes.

## One-Time PyPI Setup

Create the `genulens` project on PyPI/TestPyPI and configure trusted
publishing for this GitHub repository.

Use these GitHub environments:

- `testpypi`
- `pypi`

For PyPI, allow releases from:

- owner: `nkoshimoto`
- repository: `genulens`
- workflow: `wheels.yml`
- environment: `pypi`

For TestPyPI, use the same repository and workflow with environment
`testpypi`.

## TestPyPI

Run the workflow manually:

1. Open GitHub Actions.
2. Select `Build and publish Python wheels`.
3. Run workflow with `publish_to_testpypi = true`.
4. Install from TestPyPI in a clean environment.

```bash
python -m venv /tmp/genulens-testpypi
/tmp/genulens-testpypi/bin/python -m pip install --upgrade pip
/tmp/genulens-testpypi/bin/python -m pip install \
  --index-url https://test.pypi.org/simple/ \
  --extra-index-url https://pypi.org/simple/ \
  genulens
/tmp/genulens-testpypi/bin/python - <<'PY'
import genulens
cfg = genulens.Config(l=1.0, b=-3.9, n_simu=5, seed=1)
result = genulens.simulate(cfg)
print(result.to_numpy().shape)
PY
```

The workflow also validates the built artifacts before upload.

For Linux wheels it runs:

```bash
unzip -l wheelhouse/genulens-*.whl | grep -E 'libgsl|libgslcblas'
python -m auditwheel show wheelhouse/genulens-*.whl
docker run --rm -v "$PWD/wheelhouse:/wheelhouse:ro" python:3.11-slim ...
```

The `python:3.11-slim` smoke test installs a `cp311` wheel in a container that
does not have GSL installed and runs `import genulens` plus a minimal
`genulens.simulate(...)` call.

For macOS wheels it runs:

```bash
unzip -l wheelhouse/genulens-*.whl | grep -E 'libgsl|libgslcblas|dylib'
delocate-listdeps wheelhouse/genulens-*.whl
```

## PyPI

After TestPyPI works, publish a GitHub Release from a tag whose version matches
`pyproject.toml`.

For example:

```bash
git tag -a v2.0.0a2 -m "genulens v2.0.0 alpha 2"
git push upstream v2.0.0a2
```

Then create a GitHub Release for `v2.0.0a2`. Publishing the release triggers
the PyPI job.

## Local Source Install

Source installs still require a visible GSL installation:

```bash
GSL_ROOT=/path/to/gsl pip install .
```

This is separate from released wheels, which bundle the linked GSL shared
libraries.

## License Notices

`genulens` is distributed under GPL-3.0. Binary wheels may bundle GNU
Scientific Library shared libraries. The wheel metadata includes
`third_party_licenses/GSL.md` with the upstream GSL notice and source location.
