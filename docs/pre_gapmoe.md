# `pre_gapmoe` helper tools

The repository includes four preprocessing helper tools for
[`gapmoe`](https://github.com/NunotaKansuke/gapmoe):

- `calc_rho_profile`
- `calc_mass_dist`
- `calc_murel_dist`
- `generate_flow_samples`

They are built from `src/genulens/tools/pre_gapmoe/` and linked against the
shared `genulens_core` library.

## Current status

These tools are still implemented as command-line programs internally, but the
Python extension exposes a small `genulens.pre_gapmoe` API that runs the bundled
helpers and returns table objects with `columns`, `shape`, `command`, `stdout`,
and `to_numpy()`.

This removes the need for downstream packages to find a separate genulens
checkout when genulens is installed from a wheel. A future implementation may
move these helpers fully into library-level C++ result types, but the public
Python call shape is intended to remain stable.

## Build

From the repository root:

```bash
make pre_gapmoe
```

or through CMake:

```bash
cmake --build build --target calc_rho_profile calc_mass_dist calc_murel_dist generate_flow_samples
```

The CMake build places the tools under:

```text
build/pre_gapmoe/
```

Installed wheels also install the helpers as scripts:

```text
calc_rho_profile
calc_mass_dist
calc_murel_dist
generate_flow_samples
```

## Usage

Each tool supports `--help`:

```bash
./build/pre_gapmoe/calc_rho_profile --help
./build/pre_gapmoe/calc_mass_dist --help
./build/pre_gapmoe/calc_murel_dist --help
./build/pre_gapmoe/generate_flow_samples --help
```

Examples:

```bash
./build/pre_gapmoe/calc_rho_profile \
  l 1.0 b -3.9 Dmin 100 Dmax 16000 Dstep 100 \
  AIrc 1.5 EVIrc 1.2 Isrange 14 21

./build/pre_gapmoe/calc_mass_dist Nmass 1000

./build/pre_gapmoe/calc_murel_dist GRID 0 Dl 4000 Ds 8000

# Full event rows for source-group conditional flow training. The source is
# fixed, while lens components retain the normal Galactic mixture.
./build/pre_gapmoe/generate_flow_samples \
  SOURCEGROUP 4 l 1.0 b -3.9 Nsimu 100000 seed 42

./build/pre_gapmoe/generate_flow_samples \
  SOURCEGROUP 3 NSD 1 l 0.0 b 0.0 Nsimu 100000 seed 43
```

`generate_flow_samples` accepts `SOURCEGROUP=3` for NSD or `SOURCEGROUP=4`
for stellar-halo sources. It emits the standard 17-column `VERBOSITY=3`
weighted-event schema, so downstream flow-training parsers can combine these
rows with ordinary genulens output. Its `pi_E`, `pi_EN`, `pi_EE`, `theta_E`,
and `t_E` columns retain enough information to derive the geocentric
North/East relative proper-motion components used by the gapmoe Flow target.
`NSD 1` is supplied automatically for an NSD run when the caller does not
choose an NSD model explicitly. A source group may have zero density at a
given sightline; that invocation fails rather than producing an invalid
training sample.

## Python API

The same helpers are available from Python:

```python
import genulens

mass = genulens.pre_gapmoe.mass_distribution()
rho = genulens.pre_gapmoe.rho_profile(
    l=1.0,
    b=-3.9,
    Dmin=100,
    Dmax=16000,
    Dstep=100,
    SOURCE=1,
)
murel = genulens.pre_gapmoe.murel_distribution(
    l=1.0,
    b=-3.9,
    GRID=1,
    DLmin=0,
    DLmax=12000,
    DLstep=500,
    DSmin=0,
    DSmax=16000,
    DSstep=500,
)
flow_rows = genulens.pre_gapmoe.flow_samples(
    SOURCEGROUP=4,
    l=1.0,
    b=-3.9,
    Nsimu=100000,
    seed=42,
    REMNANT=0,
    BINARY=0,
)

rho_array = rho.to_numpy()
print(rho.columns)
```

Keyword arguments use the same option names as the command-line tools. Sequence
values are expanded, so range options can be written as `Isrange=(14, 21)`.

For development builds, the API searches for helpers in the installed script
directory, `./pre_gapmoe/`, and `./build/pre_gapmoe/`. Pass
`executable_dir="/path/to/pre_gapmoe"` to override discovery. It also passes the
installed `share/genulens/input_files` directory to the helper process; pass
`input_dir="/path/to/input_files"` to override that discovery.
