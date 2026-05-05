# pre_gapmoe helper tools

The repository includes three GAPMOE preprocessing helper tools:

- `calc_rho_profile`
- `calc_mass_dist`
- `calc_murel_dist`

They are built from `src/genulens/tools/pre_gapmoe/` and linked against the
shared `genulens_core` library.

## Current status

These tools are currently command-line programs. They write tabular output to
stdout with `printf`. Unlike `genulens.simulate()`, they do not yet have direct
Python APIs that return `numpy.ndarray` objects.

This is intentional for the first Python API milestone. The main simulator has
been moved away from stdout parsing; the preprocessing helpers can be wrapped
later once their library-level result types are designed.

## Build

From the repository root:

```bash
make pre_gapmoe
```

or through CMake:

```bash
cmake --build build --target calc_rho_profile calc_mass_dist calc_murel_dist
```

The CMake build places the tools under:

```text
build/pre_gapmoe/
```

## Usage

Each tool supports `--help`:

```bash
./build/pre_gapmoe/calc_rho_profile --help
./build/pre_gapmoe/calc_mass_dist --help
./build/pre_gapmoe/calc_murel_dist --help
```

Examples:

```bash
./build/pre_gapmoe/calc_rho_profile \
  l 1.0 b -3.9 Dmin 100 Dmax 16000 Dstep 100 \
  AIrc 1.5 EVIrc 1.2 Isrange 14 21

./build/pre_gapmoe/calc_mass_dist Nmass 1000

./build/pre_gapmoe/calc_murel_dist GRID 0 Dl 4000 Ds 8000
```

## Future Python direction

A direct Python API should not parse these stdout tables. The preferred future
shape is:

```python
rho = genulens.pre_gapmoe.rho_profile(config).to_numpy()
mass = genulens.pre_gapmoe.mass_distribution(config).to_numpy()
murel = genulens.pre_gapmoe.murel_distribution(config).to_numpy()
```

That requires extracting the current `main()` implementations into reusable
library functions with typed result containers.
