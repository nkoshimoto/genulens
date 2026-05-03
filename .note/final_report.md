# Final Report

## Summary

- Protected the starting state with commit `checkpoint: before major genulens refactor`.
- Created backup branch `backup/pre-major-refactor-20260503-084619`.
- Added a CMake-backed build while preserving root `make` entry points.
- Added `src/genulens` shared-library skeleton covering config, RNG, input path resolution, model services, density, kinematics, mass function, extinction, likelihood, simulation, output, and tool helpers.
- Added pybind11 binding source with `Config`, `simulate()`, Python callable likelihood support, `SimulationResult.columns`, and `SimulationResult.to_numpy()`.
- Added smoke/unit/Python tests and regression comparison against the pre-refactor binary backup.
- Continued object extraction beyond the initial wrapper:
  - `ScientificEngine` owns grouped runtime/scientific state.
  - `RandomEngine` owns GSL RNG via RAII.
  - `ObservationLikelihood` owns Gaussian/uniform observation likelihood dispatch.
  - `math::Interpolation` owns interpolation routines.
  - `math::NewtonCotes` owns integration coefficients.
  - `model::CoordinateTransformer` owns coordinate and position-angle transforms.
  - `model::BrokenPowerLawIMF` owns broken-power-law IMF grid generation.

## Test Results

- `make`: passed.
- `make pre_gapmoe`: passed.
- `make -C pre_gapmoe`: passed.
- `make python`: passed.
- `make test`: passed all configured tests:
  - `unit_core`
  - `smoke_cli`
  - `regression_cli`
  - `python_binding`
- Manual CLI check passed:
  `./genulens l 0.5 b 0.2 Nsimu 3 seed 1234`
- Manual Python check passed:
  `Config`, `simulate()`, callable likelihood, `to_numpy()`, and `ruc(...)`.

## Current Architecture

- `./genulens` is a thin wrapper over the shared core entry point.
- The production scientific `genulens` implementation lives under `src/genulens/simulation/scientific_engine.cpp` and is invoked through `genulens::run_scientific_cli()`.
- `EventSimulator` uses `ScientificSimulationBackend`, which runs the scientific backend with event-level verbosity, parses the generated physical events into `SimulationResult`, and applies C++/Python likelihood callables.
- The Python binding calls the same `EventSimulator` path as C++.
- The pre-gapmoe implementation sources have been moved under `src/genulens/tools/pre_gapmoe/` and are built from there.
- Input data resolution no longer depends on the caller's working directory. It checks direct path, cwd `input_files/`, `GENULENS_INPUT_DIR`, source-tree data, and installed shared-data locations.

## Known Limitations

- The scientific engine is now callable from the shared core and Python, and ownership of the old global state has moved into `ScientificEngine`/`ScientificState`. A compatibility bridge still exists inside `scientific_engine.cpp` so the remaining C-style subroutines can be extracted incrementally without changing numerical behavior.
- Population normalization, luminosity-function construction, density formulas, and kinematic sampling still need further class extraction. The first real cuts are in place for RNG, likelihood, interpolation, quadrature, coordinates, and IMF grid generation.
- Python custom likelihood is applied to parsed generated events as an additional weight. It does not yet alter the scientific backend importance-sampling proposal distribution before event generation.

## Python Example

```python
import genulens

cfg = genulens.Config(l=1.0, b=-3.9, n_simu=10000, seed=1234)
result = genulens.simulate(cfg)
arr = result.to_numpy()

def my_like(event):
    return 1.0 if event.tE > 10 else 0.0

custom = genulens.simulate(cfg, likelihood=my_like)

quick = genulens.ruc(l=0.0, b=0.0, n_simu=1000, seed=1)
```

## Human Review Points

- Confirm that the staged CMake fallback GSL search paths are acceptable for the target deployment environments.
- Continue extracting the remaining `scientific_engine.cpp` subroutines into:
  `StellarPopulationModel`, `DensityModel`, `KinematicsSampler`, `LuminosityFunction`, and `EventSampler`.
- Decide whether Python custom likelihood should affect proposal sampling, not just final event weights.
