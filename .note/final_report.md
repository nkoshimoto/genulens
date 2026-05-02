# Final Report

## Summary

- Protected the starting state with commit `checkpoint: before major genulens refactor`.
- Created backup branch `backup/pre-major-refactor-20260503-084619`.
- Added a CMake-backed build while preserving root `make` entry points.
- Added `src/genulens` shared-library skeleton covering config, RNG, input path resolution, model services, density, kinematics, mass function, extinction, likelihood, simulation, output, and tool helpers.
- Added pybind11 binding source with `Config`, `simulate()`, Python callable likelihood support, `SimulationResult.columns`, and `SimulationResult.to_numpy()`.
- Added smoke/unit/Python tests and regression comparison against the pre-refactor binary backup.

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

## Known Limitations

- The new `src/genulens` core currently provides a reusable API boundary and testable implementations, but the full scientific Monte Carlo logic from `genulens.cpp` has not yet been completely decomposed into the core.
- The production `./genulens` binary still compiles from the existing `genulens.cpp` to preserve CLI/output compatibility.
- The production `pre_gapmoe` binaries still compile from the existing `pre_gapmoe/*.cpp` sources through CMake. The new `src/genulens/tools/*` helpers are present as the target abstraction layer but are not yet the sole implementation for those tools.
- Full removal of duplicated scientific logic remains the next refactor step and should be done with narrower, function-by-function behavioral tests.

## Python Example

```python
import genulens

cfg = genulens.Config(l=1.0, b=-3.9, n_simu=10000, seed=1234)
result = genulens.simulate(cfg)
arr = result.to_numpy()

def my_like(event):
    return 1.0 if event.tE > 10 else 0.0

custom = genulens.simulate(cfg, likelihood=my_like)
```

## Human Review Points

- Confirm that the staged CMake fallback GSL search paths are acceptable for the target deployment environments.
- Decide the next extraction order for the large `genulens.cpp` global state.
- Review whether the temporary approximate implementations in `src/genulens/model/*` should be replaced first by direct moves from `pre_gapmoe/galactic_model.cpp` or by wrapping that code behind a stateful adapter.

