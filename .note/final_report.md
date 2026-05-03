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

## Current Architecture

- `./genulens` is now a thin wrapper.
- The production scientific `genulens` implementation lives under `src/genulens/simulation/scientific_engine.cpp` and is invoked through `genulens::run_scientific_cli()`.
- `EventSimulator` uses `ScientificSimulationBackend`, which runs the scientific backend with event-level verbosity, parses the generated physical events into `SimulationResult`, and applies C++/Python likelihood callables.
- The Python binding calls the same `EventSimulator` path as C++.
- The pre-gapmoe implementation sources have been moved under `src/genulens/tools/pre_gapmoe/` and are built from there.

## Known Limitations

- The scientific engine is now callable from the shared core and Python, but some internal state still remains in the migrated scientific engine file. Further cleanup can split that file into the existing `model/`, `io/`, and `simulation/` modules without changing the public CLI/Python surface.
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
```

## Human Review Points

- Confirm that the staged CMake fallback GSL search paths are acceptable for the target deployment environments.
- Decide the next extraction order for the large `genulens.cpp` global state.
- Review whether the temporary approximate implementations in `src/genulens/model/*` should be replaced first by direct moves from `pre_gapmoe/galactic_model.cpp` or by wrapping that code behind a stateful adapter.
