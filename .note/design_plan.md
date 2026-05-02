# Design Plan

## Goal

Move duplicated galaxy model, RNG, density, mass function, kinematics, likelihood, input, output, and simulation logic into a reusable C++ core under `src/genulens/`, then use that core from:

- the existing `./genulens` CLI,
- the `pre_gapmoe` tools,
- a pybind11 Python extension.

## Target Shape

- `GenulensConfig` owns CLI/Python shared configuration: seed, `l`, `b`, `Nsimu`, likelihood parameters, input path settings, and output settings.
- `RandomEngine` wraps GSL RNG with RAII, copy disabled, move enabled, and `uniform()` / `gaussian()` helpers.
- `InputDataRepository` centralizes `input_files/` lookup and loading. Search order: current directory, `GENULENS_INPUT_DIR`, installed/shared data path.
- `GalacticModel` owns coordinate conversion, density, mass function, and kinematics services used by both main simulation and `pre_gapmoe`.
- `EventSimulator` performs Monte Carlo simulation and returns `SimulationResult` / `std::vector<Event>`.
- `Likelihood` defaults to current Gaussian behavior and also supports custom C++ callables.
- Python binding exposes `Config`, `simulate()`, `SimulationResult.to_numpy()`, `SimulationResult.columns`, and Python callable likelihood dispatch.

## Compatibility Rules

- Keep existing CLI argument compatibility for `./genulens`.
- Keep executable names for `pre_gapmoe/calc_rho_profile`, `pre_gapmoe/calc_mass_dist`, and `pre_gapmoe/calc_murel_dist`.
- Avoid changing data files, `Usage.pdf`, and `genulens_samples.ipynb`.
- Preserve existing output formats unless a change is explicitly documented in `.note/breaking_changes.md`.
- Do not commit build products, caches, or generated bulk artifacts.

## Phases

1. Phase 0: protect current state, create checkpoint and backup branch, capture current build behavior.
2. Phase 1: create smoke/regression tests around representative current CLI and pre-gapmoe behavior.
3. Phase 2: add `src/genulens` library skeleton and move logic with minimal behavioral change.
4. Phase 3: replace duplicated CLI/pre-gapmoe code with shared-library calls.
5. Phase 4: add pybind11 binding and Python tests.
6. Phase 5: add CMake, keep Makefile compatibility, add pyproject if practical.
7. Phase 6: run tests, update docs, write final report, commit final refactor.

## Risk Controls

- Snapshot existing behavior before refactoring.
- Prefer mechanical moves over scientific rewrites.
- Record every failed build/test and unresolved issue in `.note/progress_log.md`.
- Use narrow commits/checkpoints so recovery is practical.

