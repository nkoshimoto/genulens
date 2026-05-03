# Progress Log

## 2026-05-03 08:46 JST - Phase 0 Start

- Checked initial branch: `nunota-dev`.
- Initial tracked tree was clean, but untracked `.gaplan/` and `gapmoe_src/` existed.
- Recorded untracked paths in `.note/pre_refactor_uncommitted_status.md`.
- Created `.note/design_plan.md` and this progress log.
- Decision: do not overwrite or delete untracked material; do not commit `gapmoe_src/` because it includes binaries and object files.

### Status Before Phase 0 Checkpoint

```text
?? .gaplan/
?? gapmoe_src/
?? .note/
```

## 2026-05-03 08:48 JST - Phase 0 Checkpoint

- Committed `.note` protection records with:
  `checkpoint: before major genulens refactor`
- Created backup branch:
  `backup/pre-major-refactor-20260503-084619`
- Existing untracked `.gaplan/` and `gapmoe_src/` were left untouched.

## 2026-05-03 08:50 JST - Phase 0 Build Baseline

- `make`: succeeded; existing root binary was already up to date.
- `./genulens --help`: succeeded and produced the standard model summary plus Monte Carlo output.
- `make -C pre_gapmoe`: succeeded; existing helper binaries were already up to date.
- `pre_gapmoe` helper `--help` checks succeeded for:
  - `calc_rho_profile`
  - `calc_mass_dist`
  - `calc_murel_dist`
- Copied current binaries to `backup/pre-refactor-bin/` for regression comparison. This directory is a local generated backup and should not be committed.

## 2026-05-03 09:05 JST - Build System Refactor

- Added initial `src/genulens` core API, apps directory, Python binding source, tests, CMake, and a CMake-backed root Makefile.
- `make`: succeeded and rebuilt `./genulens` through CMake.
- `make pre_gapmoe`: initially failed while linking `calc_rho_profile` because `galactic_model.cpp` references `store_cumuP_Shu`, which lives in `galactic_kinematics.cpp`.
- Fix: include `pre_gapmoe/galactic_kinematics.cpp` in every CMake `pre_gapmoe` target, matching the old shared-object set.

## 2026-05-03 09:09 JST - First Test Run

- `make test` built the CLI and pre-gapmoe targets.
- Passed:
  - `smoke_cli`
  - `regression_cli`
- Failed:
  - `unit_core`: `InputDataRepository` could not resolve `Minidie.dat` because CTest ran from `build/`.
  - `python_binding`: Python extension target was configured but not built before CTest.
- Fix in progress:
  - Set `unit_core` and `python_binding` working directory to the source root.
  - Make root `make test` attempt to build `genulens_python` before running CTest.

## 2026-05-03 09:12 JST - Tests Passing

- `make test`: passed all configured tests.
  - `unit_core`
  - `smoke_cli`
  - `regression_cli`
  - `python_binding`
- `make python`: succeeded and built the pybind11 extension.
- Updated `pre_gapmoe/Makefile` to be a compatibility wrapper around the root CMake-backed `make pre_gapmoe`, avoiding the old `gsl-config`-only path.
- Added `build/` and `backup/` to `.gitignore` to avoid committing generated binaries and local backups.

## 2026-05-03 09:18 JST - Finalization

- Verified Python usage manually with `PYTHONPATH=build python3`.
- Updated README build instructions for CMake-backed Makefile and Python target.
- Added `.note/breaking_changes.md`; no intentional breaking changes.
- Added `.note/final_report.md` with test results and known limitations.

## 2026-05-03 09:35 JST - Genulens-first Refactor Direction

- User clarified that the priority is not pre-gapmoe compatibility first; the priority is proper `genulens` modularization and pybind exposure of the real simulation logic.
- Moved the legacy scientific `genulens.cpp` implementation into `src/genulens/simulation/legacy_genulens.cpp`.
- Added thin root `genulens.cpp` wrapper that calls `genulens::run_legacy_cli()`.
- Added `LegacySimulationBackend`, which invokes the legacy scientific simulation with `VERBOSITY=3`, captures output, parses real event rows into `SimulationResult`, and applies C++/Python likelihood callables to event weights.
- Updated CMake so `genulens_core` owns the legacy scientific backend and Python calls the same backend path.
- Moved pre-gapmoe implementation sources under `src/genulens/tools/pre_gapmoe/`; this is a mechanical relocation, not yet a semantic rewrite.
- `make test`: passed after the genulens-first backend change.
- Manual check passed: `./genulens l 0.5 b 0.2 Nsimu 10 seed 1234`.
- Manual Python check passed: `genulens.Config(...)`, `genulens.simulate(cfg)`, `to_numpy()`, and Python callable likelihood.
