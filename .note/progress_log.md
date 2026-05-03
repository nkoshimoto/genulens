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
- Moved the scientific `genulens.cpp` implementation into `src/genulens/simulation/scientific_engine.cpp`.
- Added thin root `genulens.cpp` wrapper that calls `genulens::run_scientific_cli()`.
- Added `ScientificSimulationBackend`, which invokes the scientific simulation with `VERBOSITY=3`, captures output, parses real event rows into `SimulationResult`, and applies C++/Python likelihood callables to event weights.
- Updated CMake so `genulens_core` owns the scientific backend and Python calls the same backend path.
- Moved pre-gapmoe implementation sources under `src/genulens/tools/pre_gapmoe/`; this is a mechanical relocation, not yet a semantic rewrite.
- `make test`: passed after the genulens-first backend change.
- Manual check passed: `./genulens l 0.5 b 0.2 Nsimu 10 seed 1234`.
- Manual Python check passed: `genulens.Config(...)`, `genulens.simulate(cfg)`, `to_numpy()`, and Python callable likelihood.

## 2026-05-03 09:48 JST - Robust Input File Resolution

- Confirmed the migrated scientific backend still had direct `input_files/...` `fopen()` calls.
- Added C++ input file compatibility helpers:
  - direct path first,
  - current working directory `input_files/`,
  - `GENULENS_INPUT_DIR`,
  - build/source-tree `input_files/`,
  - installed shared data locations.
- Wired the scientific genulens backend and moved pre-gapmoe sources through `genulens::open_input_file()`.

## 2026-05-03 10:05 JST - Naming and Parameter Cleanup

- Removed `legacy` naming from the active genulens source layout.
- Renamed the core simulation entry points to `scientific_engine`, `scientific_cli`, and `ScientificSimulationBackend`.
- Added `model::IMFParameters` and `model::ModelParameters` so IMF break masses and slopes are centralized.
- Exposed model/IMF parameters through `GenulensConfig` and pybind.

## 2026-05-03 10:20 JST - First Real Object Extraction

- User clarified that the previous backend wrapper was still not true objectification.
- Started extracting actual scientific-engine responsibilities:
  - replaced direct global GSL RNG access with `RandomEngine`;
  - added `ObservationLikelihood` and `ObservationConstraint`;
  - routed the scientific engine's `like_obs()` through the likelihood object.
- This is the first incremental cut; remaining global model, density, luminosity-function, and kinematic state still need to be moved into owned objects.

## 2026-05-03 10:35 JST - Scientific State Ownership

- Added `ScientificState` with explicit grouped state objects:
  - runtime/RNG,
  - stellar population and IMF tables,
  - density and structural parameters,
  - luminosity-function tables,
  - kinematics tables,
  - NSD moments.
- Added `ScientificEngine` as the object that owns `ScientificState`.
- Converted the scientific CLI to instantiate `ScientificEngine` and call `ScientificEngine::run()`.
- Removed the large file-scope global variable definitions from `scientific_engine.cpp`; the remaining C-style subroutines now access state through a transitional active-state bridge.
- This is not the final decomposition, but ownership has moved from process-global variables into an engine object, enabling the next cuts into `DensityModel`, `KinematicsModel`, and `StellarPopulationModel`.

## 2026-05-03 10:50 JST - Math Helper Extraction

- Added `math::Interpolation` under `src/genulens/math/`.
- Moved the substantive implementations of:
  - linear interpolation,
  - indexed interpolation,
  - uniform-grid interpolation,
  - bilinear interpolation,
  - bilinear interpolation coefficients,
  - inverse-CDF lookup with linearly interpolated density.
- Left thin compatibility wrappers in `scientific_engine.cpp` so the scientific path remains stable while the rest of the engine is split.
- `make test`: passed.

## 2026-05-03 11:05 JST - Coordinate Transformation Extraction

- Added `model::CoordinateTransformer` and `model::PositionAngle`.
- Moved the substantive line-of-sight Cartesian conversion and position-angle calculations out of `scientific_engine.cpp`.
- Left compatibility functions `Dlb2xyz()` and `calc_PA()` in the scientific engine as thin delegates to the coordinate object.
- Added a unit smoke assertion that position-angle calculation returns finite output.
- `make test`: passed.

## 2026-05-03 11:15 JST - Quadrature Extraction

- Added `math::NewtonCotes` for the Newton-Cotes integration locations and weights used by density and optical-depth integrations.
- Replaced the old `get_p_integral()` body with a thin delegate to `math::NewtonCotes::coefficients()`.
- Added unit checks for the order-4 coefficients to pin the table values.
- `make test`: passed.

## 2026-05-03 11:30 JST - IMF Grid Extraction

- Added `model::BrokenPowerLawIMF` and `model::MassFunctionGrid`.
- Moved broken-power-law IMF grid construction, cumulative number normalization, cumulative mass normalization, and percentile-index initialization out of `store_IMF_nBs()`.
- Kept population normalization and remnant-mass accounting in `store_IMF_nBs()` for now; it still depends on age tables and density normalization state.
- First `make test` failed because `mass_function.hpp` forward-declared `IMFParameters` while storing it by value.
- Second `make test` failed because vector-backed cumulative mass arrays needed `.data()` when passed to existing interpolation wrappers.
- Fixed both issues.
- `make test`: passed.

## 2026-05-03 11:45 JST - Python `ruc` Alias

- Manual Python check found `genulens.simulate()` worked, but `genulens.ruc(...)` was missing.
- Added `genulens.ruc(l=..., b=..., n_simu=..., seed=..., likelihood=...)` as a convenience API that builds `Config` and calls the same simulation path.
- Added smoke coverage for the alias.
- `make test`: passed.
- Manual Python check passed for `simulate()` with Python callable likelihood and `ruc()`.

## 2026-05-03 12:05 JST - Engine Naming and Init/Sampler Pipeline

- Removed the `ScientificEngine` naming from the active C++ API and source layout.
- Renamed the run state to `GenulensRunContext`.
- Added `GenulensInitializer`, currently responsible for creating the run context.
- Added `GenulensSampler`, currently responsible for executing the existing CLI sampling path against a run context.
- Changed `./genulens` and `apps/genulens_main.cpp` to call `run_genulens_cli()`, which now composes initializer + sampler.
- Renamed `ScientificSimulationBackend` to `GenulensSimulationBackend`.
- Confirmed no remaining active references to `ScientificEngine`, `ScientificState`, `scientific_engine`, `scientific_cli`, or `scientific_backend`.
- `make test`: passed.

## 2026-05-03 12:25 JST - Event-Loop Submodel Extraction

- Added `model::RemnantMassModel` for the Lam/Raithel initial-final mass and remnant-type calculation.
- Replaced `Mini2Mrem()` with a thin delegate to `RemnantMassModel`.
- Added `model::BinaryLensSampler` for the Koshimoto+20 projected binary separation sampling.
- Replaced `getaproj()` with a thin delegate to `BinaryLensSampler`.
- Added unit coverage for deterministic WD remnant mass and binary separation sanity.
- `make test`: passed.

## 2026-05-03 12:40 JST - Monte Carlo Statistics State

- Added `MonteCarloStats` inside the sampler translation unit.
- Moved event-loop counters, likelihood-weight totals, remnant/binary counts, and tE histogram storage into that state object.
- Kept temporary macro aliases so the numerical loop can be split further without changing expressions in the same patch.
- `make test`: passed.

## 2026-05-03 12:55 JST - Sampling Option Initialization

- Added `SamplingOptions` to `GenulensRunContext`.
- Moved parsing of simulation controls (`Nsimu`, `NlikeMIN`, Earth velocity, importance-sampling weights, verbosity, binary/remnant flags, and prior flags) into `GenulensInitializer::read_sampling_options()`.
- `GenulensSampler` now consumes those prepared options from context.
- `make test`: passed.

## 2026-05-03 13:10 JST - Extinction Object Extraction

- Added `model::ExponentialDustExtinction`.
- Moved dust-scale normalization for `AI0`, `AK0`, and `EVI0` into the extinction object.
- Replaced repeated distance-extinction expressions in the density-grid preparation and lens I/K brightness calculations with `ExponentialDustExtinction` calls.
- `make test`: passed.
- Added unit/smoke coverage for resolving files from a different current working directory.
