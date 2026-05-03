# Current genulens Refactor Structure

This document records the current implementation structure after the first real sampler split. It is intentionally descriptive: it says which file owns which responsibility, what it reads, and how data flows through CLI and Python.

## Top-Level Execution Paths

### CLI: `./genulens`

Entry point:

- `genulens.cpp`
  - Includes `genulens/simulation/cli.hpp`.
  - Calls `genulens::run_cli(argc, argv)`.

CLI dispatcher:

- `src/genulens/simulation/cli.cpp`
  - Handles `--help` and `-h` directly.
  - Creates an `Initializer`.
  - Calls `Initializer::create_context()`.
  - Calls `Initializer::initialize_rng(context, argc, argv)`.
  - Creates `Sampler`.
  - Calls `Sampler::run_cli(context, argc, argv)`.

Current CLI flow:

```text
genulens.cpp
  -> simulation/cli.cpp::run_cli()
      -> Initializer::create_context()
      -> Initializer::initialize_rng()
      -> Sampler::run_cli()
          -> Initializer::read_model_options()
          -> Initializer::finalize_spatial_model()
          -> PopulationRuntime
          -> KinematicRuntimeTables
          -> NsdMomentRuntime
          -> density/kinematics/stellar/likelihood helper functions
          -> Monte Carlo loop
          -> stdout output
```

### Python: `import genulens`

Python binding:

- `python/genulens_pybind.cpp`
  - Exposes `genulens.Config`.
  - Exposes `genulens.Event`.
  - Exposes `genulens.SimulationResult`.
  - Exposes `genulens.simulate(cfg, likelihood=None)`.
  - Exposes `genulens.ruc(l=..., b=..., n_simu=..., seed=..., likelihood=None)`.
  - Converts `SimulationResult` to `numpy.ndarray` via `to_numpy()`.
  - Accepts a Python callable likelihood. The callable receives a bound `Event` object and returns a double multiplier.

Python simulation backend:

- `src/genulens/simulation/simulator.cpp`
  - `EventSimulator` is the public C++ interface.
  - If no likelihood is provided, installs default Gaussian likelihood.
  - Calls `SimulationBackend`.

- `src/genulens/simulation/backend.cpp`
  - Converts `GenulensConfig` into CLI-style argument tokens.
  - Forces `VERBOSITY 3` for parseable event rows.
  - Captures stdout from `run_cli()`.
  - Parses non-comment output rows into `Event`.
  - Applies C++ or Python likelihood to each parsed event weight.

Current Python flow:

```text
python/genulens_pybind.cpp
  -> genulens::simulate(cfg, likelihood)
      -> EventSimulator::simulate()
          -> SimulationBackend::build_event_args()
          -> SimulationBackend::run_cli_capture()
              -> simulation/cli.cpp::run_cli()
                  -> Sampler::run_cli()
          -> SimulationBackend::parse_verbosity3_events()
              -> optional likelihood(event)
              -> SimulationResult
```

Important current limitation:

- Python still reaches the scientific engine through captured CLI stdout.
- That is functional and shared, but not the final clean design. The next target is to make the event loop return `SimulationResult` directly instead of round-tripping through text.

## Shared Runtime State

### `RunContext`

File:

- `src/genulens/simulation/run_context.hpp`

Purpose:

- Owns the mutable simulation state that used to be scattered as globals.
- Groups related state into:
  - RNG/runtime
  - sampling options
  - IMF options
  - spatial options
  - stellar population state
  - density state
  - luminosity state
  - kinematic state
  - NSD moment state

Current bridge:

- `src/genulens/simulation/internal/runtime.hpp` declares:

```cpp
extern RunContext *active_state;
```

- Runtime helper files still access current state through `active_state` and compatibility macros.
- This is transitional. It allowed moving functions out of `sampler.cpp` without changing scientific formulas at the same time.
- Final direction should replace these macros with methods/classes that take explicit references to the relevant model objects.

## Initialization

### `Initializer`

Files:

- `src/genulens/simulation/initialize.hpp`
- `src/genulens/simulation/initialize.cpp`

Responsibilities:

- `create_context()`
  - Creates a default `RunContext`.

- `initialize_rng(context, argc, argv)`
  - Reads `seed`.
  - Creates `RandomEngine`.

- `read_model_options(context, argc, argv)`
  - Reads IMF options:
    - `M0`, `M1`, `M2`, `M3`, `Ml`, `Mu`
    - `alpha0` through `alpha4`
  - Reads disk/bar/bulge/X-shape density parameters.
  - Reads disk/bar kinematic parameters.
  - Applies named legacy presets:
    - `B14disk`
    - `B14bar`
    - `E_fg0`
    - `G_fg0`
    - `EXE_fg0`
    - `GXG_fg0`

- `finalize_spatial_model(context, argc, argv)`
  - Computes `costheta` and `sintheta`.
  - Reads `CenSgrA`.
  - Computes Sgr A* coordinate offset using `CoordinateTransformer`.
  - Reads stellar halo options:
    - `SH`
    - `rho0SHMS`
    - `sigU_SH`
    - `sigV_SH`
    - `sigW_SH`

- `read_sampling_options(context, argc, argv, cos_pa, sin_pa)`
  - Reads Monte Carlo controls:
    - `Nsimu`
    - `NlikeMIN`
    - `VERBOSITY`
    - `REMNANT`
    - `BINARY`
    - `UNIFORM`
    - importance-sampling flags
    - Earth velocity projection

Current status:

- Model option initialization is no longer embedded directly in `sampler.cpp`.
- Observation-specific scalar options such as `tE`, `thetaE`, `piE`, `IL`, `KL`, `Isrange`, extinction options, and importance ranges are still read in `Sampler::run_cli()`.

## Random Numbers

Files:

- `src/genulens/rng.hpp`
- `src/genulens/rng.cpp`

Responsibilities:

- RAII wrapper around GSL RNG.
- Copy is disabled.
- Move is supported.
- Provides:
  - `uniform()`
  - `gaussian()`

Compatibility bridge:

- `sampler.cpp` defines:
  - `ran1()`
  - `gasdev()`

These call `active_state->runtime.rng`.

## Input Data Resolution

Files:

- `src/genulens/io/input_data.hpp`
- `src/genulens/io/input_data.cpp`

Responsibilities:

- `InputDataRepository::resolve(filename)` resolves data files robustly.
- `open_input_file(filename, mode)` is used by legacy-style `fopen()` calls through a macro.

Search order:

1. Direct path as provided.
2. `./input_files` from current working directory.
3. `GENULENS_INPUT_DIR`.
4. Source tree `input_files` using `GENULENS_SOURCE_DIR`.
5. Installed shared data:
   - `/usr/local/share/genulens/input_files`
   - `/usr/share/genulens/input_files`

Files currently read through this mechanism include:

- `input_files/Minidie.dat`
- `input_files/MLemp.dat`
- `input_files/NbleNall_bin.dat`
- `input_files/Rotcurve_BG16.dat`
- `input_files/NSD_moments.dat`
- stellar evolution tables used by luminosity/color distribution routines.

## Runtime Helper Split

These files are the current split of the old scientific helper functions. They still share `RunContext` through `internal/runtime.hpp`, but the function bodies are no longer all inside `sampler.cpp`.

### `src/genulens/simulation/internal/runtime.hpp`

Purpose:

- Transitional internal header.
- Declares `active_state`.
- Declares compatibility functions used across runtime translation units.
- Defines constants and state-access macros that were previously local to `sampler.cpp`.
- Defines small runtime-owner structs:
  - `PopulationRuntime`
  - `KinematicRuntimeTables`
  - `NsdMomentRuntime`

This file should shrink over time. It is currently the main compatibility layer keeping behavior stable while functions are moved.

### `src/genulens/simulation/stellar_population_runtime.cpp`

Owns:

- `PopulationRuntime`
  - Allocates and initializes IMF mass grid.
  - Allocates luminosity functions when `Isrange`/`AIrc` are active.
  - Allocates color-magnitude source selection when `VIsrange`/`EVIrc` are active.
  - Reads empirical mass-luminosity table `MLemp.dat`.
  - Releases luminosity and population tables.

Scientific helper functions:

- `store_IMF_nBs`
  - Builds broken-power-law IMF grid using `model::BrokenPowerLawIMF`.
  - Updates density normalization factors.
  - Reads `Minidie.dat`.

- `Mini2Mrem`
  - Evolves initial mass to remnant mass using `model::RemnantMassModel`.

- `fLF_detect`
  - Computes source fraction in an I-band selection.

- `fIVI_detect`
  - Computes source fraction in I and V-I selection.

- `read_MLemp`
  - Reads empirical mass-luminosity table.

- `make_LFs`
  - Builds luminosity functions.

- `store_VI_MI`
  - Builds V-I versus I source distribution.

### `src/genulens/simulation/kinematics_runtime.cpp`

Owns:

- `KinematicRuntimeTables`
  - Allocates Shu DF tables.
  - Reads `Rotcurve_BG16.dat`.
  - Calls `store_cumuP_Shu`.
  - Releases Shu DF tables.

- `NsdMomentRuntime`
  - Allocates NSD moment tables when `ND == 3`.
  - Reads `NSD_moments.dat`.
  - Releases NSD moment tables.

Scientific helper functions:

- `store_cumuP_Shu`
  - Builds cumulative probability distribution for `fg = Rg/R`.

- `get_PRRGmax2`
  - Finds bounds and peak for the Shu radial-guiding-radius distribution.

- `calc_dpdfg`
  - Numerical derivative for the Shu distribution.

- `get_vxyz_ran`
  - Samples 3D velocities for disk, bulge, NSD, and halo components.

- `getaproj`
  - Samples binary projected separation through `model::BinaryLensSampler`.

- `getcumu2xist`
  - Inverse cumulative helper.

- `calc_PRRg`, `calc_gc`, `calc_SigRg`, `calc_faca`
  - Shu distribution pieces.

- `calc_sigvb`
  - Bar velocity dispersion model.

### `src/genulens/simulation/density_runtime.cpp`

Owns:

- `crude_integrate`
  - Integrates bulge/bar density for normalization.

- `calc_rho_n`
  - Computes total mass density and number density at distance.

- `calc_rho_each`
  - Computes density contribution for each Galactic component:
    - thin disk bins
    - thick disk
    - bar/bulge
    - NSD
    - stellar halo

- `calc_rhoB`
  - Bar/bulge/X-shape density law.

### `src/genulens/simulation/observation_runtime.cpp`

Owns:

- `calc_PA`
  - Computes position angle through `model::CoordinateTransformer`.

- `calc_opticaldepth`
  - Computes optical depth and source number density when `CALCTAU` is enabled.

- `store_NSDmoments`
  - Reads NSD moments into runtime arrays.

- `like_obs`
  - Dispatches to `ObservationLikelihood`.
  - Supports Gaussian or uniform likelihood behavior based on `UNIFORM`.

### `src/genulens/simulation/math_runtime.cpp`

Owns:

- Vector helpers:
  - `cross`
  - `dot`
  - `norm_vec`

- Coordinate wrapper:
  - `Dlb2xyz`
  - Uses `model::CoordinateTransformer`.

- Numerical wrappers:
  - `get_p_integral`
  - `getx2y`
  - `getx2y_ist`
  - `getx2y_khi`
  - `interp_x`
  - `interp_xy`
  - `interp_xy_coeff`

These wrappers call the newer implementations in:

- `src/genulens/math/interpolation.cpp`
- `src/genulens/math/quadrature.cpp`

## Model Layer

### `src/genulens/model/parameters.*`

Owns:

- Default model parameters.
- IMF default values.
- Disk, density, stellar, kinematic defaults embedded in model structs.

### `src/genulens/model/mass_function.*`

Owns:

- `BrokenPowerLawIMF`
  - Builds mass grids and cumulative distributions.

- `RemnantMassModel`
  - Maps initial mass to remnant type/mass.

- `BinaryLensSampler`
  - Samples binary projected separation.

These are already real model objects, not just moved C functions.

### `src/genulens/model/coordinates.*`

Owns:

- Galactic coordinate conversion.
- Position angle calculation.

Used by:

- `Initializer::finalize_spatial_model`
- `math_runtime.cpp::Dlb2xyz`
- `observation_runtime.cpp::calc_PA`

### `src/genulens/model/extinction.*`

Owns:

- `ExponentialDustExtinction`.
- Computes:
  - `AI0`
  - `AK0`
  - `EVI0`
  - extinction at distance
  - distance modulus term

Used by:

- `Sampler::run_cli()` for source/lens brightness calculations and density-grid weighting.

### `src/genulens/model/density.*`, `galactic_model.*`, `kinematics.*`

Current status:

- Present as model-layer files.
- Some shared/simple model pieces live there.
- The full scientific density and velocity runtime is still in `density_runtime.cpp` and `kinematics_runtime.cpp`.
- The next refactor should move logic from runtime functions into explicit classes in these model-layer files.

## Likelihood Layer

Files:

- `src/genulens/simulation/likelihood.hpp`
- `src/genulens/simulation/likelihood.cpp`
- `src/genulens/simulation/observation_likelihood.hpp`
- `src/genulens/simulation/observation_likelihood.cpp`

Responsibilities:

- `LikelihoodFunction`
  - C++ `std::function<double(const Event&)>`.

- `default_likelihood`
  - Default Gaussian likelihood for Python-facing simulation.

- `ObservationLikelihood`
  - CLI-style observational constraint likelihood.
  - Supports detection mode:
    - detection
    - upper limit
    - lower limit
  - Supports Gaussian or uniform behavior.

Current split:

- CLI event loop still uses `like_obs()` from `observation_runtime.cpp`.
- Python-facing `simulate()` applies `LikelihoodFunction` after parsing CLI rows.

## Result/Event Layer

Files:

- `src/genulens/simulation/event.hpp`
- `src/genulens/types.hpp`
- `src/genulens/types.cpp`

Responsibilities:

- `Event`
  - Stores parsed physical event quantities:
    - weight
    - `tE`
    - `thetaE`
    - `piE`
    - lens/source distances
    - lens mass
    - relative proper motion
    - component IDs

- `SimulationResult`
  - Holds `std::vector<Event>`.
  - Provides `columns`.
  - Provides flattened rows for Python `to_numpy()`.

## Current `Sampler::run_cli()` Responsibilities

File:

- `src/genulens/simulation/sampler.cpp`

Still owns:

1. CLI scalar option read for observation constraints.
2. Runtime context activation through `active_state`.
3. Population table setup orchestration through `PopulationRuntime`.
4. Kinematic table setup orchestration through `KinematicRuntimeTables`.
5. NSD table setup orchestration through `NsdMomentRuntime`.
6. Header printing.
7. Bulge and NSD normalization orchestration.
8. Extinction setup.
9. Line-of-sight density grid construction:
   - `D`
   - `rhoD_S`
   - `rhoD_L`
   - `cumu_rho_S`
   - `cumu_rho_L`
   - percentile tables
   - BH scale-height correction table `fBH`
10. Optional `CALCTAU`.
11. Optional `CheckD`.
12. Monte Carlo event loop:
   - sample source component and distance
   - sample lens component and distance
   - sample velocities
   - sample lens mass/remnant/binary
   - compute microlensing observables
   - apply likelihood constraints
   - update summary statistics
   - print event rows by `VERBOSITY`
13. Summary printing.
14. Cleanup of the remaining line-of-sight arrays.

This is better than before, but still too much. The main remaining design problem is that the event loop and line-of-sight grid are still procedural inside `Sampler::run_cli()`.

## Pre-Gapmoe Tools

Files:

- `src/genulens/tools/pre_gapmoe/calc_rho_profile.cpp`
- `src/genulens/tools/pre_gapmoe/calc_mass_dist.cpp`
- `src/genulens/tools/pre_gapmoe/calc_murel_dist.cpp`
- `src/genulens/tools/pre_gapmoe/galactic_model.cpp`
- `src/genulens/tools/pre_gapmoe/galactic_kinematics.cpp`
- `src/genulens/tools/pre_gapmoe/option.cpp`

Current status:

- Built through CMake targets:
  - `calc_rho_profile`
  - `calc_mass_dist`
  - `calc_murel_dist`
- `make pre_gapmoe` copies built binaries into `pre_gapmoe/`.
- These tools still have their own compatibility code.
- They link against `genulens_core`, but are not yet fully rewritten to use the same `GalacticModel`/`EventSimulator` path as `./genulens`.

Refactor direction:

- Finish `genulens` objectization first.
- Then replace tool-specific logic with shared model/runtime classes.

## Build System

Files:

- `CMakeLists.txt`
- `Makefile`
- `pyproject.toml`

Current behavior:

- `make`
  - Configures CMake.
  - Builds `genulens`.
  - Copies `build/genulens` to `./genulens`.

- `make pre_gapmoe`
  - Builds pre-gapmoe tools.
  - Copies binaries to `pre_gapmoe/`.

- `make python`
  - Builds pybind11 module when Python/pybind11 are available.

- `make test`
  - Builds CLI, tools, tests, and Python module.
  - Runs CTest:
    - `unit_core`
    - `smoke_cli`
    - `regression_cli`
    - `python_binding`

## Current Dependencies by Runtime File

```text
sampler.cpp
  reads CLI scalar options
  calls Initializer
  calls PopulationRuntime
  calls KinematicRuntimeTables
  calls NsdMomentRuntime
  calls density_runtime functions
  calls kinematics_runtime functions
  calls observation_runtime likelihood/optical-depth functions
  calls model::ExponentialDustExtinction

stellar_population_runtime.cpp
  reads Minidie.dat
  reads NbleNall_bin.dat
  reads stellar color/magnitude tables
  reads MLemp.dat
  uses model::BrokenPowerLawIMF
  uses model::RemnantMassModel

kinematics_runtime.cpp
  reads Rotcurve_BG16.dat
  reads NSD_moments.dat
  uses RNG through ran1/gasdev
  uses interpolation wrappers

density_runtime.cpp
  uses coordinate wrapper Dlb2xyz
  uses interpolation for NSD density
  uses quadrature wrapper for crude integration

observation_runtime.cpp
  uses CoordinateTransformer for PA
  uses ObservationLikelihood
  uses density and LF functions for optical-depth calculation

math_runtime.cpp
  uses CoordinateTransformer
  uses math::Interpolation
  uses math::NewtonCotes

backend.cpp
  calls run_cli()
  captures stdout
  parses VERBOSITY=3 rows into Event
```

## Remaining Refactor Plan

The next splits should be:

1. `LineOfSightDensityGrid`
   - Move arrays `D`, `rhoD_S`, `rhoD_L`, `cumu_rho_*`, percentile tables, and `fBH` ownership out of `Sampler::run_cli()`.
   - Give it methods:
     - `build()`
     - `sample_source()`
     - `sample_lens()`
     - `bh_correction(component, distance)`
     - `release()` or RAII destructor.

2. `EventLoop` or `EventSampler`
   - Move the `for (long j=0; j<NSIMU; j++)` loop out of `Sampler::run_cli()`.
   - Inputs:
     - `RunContext`
     - `LineOfSightDensityGrid`
     - `PopulationRuntime`
     - observational constraints
     - output writer or result collector.
   - Output:
     - `SimulationResult` for Python.
     - existing stdout rows for CLI through an output adapter.

3. `ObservationConfig`
   - Move observed quantities and likelihood settings from local variables into a struct.
   - Should include:
     - `tE`, `thetaE`, `piE`, `piEN`, `piEE`
     - source proper motion constraints
     - heliocentric proper motion constraints
     - lens brightness constraints
     - detection modes
     - importance-sampling ranges.

4. Replace `active_state` macros
   - After the above split, convert runtime functions into classes:
     - `DensityModel`
     - `KinematicModel`
     - `StellarPopulation`
     - `ObservationModel`
   - Each should receive explicit references to `RunContext` or narrower state structs.

5. Python direct path
   - Stop capturing stdout in `SimulationBackend`.
   - Have `EventSampler` return `SimulationResult` directly.
   - Keep CLI output compatibility by formatting the same result/event stream.

## Practical Mental Model

At the moment:

- `Initializer` owns option-to-context setup.
- `RunContext` owns mutable model state.
- `PopulationRuntime` owns mass and luminosity tables.
- `KinematicRuntimeTables` owns Shu DF velocity tables.
- `NsdMomentRuntime` owns NSD moment tables.
- Runtime `.cpp` files own extracted scientific helper functions.
- `Sampler::run_cli()` still owns orchestration, line-of-sight grid construction, event loop, and text output.
- Python uses the same CLI engine, then parses event rows back into objects.

The next meaningful step is not adding more files. It is cutting `LineOfSightDensityGrid` and `EventSampler` out of `sampler.cpp`.
