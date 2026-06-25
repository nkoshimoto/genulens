# Architecture notes

This document describes the current refactored structure for collaborators. It
is intentionally practical: it records the boundaries that exist now and the
parts that are still transitional.

## Goals

The refactor is organized around four goals:

1. Provide a direct Python API that returns typed results rather than parsing
   command-line stdout.
2. Keep the command-line program compatible with existing workflows.
3. Move the simulation toward explicit objects for configuration, prepared
   runtime state, and event generation.
4. Eventually share the same model objects with the `pre_gapmoe` helper tools.

The v2 public branch focuses on goals 1 and 2 while exposing the source-forward
and event-rate APIs needed for physically motivated prior construction.

## Current simulation flow

CLI path:

```text
genulens.cpp
  -> simulation/cli.cpp::run_cli()
      -> Initializer::create_context()
      -> Initializer::initialize_rng()
      -> Sampler::run_cli()
          -> read/finalize model options
          -> prepare population, kinematic, NSD, density-grid state
          -> PreparedSimulation
          -> EventSampler::run_cli()
              -> EventSampler::run()
              -> CliEventReporter for Monte Carlo rows and summary output
```

Python path:

```text
genulens.simulate(Config, likelihood=None)
  -> SimulationBackend::simulate()
      -> Sampler::simulate()
          -> PreparedSimulation
          -> EventSampler::run(..., custom_likelihood, event_sink, emit_cli_output=false)
              -> SimulationResult.events
```

The Python path does not invoke `./genulens` and does not parse stdout.

## Main objects

| Object | Role |
|---|---|
| `GenulensConfig` | Public typed configuration for Python/C++ API use |
| `RunContext` | Mutable runtime state used by the legacy scientific kernels |
| `Initializer` | Reads legacy options and initializes runtime state |
| `PreparedSimulation` | Owns prepared population tables, kinematic tables, NSD moments, density grid, mass function, event config, and cleanup |
| `LineOfSightDensityGrid` | Builds and samples the line-of-sight density grid |
| `MassFunction` | Wraps lens mass sampling and magnitude conversion data |
| `EventSampler` | Monte Carlo event loop, observation likelihood, custom likelihood, and event sink notification |
| `CliEventReporter` | CLI Monte Carlo row formatting and summary output |
| `SimulationResult` | Typed event collection exposed to Python and convertible to ndarray |

## Configuration boundary

The Python API exposes nested typed config objects:

- `cfg.observation`
- `cfg.source`
- `cfg.sampling`
- `cfg.runtime`
- `cfg.model.imf`
- `cfg.model.density`
- `cfg.model.kinematics`
- `cfg.model.nsd`
- `cfg.model.bh_kick`

`raw_cli_args` remains as a compatibility escape hatch for legacy options that
have not yet been promoted to typed fields. New public Python examples should
prefer typed fields.

## Custom likelihood boundary

`LikelihoodFunction` is evaluated inside `EventSampler::run()` after an `Event`
has been assembled and before the accepted event is sent to the sink. This means
Python custom likelihoods participate in the sampling loop rather than being
applied after parsing text output.

The event sink is used differently by each frontend:

- Python stores accepted events into `SimulationResult`.
- CLI prints accepted events through `CliEventReporter`.

## Transitional pieces

The following parts are intentionally not fully cleaned up yet:

- `sampler.cpp` still contains a large orchestration function.
- Some CLI header and input-parameter printing still lives in `sampler.cpp`.
- `active_state` and related compatibility macros remain in the runtime layer.
- `raw_cli_args` still exists for legacy options.
- `pre_gapmoe` tools are still command-line stdout tools, not Python direct APIs.

These are not blockers for the v2 public branch. They should be handled
incrementally to avoid changing too much scientific logic at once.

## Recommended next steps

1. Keep the Python API examples and docs focused on typed config.
2. Promote frequently used `raw_cli_args` options to typed config fields as
   users need them.
3. Move remaining CLI header printing out of `sampler.cpp` when it becomes
   painful, not preemptively.
4. Treat `active_state` removal as a larger scientific-kernel refactor requiring
   careful regression tests.
5. Add direct Python APIs for `pre_gapmoe` only after the CLI helper behavior is
   stable and covered by tests.
