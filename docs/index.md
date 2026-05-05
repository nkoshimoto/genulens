# genulens documentation

This directory collects user-facing and developer-facing documentation for the
refactored `genulens` codebase.

## User guides

- [Python API](python_api.md): how to call the simulator from Python, work with
  `SimulationResult`, configure typed options, and provide a custom likelihood.
- [pre_gapmoe tools](pre_gapmoe.md): current status and command-line usage of
  the GAPMOE preprocessing helper tools.

## Developer notes

- [Architecture](architecture.md): the current object boundaries, simulation
  flow, and remaining refactor targets.

## Examples

- [`examples/python_binding.ipynb`](../examples/python_binding.ipynb): public
  notebook for the direct Python API.
- [`examples/cli_examples.ipynb`](../examples/cli_examples.ipynb): examples for
  the command-line program.
