# Breaking Changes

No intentional user-facing breaking changes were introduced.

Compatibility notes:

- `./genulens` is still built by `make`.
- `pre_gapmoe/calc_rho_profile`, `pre_gapmoe/calc_mass_dist`, and `pre_gapmoe/calc_murel_dist` keep their executable names.
- The root Makefile is now a CMake-backed compatibility wrapper.
- `pre_gapmoe/Makefile` delegates to the root Makefile.
- The implementation sources for `genulens` and pre-gapmoe have moved under `src/genulens/`; the executable names and command-line option style are preserved.
