# Pre-refactor Uncommitted Status

Timestamp: 2026-05-03 08:46:19 JST

Initial command:

```sh
git status --porcelain=v1 -uall
```

Output:

```text
?? .gaplan/README.md
?? .gaplan/gapmoe_rewrite.md
?? gapmoe_src/calc_mass_distribution_each
?? gapmoe_src/calc_mass_distribution_each.c
?? gapmoe_src/calc_rhon_at_Dlb
?? gapmoe_src/calc_rhon_at_Dlb.c
?? gapmoe_src/compile.sh
?? gapmoe_src/genulens_helio
?? gapmoe_src/genulens_helio.c
?? gapmoe_src/murel_sampling
?? gapmoe_src/murel_sampling.c
?? gapmoe_src/murel_sampling_2d
?? gapmoe_src/murel_sampling_2d.c
?? gapmoe_src/option.c
?? gapmoe_src/option.h
?? gapmoe_src/random.c
?? gapmoe_src/random.h
```

Notes:

- `.gaplan/` appears to contain planning material.
- `gapmoe_src/` appears to contain external or legacy GAPMOE source plus generated binaries and object files.
- These paths were not overwritten or removed.
- Because `gapmoe_src/` contains build products, it is not included in the checkpoint commit unless explicitly requested later.

