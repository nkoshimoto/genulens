[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-yellow.svg)](http://www.gnu.org/licenses/gpl-3.0)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4784948.svg)](https://doi.org/10.5281/zenodo.4784948)
[![arXiv](http://img.shields.io/badge/arXiv-2104.03306-orange.svg?style=flat)](https://arxiv.org/abs/2104.03306)



# genulens
`genulens`, which stands for "generate microlensing", is a tool to simulate microlensing events using Monte Carlo simulation of the Galactic model developed by [Koshimoto, Baba & Bennett (2021), ApJ, 917, 78](https://ui.adsabs.harvard.edu/abs/2021ApJ...917...78K/abstract).
The Galactic model is optimized for the bulge direction, and it is best to be used for analyzing microlensing events in |l| < 10 deg. and |b| < 7 deg.  
The code itself has been published as [Koshimoto & Ranc (2021), Zenodo.4784948](http://doi.org/10.5281/zenodo.4784948).
Please cite the papers if you find this code useful in your research.

A simulator for stars, [`genstars`](https://github.com/nkoshimoto/genstars), is also available.

The copyright of an included supplementary code, "option.cpp", belongs to Ian A. Bond and Takahiro Sumi.


## Major updates
### v1.2 and star generator (June-July 2022)
- A much more efficient simulation is available thanks to a newly introduced importance sampling feature that only simulates events with parameters that are close to the input parameters. This allows most calculations to be performed 5 to 15 times faster than before.  
- Although [Koshimoto, Baba & Bennett (2021), ApJ, 917, 78](https://ui.adsabs.harvard.edu/abs/2021ApJ...917...78K/abstract) did not developed their Galactic model for the Galactic Central region (|b| < 1 deg.), the nuclear stellar disk (NSD) component was added to `genulens` due to a demand, which affects |b| < ~0.5 deg. Although it was not fitted to the data for the Galactic Central region, we have confirmed that there is no major discrepancy with the star count data by the GALACTICNUCLEUS survey (Koshimoto et al., in preparation).  
- Accordingly, the position of the Galactic Center in the model was changed from (l, b) = (0, 0) to (l, b) = (-0.056, -0.046) deg., which is the position of Sgr A*.  
- The directory hierarchy was also changed. The former nkoshimoto/genulens/genulens was deleted and its contents were brought under nkoshimoto/genulens/. Therefore, if you have been using a wrapper program that assumes the previous directory structure, please change it.
- I have made a significant revision on the [Usage.pdf](https://github.com/nkoshimoto/genulens/blob/main/Usage.pdf). Now it explains many options available in `genulens`.  
- A jupyter-notebook, [genulens_samples.ipynb](https://github.com/nkoshimoto/genulens/blob/main/genulens_samples.ipynb), is added to the repository to introduce practical analysis of some events.  
- `genstars`, a version to simulate stars rather than microlensing events, is now available [here](https://github.com/nkoshimoto/genstars)


### v1.1 (June 2021)
On June 3 2021, we released version 1.1 using a random number generator from GSL (GNU Scientific Library).
This is because we realized that the random number generator from the C++ standard library, which was used in v1.0, does not have very good randomness.
The statistics (the median, 1 sigma, 2 sigma values) were probably OK, but the distributions were jagged compared to the GSL one with a same number of simulation.
Note that this version requires that you have GSL in your environment.

### v1.0 (May 2021)
First release on May 24 2021.

## Before installation
Please ensure that GSL (GNU Scientific Library) is installed in your environment.
If you do not have GSL, please download it from the link provided on the [GSL page](https://www.gnu.org/software/gsl/), and install it following README or INSTALL in the downloaded directory.

You can check where you have GSL by
```
gsl-config --libs
```
If the command gsl-config does not work in the terminal, it probably means that the GSL lib is not installed, or unknown to the OS.


## Installation
If you have `git`, you can download the package by
``` 
git clone https://github.com/nkoshimoto/genulens.git
```
This is recommended because that way you can track any future updates with `git`.

If you do not have `git`, you can simply download the package by clicking the green button "Code" on the upper right in [the repository page](https://github.com/nkoshimoto/genulens), selecting "Download ZIP", and then unzipping it.

You need a C++ compiler to `make` using Makefile to compile the program genulens.cpp in the downloaded directory.  
The default compiler is `clang++`, which is available in macOS.  
Please replace the first uncommented line in Makefile with `g++` or any other C++ compiler if you prefer.  
If you are not sure if you have `clang++` or `g++`, you can check it by
```
which clang++ (or g++)
```
If your terminal returns a full path for it, then you have the compiler.

You might also need to change the paths for GSL specified in Makefile.  
There are two lines to specify paths for GSL in the file;
> INCLUDE = -I/opt/local/include  
> LINK = -L/opt/local/lib

Please replace the path in INCLUDE with the path returned by
```
gsl-config --cflags
```
and replace the path in LINK with the path returned by
```
gsl-config --libs
```



After making sure that you specify your C++ compiler and the paths for GSL in Makefile, you can compile the program by
```
make
```

If
```
./genulens
```
returns output lines that starts from
> \#   Output of "./genulens "

and ends with
> \# Nlike/N= 100000 / 100000      wtlike/wtlike_tE= 112448 / 112448 = 1.000000

you are ready to use `genulens`. Note that the exact numbers of the end line might depend on your environment because the calculation uses random numbers.
Please make sure that the input\_files/ directory is in the same directory as where you run `genulens`.


## Usage
See [Usage.pdf](https://github.com/nkoshimoto/genulens/blob/main/Usage.pdf).  
[genulens_samples.ipynb](https://github.com/nkoshimoto/genulens/blob/main/genulens_samples.ipynb) is also available and you can follow some practical event analysis using `genulens`.  


## GAPMOE preprocessing helpers

This repository also includes a small set of standalone preprocessing tools in
`pre_gapmoe/`.  They are C++ rewrites of the older GAPMOE helper programs and
share the Galactic model implementation used by `genulens`.  The tools are
intended for generating tabulated inputs such as line-of-sight density profiles,
present-day mass functions, and relative proper-motion distributions.

### Building the helper tools

The helper tools require GSL and a C++11 compiler.  If `gsl-config` is available
on your `PATH`, build them with

```
make -C pre_gapmoe
```

If GSL is installed in a non-standard location, pass the GSL flags explicitly:

```
make -C pre_gapmoe \
  GSL_CFLAGS="-I/path/to/gsl/include" \
  GSL_LIBS="-L/path/to/gsl/lib -lgsl -lgslcblas -lm"
```

For a non-standard shared-library path, also set `LD_LIBRARY_PATH` when running
the tools, or add an rpath through `LDFLAGS`, for example
`LDFLAGS="-Wl,-rpath,/path/to/gsl/lib"`.

The build creates three executables:

```
pre_gapmoe/calc_rho_profile
pre_gapmoe/calc_mass_dist
pre_gapmoe/calc_murel_dist
```

Each program supports `--help`:

```
./pre_gapmoe/calc_rho_profile --help
./pre_gapmoe/calc_mass_dist --help
./pre_gapmoe/calc_murel_dist --help
```

### `calc_rho_profile`

`calc_rho_profile` outputs number-density profiles along a line of sight.
It can also apply source-selection weights using extinction, reddening, red-clump
distance modulus, and apparent source magnitude/color cuts.

Example:

```
./pre_gapmoe/calc_rho_profile \
  l 1.0 b -3.9 Dmax 16000 Dstep 100 \
  AIrc 1.5 EVIrc 1.2 Isrange 14 21
```

Important source-selection options include `AIrc` or `IsAIrc`, `EVIrc`, `DMrc`,
`Isrange`, `Is`, `Iserr`, `VIsrange`, `VIs`, `VIserr`, `hdust`, and `gammaDs`.
When source selection is active, the output appends `rhoD_S[0..10]` and
`rhoD_S_tot` columns.

### `calc_mass_dist`

`calc_mass_dist` outputs the present-day mass function for each Galactic
component.

Example:

```
./pre_gapmoe/calc_mass_dist
```

The output columns are `logM[Msun]`, `dN/dlogM[0..10]`, and `dN/dlogM_tot`.

### `calc_murel_dist`

`calc_murel_dist` outputs the relative proper-motion distribution either on a
`Dl x Ds` grid or for a single `(Dl, Ds)` pair.  It always reports both the
proper-motion amplitude `murel` and the direction angle `phi`.

By default, automatic precision control is enabled:

```
AUTOERR 1
ERR_TARGET 0.03
ERR_CHECK 100000
Nsimu 10000000
```

`Nsimu` is treated as the maximum number of accepted Monte Carlo draws.  The run
stops early only when both the `murel` and `phi` histograms satisfy the requested
relative Poisson precision.  Precision diagnostics are printed in comment lines,
including `mu_overflow_frac`, which reports the fraction of draws outside the
configured `murel` histogram range.

Examples:

```
./pre_gapmoe/calc_murel_dist

./pre_gapmoe/calc_murel_dist GRID 0 Dl 4000 Ds 8000

./pre_gapmoe/calc_murel_dist \
  AUTOERR 1 ERR_TARGET 0.03 ERR_CHECK 100000 \
  mumax 1000 dmu 0.5
```

The default `mumax` is `300 mas/yr`, which is sufficient for many sight lines.
Very nearby lens bins can have larger relative proper motions; in those cases,
increase `mumax` until `mu_overflow_frac` is acceptably small.

In grid mode the main output columns are

```
DS[pc]  DL[pc]  murel[mas/yr]  phi[rad]  dP/dmurel  dP/dphi
```

With `AUTOERR 1`, the following diagnostic columns are appended:

```
mu_count  mu_relerr  phi_count  phi_relerr  Ndraw
```

`phi` is defined as `atan2(mu_E, mu_N)` after rotating the Galactic
proper-motion components into north/east coordinates, matching the convention
used by the legacy `murel_sampling.c` helper.
