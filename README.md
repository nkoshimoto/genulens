[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-yellow.svg)](http://www.gnu.org/licenses/gpl-3.0)
[![DOI](https://zenodo.org/badge/369252917.svg)](https://zenodo.org/badge/latestdoi/369252917)
[![arXiv](http://img.shields.io/badge/arXiv-2104.03306-orange.svg?style=flat)](https://arxiv.org/abs/2104.03306)



# genulens
`genulens`, which stands for "generate microlensing", is a tool to simulate microlensing events using Monte Carlo simulation of the Galactic model developed by [Koshimoto, Baba & Bennett (2021), arXiv:2104.03306](https://arxiv.org/abs/2104.03306).  
The code itself has been published as [Koshimoto & Ranc (2021), Zenodo.4784949](http://doi.org/10.5281/zenodo.4784949).   
Please cite the papers if you find this code useful in your research. 

The copyright of an included supplementary code, "option.cpp", belongs to Ian A. Bond and Takahiro Sumi.


## Warning
Since June 2 2021, the source code, genulens.cpp, has been replaced with a version with a random number generator from GSL (GNU Scientific Library).  
This is because we realized that the random number generator from the C++ standard library, which was used before, does not have a very good randomness.  
The statistics (the median, 1 sigma, 2 sigma values) were probably OK, but the distributions were jagged compared to the GSL one with a same number of simulation.  
Please install the new version with the GSL random number functions if you are using the version before June 2 2021.  
Note that this new version requires that you have GSL in your environment.  
 
## Before Installation
Please ensure that GSL (GNU Scientific Library) is installed in your environment.  
If you do not have GSL, please download it from the link provided on the [GSL page](https://www.gnu.org/software/gsl/), and install it following README or INSTALL in the downloaded derectory.

## Installation
If you have `git`, you can download the package by 
``` 
git clone https://github.com/nkoshimoto/genulens.git
```
This is recommended because that way you can track any future updates with `git`.

If you do not have `git`, you can simply download the package by clicking the green button "Code" on the upper right in [the repository page](https://github.com/nkoshimoto/genulens), selecting "Download ZIP", and then unzipping it.

In the downloaded directory (genulens), move into another genulens directory  
``` 
cd genulens  
```
There is Makefile in the directory, and you need a C++ compiler to `make`.  
The default compiler is `clang++`, which is available in macOS.  
Please replace the first line in Makefile with `g++` or any other C++ compiler if you prefer.  
If you are not sure if you have `clang++` or `g++`, you can check it by  
```
which clang++ (or g++)
```
If your terminal returns a full path for it, then you have the compiler.

After making sure that you have a C++ compiler and specify it in Makefile, you can compile the program by  
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
> \# Nlike/N= 100000 / 100000      wtlike/wtlike\_tE= 112366 / 112366 = 1.000000   

you are ready to use `genulens`. Note that the exact numbers of the end line might depend on your environment because the calculation uses random numbers.  
Please make sure that the input\_files/ directory in the same directory where you run `genulens`.


## Usage
See [Usage.pdf](https://github.com/nkoshimoto/genulens/blob/main/Usage.pdf).
This document is still incomplete and will be revised in the future.


