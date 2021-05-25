[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/369252917.svg)](https://zenodo.org/badge/latestdoi/369252917)
[![arXiv](http://img.shields.io/badge/arXiv-2104.03306-orange.svg?style=flat)](https://arxiv.org/abs/2104.03306)


# genulens
`genulens`, which stands for "generate microlensing", is a tool to simulate microlensing events using Monte Carlo simulation of the Galactic model developed by [Koshimoto, Baba & Bennett (2021), arXiv:2104.03306](https://arxiv.org/abs/2104.03306).  
The code itself has been published as [Koshimoto & Ranc (2021), Zenodo.4784949](http://doi.org/10.5281/zenodo.4784949).   
Please cite the papers if you find this code useful in your research. 

Data values of files in input\_files/ directory mostly come from other people's studies.  
See the header of each file for references.  
Please cite those references as well if they are specially relevant to your work.  

The copyright of an included supplementary code, "option.cpp", belongs to Ian A. Bond and Takahiro Sumi.
 

## Installation
If you have `git`, you can download the package by 
``` 
$ git clone https://github.com/nkoshimoto/genulens.git
```
This is recommended because that way you can track any future updates with `git`.

If you do not have `git`, you can simply download the package by clicking the green button "Code" on the upper right in [the repository page](https://github.com/nkoshimoto/genulens), selecting "Download ZIP", and then unzipping it.

In the downloaded directory (genulens), move into another genulens directory  
``` 
$ cd genulens  
```
There is Makefile in the directory, and you need a C++ compiler to `make`.  
The default compiler is `clang++`, which is available in macOS.  
Please replace the first line in Makefile with `g++` or any other C++ compiler if you prefer.  
If you are not sure if you have `clang++` or `g++`, you can check it by  
```
$ which clang++ (or g++)
```
If your terminal returns a full path for it, then you have the compiler.

After making sure that you have a C++ compiler and specify it in Makefile, you can compile the program by  
```
$ make
```

If  
```
$ ./genulens  
```
returns output lines that starts from   
> \#   Output of genulens   

and ends with  
> \# Nlike/N= 100000 / 100000      wtlike/wtlike\_tE= 112596 / 112596 = 1.000000  

you are ready to use `genulens`. Note that the exact numbers of the end line might depend on your environment because the calculation uses random numbers.  
Please make sure that the input\_files/ directory in the same directory where you run `genulens`.


## Usage
See [Usage.pdf](https://github.com/nkoshimoto/genulens/blob/main/Usage.pdf).
This document is still incomplete and will be revised in the future.


