# genulens
`genulens`, which stands for "generate microlensing", is a tool to simulate microlensing events using Monte Carlo simulation of the Galactic model developed by [Koshimoto, Baba & Bennett (2021), arXiv:2104.03306](https://arxiv.org/abs/2104.03306).  
Please cite the paper if you find this code (or the Galactic model itself) useful in your research.  

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
you are ready to use `genulens`. Note that the exact numbers of the end line might depend on your environment because it uses random numbers.  
Please make sure that the "input\_files" directory in the same directory where you run `genulens`.



