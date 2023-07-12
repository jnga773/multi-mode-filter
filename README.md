# Multi-Mode Filter Array
This is part of the PhD project of Jacob Peter Kia Ngaha. Contact email: j.ngaha@auckland.ac.nz

We simulate a two- or three-level atom, emitting fluorescence into a filter consisting of *2N+1* single-mode cavities, equally spaced in frequency with spacing *δω*. We couple the atom into the cavity by using *[open quantum cascaded system theory](https://link.aps.org/doi/10.1103/PhysRevLett.70.2273 "Phys. Rev. Lett. 70, 2273 (1993)")*.

This repository is split into three main directories: **two_level**, **three_level**, and **plots_for_thesis**. The folders **two_level** and **three_level** each contain a pdf of the relevant moment equations, as well as the LaTeX files for that pdf, and programs that:
- numerically integrate the operator moments, 
- solve for the steady state moments, 
- solve for the first-order correlation function,
- solve for the second-order correlation function, and
- solve for the second-order cross-correlation of two different filters.

The folder **plots_for_thesis** contain all of the simulation code, data, and plotting files used for each of the figures presented in the PhD Thesis.

The programs are written in Fortran90 and and can be compiled using Intel's oneAPI `ifort` compiler, with the command in Linux
```
ifort -O3 -mkl -o [executable_name] [module_file.f90] [filename.f90]
```
```
ifort /O3 /Qmkl /o [executable_name] [module_file.f90] [filename.f90]
```
where the `-mkl`(`/Qmkl`) flag links to Intel's MKL libraries, which are included with the compiler, and are needed to access the `LAPACK` routines.

Similarly, you can compile them with `gfortran` using the command
```
gfortran -O3 -o [executable_name] [module_file.f90] [filename.f90] -I/path/to/lapack -L/path/to/lapack -lblas -llapack
```
You will need to compile and install the necessary lapack routine files ([BLAS](https://www.netlib.org/blas/), [LAPACK](https://www.netlib.org/lapack), and [LAPACK95](https://www.netlib.org/lapack95)).

The programs take the necessary parameters from the `ParamList.nml` file, which is included in each directory, hence the code only needs to be compiled once.

The code is structure in a way such that all of the necessary subroutines are contained in the module files `MODULE_single_filter.f90` or `MODULE_two_filters.f90`, for single- and two-filter correlations, respectively. 