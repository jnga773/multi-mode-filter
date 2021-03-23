# Multi-Mode Filter Array
This is part of the PhD project of Jacob Peter Kia Ngaha. Contact email: j.ngaha@auckland.ac.nz

We simulate a two- or three-level atom, emitting fluorescence into a filter consisting of *2N+1* single-mode cavities, equally spaced in frequency with spacing *δω*. We couple the atom into the cavity by using *[open quantum cascaded system theory](https://link.aps.org/doi/10.1103/PhysRevLett.70.2273 "Phys. Rev. Lett. 70, 2273 (1993)")*.

This repository is split into two directory: **two_level** and **three_level**. Each folder contains a pdf of the relevant moment equations, as well as the LaTeX files for that pdf, and programs that:
- numerically integrate the operator moments, 
- solve for the steady state moments, 
- solve for the first-order correlation function,
- solve for the second-order correlation function, and
- solve for the second-order cross-correlation of two different filters.

The programs are written in Fortran90 and I would recommend compiling them with Intel's Parallel Studio `ifort` compiler, with the command in Linux
```
ifort -O3 -mkl -o [executable_name] [filename.f90]
```
or, in Windows,
```
ifort /O3 /Qmkl /o [executable_name] [filename.f90]
```
where the `-mkl`(`/Qmkl`) flag links to Intel's MKL libraries, which are included with the compiler, which is needed to access the `LAPACK` routines. The programs take the necessary parameters from the `ParamList.nml` file, which is included in each directory, hence the code only needs to be compiled once.
