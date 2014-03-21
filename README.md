PolarDecomp
===========

A collection of decomposition routines for polarimetric synthetic aperture radar (PolSAR) applications written in Fortran. Available decomposition methods include Pauli, eigenvalue/eigenvector, Cloude-Pottier H/A/&#945;, and Shannon.

requirements
------------
- Fortran 90 compiler (only tested on gfortran)

installation
------------
   $ make 

optionally:

   $ make install [prefix=/install/directory]

usage
-----
Copy PolarDecomp.cmd to a directory and fill out the necessary blanks. Then do:

   $ PolarDecomp PolarDecomp.cmd


