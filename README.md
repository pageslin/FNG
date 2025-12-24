------------
Description:
------------
This C program generates uncorrelated or Gaussian correlated random force fields on a regular grid. It first generates random amplitudes for Fourier modes and perform an fft to obtain the real space force field. The resulting field as average 0 and standard deviation 1 (on average). Compilation requires the FFTW3 library. The output is a vtk binary file containing the values of the field on a regular grid.

-------
Inputs:
-------
The arguments of the program are specified as line command as follows: \
#1: Nx, number of grid points along x direction [int] \
#2: Ny, number of grid points along y direction [int] \
#3: dx, grid spacing in x direction [AA] \
#4: dy, grid spacing in y direction [AA] \
#5: sigma, Gaussian spreading distance [AA]. If sigma<0, uncorrelated noise is generated \
#6: E_per, this flag is 1 to indicate that the energy profile corresponding to the force field is periodic along the y direction; E_per=0 otherwise [0 or 1]. \
#7: random seed [int] \
#8: file name for the output random force field. \

-------
Ouputs:
-------
The outputs is a vtk file with name given in the command line containing the noise on a 2D regular grid of dimension Nx x Ny. Note that the file is binary and each float is written with the Big Endian convention, which is traditional for vtk files. Also, the values of the field f[i,j] are listed in row-major order: f[i,j],f[i,j+1],...,f[i+1,j],f[i+1,j+1]... 

---------
Examples:
---------

Generate an uncorrelated force field without periodic energy profile (sigma<=0 and E_per=0): \
./noise_generator.x 256 256 1.0 1.0 -1 0 1234 noise.vtk

Generate an uncorrelated force field with periodic energy profile (sigma<=0 and E_per=1): \
./noise_generator.x 256 256 1.0 1.0 -1 1 1234 noise.vtk

Generate a Gaussian correlated force field with periodic energy profile (sigma=2.0 and E_per=1): \
./noise_generator.x 256 256 1.0 1.0 2.0 1 1234 noise.vtk
