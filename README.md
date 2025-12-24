This project contains two programs enabling to generate uncorrelated and correlated random force field. 

--------------------------------------------------------------------------
--------------------------------------------------------------------------
FNG_iso: generates uncorrelated and Gaussian correlated random force field
--------------------------------------------------------------------------
--------------------------------------------------------------------------

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
#8: file name for the output random force field.

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

-----------------------------------------------------------------------------------------------------------------------------------------------------
-----------------------------------------------------------------------------------------------------------------------------------------------------
FNG_elastic: generates correlated random force field emerging from elastic interactions between a random solid solution and a dissociated dislocation
-----------------------------------------------------------------------------------------------------------------------------------------------------
-----------------------------------------------------------------------------------------------------------------------------------------------------

------------
Description:
------------
This C program generates a correlated random field emerging from a random alloy assuming elastic interactions between the dislocation and the solute field. The dislocation line is assumed to be aligned along the x direction, and propagating in the y direction, remaining on the plane of normal z. The resulting force field incoporates the spatial correlations of the stress field discussed in Refs. [1,2] acting on the partials of the dissociated dislocation. We assume the spacing between partials to remain constant.

[1] Geslin, Pierre-Antoine, and David Rodney. "Microelasticity model of random alloys. Part I: mean square displacements and stresses." Journal of the Mechanics and Physics of Solids 153 (2021): 104479.\
[2] Geslin, Pierre-Antoine, Ali Rida, and David Rodney. "Microelasticity model of random alloys. Part II: displacement and stress correlations." Journal of the Mechanics and Physics of Solids 153 (2021): 104480.

-------
Inputs:
-------
The arguments of the program are specified as line command as follows: \
#1: Nx, number of grid points along x direction [int] \
#2: Ny, number of grid points along y direction [int] \
#3: dx, grid spacing in x direction [AA] \
#4: dy, grid spacing in y direction [AA] \
#5: a, spreading parameters partial core structure [AA] \
#6: diss, dissociation distance [AA] \
#7: phi, dislocation character in radian ranging from 0 (screw) to 90 (edge) [rad] \
#8: E_per, this flag is 1 to indicate that the energy profile corresponding to the force field is periodic along the y direction; E_per=0 otherwise [0 or 1]. \
#9: random seed [int] \
#10: output file name 

-------
Ouputs:
-------

The outputs is a vtk file with name given in the command line containing the noise on a 2D regular grid of dimension Nx x Ny. Note that the file is binary and each float is written with the Big Endian convention, which is traditional for vtk files. Also, the values of the field f[i,j] are listed in row-major order: f[i,j],f[i,j+1],...,f[i+1,j],f[i+1,j+1]...

--------
Example:
--------

Generate a noise with a spreading parameter of 1AA for an undissociated edge dislocations, while imposing the periodicity of the energy landscape (a=3AA, E_per=1):
./noise_generator.x 256 256 1.0 1.0 3.0 0 90 1 1234 noise.vtk



