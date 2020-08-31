# VMC-1D (Under Construction)
Variational Monte Carlo implemented for the 1D Heisenberg Model and the Haldane-Shastry Model using a Gutzwiller projected wave function as the initial ansatz. (Fortran90) 


 ---------------------- Variational Monte Carlo 1D --------------------

 Author: João Augusto Sobral da Silva / Github: @joaosds

 The code is commented, and each subroutine has a brief description (B.D.)
 of its functionality. For any discussions you can contact me at Github!

 History:

  v1.0 (08/27/20) - First implementation in f90
  
 Future intended changes:
  - implementation of modern fortran (>=2008);
  - python script for plots;
  - python interface to link both the fortran file and the plots automatically;
  - spin spin correlation for long range haldane-shastry;

 ---------------------- Informations ---------------------------------

 File 'random.f90' - All codes responsible for random generating numbers
 and related topics/ we use the implementation 'rkiss05.f90' by Thomas Vojta,
 vojta@mst.edu

 File 'vmcvariables.f90' - Global variables

 File 'modelvmc.f90' - retains the lattice definition, energy and spin-spin
 calculation and initial wavefunction.
 
 File 'inputpar.dat': Input parameters as following:
 N, seed -> Chain Length and random seed;
 nz, sweep -> Number of next neighbors considered in the model / values
 used after the thermalization for the mean calculation;
 thermalization -> steps for themalization of the MC algorithm;
 Nbin - > Bin number for calculation of the means.

 File 'meanene.dat': saves the energy obtained from nonlocal+local method
 and only local methods (best suitable for large N).

 File 'spincor.dat': saves the z component of the spin-spin correlations;

 Note: One can change line 116 and 117 to inversesvd if the input matrices
 are ill-conditioned. This may improve a little the accuracy of the algorithm.
