# cz-index-matlab

MATLAB functions for computing Conley-Zehnder (CZ) indices of CR3BP periodic orbits
The bulk of the code contained in this library are helper functions translated
to MATLAB from the original Python library "cz-index" written by Otto van Koert 
(Github: ovkoert). The original Python library can be found at 
https://github.com/ovkoert/cz-index.

Two new functions "get_cz_index" and "get_split_cz_index" (not translated from Python) 
are provided which provide a simple 1-line-of-code way to compute the CZ index of a 
CR3BP periodic orbit. A few other corrections and improvements to the Python library 
were also made in this MATLAB functions. 

If you use this code, please cite the below two papers:

For this MATLAB code:
B. Kumar and A. Moreno, “Networks of Periodic Orbits in the Earth-Moon System Through a Regularized
and Symplectic Lens,” AAS/AIAA Astrodynamics Specialist Conference, Aug 2025. Paper AAS-25-677. 

For the methodology and original Python library by Otto van Koert: 
A. Moreno, C. Aydin, O. v. Koert, U. Frauenfelder, and D. Koh, “Bifurcation Graphs for the CR3BP via
Symplectic Methods,” J. Astronaut. Sci., Vol. 71, No. 6, 2024, p. 51.
