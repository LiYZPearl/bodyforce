How to compile the bodyForce for myPimpleFoam solver / Yuzhu (Pearl) Li

- Open a terminal 

- Load OpenFOAM v1912 (for example when using DTU cluster) and create user directory by 

      module load OpenFoam/v1912/gcc-9.2.0-openmpi-3.1.5
      cd $WM_PROJECT_USER_DIR
      mkdir -p applications/solvers
      mkdir applications/utilities

- Copy bodyForce folder to applications/utilities
  
- Copy myPimpleFoam folder to applications/solvers

- Compile bodyForce library by:

      cd $WM_PROJECT_USER_DIR/applications/utilities/bodyForce
      wclean
      wmake libso


- Then compile myPimpleFoam solver by:

      cd $WM_PROJECT_USER_DIR/applications/solvers/myPimpleFoam
      wclean
      wmake 

- Now you can use myPimpleFoam with specifying bodyForce in a case


Here is a potential error could be encountered in OpenFOAM newer version:

Error: 

correctPhi.H:2:1: error: expected constructor, destructor, or type conversion before ‘(’ token
(
^
 
 
Solution:

Modify the main file(myPimpleFoam.C) by changing the

#include "CorrectPhi.H" line 

with:

#include "path/where/CorrectPhi/is/located", (e.g. in my case /usr/lib/openfoam/openfoam2006/src/finiteVolume/lnInclude/CorrectPhi.H)
