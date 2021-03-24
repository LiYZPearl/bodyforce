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
