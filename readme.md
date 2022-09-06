# Install AMREX lib
1. git clone git@github.com:yuxi-chen/amrex_practice.git
2. cd amrex_practice; gitlabclone AMREX
3. Configure AMREX:
   * cd AMREX; 
   * ./configure --prefix InstallDir #Choose InstallDir as the install directory
4. Install AMREX:
   * cd AMREX
   * make -j 8; make install  # then amrex lib will be installed in amrex_practice/AMREX/InstallDir

5. Explor the InstallDir directory. It contains two directories: include and lib. include contains all the headers and lib contains the library libamrex.a. 

# Get familar with basic AMREX data types. 
Finish the Task 0-9 described in amrex_practice/examples/ex.cpp

# Particle
## Read amrex online document 'Particles' chapter. 
Skip the 'Adding particle components at runtime', 'Passing particle data into Fortran routines', 'Short Range Forces', 'Short Range Forces' and 'Inputs parameters' parts. Think about the following questions:
1. The memory layout of the Particle and ParticleContainer classes. 
2. What is the difference between Arrays-of-Structs and Structs-of-Arrays?
## Coding tasks. Finish Task 10-13 in ex.cpp
Requirement: 
* Do NOT copy and paste any code from amrex_tutorials or online document. Do NOT even use the same variable names as the tutorials.
* Do NOT use lambdas, such as AMREX_FOR_1D or ParallelFor. Explicitly use for loops, instead.
* Only keep the code that is necessary and push it to github after finished.  

## Read the mover method in FLEKS/src/Particles.cpp. Ignore the numerical details The goal is to understand the particle-mesh interaction operations. 

