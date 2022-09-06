#include "../AMREX/InstallDir/include/AMReX.H"
#include "../AMREX/InstallDir/include/AMRex_ParallelDescriptor.H"
#include <iostream>

#include "ex.h"

int main(int argc, char* argv[]) {
  amrex::Initialize(argc, argv);

  std::cout << "Number of processors = " << amrex::ParallelDescriptor::NProcs()
            << " rank = " << amrex::ParallelDescriptor::MyProc() << std::endl;

  // Best practice: you are encouraged to find useful information from AMREX
  // online document, amrex-tutorials and the AMREX source code. But do not
  // directly copy and paste ANY source code! It will be much useful to TYPE,
  // compile and debug the code.

  // Task 0: Compile and run this file. See Makefile

  // Task 1: Create an amrex Vector v1.

  // Task 2: Create a cell based Box bc1, a node based Box bn1. Convert bc1 to a
  // node based box bn2.

  // Task 3: Create a RealBox rb1

  // Task 4: Create a 3D Gemetry object geom with 16x16x4 cells,
  // xyzmin=[-4,-4,-1], xyzmax=[4,4,1], cartesian grid and periodic in all
  // directions

  // Task 5: Create a BoxArray ba with 16x16x4 cells. Chop into boxes with 8x8x4
  // cells. Print ba to the screen.

  // Task 6: Create a DistributionMapping dm with Task 5 ba.

  // Task 7: Create 3 FArrayBox object fb1, fb2 and fb3 with 16x16x4 cells and 2
  // variables per cell. a) Loop through all cells of fb1 and fb2, set their
  // values. b) calculate fb3 = fb1 + fb2. c) loop through all cells of fb1, fb2
  // and fb3, print their cell indices and values to the screen.

  // Task 8: a)Create a MultiFab mf with Task 5 ba, Task 6 dm, 2 variables per
  // cell, and 1 ghost cells. b) Loop through all physical cells of mf,
  // with Task 4 geometry, set the first mf variable to be x, and the second
  // variable to y. c) Now, print the value of all cells, including the
  // ghost cells, to screen. You should see the ghost cell values are either 0
  // or not well defined. d) Find a proper MultiFab method to fill in the ghost
  // cells.

  // Task 9 [optional]: redo Task 8 with multiple processors.

  //=============== Particles =====================
  // Task 10: Initialize a Particle object. Set particle initial velocity, location and ID. 

  // Task 11: Do NOT use ParticleContainer. Creat a list of Particles with std::vector or std::array or std::list. What is the difference between these three containers? Which do you think is the best for our purpose? 

  // Task 12: Create a ParticleContainer to store particles. Each particle store x,y,z,vx,vy,vz,w and ID, where w is the weight or 'mass' of a particle and it is a real number. Loop through all particles and set initial conditions. 

  // Task 13: update particle locations with the corresponding velocities. 1st order accuracy is good enough. 


  amrex::Finalize();
}
