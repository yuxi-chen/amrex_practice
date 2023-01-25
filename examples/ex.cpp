#include "ex.h"
#include "../AMREX/InstallDir/include/AMReX.H"
#include "../AMREX/InstallDir/include/AMRex_ParallelDescriptor.H"
#include "../AMREX/InstallDir/include/AMRex_ParticleContainer.H"
#include "../AMREX/InstallDir/include/AMRex_Particles.H"
#include "../AMREX/InstallDir/include/AMRex_Vector.H"
#include "AmrCoreAdv.H"
#include "MyVector.h"
#include <iostream>
using namespace amrex;
int main(int argc, char *argv[]) {
  amrex::Initialize(argc, argv);
  std::cout << "Number of processors = " << amrex::ParallelDescriptor::NProcs()
            << " rank = " << amrex::ParallelDescriptor::MyProc() << std::endl;
  /////////////////// Initialize ////////////////////////////
  Real time = 0.0;
  AmrCoreAdv core; // First calls geometry::seutp - sets up level 0 geometry
                   // lo/hi bound and periodicity
  // and then initamrmesh sets levels of geom dmap grids  nerrorbuff
  // blockingfactor maxgridsize refratio   might also be setting maxlevel
  {
    ParmParse pp;
    pp.query("time_step", core.time_step);
    pp.query("stop_time", core.stop_time);
  }
  int max_level = 0;
  {
    ParmParse pp("amr");
    pp.query("max_level", max_level);
  }
  core.t_new.resize(max_level + 1, 0.0);
  core.t_old.resize(max_level + 1, 0.0);
  core.istep.resize(max_level + 1, 0);
  core.phi_new.resize(max_level + 1);
  core.phi_old.resize(max_level + 1);
  // core.getsomething();
  int bc_lo[] = { BCType::int_dir, BCType::int_dir, BCType::int_dir };
  int bc_hi[] = { BCType::int_dir, BCType::int_dir, BCType::int_dir };
  core.bcs.resize(1);
  for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
    core.bcs[0].setLo(idim, bc_lo[idim]);
    core.bcs[0].setHi(idim, bc_hi[idim]);
  }
  core.InitFromScratch(time); // First calls MakeNewGrids which calls
                              // MakeBaseGrids(lev0) and then sets other levels
                              // grids too and then calls
  // MakeNewLevelFromScratch -   and  ErrorEst both are located in amrcoreadv
  // core.WritePlotFile();

  core.init_particles();
  core.AverageDown();
  amrex::Real dt = 0.1;
  {
    ParmParse pp;
    pp.query("dt", dt);
  }
  core.dt = dt;
  int st;
  {
    ParmParse pp;
    pp.query("steps", st);
  }
  for (int aa = 0; aa < st; ++aa) {
    core.setBC0();
    core.setBC1();
    core.Evolve();

    core.AverageDown();
  }
  core.WritePlotFile();
  //////////////////////////////////////////////////////////
  // // Print()<<aa[1];
  // //  Best practice: you are encouraged to find useful information from AMREX
  // //  online document, amrex-tutorials and the AMREX source code. But do not
  // //  directly copy and paste ANY source code! It will be much useful to
  // TYPE,
  // //  compile and debug the code.
  // // Task 0: Compile and run this file. See Makefile
  // // Task 1: Create an amrex Vector v1.
  // Vector<amrex::Real> v1;
  // // Task 2: Create a cell based Box bc1, a node based Box bn1. Convert bc1
  // to a
  // // node based box bn2.
  // // Task 3: Create a RealBox rb1 //
  // // Task 4: Create a 3D Gemetry object geom with 16x16x4 cells,
  // // xyzmin=[-4,-4,-1], xyzmax=[4,4,1], cartesian grid and periodic in all
  // // directions
  // // Task 5: Create a BoxArray ba with 16x16x4 cells. Chop into boxes with
  // 8x8x4
  // // cells. Print ba to the screen.
  // // Task 6: Create a DistributionMapping dm with Task 5 ba.
  // // Task 7: Create 3 FArrayBox object fb1, fb2 and fb3 with 16x16x4 cells
  // and 2
  // // variables per cell. a) Loop through all cells of fb1 and fb2, set their
  // // values. b) calculate fb3 = fb1 + fb2. c) loop through all cells of fb1,
  // fb2
  // // and fb3, print their cell indices and values to the screen.
  // // Task 8: a)Create a MultiFab mf with Task 5 ba, Task 6 dm, 2 variables
  // per
  // // cell, and 1 ghost cells. b) Loop through all physical cells of mf,
  // // with Task 4 geometry, set the first mf variable to be x, and the second
  // // variable to y. c) Now, print the value of all cells, including the
  // // ghost cells, to screen. You should see the ghost cell values are either
  // 0
  // // or not well defined. d) Find a proper MultiFab method to fill in the
  // ghost
  // // cells.
  // // Task 9 [optional]: redo Task 8 with multiple processors.
  // //=============== Particles =====================
  // // Task 10: Initialize a Particle object. Set particle initial velocity,
  // // location and ID.
  // Particle<1, 1> pp1;
  // // Task 11: Do NOT use ParticleContainer. Creat a list of Particles with
  // // std::vector or std::array or std::list. What is the difference between
  // // these three containers? Which do you think is the best for our purpose?
  // std::vector<amrex::Particle<1, 1> > pv1;
  // Particle<1, 1> p1;
  // Particle<1, 1> p2;
  // for (int i = 0; i <= 10; ++i)
  // {
  //   p1.pos(0) = i;
  //   p1.pos(1) = -i;
  //   pv1.push_back(p1);
  // }
  // //... I know what the difference is but cant say for sure which is better.
  // // Array wont work if particles are added or removed. List will be faster
  // for
  // // adding or removing particles but might be slower when we iterate through
  // // the particles.
  // // Task 11.1: loop through and print location of all particles with the
  // // traditional loop: for(int i = 0......)
  // for (int i = 0; i <= 10; ++i)
  // {
  //   Print() << std::endl
  //           << "POS:" << i << "  (" << pv1[i].pos(0) << "," << pv1[i].pos(1)
  //           << ")";
  // }
  // Print() << std::endl;
  // // Task 11.2: do Task 11.1 with vector iterator
  // for (auto i : pv1) {
  //   Print() << std::endl
  //           << "POS:"
  //           << "  (" << i.pos(0) << "," << i.pos(1) << ")";
  // }
  // // Task 11.3: We do not directly use std::vector too much, but it can help
  // us
  // // to get familar with c++ containers. The best way to lean vector is to
  // // design and implement an simple vector template. Here is a simple
  // //
  // implementation:https://www.geeksforgeeks.org/how-to-implement-our-own-vector-class-in-c/
  // // Think about the questions: 1. What is the difference between
  // // 'size' and 'capacity'? 2. What will happen if a vector is already
  // 'full',
  // // but push() is still called to add more more elements?
  // // std::vector is never full until memory runs out.
  // // The vector I implemented does the same so I do not know what does vector
  // // being "full" mean
  // // Task 11.4: Implement you own vector class MyClass. See instructions in
  // // MyVector.h
  // vectorClass<int> aa;
  // Print() << aa.capacity;
  // Print() << aa.current;
  // aa.push_back(1);
  // Print() << aa.capacity;
  // Print() << aa.current;
  // aa.push_back(2);
  // Print() << aa.capacity;
  // Print() << aa.current;
  // aa.push_back(3);
  // Print() << aa.capacity;
  // Print() << aa.current;
  // aa.push_back(4);
  // Print() << aa.capacity;
  // Print() << aa.current;
  // aa.push_back(5);
  // Print() << aa.capacity;
  // Print() << aa.current;
  // aa.print();
  // // Task 12: Create a ParticleContainer to store particles. Each particle
  // store
  // // x,y,z,vx,vy,vz,w and ID, where w is the weight or 'mass' of a particle
  // and
  // // it is a real number. Loop through all particles and set initial
  // conditions. int n = 4; Box domain({ 1, 1, 1 }, { n, n, n }); RealBox
  // real_box({ 0.0, 0.0, 0.0 }, { 4.0, 4.0, 4.0 }); Geometry geom(domain,
  // real_box, 0, { 1, 1, 1 }); BoxArray ba(domain); ba.maxSize(2);
  // DistributionMapping dm{ ba };
  // MultiFab mf(ba, dm, 1, 1);
  // int lev = 0;
  // ParticleContainer<1, 0> PC1(geom, dm, ba);
  // auto& particle_tiles1 = PC1.GetParticles(
  //     lev); //[std::make_pair(grid_id,tile_id)];  - This is as Vec
  //     ParticleLevel
  // int tnp = 0; // - This is a ParticleTile
  // for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
  //   auto& particle_tile1 =
  //       particle_tiles1[std::make_pair(mfi.index(), mfi.LocalTileIndex())];
  //   const Box& bx = mfi.validbox();
  //   FArrayBox& fab = mf[mfi];
  //   Array4<Real> const& a = fab.array();
  //   const auto lo = lbound(bx);
  //   const auto hi = ubound(bx);
  //   for (int k = lo.z; k <= hi.z; ++k) {
  //     for (int j = lo.x; j <= hi.x; ++j) {
  //       for (int i = lo.x; i <= hi.x; ++i) {
  //         tnp = tnp + 1;
  //         Particle<1, 0> p;
  //         p.id() = tnp;
  //         p.cpu() = ParallelDescriptor::MyProc();
  //         p.pos(0) = i - 0.5;
  //         p.pos(1) = j - 0.5;
  //         p.pos(2) = k - 0.5;
  //         particle_tile1.push_back(p);
  //         Print() << "particle = " << p << "\n";
  //       }
  //     }
  //   }
  // }
  // PC1.Redistribute();
  // for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
  //   auto& particle_tile1 = PC1.GetParticles(
  //       lev)[std::make_pair(mfi.index(), mfi.LocalTileIndex())];
  //   auto& particles1 = particle_tile1.GetArrayOfStructs();
  //   const Box& bx = mfi.validbox();
  //   FArrayBox& fab = mf[mfi];
  //   Array4<Real> const& a = fab.array();
  //   const auto lo = lbound(bx);
  //   const auto hi = ubound(bx);
  //   {
  //     Print() << std::endl;
  //     Print() << particles1.numParticles();
  //     Print() << std::endl;
  //   }
  // }
  // // ParticleVector is PODvector??
  // // PODVector<ParticleType, Allocator<ParticleType> >;
  // // Print()<<particle_tile[1];
  // // Print()<<
  // // Task 13: update particle locations with the corresponding velocities.
  // 1st
  // // order accuracy is good enough.
  amrex::Finalize();
}
