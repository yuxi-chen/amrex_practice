

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

## BoxArray 
BoxArray: BoxArray is a class in AMReX_BoxArray.H for storing a collection of Boxes on a single AMR level. In AMReX, BoxArray is a global data structure. It holds all the Boxes in a collection, even though a single process in a parallel run only owns some of the Boxes via domain decomposition. BA itself only contains a shared pointer to the BARef, which contains a Vector<Box>.  So the shallow copy of BA is fast. 
```cpp
class BoxArray{
    //! The data -- a reference-counted pointer to a Ref.
    std::shared_ptr<BARef> m_ref;
}

struct BARef{
    //! The data.
    Vector<Box> m_abox;
}

```

## DistributionMapping
```cpp
class DistributionMapping
{
    struct Ref{        
        Vector<int> m_pmap; //!< index array for all boxes
        Vector<int> m_index_array;  //!< index array for local boxes owned by the team
        std::vector<bool> m_ownership; //!< true ownership
    };
    
    //! The data -- a reference-counted pointer to a Ref.
    std::shared_ptr<Ref> m_ref;
}
```

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

## AMREX particle data structure 
```cpp
class ArrayOfStructs {
  ParticleVector m_data;
}

class StructOfArrays {
private:
  std::array<RealVector, NReal> m_rdata;
  std::array<IntVector, NInt> m_idata;
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt,
          template <class> class Allocator = DefaultAllocator>
struct ParticleTile {
  // In FLEKS, we assume tile_size = 1 so far. See FLEKS::Particles.cpp::line31
  AoS m_aos_tile;
  SoA m_soa_tile;
}

class ParticleContainer : public ParticleContainerBase {

  using ParticleTileType =
      ParticleTile<NStructReal, NStructInt, NArrayReal, NArrayInt, Allocator>;

  //! A single level worth of particles is indexed (grid id, tile id)
  //! for both SoA and AoS data.
  using ParticleLevel = std::map<std::pair<int, int>, ParticleTileType>;

private:
  Vector<ParticleLevel> m_particles;
}
```
## Read the mover method in FLEKS/src/Particles.cpp. Ignore the numerical details The goal is to understand the particle-mesh interaction operations. 

