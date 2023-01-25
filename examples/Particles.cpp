#include "Particles.h"

using namespace amrex;

template <int NStructReal, int NStructInt>
void Particles<NStructReal, NStructInt>::init_particles() {
  printf("init_particles\n");

  for (int lev = 0; lev <= maxLevel(); lev++) {
    printf("lev=%d\n", lev);
    for (MFIter mfi = MakeMFIter(lev, false); mfi.isValid(); ++mfi) {
      const Box& tile_box = mfi.validbox();
      Print() << "box = " << tile_box << std::endl;
      const auto lo = amrex::lbound(tile_box);
      const auto hi = amrex::ubound(tile_box);

      int iMax = hi.x, jMax = hi.y, kMax = hi.z;
      int iMin = lo.x, jMin = lo.y, kMin = lo.z;

      for (int i = iMin; i <= iMax; ++i)
        for (int j = jMin; j <= jMax; ++j)
          for (int k = kMin; k <= kMax; ++k) {
            // printf("i=%d,j=%d,k=%d\n",i,j,k);
          }
    }
  }
}

template class Particles<nPicPartReal>;