#ifndef AmrCoreAdv_H_
#define AmrCoreAdv_H_

#include <limits>
#include <memory>
#include <string>

#ifdef AMREX_USE_OMP
#include <omp.h>
#endif

#include <AMReX_AmrCore.H>
#include <AMReX_BCRec.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_FluxRegister.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_Print.H>
#include <AMReX_RealBox.H>
#include <AMReX_Vector.H>

#include "Particles.h"

class AmrCoreAdv : public amrex::AmrCore {
public:
  AmrCoreAdv();

  virtual ~AmrCoreAdv();

  ////////////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  virtual void MakeNewLevelFromCoarse(
      int lev, amrex::Real time, const amrex::BoxArray &ba,
      const amrex::DistributionMapping &dm) override;

  virtual void RemakeLevel(int lev, amrex::Real time, const amrex::BoxArray &ba,
                           const amrex::DistributionMapping &dm) override;

  virtual void ClearLevel(int lev) override;

  virtual void MakeNewLevelFromScratch(
      int lev, amrex::Real time, const amrex::BoxArray &ba,
      const amrex::DistributionMapping &dm) override;

  virtual void ErrorEst(int lev, amrex::TagBoxArray &tags, amrex::Real time,
                        int ngrow) override;

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  void getsomething() { amrex::Print() << max_level; }
  void setBC0() { phi_new[0].FillBoundary(geom[0].periodicity()); }
  void setBC1();
  void FillPatch(int lev, amrex::Real time, amrex::MultiFab &mf, int icomp,
                 int ncomp);

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  int nSpecies;
  amrex::Vector<Particles<>> pts;

  amrex::Real time_step;
  amrex::Real stop_time;
  amrex::Vector<amrex::Real> ts;
  amrex::Vector<amrex::Real> t_new;
  amrex::Vector<amrex::Real> t_old;
  amrex::Real dt;

  // array of multifabs to store the solution at each level of refinement
  // after advancing a level we use "swap".
  amrex::Vector<amrex::MultiFab> phi_new;
  amrex::Vector<amrex::MultiFab> phi_old;
  amrex::Vector<amrex::BCRec> bcs;

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // get plotfile name
  std::string PlotFileName(int lev) const;

  // put together an array of multifabs for writing
  amrex::Vector<const amrex::MultiFab *> PlotFileMF() const;

  // set plotfile variables names
  amrex::Vector<std::string> PlotFileVarNames() const;

  // write plotfile to disk
  void WritePlotFile() const;
  void AverageDown();
  void Evolve();

  amrex::Vector<int> istep;

private:
  //////////////////////////////////////////////

  //////////////////////////////////////////////////

  //////////////////////////////////////////////////

  ///////////////////////////////////////////////////
};

#endif
