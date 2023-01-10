#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>
#include <AMReX_PhysBCFunct.H>
#include "Tagging.H"
#include "AmrCoreAdv.H"
using namespace amrex;
AmrCoreAdv::AmrCoreAdv()
{
}
AmrCoreAdv::~AmrCoreAdv()
{
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////
void AmrCoreAdv::MakeNewLevelFromScratch(int lev, Real time, const BoxArray &ba,
                                         const DistributionMapping &dm)
{
    const int ncomp = 1;
    const int nghost = 1;
    phi_new[lev].define(ba, dm, ncomp, nghost);
    phi_old[lev].define(ba, dm, ncomp, nghost);
    t_new[lev] = time;
    t_old[lev] = time - 1.e200;
    MultiFab &state = phi_new[lev];
    const auto problo = Geom(lev).ProbLoArray();
    const auto dx = Geom(lev).CellSizeArray();
    for (MFIter mfi(state); mfi.isValid(); ++mfi)
    {
        Array4<Real> fab = state[mfi].array();
        const Box &box = mfi.tilebox();
        const auto lo = lbound(box);
        const auto hi = ubound(box);
        for (int k = lo.z; k <= hi.z; ++k)
        {
            for (int j = lo.y; j <= hi.y; ++j)
            {
                Real y = problo[1] + (0.5 + j) * dx[1];
                for (int i = lo.x; i <= hi.x; ++i)
                {
                    Real x = problo[0] + (0.5 + i) * dx[0];
                    Real r2 = (std::pow(x - 0.75, 2) + std::pow((y - 0.5), 2)) / 0.01;
                    fab(i, j, k) = 1.0 + std::exp(-r2);
                }
            }
        }
    }
}
void AmrCoreAdv::ErrorEst(int lev, amrex::TagBoxArray &tags, amrex::Real time, int ngrow)
{
    amrex::Vector<Real> phierr;
    ParmParse pp;
    pp.getarr("phierr", phierr);
    const int tagval = TagBox::SET;
    const MultiFab &state = phi_new[lev];
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
        for (MFIter mfi(state); mfi.isValid(); ++mfi)
        {
            const Box &bx = mfi.tilebox();
            const auto statefab = state.array(mfi);
            const auto tagfab = tags.array(mfi);
            Real phierror = phierr[lev];

            const auto lo = lbound(bx);
            const auto hi = ubound(bx);

            for (int k = lo.z; k <= hi.z; ++k)
            {
                for (int j = lo.y; j <= hi.y; ++j)
                {

                    for (int i = lo.x; i <= hi.x; ++i)
                    {

                        state_error(i, j, k, tagfab, statefab, phierror, tagval);
                    }
                }
            }
        }
    }
}
void AmrCoreAdv::MakeNewLevelFromCoarse(int lev, amrex::Real time, const amrex::BoxArray &ba,
                                        const amrex::DistributionMapping &dm)
{
}
void AmrCoreAdv::RemakeLevel(int lev, amrex::Real time, const amrex::BoxArray &ba,
                             const amrex::DistributionMapping &dm)
{
}
void AmrCoreAdv::ClearLevel(int lev)
{
}
////////////////////////////////////////////////////////////////////////////////////////////////////////
void AmrCoreAdv::WritePlotFile() const
{
    const std::string &plotfilename = PlotFileName(istep[0]);
    const auto &mf = PlotFileMF();
    const auto &varnames = PlotFileVarNames();
    amrex::Print() << "Writing plotfile " << plotfilename << "\n";
    amrex::WriteMultiLevelPlotfile(plotfilename, finest_level + 1, mf, varnames,
                                   Geom(), t_new[0], istep, refRatio());
}
std::string
AmrCoreAdv::PlotFileName(int lev) const
{
    return amrex::Concatenate("plt", lev, 5);
}
Vector<const MultiFab *>
AmrCoreAdv::PlotFileMF() const
{
    Vector<const MultiFab *> r;
    for (int i = 0; i <= finest_level; ++i)
    {
        r.push_back(&phi_new[i]);
    }
    return r;
}
Vector<std::string>
AmrCoreAdv::PlotFileVarNames() const
{
    return {"phi"};
}
void AmrCoreAdv::AverageDown()
{
    for (int lev = finest_level - 1; lev >= 0; --lev)
    {
        amrex::average_down(phi_new[lev + 1], phi_new[lev],
                            geom[lev + 1], geom[lev],
                            0, phi_new[lev].nComp(), refRatio(lev));
    }
}
void AmrCoreAdv::Evolve() 
{
    
    for (int ii=0; ii<2; ++ii)

    {

    MultiFab &state = phi_new[ii];
    MultiFab &state2 = phi_old[ii];
    const auto problo = Geom(ii).ProbLoArray();
    const auto dx = Geom(ii).CellSizeArray();
    for (MFIter mfi(state); mfi.isValid(); ++mfi)
    {
        Array4<Real> fab = state[mfi].array();
        Array4<Real> fab2 = state2[mfi].array();
        const Box &box = mfi.tilebox();
        const auto lo = lbound(box);
        const auto hi = ubound(box);
        for (int k = lo.z; k <= hi.z; ++k)
        {
            for (int j = lo.y; j <= hi.y; ++j)
            {
                for (int i = lo.x; i <= hi.x; ++i)
                {
                    fab2(i, j, k) = dt * (fab(i - 1, j, k) - 2 * fab(i, j, k) + fab(i + 1, j, k) + fab(i, j - 1, k) - 2 * fab(i, j, k) + fab(i, j + 1, k)) / (dx[0] * dx[0]) + fab(i, j, k);
                }
            }
        }
        for (int k = lo.z; k <= hi.z; ++k)
        {
            for (int j = lo.y; j <= hi.y; ++j)
            {
                for (int i = lo.x; i <= hi.x; ++i)
                {
                    fab(i, j, k) = fab2(i, j, k);
                }
            }
        }
    }


    }
}


void AmrCoreAdv::setBC1()
{
    MultiFab Sborder(grids[1], dmap[1], phi_new[1].nComp(), 1); 
    FillPatch(1, 0.0, phi_new[1], 0, phi_new[1].nComp());
    phi_new[1].FillBoundary();
}
void AmrCoreAdv::FillPatch(int lev, Real time, MultiFab &mf, int icomp, int ncomp)
{
    // if (lev == 0)
    // {
    //     Vector<MultiFab *> smf;
    //     Vector<Real> stime;
    //     GetData(0, time, smf, stime);
    //     CpuBndryFuncFab bndry_func(nullptr); // Without EXT_DIR, we can pass a nullptr.
    //     PhysBCFunct<CpuBndryFuncFab> physbc(geom[lev], bcs, bndry_func);
    //     amrex::FillPatchSingleLevel(mf, time, smf, stime, 0, icomp, ncomp,
    //                                 geom[lev], physbc, 0);
    // }
    // else
    {

        
        Vector<MultiFab *> cmf, fmf;
        Vector<Real> ctime, ftime,cc,ff;
        cmf.push_back(&phi_new[0]);
        fmf.push_back(&phi_new[1]);
        
        Interpolater *mapper = &cell_cons_interp;
        CpuBndryFuncFab bndry_func(nullptr); // Without EXT_DIR, we can pass a nullptr.
        PhysBCFunct<CpuBndryFuncFab> cphysbc(geom[lev - 1], bcs, bndry_func);
        PhysBCFunct<CpuBndryFuncFab> fphysbc(geom[lev], bcs, bndry_func);
        
    
        amrex::FillPatchTwoLevels(mf, time, cmf, cc, fmf, ff,
                                  0, icomp, ncomp, geom[lev - 1], geom[lev],
                                  cphysbc, 0, fphysbc, 0, refRatio(lev - 1),
                                  mapper, bcs, 0);
    }
}
