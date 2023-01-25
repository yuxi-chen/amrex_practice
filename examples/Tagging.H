#ifndef TAGGING_H
#define TAGGING_H

#include <AMReX_Array4.H>
#include <AMReX_ParmParse.H>


AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
state_error (int i, int j, int k,
             amrex::Array4<char> const& tag,
             amrex::Array4<amrex::Real const> const& state,
             amrex::Real phierr, char tagval)
{
    //if (state(i,j,k) > phierr)

    // {
    //     int tagn;
    //     amrex::ParmParse pp; 
    //     pp.query("tag", tagn);
    //     tag(i,j,k) = tagn;
    // }

    int tagn;
    amrex::ParmParse pp; 
   pp.query("tag", tagn);

       if(i>=3&&i<=4)
       {
        if(j>=3&&j<=4)
        {

         

           tag(i,j,k) = tagn;


        
        }

        }
}

#endif
