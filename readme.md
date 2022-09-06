Task 1: Install AMREX lib
1. git clone git@github.com:yuxi-chen/amrex_practice.git
2. cd amrex_practice
   gitlabclone AMREX
3. Configure AMREX:
   cd AMREX
   # Choose InstallDir as the install directory
   ./configure --prefix InstallDir
4. Install AMREX:
   cd AMREX
   make -j 8; make install  # then amrex lib will be installed in amrex_practice/AMREX/InstallDir

5. Explor the InstallDir directory. It contains two directories: include and lib. include contains all the headers and lib contains the library libamrex.a. 

Task 2: Finish the tasks described in amrex_practice/examples/ex.cpp



