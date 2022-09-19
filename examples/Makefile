default: 
	mpicxx  -I../AMREX/InstallDir/include -c ex.cpp -o ex.o
	mpif90 -o ex.exe ex.o -L../AMREX/InstallDir/lib -lamrex  -lstdc++ -lmpi_cxx

# run with 4 processors by default. It can be changed to any number. 
run:
	mpiexec -np 1 ./ex.exe

clean:
	rm *.exe *.o