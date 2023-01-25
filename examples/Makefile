default: 
	mpicxx  -I../AMREX/InstallDir/include -c -g -O0 ex.cpp -o ex.o
	mpicxx  -I../AMREX/InstallDir/include -c -g -O0 AmrCoreAdv.cpp -o AmrCoreAdv.o

	mpif90 -o ex.exe ex.o AmrCoreAdv.o -L../AMREX/InstallDir/lib -lamrex  -lstdc++ -lmpi_cxx

# run with 4 processors by default. It can be changed to any number. 
run:
	mpiexec -np 1 ./ex.exe inputs

clean:
	rm *.exe *.o