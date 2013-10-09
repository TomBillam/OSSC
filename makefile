SOURCEFILES = constants.f90 fftw_3_3_2.f90 types.f90 output.f90 evolution.f90 main.f90
LIBS = -L. -lhdf5_fortran -lhdf5 -lz -lessl
BINARY = prog.out
OPTFLAGS = -qsmp=omp -q64 -qarch=pwr7 -qtune=pwr7 -O3 -qhot

program : $(SOURCEFILES)
	$(FC) $(OPTFLAGS) $(SOURCEFILES) $(FFLAGS) -o $(BINARY) $(LDFLAGS) $(LIBS)

clean :
	rm -rf *.mod *.so *.out *.o *.exe

