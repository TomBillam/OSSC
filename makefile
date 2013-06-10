SOURCEFILES = constants.f90 fftw_3_3_2.f90 types.f90 output.f90 evolution.f90 main.f90
LIBS = -L. -lfftw3_omp -lfftw3 -lhdf5_fortran -lhdf5 -lz
BINARY = prog.out

program : $(SOURCEFILES)
	$(FC) $(OPTFLAGS) $(SOURCEFILES) $(FFLAGS) -o $(BINARY) $(LDFLAGS) $(LIBS)

clean :
	rm -rf *.mod *.so *.out *.o *.exe

