!***************************************************************************************
!** Otago Superfluid Simulation Code
!** Two-dimensional Gross-Pitaevskii simulation using Fourier pseudospectral method in a
!** periodic, uniform domain, with 4/5 order adaptive Runge_Kutta timestepping.
!**
!** 2013  Thomas P. Billam, University of Otago
!***************************************************************************************

program Program2DGrossPitaevskiiHPC
use constants
use types
use output
use evolution
use, intrinsic :: iso_c_binding
use fftw_3_3_2
use omp_lib
use hdf5
implicit none

type (vector):: length, dk_sq
type (grid):: gd

double precision:: g, mu, gam, output_delta_t, goal, time
integer:: output_number, error, nthreads, input_file_number
integer(C_INT) :: ret

complex(C_DOUBLE_COMPLEX), pointer :: phi(:,:)
type(C_PTR) :: p

type (transform_pair):: phi_transforms

character(len=30):: filename
character(len=255,kind=C_CHAR):: wisdom_file_name
!** Command line character vars
character(len=255):: cmd_string

double precision:: prog_run_time, prog_start_time

!** PARSE COMMAND LINE
if (command_argument_count() .lt. 8) then
  write(*,*) 'Usage: program nthreads tol num_steps output_delta_t wisdom_filename input_file_number prog_run_time gamma'
  call get_command(cmd_string)
  write(*,*) cmd_string
  stop
end if
call get_command_argument(1,cmd_string)
 cmd_string = trim(cmd_string)
  read(cmd_string,*) nthreads
call get_command_argument(2,cmd_string)
 cmd_string = trim(cmd_string)
  read(cmd_string,*) goal
call get_command_argument(3,cmd_string)
 cmd_string = trim(cmd_string)
  read(cmd_string,*) output_number
call get_command_argument(4,cmd_string)
 cmd_string = trim(cmd_string)
  read(cmd_string,*) output_delta_t
call get_command_argument(5,cmd_string)
 cmd_string = trim(cmd_string)
  read(cmd_string,*) wisdom_file_name
call get_command_argument(6,cmd_string)
 cmd_string = trim(cmd_string)
  read(cmd_string,*) input_file_number
call get_command_argument(7,cmd_string)
 cmd_string = trim(cmd_string)
  read(cmd_string,*) prog_run_time
call get_command_argument(8,cmd_string)
 cmd_string = trim(cmd_string)
  read(cmd_string,*) gam

prog_start_time = omp_get_wtime()

!** Initialize hdf5
call h5open_f(error)

!** Make input filename
write(filename,*) input_file_number
filename = trim(adjustl(filename))//'.h5'
write(*,*) filename

!** Get grid params from input file attributes
call get_grid_from_file(filename, gd)

length%x = dble(gd%n%x)*gd%delta%x
length%y = dble(gd%n%y)*gd%delta%y
dk_sq=vector((2.0d0*pi/length%x)**2,(2.0d0*pi/length%y)**2)

call omp_set_num_threads(nthreads)
!$OMP PARALLEL
nthreads = omp_get_num_threads()
!$OMP END PARALLEL
write(*,*) 'Running on ', nthreads, 'cores'

!** Use openmp on ffts
ret = fftw_init_threads()
call fftw_plan_with_nthreads(nthreads)

!** Allocate aligned array for fftw
p = fftw_alloc_complex(int(gd%n%x * gd%n%y, C_SIZE_T))
       call c_f_pointer(p, phi, [gd%n%x,gd%n%y])

!** Try to find wisdom file
wisdom_file_name = trim(adjustl(wisdom_file_name))
ret = fftw_import_wisdom_from_filename( trim(wisdom_file_name)//C_NULL_CHAR)
if (ret .eq. 0) then
  write(*,*) 'error importing wisdom from specified file'
end if
write(*,*) 'Planning fast Fourier transforms...'
phi_transforms%forward = fftw_plan_dft_2d(gd%n%y,gd%n%x, phi, phi, FFTW_FORWARD, FFTW_PATIENT)
phi_transforms%backward = fftw_plan_dft_2d(gd%n%y,gd%n%x, phi, phi, FFTW_BACKWARD, FFTW_PATIENT)
write(*,*) '...done.'
ret = fftw_export_wisdom_to_filename( trim(wisdom_file_name)//C_NULL_CHAR)
       if (ret .eq. 0) stop 'error exporting wisdom to file'

!** LOAD INPUT FILE
call read_hdf5_named(filename,phi,gd,time)
write(*,*) 'Loaded input file'


!***** REAL TIME EVOLUTION *****
!** Set system parameters for real time evolution
g = 1.0d0
mu = 1.0d0

!** Evolve damped GPE
call evolve_gpe_adaptive(phi, gd, phi_transforms, dk_sq, output_number, output_delta_t, &
            write_hdf5, goal, g, mu, gam, input_file_number, time, prog_start_time, prog_run_time)

!** Tidy up **
call fftw_destroy_plan(phi_transforms%forward)
call fftw_destroy_plan(phi_transforms%backward)
call fftw_free(p)

! Close FORTRAN interface.
call h5close_f(error)

end program 
