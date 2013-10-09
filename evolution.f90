!***************************************************************************************
!** Otago Superfluid Simulation Code
!** Two-dimensional Gross-Pitaevskii simulation using Fourier pseudospectral method in a
!** periodic, uniform domain, with 4/5 order adaptive Runge_Kutta timestepping.
!**
!** 2013  Thomas P. Billam, University of Otago
!***************************************************************************************

module evolution
use, intrinsic:: iso_c_binding
implicit none
save
private

public:: evolve_gpe_adaptive
double precision, allocatable:: karray(:,:)

!** ESSL FFT arrays and parameters.
!** Moved here from main since the planning overhead is so much smaller than for fftw
!** TODO: drop the fftw plans from all the function calls
double precision, allocatable:: aux1fwd(:), aux1bwd(:), aux2(:)
integer:: naux1, naux2
integer:: narr(2)
integer, parameter:: isign_fwd(2) = (/1,1/), isign_bwd(2) = (/-1,-1/)
integer, parameter:: incx(2) = (/1,1/), incy(2) = (/1,1/)

!** Also make the working arrays module based: this may make a tiny performance saving over local arrays,
!** plus it makes planning the transforms simpler
double complex, allocatable:: phi_old(:,:)
double complex, allocatable:: phi2(:,:)
double complex, allocatable:: k1(:,:)
double complex, allocatable:: k2(:,:)
double complex, allocatable:: k3(:,:)
double complex, allocatable:: k4(:,:)
double complex, allocatable:: k5(:,:)
double complex, allocatable:: k6(:,:)

!** Count Number of Iterations
integer:: iterationcount = 0

contains


subroutine evolve_gpe_adaptive(phi, gd, phi_transforms, dk_sq, output_number, output_delta_t, &
                   output_sub, goal, g, mu, gam, input_file_number, time, prog_start_time, prog_run_time, tol)
use, intrinsic:: iso_c_binding
use types
use output
use omp_lib
implicit none
type (grid), intent(in):: gd
double complex, intent(inout):: phi(0:gd%n%x-1,0:gd%n%y-1)
type (transform_pair), intent(in):: phi_transforms
type (vector), intent(in):: dk_sq
integer, intent(in):: output_number, input_file_number
double precision, intent(in):: output_delta_t, goal, g, mu, gam
double precision, intent(inout)::time
double precision, intent(in):: prog_start_time, prog_run_time
double precision, optional, intent(in):: tol
external:: output_sub
double precision:: delta_t
integer:: i, j
double precision:: step_start, step_finish, last_step_time, prog_time_remaining
character(len=30):: currentstepfilename
double precision:: evolution_start, evolution_finish
currentstepfilename = 'CurrentStep.h5'

!** Set momentum array for speedup
if (allocated(karray)) then
  deallocate(karray) !TODO: THIS IS POTENTIALLY INEFFIECIENT
end if
allocate(karray(0:gd%n%x-1,0:gd%n%y-1))
!$OMP PARALLEL DO private(i)
do j = 0, gd%n%y-1
  do i = 0, gd%n%x-1
    karray(i,j) = momentum_factor(i,j,gd%n,dk_sq,mu)
  end do
end do
!$OMP END PARALLEL DO

!** Set up ESSL transforms
naux1 = ceiling(60000.0d0*dble(2) + 28.24d0 *dble(gd%n%x+gd%n%y))
naux2 = ceiling(20000.0d0+dble(2*max(gd%n%x,gd%n%y)+256)*(64.0d0+17.12))
narr(1) = gd%n%x
narr(2) = gd%n%y
allocate(aux1fwd(naux1))
allocate(aux1bwd(naux1))
allocate(aux2(naux2))
allocate(phi_old(0:gd%n%x-1,0:gd%n%y-1))
allocate(phi2(0:gd%n%x-1,0:gd%n%y-1))
allocate(k1(0:gd%n%x-1,0:gd%n%y-1))
allocate(k2(0:gd%n%x-1,0:gd%n%y-1))
allocate(k3(0:gd%n%x-1,0:gd%n%y-1))
allocate(k4(0:gd%n%x-1,0:gd%n%y-1))
allocate(k5(0:gd%n%x-1,0:gd%n%y-1))
allocate(k6(0:gd%n%x-1,0:gd%n%y-1))
!** N.B.: As I understand it, ESSL doesn't care which arrays I plan with, as long as I pass
!** it arrays which have been allocated in the same fashion (and stick to either in-place or out-of-place as appropriate)
CALL DCFTD (1, 2, phi2, incx, 1, phi, incy, 1, narr, 1, isign_fwd, 1.0d0, aux1fwd, naux1, aux2, naux2)
CALL DCFTD (1, 2, phi, incx, 1, phi2, incy, 1, narr, 1, isign_bwd, 1.0d0, aux1bwd, naux1, aux2, naux2)


!** Set initial timestep and time
delta_t = 0.01d0*cfl_timestep(gd%delta)

if (present(tol)) then
  !** Evolve until convergence
  call step_adaptive(phi, gd, phi_transforms, time, time+10.0d0*delta_t, delta_t, goal, g, gam, tol)
  !
else
  !** Evolve for set number of steps
  evolution_start = omp_get_wtime()
  do i = 1, output_number
    step_start = omp_get_wtime()
    call step_adaptive(phi, gd, phi_transforms, time, time+output_delta_t, delta_t, goal, g, gam)
    step_finish = omp_get_wtime()
    write(*,*)'Just reached output breakpoint ', i, ' ...'
    call output_sub(input_file_number+i,phi,gd,time)
    last_step_time = step_finish - step_start
    prog_time_remaining = prog_run_time + prog_start_time - step_finish
    write(*,*) 'Last step took', last_step_time, 'secs, and there are', prog_time_remaining, 'secs remaining...'
    if (last_step_time*1.5 .gt. prog_time_remaining) then
      !** Stop already if we can't fit another step in...
      write(*,*) 'Not sure there is enough time for another step: exiting gracefully now'
      exit
    end if
  end do
  evolution_finish = omp_get_wtime()
  write(*,*) 'Did', iterationcount, 'iterations in', evolution_finish-evolution_start, 'seconds'
end if

deallocate(phi_old, k1, k2, k3, k4, k5, k6, phi2)

end subroutine


subroutine step_adaptive(phi, gd, pt, time, end_time, delta_t, goal, g, gam, tol)
use, intrinsic:: iso_c_binding
use types
use output
implicit none
type (grid), intent(in):: gd
double complex, intent(inout):: phi(0:gd%n%x-1,0:gd%n%y-1)
type (transform_pair), intent(in):: pt
double precision, intent(inout):: time, delta_t
double precision, intent(in):: end_time, goal, g, gam
double precision, optional, intent(in):: tol
logical:: can_has_exit !** Keep track of if we might make a successful final timestep
double precision:: error, start_time, conv
integer:: i,j
double precision:: norm, old_norm, newmax, oldmax
double precision, parameter:: error_cut_off = 1.0d-3
double precision, parameter:: b21 = 1.0d0/5.0d0
double precision, parameter:: b31 = 3.0d0/40.0d0
double precision, parameter:: b32 = 9.0d0/40.0d0
double precision, parameter:: b41 = 3.0d0/10.0d0
double precision, parameter:: b42 = -9.0d0/10.0d0
double precision, parameter:: b43 = 6.0d0/5.0d0
double precision, parameter:: b51 = -11.0d0/54.0d0
double precision, parameter:: b52 = 5.0d0/2.0d0
double precision, parameter:: b53 = -70.0d0/27.0d0
double precision, parameter:: b54 = 35.0d0/27.0d0
double precision, parameter:: b61 = 1631.0d0/55296.0d0
double precision, parameter:: b62 = 175.0d0/512.0d0
double precision, parameter:: b63 = 575.0d0/13824.0d0
double precision, parameter:: b64 = 44275.0d0/110592.0d0
double precision, parameter:: b65 = 253.0d0/4096.0d0
double precision, parameter:: c1 = 37.0d0/378.0d0
double precision, parameter:: c2 = 0.0d0
double precision, parameter:: c3 = 250.0d0/621.0d0
double precision, parameter:: c4 = 125.0d0/594.0d0
double precision, parameter:: c5 = 0.0d0
double precision, parameter:: c6 = 512.0d0/1771.0d0
double precision, parameter:: d1 = 2825.0d0/27648.0d0
double precision, parameter:: d2 = 0.0d0
double precision, parameter:: d3 = 18575.0d0/48384.0d0
double precision, parameter:: d4 = 13525.0d0/55296.0d0
double precision, parameter:: d5 = 277.0d0/14336.0d0
double precision, parameter:: d6 = 1.0d0/4.0d0

!** Increment iteration count
iterationcount = iterationcount + 1

old_norm = 0.0d0
oldmax = 0.0d0
!$OMP PARALLEL DO private(i) reduction(+:old_norm) reduction(MAX:oldmax)
do j = 0, gd%n%y-1
  do i = 0, gd%n%x-1
    old_norm = old_norm + abs(phi(i,j))**2
    if(abs(phi(i,j)) .gt. oldmax) then
      oldmax = abs(phi(i,j))
    end if
  end do
end do
!$OMP END PARALLEL DO
old_norm = old_norm*gd%delta%x*gd%delta%y
start_time = time

do
  if (.not. present(tol)) then
    !** In real time, check if we might be able to make a successful final timestep
    if (end_time-time .le. delta_t) then
      can_has_exit = .true.
    else
      can_has_exit = .false.
    end if
  end if
  
  !$OMP PARALLEL DO private(i)
  do j = 0, gd%n%y-1
    do i =0, gd%n%x-1
     !** Save the old wavefunction in case we have to go back
      phi_old(i,j) = phi(i,j)
    end do
  end do
  !$OMP END PARALLEL DO
  
  !** Obtain k1 coefficients 
  CALL DCFTD (0, 2, phi_old, incx, 1, phi, incy, 1, narr, 1, isign_fwd, 1.0d0, aux1fwd, naux1, aux2, naux2)
  !$OMP PARALLEL DO private(i)
  do j = 0, gd%n%y-1
    do i = 0, gd%n%x-1
      phi(i,j) = karray(i,j) * phi(i,j) / (gd%n%x*gd%n%y)
    end do
  end do
  !$OMP END PARALLEL DO
  CALL DCFTD (0, 2, phi, incx, 1, phi2, incy, 1, narr, 1, isign_bwd, 1.0d0, aux1bwd, naux1, aux2, naux2)
  !$OMP PARALLEL DO private(i)
  do j = 0, gd%n%y-1
    do i = 0, gd%n%x-1
      k1(i,j) = -(0.0d0,1.0d0)*dcmplx(1.0d0,-gam)*delta_t * (phi2(i,j) &
                + (g*abs(phi_old(i,j))**2)*phi_old(i,j))
      k2(i,j) = phi_old(i,j) + b21*k1(i,j)
    end do
  end do
  !$OMP END PARALLEL DO

  !** Obtain k2 coefficients
  CALL DCFTD (0, 2, k2, incx, 1, phi, incy, 1, narr, 1, isign_fwd, 1.0d0, aux1fwd, naux1, aux2, naux2)
  !$OMP PARALLEL DO private(i)
  do j = 0, gd%n%y-1
    do i = 0, gd%n%x-1
      phi(i,j) = karray(i,j) * phi(i,j) / (gd%n%x*gd%n%y)
    end do
  end do
  !$OMP END PARALLEL DO
  CALL DCFTD (0, 2, phi, incx, 1, phi2, incy, 1, narr, 1, isign_bwd, 1.0d0, aux1bwd, naux1, aux2, naux2)
  !$OMP PARALLEL DO private(i)
  do j = 0, gd%n%y-1
    do i = 0, gd%n%x-1
      k2(i,j) = -(0.0d0,1.0d0)*dcmplx(1.0d0,-gam)*delta_t * (phi2(i,j) &
                + (g*abs(k2(i,j))**2)*k2(i,j))
      k3(i,j) = phi_old(i,j) + b31*k1(i,j) + b32*k2(i,j)
    end do
  end do
  !$OMP END PARALLEL DO

  !** Obtain k3 coefficients
  CALL DCFTD (0, 2, k3, incx, 1, phi, incy, 1, narr, 1, isign_fwd, 1.0d0, aux1fwd, naux1, aux2, naux2)
  !$OMP PARALLEL DO private(i)
  do j = 0, gd%n%y-1
    do i = 0, gd%n%x-1
      phi(i,j) = karray(i,j) * phi(i,j) / (gd%n%x*gd%n%y)
    end do
  end do
  !$OMP END PARALLEL DO
  CALL DCFTD (0, 2, phi, incx, 1, phi2, incy, 1, narr, 1, isign_bwd, 1.0d0, aux1bwd, naux1, aux2, naux2)
  !$OMP PARALLEL DO private(i)
  do j = 0, gd%n%y-1
    do i = 0, gd%n%x-1
      k3(i,j) = -(0.0d0,1.0d0)*dcmplx(1.0d0,-gam)*delta_t * (phi2(i,j) &
                + (g*abs(k3(i,j))**2)*k3(i,j))
      k4(i,j) = phi_old(i,j) + b41*k1(i,j) + b42*k2(i,j) + b43*k3(i,j)
    end do
  end do
  !$OMP END PARALLEL DO

  !** Obtain k4 coefficients
  CALL DCFTD (0, 2, k4, incx, 1, phi, incy, 1, narr, 1, isign_fwd, 1.0d0, aux1fwd, naux1, aux2, naux2)
  !$OMP PARALLEL DO private(i)
  do j = 0, gd%n%y-1
    do i = 0, gd%n%x-1
      phi(i,j) = karray(i,j) * phi(i,j) / (gd%n%x*gd%n%y)
    end do
  end do
  !$OMP END PARALLEL DO
  CALL DCFTD (0, 2, phi, incx, 1, phi2, incy, 1, narr, 1, isign_bwd, 1.0d0, aux1bwd, naux1, aux2, naux2)
  !$OMP PARALLEL DO private(i)
  do j = 0, gd%n%y-1
    do i = 0, gd%n%x-1
      k4(i,j) = -(0.0d0,1.0d0)*dcmplx(1.0d0,-gam)*delta_t * (phi2(i,j) &
                + (g*abs(k4(i,j))**2)*k4(i,j))
      k5(i,j) = phi_old(i,j) + b51*k1(i,j) + b52*k2(i,j) + b53*k3(i,j) + b54*k4(i,j)
    end do
  end do
  !$OMP END PARALLEL DO

  !** Obtain k5 coefficients
  CALL DCFTD (0, 2, k5, incx, 1, phi, incy, 1, narr, 1, isign_fwd, 1.0d0, aux1fwd, naux1, aux2, naux2)
  !$OMP PARALLEL DO private(i)
  do j = 0, gd%n%y-1
    do i = 0, gd%n%x-1
      phi(i,j) = karray(i,j) * phi(i,j) / (gd%n%x*gd%n%y)
    end do
  end do
  !$OMP END PARALLEL DO
  CALL DCFTD (0, 2, phi, incx, 1, phi2, incy, 1, narr, 1, isign_bwd, 1.0d0, aux1bwd, naux1, aux2, naux2)
  !$OMP PARALLEL DO private(i)
  do j = 0, gd%n%y-1
    do i = 0, gd%n%x-1
      k5(i,j) = -(0.0d0,1.0d0)*dcmplx(1.d0,-gam)*delta_t * (phi2(i,j) &
                + (g*abs(k5(i,j))**2)*k5(i,j))
      k6(i,j) = phi_old(i,j) + b61*k1(i,j) + b62*k2(i,j) + b63*k3(i,j) + b64*k4(i,j) + b65*k5(i,j)
    end do
  end do
  !$OMP END PARALLEL DO

  !** Obtain k6 coefficients
  CALL DCFTD (0, 2, k6, incx, 1, phi, incy, 1, narr, 1, isign_fwd, 1.0d0, aux1fwd, naux1, aux2, naux2)
  !$OMP PARALLEL DO private(i)
  do j = 0, gd%n%y-1
    do i = 0, gd%n%x-1
      phi(i,j) = karray(i,j) * phi(i,j) / (gd%n%x*gd%n%y)
    end do
  end do
  !$OMP END PARALLEL DO
  CALL DCFTD (0, 2, phi, incx, 1, phi2, incy, 1, narr, 1, isign_bwd, 1.0d0, aux1bwd, naux1, aux2, naux2)
  error = 0.0d0
  newmax = 0.0d0
  !$OMP PARALLEL DO private(i) reduction(MAX:error) reduction(MAX:newmax)
  do j = 0, gd%n%y-1
    do i = 0, gd%n%x-1
      k6(i,j) = -(0.0d0,1.0d0)*dcmplx(1.0d0,-gam)*delta_t * (phi2(i,j) &
                + (g*abs(k6(i,j))**2)*k6(i,j))
      phi(i,j) = phi_old(i,j) + c1*k1(i,j) + c2*k2(i,j) + c3*k3(i,j) + c4*k4(i,j) &
                              + c5*k5(i,j) + c6*k6(i,j) !** Fifth-order accurate
      k4(i,j) = phi_old(i,j) + d1*k1(i,j) + d2*k2(i,j) + d3*k3(i,j) + d4*k4(i,j) &
                             + d5*k5(i,j) + d6*k6(i,j) !** Fourth-order accurate
      !** NB: Cheekily inherit maximum value from last timestep!
      if (abs(phi(i,j)) .gt. error_cut_off*oldmax) then
        if (abs(phi(i,j)-k4(i,j)) .gt. error) then
          error = abs(phi(i,j)-k4(i,j))
        end if
        if (abs(phi(i,j)) .gt. newmax) then
          newmax = abs(phi(i,j))
        end if
      end if
    end do
  end do
  !$OMP END PARALLEL DO
  oldmax = newmax
  !** Assess timestep
  if (error .gt. goal) then
    !write(*,*) 'Timestep rejected:', delta_t
    !$OMP PARALLEL DO private(i)
    do j = 0, gd%n%y-1
      do i = 0, gd%n%x-1
        phi(i,j) = phi_old(i,j)
      end do
    end do
    !$OMP END PARALLEL DO
    delta_t = 0.92d0 * (goal/error)**0.25d0 * delta_t
  else
    time = time + delta_t
    delta_t = 0.92d0 * (goal/error)**0.2d0 * delta_t

    if (present(tol)) then
      norm = 0.0d0
      !$OMP PARALLEL DO private(i) reduction(+:norm)
      do j = 0, gd%n%y-1
        do i = 0, gd%n%x-1
          norm = norm + abs(phi(i,j))**2
        end do
      end do
      !$OMP END PARALLEL DO
      norm = norm*gd%delta%x*gd%delta%y
      conv = abs(norm-old_norm)/norm
      write(*,*) norm, conv
      if (conv .lt. tol) then
        exit
      end if
      old_norm = norm
    else
      if (can_has_exit) then
        exit
      end if
    end if
  end if
end do

end subroutine


double precision function momentum_factor(i,j,n,dk_sq,mu)
use types
implicit none
integer, intent(in):: i,j
type (ivector), intent(in):: n
type (vector), intent(in):: dk_sq
double precision, intent(in):: mu

momentum_factor  = 0.5d0*(dk_sq%x*dble(fold(i,n%x))**2 + dk_sq%y*dble(fold(j,n%y))**2) - mu

end function


integer function fold(i,n)
implicit none
integer, intent(in):: i, n

if (i .le. n/2) then
  fold = i
else
  fold = i-n
end if

end function


double precision function cfl_timestep(delta)
use constants
use types
implicit none
type (vector), intent(in):: delta

 cfl_timestep = 2.0d0 / ( pi**2 * (1.0d0/(delta%x**2) + 1.0d0/(delta%y**2)) )

end function



end module
