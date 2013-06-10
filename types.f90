!***************************************************************************************
!** Otago Superfluid Simulation Code
!** Two-dimensional Gross-Pitaevskii simulation using Fourier pseudospectral method in a
!** periodic, uniform domain, with 4/5 order adaptive Runge_Kutta timestepping.
!**
!** 2013  Thomas P. Billam, University of Otago
!***************************************************************************************

module types
use, intrinsic:: iso_c_binding
implicit none
save

type:: ivector
  integer:: x,y
end type

type:: vector
  double precision:: x,y
end type

type:: transform_pair
  type(C_PTR)::  forward, backward
end type

type:: grid
  type (ivector):: n
  type (vector):: offset
  type (vector):: delta
end type

type:: vortex_info
  double precision:: x, y
  integer:: sig, cluster
  logical:: dipole
end type
  

end module 
