!***************************************************************************************
!** Otago Superfluid Simulation Code
!** Two-dimensional Gross-Pitaevskii simulation using Fourier pseudospectral method in a
!** periodic, uniform domain, with 4/5 order adaptive Runge_Kutta timestepping.
!**
!** 2013  Thomas P. Billam, University of Otago
!***************************************************************************************

module output
implicit none
save

contains


subroutine write_hdf5_named(filename, phi, gd, time)
use types
use hdf5
implicit none
type (grid), intent(in):: gd
double complex, intent(in):: phi(0:gd%n%x-1,0:gd%n%y-1)
character(len = 30), intent(inout):: filename
double precision, intent(in):: time
character(len=2),  parameter:: groupname = "/1"
character(len=4):: dsetname  ! dataset name
integer(hid_t) :: dataset_id    ! dataset identifier
integer(hid_t) :: dataspace_id  ! data space identifier
integer(hid_t) :: attr_id, attr_space
integer(hid_t):: file_id, group_id
integer:: error ! error flag
integer, parameter:: rank = 2 ! Datasets rank
integer(hsize_t), dimension(rank) :: dims
dims(1) = gd%n%x
dims(2) = gd%n%y
filename = trim(adjustl(filename))

!** Create a new file using default properties.
call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)

!** Create group "1" in the root group using absolute name.
call h5gcreate_f(file_id, groupname, group_id, error)

!** WRITE IMAGINARY COMPONENT

  !** Create the data space for the imaginary component
  call h5screate_simple_f(rank, dims, dataspace_id, error)
  
  !** Create the imaginary component dataset in group "1" with default properties.
  dsetname="phiI"
  call h5dcreate_f(group_id, dsetname, H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, error)
  
  !** Write the imaginary component dataset.
  call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, aimag(phi), dims, error)
  
  !** Close the dataspace for the imaginary component dataset.
  call h5sclose_f(dataspace_id, error)
  
  !** Close the imaginary component dataset.
  call h5dclose_f(dataset_id, error)

!** WRITE REAL COMPONENT

  !** Create the data space for the imaginary component
  call h5screate_simple_f(rank, dims, dataspace_id, error)
  
  !** Create the real component dataset in group "1" with default properties.
  dsetname="phiR"
  call h5dcreate_f(group_id, dsetname, H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, error)
  
  !** Write the real component dataset.
  call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, dble(phi), dims, error)
  
  !** Close the dataspace for the real component dataset.
  call h5sclose_f(dataspace_id, error)
  
  !** Close the real component dataset.
  call h5dclose_f(dataset_id, error)


!** Write attributes...
call h5screate_f(H5S_SCALAR_F, attr_id, error)
call h5acreate_f(group_id,'nx',H5T_NATIVE_INTEGER,attr_id,attr_space,error)
call h5awrite_f(attr_space, H5T_NATIVE_INTEGER, gd%n%x, dims, error)
call H5Sclose_f(attr_id,error) 
call H5Aclose_f(attr_space,error)

call h5screate_f(H5S_SCALAR_F, attr_id, error)
call h5acreate_f(group_id,'ny',H5T_NATIVE_INTEGER,attr_id,attr_space,error)
call h5awrite_f(attr_space, H5T_NATIVE_INTEGER, gd%n%y, dims, error)
call H5Sclose_f(attr_id,error) 
call H5Aclose_f(attr_space,error)

call h5screate_f(H5S_SCALAR_F, attr_id, error)
call h5acreate_f(group_id,'dx',H5T_NATIVE_DOUBLE,attr_id,attr_space,error)
call h5awrite_f(attr_space, H5T_NATIVE_DOUBLE, gd%delta%x, dims, error)
call H5Sclose_f(attr_id,error) 
call H5Aclose_f(attr_space,error)

call h5screate_f(H5S_SCALAR_F, attr_id, error)
call h5acreate_f(group_id,'dy',H5T_NATIVE_DOUBLE,attr_id,attr_space,error)
call h5awrite_f(attr_space, H5T_NATIVE_DOUBLE, gd%delta%x, dims, error)
call H5Sclose_f(attr_id,error) 
call H5Aclose_f(attr_space,error)

call h5screate_f(H5S_SCALAR_F, attr_id, error)
call h5acreate_f(group_id,'x0',H5T_NATIVE_DOUBLE,attr_id,attr_space,error)
call h5awrite_f(attr_space, H5T_NATIVE_DOUBLE, gd%offset%x, dims, error)
call H5Sclose_f(attr_id,error) 
call H5Aclose_f(attr_space,error)

call h5screate_f(H5S_SCALAR_F, attr_id, error)
call h5acreate_f(group_id,'y0',H5T_NATIVE_DOUBLE,attr_id,attr_space,error)
call h5awrite_f(attr_space, H5T_NATIVE_DOUBLE, gd%offset%y, dims, error)
call H5Sclose_f(attr_id,error) 
call H5Aclose_f(attr_space,error)

call h5screate_f(H5S_SCALAR_F, attr_id, error)
call h5acreate_f(group_id,'t',H5T_NATIVE_DOUBLE,attr_id,attr_space,error)
call h5awrite_f(attr_space, H5T_NATIVE_DOUBLE, time, dims, error)
call H5Sclose_f(attr_id,error) 
call H5Aclose_f(attr_space,error)

!** Close the group.
call h5gclose_f(group_id, error)

! Close the file.
call h5fclose_f(file_id, error)

end subroutine


subroutine get_grid_from_file(filename,gd)
use hdf5
use types
implicit none
character(len = 30), intent(inout):: filename
type(grid), intent(out):: gd
character(len=2),  parameter:: groupname = "/1"
integer(hid_t):: file_id, group_id, attr_id
integer:: error ! error flag
integer, parameter:: rank = 2 ! Datasets rank
integer(hsize_t), dimension(rank):: dims

filename = trim(adjustl(filename))

!** Open an existing file.
call h5fopen_f (filename, H5F_ACC_RDONLY_F, file_id, error)

!** Open group "1" in the root group using absolute name.
call h5gopen_f(file_id, groupname, group_id, error)

  !** READ GRID
  call H5Aopen_name_f(group_id,'nx',attr_id,error)
  call H5Aread_f(attr_id, H5T_NATIVE_INTEGER, gd%n%x, dims, error)
  call H5Aopen_name_f(group_id,'ny',attr_id,error)
  call H5Aread_f(attr_id, H5T_NATIVE_INTEGER, gd%n%y, dims, error)
  call H5Aopen_name_f(group_id,'dx',attr_id,error)
  call H5Aread_f(attr_id, H5T_NATIVE_DOUBLE, gd%delta%x, dims, error)
  call H5Aopen_name_f(group_id,'dy',attr_id,error)
  call H5Aread_f(attr_id, H5T_NATIVE_DOUBLE, gd%delta%y, dims, error)
  call H5Aopen_name_f(group_id,'x0',attr_id,error)
  call H5Aread_f(attr_id, H5T_NATIVE_DOUBLE, gd%offset%x, dims, error)
  call H5Aopen_name_f(group_id,'y0',attr_id,error)
  call H5Aread_f(attr_id, H5T_NATIVE_DOUBLE, gd%offset%y, dims, error)

!** Close the group.
call h5gclose_f(group_id, error)

! Close the file.
call h5fclose_f(file_id, error)

end subroutine


subroutine read_hdf5_named(filename, phi, gd, time)
use hdf5
use, intrinsic:: iso_c_binding
use types
implicit none
type (grid), intent(in):: gd
double complex, target, intent(inout):: phi(0:gd%n%x-1,0:gd%n%y-1)
character(len = 30), intent(inout):: filename
double precision, intent(out):: time
!** Pointer gymnastics to save memory...
type(C_PTR):: c
double precision, pointer:: p(:)
!** HDF5 vars
character(len=2),  parameter:: groupname = "/1"
character(len=4):: dsetname  ! dataset name
integer(hid_t) :: dataset_id    ! dataset identifier
integer(hid_t) :: memspace_id  ! data space identifier
integer(hid_t):: file_id, group_id, attr_id
integer:: error ! error flag
integer, parameter:: rank = 2 ! Datasets rank
integer(hsize_t), dimension(rank):: dims
integer, parameter:: mem_rank = 1 ! Datasets rank
integer(hsize_t), dimension(mem_rank):: mem_dims
integer(hsize_t), dimension(mem_rank):: stride = (/2/), off_real=(/0/), off_imag=(/1/), length

!** Pointer gym...
c = c_loc(phi)
call c_f_pointer(c,p,[gd%n%x*gd%n%y*2])

dims(1) = gd%n%x
dims(2) = gd%n%y
mem_dims(1) = gd%n%x*gd%n%y*2
length(1) = gd%n%x*gd%n%y
filename = trim(adjustl(filename))

!** Open an existing file.
call h5fopen_f (filename, H5F_ACC_RDONLY_F, file_id, error)

!** Open group "1" in the root group using absolute name.
call h5gopen_f(file_id, groupname, group_id, error)


!** Make a memory dataspace...
call h5screate_simple_f(mem_rank, mem_dims, memspace_id, error)


!** READ IMAGINARY COMPONENT
  !** Open imaginary component dataset
  dsetname="phiI"
  call h5dopen_f(group_id, dsetname, dataset_id, error)

  !** Select memory hyperslab
  CALL h5sselect_hyperslab_f(memspace_id, H5S_SELECT_SET_F, off_imag, length, error, stride)

  !** Read the imaginary component dataset.
  call h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, p, mem_dims, error,memspace_id)
  
  !** Close the imaginary component dataset.
  call h5dclose_f(dataset_id, error)

!** READ REAL COMPONENT
  !** Open real component dataset
  dsetname="phiR"
  call h5dopen_f(group_id, dsetname, dataset_id, error)
  
  !** Select memory hyperslab
  CALL h5sselect_hyperslab_f(memspace_id, H5S_SELECT_SET_F, off_real, length, error, stride)

  !** Read the real component dataset.
  call h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, p, mem_dims, error,memspace_id)
  
  !** Close the imaginary component dataset.
  call h5dclose_f(dataset_id, error)

  !** READ TIME 
  call H5Aopen_name_f(group_id,'t',attr_id,error)
  call H5Aread_f(attr_id, H5T_NATIVE_DOUBLE, time, dims, error)

!** Close the group.
call h5gclose_f(group_id, error)

! Close the file.
call h5fclose_f(file_id, error)

end subroutine


subroutine write_hdf5(filenum, phi, gd, time)
use types
use hdf5
implicit none
type (grid), intent(in):: gd
double complex, intent(in):: phi(0:gd%n%x-1,0:gd%n%y-1)
integer, intent(in):: filenum
double precision, intent(in):: time
character(len = 30):: filename

write(filename,*) filenum
filename = trim(adjustl(filename))//'.h5'
call write_hdf5_named(filename,phi,gd,time)

!** Record this as the highest good output
open(417,file='LastGoodSave.adat')
write(417,*) filenum
flush(417)
close(417)

end subroutine


end module 
