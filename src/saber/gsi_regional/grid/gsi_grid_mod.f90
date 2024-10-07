! (C) Copyright 2022 United States Government as represented by the Administrator of the National
!     Aeronautics and Space Administration
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module gsi_grid_mod

! netcdf
use netcdf

! atlas
use atlas_module,                   only: atlas_field, atlas_fieldset, atlas_real

! fckit
use fckit_mpi_module,               only: fckit_mpi_comm
use fckit_configuration_module,     only: fckit_configuration

! oops
use kinds,                          only: kind_real

! saber
use gsi_utils_mod,                  only: nccheck

! gsifull
use m_gsimod,                         only: gsimain_gridopts
use m_gsi,                            only: gsi_get_grid
use m_gsi,                            only: gsi_set_grid

implicit none
private
public gsi_grid

! Fortran class header
type :: gsi_grid
  type(fckit_mpi_comm) :: comm
  character(len=2055) :: filename
  integer :: npx, npy, npz          ! Grid points in global grid
  integer :: layout(2)              ! Number of processors in x (index 1) and y (index 2) directions
  integer :: lat2,lon2
  integer :: isc, iec, jsc, jec     ! Start and ending grid points for each processor
  logical :: vflip                  ! Flip vertical grid (gsi k=1=top)
  logical :: noGSI
  real(kind=kind_real), allocatable :: lats(:,:), lons(:,:)
  real(kind=kind_real), allocatable :: grid_lats(:,:), grid_lons(:,:)
  integer :: ngrid ! Number of grid points for each processor
  logical :: debug
  contains
    procedure, public :: create
    procedure, public :: delete
    procedure, public :: print
    procedure, public :: get_levels
    procedure, public :: set_atlas_lonlat
end type gsi_grid

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

subroutine create(self, conf, comm)

! Arguments
class(gsi_grid),           intent(inout) :: self
type(fckit_configuration), intent(in)    :: conf
type(fckit_mpi_comm),      intent(in)    :: comm

! Locals
integer :: ncid, dimid(3), varid(2)
character(len=:), allocatable :: str
character(len=:), allocatable :: nml, vgrdfn
integer :: posx, posy, i, j
logical :: verbose

! Create copy of comm
! -------------------
self%comm = comm
verbose = comm%rank()==0

! Debug mode
! ----------
call conf%get_or_die("debugging mode", self%debug)
call conf%get_or_die("debugging bypass gsi", self%noGSI)

! Domain decomposition
! --------------------
if (conf%has("processor layout x direction").and.conf%has("processor layout y direction")) then
  call conf%get_or_die("processor layout x direction", self%layout(1))
  call conf%get_or_die("processor layout y direction", self%layout(2))
else
  self%layout(1) = floor(sqrt(real(comm%size(),kind_real)))
  self%layout(2) = comm%size()/self%layout(1)
end if

! Handle vertical grid opt
! ------------------------
call conf%get_or_die("flip vertical grid", self%vflip)

! Open file with GSI grid info (for now here)
! ----------------------------
if (comm%rank() == 0) then

  ! Get filename
  call conf%get_or_die("gsi error covariance file", str)
  self%filename = str

endif

call wGSI()

! Create arrays of lon/lat to be compatible with interpolation
if(.not.allocated(self%grid_lons)) allocate(self%grid_lons(self%isc:self%iec, self%jsc:self%jec))
if(.not.allocated(self%grid_lats)) allocate(self%grid_lats(self%isc:self%iec, self%jsc:self%jec))

do i = self%isc, self%iec
  do j = self%jsc, self%jec
    self%grid_lons(i,j) = self%lons(i,j)
    self%grid_lats(i,j) = self%lats(i,j)
  enddo
enddo

if ( self%debug ) then
    do j=1,self%layout(1)*self%layout(2)
       if(self%comm%rank() == j-1) then
       write(6,'(a,6(i5,1x))') 'grid dist indexes: task, is,ie, js,je ', j-1, &
                                self%isc, self%iec, &
                                self%jsc, self%jec, &
                                self%ngrid
       endif
    enddo
endif

contains
! ------------------
! Actual hook to GSI (lat/lon to come from GSI)
! ------------------
  subroutine wGSI

  integer :: npe,igdim
  logical :: eqspace

  npe = self%layout(1)*self%layout(2)

  ! Check that user choices match comm size
  if (.not. self%layout(1)*self%layout(2) == comm%size()) &
    call abor1_ftn("GSI grid: number of processor in layout does not match number in communicator")

  ! Get required name of resources for GSI B error
  ! ----------------------------------------------
  !call conf%get_or_die("gsi berror namelist file2",  nml)
  nml="Data/gsibec/gsiparm.anl"
  call conf%get_or_die("gsi akbk",  vgrdfn)

  ! Initialize GSIbec grid
  ! ----------------------
  call gsimain_gridopts (nml,comm%rank(),self%layout(1),self%layout(2),&
                         self%npy,self%npx,self%npz,eqspace,&
                         self%lon2,self%lat2,&
                         self%isc,self%iec,self%jsc,self%jec,igdim)

  ! Allocate the lat/lon arrays
  ! ---------------------------
  if(.not.allocated(self%lons)) allocate(self%lons(self%npx,self%npy))
  if(.not.allocated(self%lats)) allocate(self%lats(self%npx,self%npy))

  ! Read the latitudes and longitudes per GSIbec
  ! --------------------------------------------
  call gsi_get_grid ('degree',self%lats,self%lons)
  call gsi_set_grid (comm%rank(),vgrdfn)

  self%ngrid = (self%iec-self%isc+1)*(self%jec-self%jsc+1)
  if (self%ngrid /= igdim) then
    call abor1_ftn("gsi_grid_mod: inconsistent distribution")
  endif

  end subroutine wGSI

end subroutine create

! --------------------------------------------------------------------------------------------------

subroutine delete(self)

! Arguments
class(gsi_grid), intent(inout) :: self

! Deallocate arrays
deallocate(self%lons)
deallocate(self%lats)
deallocate(self%grid_lons)
deallocate(self%grid_lats)

! Set grid to zero
self%npx = 0
self%npy = 0
self%npz = 0
self%layout = 0
self%isc = 0
self%iec = 0
self%jsc = 0
self%jec = 0

end subroutine delete

! --------------------------------------------------------------------------------------------------

subroutine print(self)

! Arguments
class(gsi_grid), intent(in) :: self

! Root PE prints grid info
if (self%comm%rank() == 0) then

  write(*,'(A28)')      "ErrorCovarianceGSI GSI grid:"
  write(*,'(A38, I5)')  "  Number of longitudinal grid points: ", self%npx
  write(*,'(A37, I5)')  "  Number of latitudinal grid points: ", self%npy
  write(*,'(A34, I5)')  "  Number of vertical grid points: ", self%npz
  write(*,'(A1)')       " "
  write(*,'(A43, I5)')  "  Number of processors in the x direction: ", self%layout(1)
  write(*,'(A43, I5)')  "  Number of processors in the y direction: ", self%layout(2)

endif

if (self%debug) then
  ! Print index ranges
  write(*,'(A7, I6, A7, I6, A7, I6, A7, I6, A7, I6)')  "  Proc ", self%comm%rank(), &
                        ' isc = ', self%isc, ' iec = ', self%iec, &
                        ' jsc = ', self%jsc, ' jec = ', self%jec

  ! Print latlon
  write(*,'(A10, F10.3, A10, F10.3, A10, F10.3, A10, F10.3)')  &
        "  Lat min ", minval(self%grid_lats), &
        "  Lat max ", maxval(self%grid_lats), &
        "  Lon min ", minval(self%grid_lons), &
        "  Lon max ", maxval(self%grid_lons)
endif

end subroutine print

! --------------------------------------------------------------------------------------------------

subroutine get_levels(self, levels)

! Arguments
class(gsi_grid), intent(in)    :: self
integer,         intent(inout) :: levels

! Get number of levels
! --------------------
levels = self%npz

end subroutine get_levels

! --------------------------------------------------------------------------------------------------

subroutine set_atlas_lonlat(self, grid_fieldset)

!Arguments
class(gsi_grid),      intent(inout) :: self
type(atlas_fieldset), intent(inout) :: grid_fieldset

!Locals
real(kind_real), pointer :: real_ptr(:,:)
type(atlas_field) :: lonlat_field

! Create lon/lat field
lonlat_field = atlas_field(name="lonlat", kind=atlas_real(kind_real), shape=(/2,self%ngrid/))

! Get pointer to the data
call lonlat_field%data(real_ptr)

! Fill lon/lat
real_ptr(1,:) = reshape(self%grid_lons(self%isc:self%iec, &
                                       self%jsc:self%jec), (/self%ngrid/))
real_ptr(2,:) = reshape(self%grid_lats(self%isc:self%iec, &
                                       self%jsc:self%jec), (/self%ngrid/))

! Add field to fieldset
call grid_fieldset%add(lonlat_field)

end subroutine set_atlas_lonlat

! --------------------------------------------------------------------------------------------------

end module gsi_grid_mod
