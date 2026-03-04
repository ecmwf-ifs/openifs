! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

Interface

SUBROUTINE RTTOV_EC_ALLOC( &
     & errorstatus, ksat, nchannels, nprofiles, nlevels,  &
     & ec_opts, &
     & asw)

  Use parkind1, Only : jpim
  Use rttov_ec_mod, only : rttov_ec_opts

  IMPLICIT NONE

 ! Subroutine arguments
  ! Scalar arguments with intent(in):
  Integer(Kind=jpim), INTENT(in) :: ksat        ! Satellite index (see rttvi)
  Integer(Kind=jpim), INTENT(in) :: nchannels   ! Number of radiances
  Integer(Kind=jpim), INTENT(in) :: nprofiles   ! Number of profiles
  Integer(Kind=jpim), INTENT(in) :: nlevels     ! Number of input levels
  Integer(Kind=jpim), INTENT(in) :: asw         ! Allocate/Deallocate flag

  Type(rttov_ec_opts), INTENT(in) :: ec_opts ! Options

  ! Array  arguments with intent(out):
  Integer(Kind=jpim) , INTENT(out):: errorstatus  !  return flag 

END SUBROUTINE RTTOV_EC_ALLOC

End Interface
