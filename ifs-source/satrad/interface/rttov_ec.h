! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

Interface
SUBROUTINE RTTOV_EC( &
     & errorstatus, ksat, nchannels, nprofiles, nlevels, ntoplevels, nav, &
     & channels, lprofiles,                                          &
     & ppres, pav, pt2m, pq2m, ppsurf, pu10m, pv10m,                 &
     & pctp, pcfrac, ptskin, ksurf, emissivity,                      &
     & pzenith, pazimuth, psolzenith, psolazimuth,                   &
     & plat, plon, porog,                                            &
     & ec_opts,                                      &
     & ptbclr, pradclr, ptbtot, pradtot, pradcld, pradcldlev,        &
     & tau, tausfc, lretr_emis_tskin, pobs, pvarout, pcape, specularity)

  Use parkind1, Only : jpim     ,jprb
  Use rttov_ec_mod, only : rttov_ec_opts
  
  IMPLICIT NONE


  ! Subroutine arguments
  ! Scalar arguments with intent(in):
  Integer(Kind=jpim) , INTENT(in) :: ksat        ! Satellite index (see rttvi)
  Integer(Kind=jpim), INTENT(in)  :: nchannels   ! Number of radiances
  Integer(Kind=jpim) , INTENT(in) :: nprofiles   ! Number of profiles
  Integer(Kind=jpim) , INTENT(in) :: nlevels     ! Number of input levels
  Integer(Kind=jpim) , INTENT(in) :: ntoplevels  ! Number of top RTTOV levels
                                                 ! required for extension
  Integer(Kind=jpim) , INTENT(in) :: nav

  ! Array  arguments with intent(in):
  Integer(Kind=jpim) , INTENT(in) :: channels(nchannels) ! Channel indices
  Integer(Kind=jpim) , INTENT(in) :: lprofiles(nchannels)! Profile indices
  Real(Kind=jprb)    , INTENT(in) :: ppres(nlevels,nprofiles)! Pressure levels (hpa) of 
                                         !   atmospheric profile vectors

  Real(Kind=jprb) , INTENT(in)    :: pav(nlevels,nav,nprofiles)! Atmosp. profile variables
  Real(Kind=jprb) , INTENT(in)    :: pt2m(nprofiles)! 2m-temperature
  Real(Kind=jprb) , INTENT(in)    :: pq2m(nprofiles)! 2m-humidity
  Real(Kind=jprb) , INTENT(in)    :: ppsurf(nprofiles)! Surface pressure [hPa]
  Real(Kind=jprb) , INTENT(in)    :: pu10m(nprofiles)! 10m-wind (u-component)
  Real(Kind=jprb) , INTENT(in)    :: pv10m(nprofiles)! 10m-wind (v-component)
  Real(Kind=jprb) , INTENT(in)    :: pctp (nprofiles) ! Cloud top pressure [hPa]
  Real(Kind=jprb) , INTENT(in)    :: pcfrac(nprofiles)! Cloud fraction
  Real(Kind=jprb) , INTENT(in)    :: ptskin(nprofiles)! Surface skin temperature
  Integer(Kind=jpim) , INTENT(in) :: ksurf(nprofiles)   ! Surface type index
  Real(Kind=jprb) , INTENT(in)    :: pzenith(nprofiles)! Local satellite zenith angle (deg)
  Real(Kind=jprb) , INTENT(in)    :: pazimuth(nprofiles)! Local satellite azimuth angle (deg)
  Real(Kind=jprb) , INTENT(in)    :: psolzenith(nprofiles)! Solar zenith angle (deg)
  Real(Kind=jprb) , INTENT(in)    :: psolazimuth(nprofiles)! Solar azimuth angle (deg)

  Real(Kind=jprb) , INTENT(in)    :: plat(nprofiles)   ! Latitude [deg]
  Real(Kind=jprb) , INTENT(in)    :: plon(nprofiles)   ! Longitude [deg]
  Real(Kind=jprb) , INTENT(in)    :: porog(nprofiles)  ! Orography [m]

  Type(rttov_ec_opts), INTENT(in) :: ec_opts ! Options

  ! Array arguments with intent(inout):
  Real(Kind=jprb) , INTENT(inout) :: emissivity(nchannels)!  surface emissivities

  ! Array  arguments with intent(out):
  Integer(Kind=jpim) , INTENT(out):: errorstatus  !  return flag 

  Real(Kind=jprb) , INTENT(out)   :: ptbclr(nchannels)! clear brightness temperatures (K)
  Real(Kind=jprb) , INTENT(out)   :: pradclr(nchannels)! clear radiances (mw/cm-1/ster/sq.m)
  Real(Kind=jprb) , INTENT(out)   :: ptbtot(nchannels)! total brightness temperatures (K)
  Real(Kind=jprb) , INTENT(out)   :: pradtot(nchannels)!total radiance mw/m2/sr/cm-1
  Real(Kind=jprb) , INTENT(out)   :: pradcld(nchannels) ! 100% cloudy radiance at given
                                                !   cloud top in mw/m2/sr/cm-1
  REAL(Kind=jprb) , INTENT(out)   :: pradcldlev(nlevels,nchannels) 
    ! 100% cloudy radiance for each level

  REAL(Kind=jprb) , INTENT(out)    :: tau(nlevels,nchannels)  ! transmittance from each 
                                                !   standard pressure level
  Real(Kind=jprb) , INTENT(out)    :: tausfc(nchannels)      ! transmittance from surface
  Real(Kind=jprb) , INTENT(in)     :: pobs(nchannels) !   Observed TB
  Real(Kind=jprb) , INTENT(out)    :: pvarout(4,nchannels)
  Logical, INTENT(in)   :: lretr_emis_tskin  ! switch for emis/tskin retrievals

  Real(Kind=jprb) , OPTIONAL, INTENT(in) :: pcape(nprofiles)
  Real(Kind=jprb) , OPTIONAL, INTENT(in) :: specularity(nchannels)

END SUBROUTINE RTTOV_EC

End Interface
