! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction
! 
! (C) Copyright 1989- Meteo-France.
! 

MODULE GOM2D_SUPPORT

! Module to support the use of 2d GOMs: 
!   - Includes routines how the locations of the GOM profiles should
!     be calulated (currently limb_plane only)
!   - Defines IDs describing what further processing should be done
!     with the 2D GOMs. 
! Currently supports NSATEM and NALLSKY, and partly GPSRO, but could 
! be extended to other observation types. 2D-GOMs are activated by
! specifying the number of profiles used per observation type (NOBSPROF) 
! in the namelist NAMNPROF.
!
! Niels Bormann, 5 August 2020

USE PARKIND1   , ONLY : JPIM, JPRB
USE YOMCOCTP , ONLY : NGTHRB, NSSMI, NGPSRO, NSTRESAT, NSMOS
USE RTTOV_CONST , ONLY : INST_ID_AMSUA, INST_ID_ATMS, INST_ID_MHS, &
                       & INST_ID_MWHS2, INST_ID_MWTS3, INST_ID_TROPICS
USE YOMHOOK   , ONLY : LHOOK, DR_HOOK, JPHOOK

IMPLICIT NONE
SAVE

! IDs referring to specific treatments for 2D GOMs (if 2d GOMS are activated through
! NOBSPROF in NAMNPROF)
INTEGER(KIND=JPIM), PARAMETER :: JP_TREAT_2DGOM_UNSPEC = -1    ! Nothing specified
INTEGER(KIND=JPIM), PARAMETER :: JP_TREAT_2DGOM_SINGLE = 0     ! Use single profile at obs location instead
INTEGER(KIND=JPIM), PARAMETER :: JP_TREAT_2DGOM_LIMB_PLANE = 1 ! Use limb-plane centred at obs location
INTEGER(KIND=JPIM), PARAMETER :: JP_TREAT_2DGOM_SLANT_PATH = 2 ! Do slant-path

REAL(KIND=JPRB), PARAMETER:: RDXMAX_SLANT = 160000._JPRB ! Maximum extent [m] of slant-path plane

CONTAINS
 
! ---------------------------------------------------------------------

FUNCTION GET_2DGOM_TREATMENT(KCDTYPE, KSENSOR)

! Function to define what kind of treatment should be used for the for the 2d GOMs.

INTEGER (KIND=JPIM), INTENT(IN) :: KCDTYPE
INTEGER (KIND=JPIM), INTENT(IN) :: KSENSOR

INTEGER (KIND=JPIM) :: GET_2DGOM_TREATMENT

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
! ---------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GOM2D_SUPPORT:GET_2DGOM_TREATMENT',0,ZHOOK_HANDLE)

IF(KCDTYPE == NGTHRB)THEN
  ! Slant-path treatment for all radiances going through the clear-sky route.
  GET_2DGOM_TREATMENT = JP_TREAT_2DGOM_SLANT_PATH

ELSEIF(KCDTYPE == NSSMI)THEN
  ! Slant-path treatment for certain sensors for all-sky radiances
  IF(ANY(KSENSOR == (/INST_ID_AMSUA, INST_ID_ATMS, INST_ID_MHS, & 
   &                  INST_ID_MWHS2, INST_ID_MWTS3, INST_ID_TROPICS/) )) THEN
    GET_2DGOM_TREATMENT = JP_TREAT_2DGOM_SLANT_PATH
  ELSE
  ! Others remain single-profile at obs location
    GET_2DGOM_TREATMENT = JP_TREAT_2DGOM_SINGLE
  ENDIF 

ELSEIF(KCDTYPE == NGPSRO)THEN   ! This is currently never used later
  ! Provide a limb-plane for GPRSO observations
  GET_2DGOM_TREATMENT = JP_TREAT_2DGOM_LIMB_PLANE

ELSEIF(KCDTYPE == NSTRESAT .OR. KCDTYPE == NSMOS)THEN
  ! Single-profile treatment for retrievals or SMOS data
  GET_2DGOM_TREATMENT = JP_TREAT_2DGOM_SINGLE

ELSE
  GET_2DGOM_TREATMENT = JP_TREAT_2DGOM_UNSPEC
ENDIF

IF (LHOOK) CALL DR_HOOK('GOM2D_SUPPORT:GET_2DGOM_TREATMENT',1,ZHOOK_HANDLE)

END FUNCTION GET_2DGOM_TREATMENT


!****************************************************************
SUBROUTINE LIMB_PLANE(PLAT,PLON,PAZIM1,KPROF,LDCENTRED,PLAT_2D,PLON_2D,KERROR)

! Routine to calculate a series of points to define the viewing
! plane for a limb instrument or for slant-path radiative transfer, 
! given the latitude and longitude of the sub-tangent point (or the 
! surface location for slant-path RT), the angle pazim1 between the plane and
! a longitude line. The definition of the viewing 
! plane is invalid if the sub-tangent point is at the pole, so the
! routine returns error=1.
! Formulae derived by considering that
!  cos(theta) = DOT_PRODUCT(sub tangent vector, vector to calculate)
! and
!  cos(pazim1) = DOT_PRODUCT(normalised vector orthogonal to longitude plane,
!                            normalised vector orthogonal to viewing plane)

! Conventions: pazim1 clockwise (ie eastward is positive)
!              theta  positive is northward if pazim1 is between -90 and 90 deg

! Niels Bormann, 11 March 2004 
!
! Modifications
! -------------
!  1-Feb-2016 Niels Bormann     Non-centred points for slant-path calculations
!-------------------------------------------------------------

USE YOMLIMB   ,ONLY : DTHETA, INUM_TRAJ
USE YOMCST   , ONLY : RA, RPI
USE YOMANCS  , ONLY : RMDI

IMPLICIT NONE

REAL(KIND=JPRB),    INTENT(IN)    :: PLAT    ! Lat of sub-tangent pt [deg]
REAL(KIND=JPRB),    INTENT(IN)    :: PLON    ! Lon of sub-tangent pt [deg]
REAL(KIND=JPRB),    INTENT(IN)    :: PAZIM1  ! Angle [deg to north] of vieving plane,
                                             ! positive is eastward
INTEGER(KIND=JPIM), INTENT(IN)    :: KPROF          ! Number of points to calculate
LOGICAL,            INTENT(IN)    :: LDCENTRED      ! T: Plane centred around plat/plon 
                                                    ! F: Plane in direction of pazim1
REAL(KIND=JPRB),    INTENT(OUT)   :: PLAT_2D(KPROF) ! Lats of points [deg] 
REAL(KIND=JPRB),    INTENT(OUT)   :: PLON_2D(KPROF) ! Lons of points [deg]
INTEGER(KIND=JPIM), INTENT(OUT)   :: KERROR         ! Error code

! local

REAL(KIND=JPRB)          :: ZTHETA(KPROF)
REAL(KIND=JPRB)          :: ZCOS_PAZIM
REAL(KIND=JPRB)          :: ZCOS_PLON
REAL(KIND=JPRB)          :: ZCOS_PLAT
REAL(KIND=JPRB)          :: ZSIN_PLON
REAL(KIND=JPRB)          :: ZSIN_PLAT
REAL(KIND=JPRB)          :: ZSIN_THETA
REAL(KIND=JPRB)          :: ZCOS_THETA
REAL(KIND=JPRB)          :: ZSIN_PLAT_2D
REAL(KIND=JPRB)          :: ZCOS_DELTA
REAL(KIND=JPRB)          :: ZRPLAT
REAL(KIND=JPRB)          :: ZAZIM,ZDTHETA_USE

INTEGER(KIND=JPIM)       :: I,IMID

REAL(KIND=JPHOOK)        :: ZHOOK_HANDLE

!-------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GOM2D_SUPPORT:LIMB_PLANE',0,ZHOOK_HANDLE)

KERROR = 0

! Only proceed if azimuth angle is not missing
IF(PAZIM1 /= RMDI) THEN

  ! scale the separation for inner loop and decide on profile order

  IF(LDCENTRED)THEN
    ! This applies to real limb data, e.g., GPSRO
    ZDTHETA_USE = REAL(INUM_TRAJ-1)/REAL(KPROF-1)*DTHETA
    IMID = KPROF/2 + 1  ! THE CENTRAL PROFILE
  ELSE
    ! This one to the slant-path RT calculations
    ! This hard-wiring needs to be consistent with what
    ! is later used in the slant path calculations
    ! in create_gom_plus.
    ! There should be a better way of doing this...
    ZDTHETA_USE = RDXMAX_SLANT / REAL(KPROF-1) / RA
    IMID = 1
  ENDIF

  ! calculate the theta values in the plane

  DO I = 1, KPROF
     ZTHETA(I) = REAL(I-IMID)*ZDTHETA_USE  ! dtheta in radians
  ENDDO   

  IF ( ABS(PLAT) /= 0.5_JPRB*RPI ) THEN

    ZAZIM = PAZIM1
    IF (PAZIM1 > RPI ) ZAZIM =  PAZIM1 - 2.0_JPRB*RPI  

    ZCOS_PAZIM = COS(ZAZIM)
    ZCOS_PLON  = COS(PLON)
    ZCOS_PLAT  = COS(PLAT)
    ZSIN_PLON  = SIN(PLON)
    ZSIN_PLAT  = SIN(PLAT)

    DO I = 1, KPROF

      ZSIN_THETA = SIN(ZTHETA(I))  ! dtheta in radians
      ZCOS_THETA = COS(ZTHETA(I))

      ZSIN_PLAT_2D = ZSIN_THETA * ZCOS_PAZIM * ZCOS_PLAT + ZSIN_PLAT * ZCOS_THETA

      ZRPLAT = ASIN(ZSIN_PLAT_2D)

      IF( ABS(ZRPLAT) == RPI/2.0_jprb ) THEN
        ZCOS_DELTA = 1.0_JPRB
      ELSE
        ZCOS_DELTA = ( ZCOS_THETA - ZSIN_PLAT * ZSIN_PLAT_2D ) / ( ZCOS_PLAT * COS(ZRPLAT) )
        IF( ZCOS_DELTA > 1.0_JPRB)  ZCOS_DELTA = 1.0_JPRB
        IF( ZCOS_DELTA < -1.0_JPRB) ZCOS_DELTA = -1.0_JPRB
      ENDIF

      IF( SIGN(1.0_JPRB,ZTHETA(I)) == SIGN(1.0_JPRB,ZAZIM) ) THEN
        PLON_2D(I) = PLON + ACOS(ZCOS_DELTA)
      ELSE
        PLON_2D(I) = PLON - ACOS(ZCOS_DELTA)
      ENDIF

      PLAT_2D(I) = ZRPLAT 

      PLON_2D(I) = MODULO(PLON_2D(I)+RPI,2.0_JPRB*RPI) - RPI

    ENDDO

    ! Avoid any chance of rounding errors...
    IF(.NOT. LDCENTRED)THEN
      PLAT_2D(1) = PLAT
      PLON_2D(1) = PLON
    ENDIF 

  ELSE

    ! Exact pole
    KERROR = 1

  ENDIF

ELSE

  ! Azimuth missing
  PLAT_2D(:) = PLAT 
  PLON_2D(:) = PLON

ENDIF

IF (LHOOK) CALL DR_HOOK('GOM2D_SUPPORT:LIMB_PLANE',1,ZHOOK_HANDLE)

END SUBROUTINE LIMB_PLANE

END MODULE GOM2D_SUPPORT
