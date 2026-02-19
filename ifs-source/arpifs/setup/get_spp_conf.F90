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

SUBROUTINE GET_SPP_CONF( YDRIP, YDCONF )
  
! Purpose :
! -------
!    Get default configuration of the SPP scheme, which represents model uncertainties
!    with stochastically perturbed parameterisations.
! 
!    Sets defaults, reads namelist variables and prints actual configuration to NULOUT

! Interface :
! ---------
!    Empty.

! External :
! --------
!    None.

! Method :
! ------
!    See Documentation.

! Reference :
! ---------

! Author :
! ------
!    M. Leutbecher (ECMWF)
!    Original : December 2014

! Modifications :
! -------------
!
!    P. Ollinaho Sep-2015: Default configuration
!    S. Lock     Jan-2016: Added option for 2 patterns for ice/liquid effective radii
!    S. Lang     Feb-2017: Added modification to cycle pattern in data assimilation
!    SJ Lock     Oct 2017: Option to reduce frequency of SPP pattern updates
!    SJ Lock     Oct-2017: Enabled options for new SPP microphysics perturbations
!    S. Lang     Mar-2019: Added option to shift seed
!    M Leutbecher Oct-2020: Abstraction allowing independent SPP configurations for ARPEGE, IFS, LAM
!
!-----------------------------------------------------------------------------
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN   , ONLY : NULOUT, NULNAM
USE YOMRIP   , ONLY : TRIP
USE SPP_MOD  , ONLY : TSPP_CONFIG, NWRMAX
USE SPP_DEF_MOD, ONLY: DEFINE_SPP_ARPEGE, DEFINE_SPP_IFS, DEFINE_SPP_LAM
USE SPP_GEN_MOD, ONLY: GET_ACTIVE_SPP_PERTS, MAP_INDICES_SPP_NML, MODIFY_SPP_PERTS, WRITE_SPP_MODEL_TABLE, JPMAXPERTS
! ------------------------------------------------------------------

IMPLICIT NONE

TYPE(TRIP),        INTENT(IN)    :: YDRIP
TYPE(TSPP_CONFIG), INTENT(OUT)   :: YDCONF
! ------------------------------------------------------------------
REAL(KIND=JPHOOK)  :: ZHOOK_HANDLE

! for namelist namspp
LOGICAL            :: LSPP
CHARACTER(LEN=32)  :: CSPP_MODEL_NAME
LOGICAL            :: LUSE_SETRAN, LNSEED_OFFS0
LOGICAL            :: LRAMIDLIMIT1
INTEGER(KIND=JPIM) :: ISEEDFAC, SHIFTSEED
REAL(KIND=JPRB)    :: XPRESS_PHR_ST

LOGICAL            :: LRDPATINIT, LWRPATTRUN, LRESETSEED, LABSTIMSEED 
INTEGER(KIND=JPIM) :: NWRPATTRUN, RESETSEEDFRQ, JPATT, NPATFR
INTEGER(KIND=JPIM), DIMENSION(NWRMAX) :: NHOUR_PATTRUN 
INTEGER(KIND=JPIM), DIMENSION(NWRMAX) :: NSTEP_PATTRUN 
CHARACTER(LEN=256) :: SPP_RDPATINIT, SPP_WRPATTRUN
!
!   for processing namelist NAM_SPP_MODIF
!
integer(kind=jpim), dimension(jpmaxperts) :: idx_map

#include "abor1.intfb.h"
#include "namspp.nam.h"
#include "posnam.intfb.h"


! ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GET_SPP_CONF',0,ZHOOK_HANDLE)


  LSPP                = .FALSE.


  LRAMIDLIMIT1= .FALSE.
  LUSE_SETRAN = .TRUE.
  LNSEED_OFFS0= .FALSE.

  ISEEDFAC   = 10000
  NPATFR = 1
  SHIFTSEED  = 0

  LRDPATINIT  = .FALSE.
  LWRPATTRUN  = .FALSE.
  LRESETSEED  = .FALSE.
  LABSTIMSEED = .FALSE.

  NWRPATTRUN   = 1
  RESETSEEDFRQ = 1
  NHOUR_PATTRUN(1:NWRMAX) = -1000
  NSTEP_PATTRUN(1:NWRMAX) = -1000
  SPP_RDPATINIT ='spp_pattern_0'
  SPP_WRPATTRUN ='spp_pattern_out'


  XPRESS_PHR_ST = 100.E2_JPRB ! 100 hPa

  CALL POSNAM(NULNAM,'NAMSPP')
  READ (NULNAM, NAMSPP)
  YDCONF%LSPP        = LSPP
  YDCONF%CSPP_MODEL_NAME = CSPP_MODEL_NAME

  YDCONF%LRAMIDLIMIT1  = LRAMIDLIMIT1
  YDCONF%LUSE_SETRAN   = LUSE_SETRAN
  YDCONF%LNSEED_OFFS0  = LNSEED_OFFS0
  
  YDCONF%ISEEDFAC      = ISEEDFAC
  YDCONF%NPATFR        = NPATFR

  YDCONF%SHIFTSEED     = SHIFTSEED
  
  YDCONF%LRDPATINIT=LRDPATINIT
  YDCONF%LWRPATTRUN=LWRPATTRUN
  YDCONF%LRESETSEED=LRESETSEED
  YDCONF%LABSTIMSEED=LABSTIMSEED

  IF (LWRPATTRUN) THEN
    IF (NWRPATTRUN>NWRMAX) THEN
      CALL ABOR1(' NWRPATTRUN>NWRMAX ')
    ENDIF
    DO JPATT=1, NWRPATTRUN
      NSTEP_PATTRUN(JPATT) = NINT( NHOUR_PATTRUN(JPATT)*3600/ YDRIP%TDT )
    ENDDO
  ENDIF

  YDCONF%NHOUR_PATTRUN(1:NWRPATTRUN)=NHOUR_PATTRUN(1:NWRPATTRUN)
  YDCONF%NSTEP_PATTRUN(1:NWRPATTRUN)=NSTEP_PATTRUN(1:NWRPATTRUN)

  YDCONF%NWRPATTRUN=NWRPATTRUN
  YDCONF%RESETSEEDFRQ=RESETSEEDFRQ

  YDCONF%SPP_RDPATINIT=SPP_RDPATINIT
  YDCONF%SPP_WRPATTRUN=SPP_WRPATTRUN

  !NPATFR: -> freq of pattern updates as integer number of timesteps
  IF (YDCONF%NPATFR<0) THEN
    YDCONF%NPATFR=NINT(-YDCONF%NPATFR*3600/YDRIP%TDT)
  ENDIF
  IF (YDCONF%NPATFR==0) THEN
    YDCONF%NPATFR=1
  ENDIF

  YDCONF%XPRESS_PHR_ST = XPRESS_PHR_ST

  ! ...

  ! write output
  WRITE(NULOUT,'(''LSPP='',L1)') YDCONF%LSPP 
  IF (YDCONF%LSPP) THEN
    ! Which SPP model is this?
    WRITE(NULOUT,'(''CSPP_MODEL_NAME='',A)') TRIM(YDCONF%CSPP_MODEL_NAME)
    SELECT CASE(YDCONF%CSPP_MODEL_NAME)
    CASE ("lam")
      CALL DEFINE_SPP_LAM( YDCONF%SM, NULOUT )
    CASE ("arpege")
      CALL DEFINE_SPP_ARPEGE( YDCONF%SM, NULOUT )
    CASE("ifs")
      CALL DEFINE_SPP_IFS( YDCONF%SM, NULOUT ) 
    CASE DEFAULT
      CALL ABOR1( "The spp model configuration for "//trim(YDCONF%CSPP_MODEL_NAME)//" is not known.")
    END SELECT
  
    CALL POSNAM( NULNAM, 'NAM_SPP_ACTIVE')
    CALL GET_ACTIVE_SPP_PERTS( YDCONF%SM, NULNAM, NULOUT)

    !
    !   First pass through NAM_SPP_MODIF to map indices (IDX_MAP)
    !
    CALL POSNAM( NULNAM, 'NAM_SPP_MODIF')
    CALL MAP_INDICES_SPP_NML( YDCONF%SM, NULNAM, NULOUT, IDX_MAP)

    !
    !   Second pass through NAM_SPP_MODIF in order to read modifications of 
    !     YDCONF%SM%PN
    !
    CALL POSNAM( NULNAM, 'NAM_SPP_MODIF')
    CALL MODIFY_SPP_PERTS( YDCONF%SM, IDX_MAP, NULNAM, NULOUT)

    WRITE( NULOUT ,*) 'The defined perturbations are '
    CALL WRITE_SPP_MODEL_TABLE( YDCONF%SM, NULOUT, ldefined=.true.)
    write( NULOUT,*) 'The active perturbations are '
    CALL WRITE_SPP_MODEL_TABLE( YDCONF%SM, NULOUT)

    !Active perturbations
    WRITE(NULOUT,'(''   ** pattern related settings **'')') 

    !***********************************************CONTINUE HERE*******************************    !

    WRITE(NULOUT,'(''XLCOR='',E10.4)') YDCONF%SM%XLCOR
    WRITE(NULOUT,'(''  TAU='',E10.4)') YDCONF%SM%TAU
    !
    WRITE(NULOUT,'(''XPRESS_PHR_ST='',E12.5)') YDCONF%XPRESS_PHR_ST
    !
    ! etc. 
    WRITE(NULOUT,'(''      NPATFR='',I6)') YDCONF%NPATFR
    !
    WRITE(NULOUT,'(''      LRDPATINIT='',L1,'', SPP_RDPATINIT='',A)') &
       &                 YDCONF%LRDPATINIT,     TRIM(YDCONF%SPP_RDPATINIT)
    WRITE(NULOUT,'(''      LWRPATTRUN='',L1,'', SPP_WRPATTRUN='',A)') &
       &                 YDCONF%LWRPATTRUN,     TRIM(YDCONF%SPP_WRPATTRUN)
    IF (YDCONF%LWRPATTRUN) THEN
      WRITE(NULOUT,'(''        NWRMAX='',I8,'', NWRPATTRUN='',I8)') &
         & NWRMAX, YDCONF%NWRPATTRUN
      DO JPATT=1,YDCONF%NWRPATTRUN
        WRITE(NULOUT,'(''        NHOUR_PATTRUN(xx)='',I8,'', NSTEP_PATTRUN(xx)='',I8)') &
           & YDCONF%NHOUR_PATTRUN(JPATT), YDCONF%NSTEP_PATTRUN(JPATT)
      ENDDO
    ENDIF
    WRITE(NULOUT,'(''LRESETSEED    ='',L1)')   YDCONF%LRESETSEED
    WRITE(NULOUT,'(''RESETSEEDFRQ  ='',I12)')  YDCONF%RESETSEEDFRQ
    WRITE(NULOUT,'(''LABSTIMSEED   ='',L1)')   YDCONF%LABSTIMSEED
    WRITE(NULOUT,'(''SHIFTSEED     ='',I12)')  YDCONF%SHIFTSEED

  ENDIF  !LSPP
  !
  IF (YDCONF%LABSTIMSEED.AND..NOT.YDCONF%LUSE_SETRAN) THEN
    CALL ABOR1('RANDOM SEED SETUP DOES NOT MAKE SENSE')
  ENDIF
  IF (LHOOK) CALL DR_HOOK('GET_SPP_CONF',1,ZHOOK_HANDLE)
END SUBROUTINE GET_SPP_CONF
