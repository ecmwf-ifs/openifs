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

MODULE SPP_MOD

! Purpose :
! -------
!    Define types for the SPP scheme, which represents model uncertainties
!    with Stochastically Perturbed Parameterisations.

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
!    Original : July 2020

! Modifications :
! -------------
!    
!  
!-----------------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE SPP_GEN_MOD, ONLY : SPP_MODEL
USE SPP_DEF_MOD, ONLY : SPP_PERT_POINTER
USE SPECTRAL_ARP_MOD, ONLY : SPECTRAL_ARP
USE GRIDPOINT_FIELDS_MIX, ONLY : GRIDPOINT_FIELD

IMPLICIT NONE
SAVE

INTEGER(KIND=JPIM), PARAMETER :: JPSPPLABLEN=16, NWRMAX=255
!=============================================================================

!      1.    TYPE DEFINITIONS
!            ----------------

!      1.1   TYPE holding the configuration of the spp scheme     
TYPE TSPP_CONFIG
  LOGICAL :: LSPP=.FALSE.        !> activates stochastically perturbed parameterisation scheme
  CHARACTER(LEN=32) :: CSPP_MODEL_NAME !> name of SPP model, selects define_spp routine
  TYPE(SPP_MODEL) :: SM          !> structure storing details of configuration
  TYPE(SPP_PERT_POINTER) :: PPTR !> perturbation pointer structure

  LOGICAL :: LRAMIDLIMIT1        !> limit perturbed value of RAMID to be <=1.

  LOGICAL :: LUSE_SETRAN=.TRUE.  !> use setran to randomize seed depending on start date, ensemble member number, ...
  LOGICAL :: LNSEED_OFFS0        !> set all seed offsets to 0 regardless of NSEED_OFFS_* given below

  ! write out pattern and reset seed
  LOGICAL            :: LRDPATINIT 
  LOGICAL            :: LWRPATTRUN
  INTEGER(KIND=JPIM) :: NWRPATTRUN
  INTEGER(KIND=JPIM), DIMENSION(NWRMAX) :: NHOUR_PATTRUN
  INTEGER(KIND=JPIM), DIMENSION(NWRMAX) :: NSTEP_PATTRUN
  LOGICAL            :: LRESETSEED 
  LOGICAL            :: LABSTIMSEED 
  INTEGER(KIND=JPIM) :: RESETSEEDFRQ
  INTEGER(KIND=JPIM) :: SHIFTSEED
  CHARACTER(LEN=256) :: SPP_RDPATINIT
  CHARACTER(LEN=256) :: SPP_WRPATTRUN

  INTEGER(JPIM) :: ISEEDFAC      !> seed factor for random numbers

  INTEGER(JPIM) :: NPATFR=1      !> frequency of pattern updates:
                                 !> >0: every NPATFR timesteps / <0: every -NPATFR hours

  REAL(KIND=JPRB) :: XPRESS_PHR_ST      !> pressure defining level of stratosphere, for heating rate perturbations


END TYPE TSPP_CONFIG

!     1.2    TYPE holding the data of the spp scheme
TYPE TSPP_DATA
  INTEGER(JPIM), DIMENSION(:), ALLOCATABLE :: IGRIBCODE !> array with grib codes of pattern

  TYPE (SPECTRAL_ARP),    DIMENSION(:), POINTER :: SP_ARP=>NULL()  !> array with patterns
  TYPE (GRIDPOINT_FIELD), DIMENSION(:), POINTER :: GP_ARP=>NULL()  !> array with grid-point fields
  TYPE (GRIDPOINT_FIELD), DIMENSION(:), POINTER :: GP_ARP0=>NULL() !> array with first  time-level (NPATFR>1)
  TYPE (GRIDPOINT_FIELD), DIMENSION(:), POINTER :: GP_ARP1=>NULL() !> array with second time-level (NPATFR>1)
  CHARACTER(LEN=JPSPPLABLEN)     , DIMENSION(:), POINTER :: LAB_ARP=>NULL() !> label to refer to ARP
END TYPE TSPP_DATA

contains

FUNCTION KGET_SEED_SPP( YDCONF, KPERT, KRF, KMEMBER, KSHIFT, LABSTIME, LDVERBOSE, PTSTEP)

USE YOMLUN   , ONLY : NULOUT
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
IMPLICIT NONE

INTEGER(KIND=JPIM)                       :: KGET_SEED_SPP
TYPE(TSPP_CONFIG), INTENT(IN)            :: YDCONF
INTEGER(KIND=JPIM), INTENT(IN)           :: KPERT
INTEGER(KIND=JPIM), INTENT(IN)           :: KRF
INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KMEMBER, KSHIFT
LOGICAL, OPTIONAL, INTENT(IN)            :: LDVERBOSE, LABSTIME
REAL(KIND=JPRB), OPTIONAL, INTENT(IN)    :: PTSTEP


LOGICAL :: LLVERBOSE=.false., LLABSTIME=.false.
INTEGER(KIND=JPIM) :: IMEMBER=0, ISHIFT=0
INTEGER(KIND=JPIM) :: ISEEDIN
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"
#include "setran.intfb.h"

IF (LHOOK) CALL DR_HOOK('SPP_MOD:KGET_SEED_SPP',0,ZHOOK_HANDLE)

IF (PRESENT(KMEMBER))   IMEMBER  =KMEMBER
IF (PRESENT(KSHIFT))     ISHIFT   =KSHIFT
IF (PRESENT(LDVERBOSE)) LLVERBOSE=LDVERBOSE
IF (PRESENT(LABSTIME))  LLABSTIME=LABSTIME

IF (LLABSTIME .AND. .NOT.PRESENT(PTSTEP)) THEN
  CALL ABOR1('SPP_MOD:KGET_SEED_SPP: LLABSTIME but PTSTEP not provided')
ENDIF

IF (YDCONF%LNSEED_OFFS0) THEN
  KGET_SEED_SPP=0
ELSE
  KGET_SEED_SPP= YDCONF%SM%PN(KPERT)%nseed_off(KRF) !ML: needs to be made generic for multiple seeds
ENDIF

KGET_SEED_SPP=KGET_SEED_SPP + YDCONF%ISEEDFAC + IMEMBER + ISHIFT

IF (YDCONF%LUSE_SETRAN) THEN
  ISEEDIN=KGET_SEED_SPP
  IF (LLABSTIME) THEN
    CALL SETRAN( ISEEDIN, KOUTSEED=KGET_SEED_SPP, LDABSTIME=LLABSTIME, PTSTEP=PTSTEP)
  ELSE
    CALL SETRAN( ISEEDIN, KOUTSEED=KGET_SEED_SPP)
  ENDIF
ELSE
  ISEEDIN=0
ENDIF

IF (LLVERBOSE) THEN
  WRITE(NULOUT,'(''KGET_SEED_SPP: '',A16,I4,2(X,I12))') YDCONF%SM%PN(KPERT)%LABEL, KRF, ISEEDIN, KGET_SEED_SPP
ENDIF

IF (LHOOK) CALL DR_HOOK('SPP_MOD:KGET_SEED_SPP',1,ZHOOK_HANDLE)


END FUNCTION KGET_SEED_SPP



END MODULE SPP_MOD
