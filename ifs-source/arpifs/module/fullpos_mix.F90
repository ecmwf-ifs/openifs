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

MODULE FULLPOS_MIX

! Purpose :
! -------
!    To define the user type "FULLPOS_TYPE" for the post-processing control :
!     - %CLNAME : the ARPEGE field name,
!     - %LLSRF  : the kind of field : surface (.T.) or upper air (.F.)
!     - %IGRIB  : the GRIB code,
!     - %IBITS  : the number of bits for packing before writing out to file,
!     - %INTER  : the kind of horizontal interpolation
!                 (quadratic, bilinear, nearest point)
!     - %IORDR  : (INTERNAL) the horizontal derivative order of the field 
!                  0=scalar, 
!                  1=vector/U component,
!                 -1=vector/V component,
!                  2=horizontal derivative on vector quantity
!     - %CLPAIR : keyworkd for a pair of vector components or for a pair of
!                 any fields which cannot be computed independantly
!     - %LLGP   : the horizontal representation of the field :
!                  gridpoint (T) or spectral (F)
!     - %ISF    : Indicator for spectral computation :
!                 1 = spectral fit only
!                 2 = spectral fit + gaussian filter in the model spectral
!                     space
!                 3 = spectral fit + gaussian filter in the spectral space
!                     of homogenous resolution (stretched model only)
!     - %ZFK    : tuning coefficient for the gaussian filter in the model
!                 spectral space
!     - %ICOD   : an INTERNAL code number to localize a specific field among
!                 all fields
!     - %CLNIL  : the INTERNAL list of STEPO(5:5) configurations not working
!                 for this field
!     - %LLBIP  : (INTERNAL) tells if the field needs to be made biperiodic
!                 between the vertical and the horizontal interpolations.
!                 This is the case for such fields as moisture convergence
!                 which is build by the spectral transforms package as a field of divergence
         
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
!    Fullpos technical & users guide.

! Author :
! ------
!    Ryad El Khatib *METEO-FRANCE* thanks to Mike Fisher *ECMWF*

! Modifications :
! -------------
! Original : 2000-08-18
! R. El Khatib : 01-03-28 Redefinition of QFPTYPE%ILED
! R. El Khatib : 03-02-05 Split module into fullpos_mix + type_fprqdyns
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
! P.Marguinaud : 10-10-13 Add structure members for reduced coupling fields
!                         Fix INTENT in CTRL_TFP
!-----------------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

IMPLICIT NONE
SAVE

PRIVATE
PUBLIC FULLPOS_TYPE, TFP_SU, CTRL_TFP

TYPE FULLPOS_TYPE 

INTEGER(KIND=JPIM) :: IGRIB = -99999
INTEGER(KIND=JPIM) :: IBITS = 99
INTEGER(KIND=JPIM) :: INTER = -9
INTEGER(KIND=JPIM) :: IORDR = 0
INTEGER(KIND=JPIM) :: ISF = 0
INTEGER(KIND=JPIM) :: ICOD = 0
REAL(KIND=JPRB)    :: ZFK = 0._JPRB
LOGICAL            :: LLGP = .TRUE.
LOGICAL            :: LLSRF
LOGICAL            :: LLBIP =.FALSE.
CHARACTER(LEN=16) :: CLNIL= ' '
CHARACTER(LEN=16) :: CLNAME = ' '
CHARACTER(LEN=8)  :: CLPAIR = ' '

! Levels with hollow fields
INTEGER(KIND=JPIM) :: ILEVHOLI(2) = (/ -1_JPIM, -1_JPIM /) ! First levels
INTEGER(KIND=JPIM) :: ILEVHOLF(2) = (/ -1_JPIM, -1_JPIM /) ! End levels
! For hollow fields, compaction parameters
INTEGER(KIND=JPIM) :: ICPLSIZE = 10_JPIM
INTEGER(KIND=JPIM) :: ICPLBITS = -1_JPIM

END TYPE FULLPOS_TYPE

CONTAINS

!-----------------------------------------------------------------------------

TYPE(FULLPOS_TYPE) FUNCTION TFP_SU &
 & (LDSRF,CDNAME,KGRIB,KBITS,KBITG,KNTER,LDGP,KSF,PFK,KORDR,CDPAIR,LDBIP)

! Purpose :
! -------
!    To set default values to the type

LOGICAL,            INTENT(IN) :: LDSRF
CHARACTER(LEN=*),   INTENT(IN) :: CDNAME
INTEGER(KIND=JPIM), INTENT(IN) :: KGRIB
INTEGER(KIND=JPIM), INTENT(IN) :: KBITS
INTEGER(KIND=JPIM), INTENT(IN) :: KBITG
INTEGER(KIND=JPIM), INTENT(IN) :: KNTER
LOGICAL,            INTENT(IN) :: LDGP
INTEGER(KIND=JPIM), INTENT(IN) :: KSF
REAL(KIND=JPRB),    INTENT(IN) :: PFK
INTEGER(KIND=JPIM), INTENT(IN) :: KORDR
CHARACTER(LEN=*),   INTENT(IN) :: CDPAIR
LOGICAL,            INTENT(IN), OPTIONAL :: LDBIP

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('FULLPOS_MIX:TFP_SU',0,ZHOOK_HANDLE)
TFP_SU%CLNAME=CDNAME
TFP_SU%CLPAIR=CDPAIR
TFP_SU%IGRIB=KGRIB
TFP_SU%INTER=KNTER
TFP_SU%IORDR=KORDR
TFP_SU%LLGP=LDGP
TFP_SU%LLSRF=LDSRF
TFP_SU%ICOD=0
TFP_SU%CLNIL=' '
SELECT CASE (LDGP)
CASE (.TRUE.)
  TFP_SU%ISF=0
  TFP_SU%ZFK=0.0_JPRB
  TFP_SU%IBITS=KBITG
CASE DEFAULT
  TFP_SU%ISF=MAX(1,MIN(3,KSF))
  TFP_SU%ZFK=MAX(0.0_JPRB,PFK)
  TFP_SU%IBITS=KBITS
END SELECT 
IF (PRESENT(LDBIP)) THEN
  TFP_SU%LLBIP=LDBIP
ENDIF
IF (LHOOK) CALL DR_HOOK('FULLPOS_MIX:TFP_SU',1,ZHOOK_HANDLE)

END FUNCTION TFP_SU

!-----------------------------------------------------------------------------

SUBROUTINE CTRL_TFP(YDFP_X,KCOD,CDNIL,LDSF,YDFP_DYNDS,KMAX)

! LDSF : .FALSE. to disable spectral filters.

INTEGER(KIND=JPIM),          INTENT(IN)    :: KMAX
INTEGER(KIND=JPIM),          INTENT(INOUT) :: KCOD
TYPE(FULLPOS_TYPE), INTENT(INOUT) :: YDFP_X
TYPE(FULLPOS_TYPE), INTENT(INOUT) :: YDFP_DYNDS(KMAX) ! A single element is modified
LOGICAL,            INTENT(IN)    :: LDSF
CHARACTER(LEN=*),   INTENT(IN)    :: CDNIL
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('FULLPOS_MIX:CTRL_TFP',0,ZHOOK_HANDLE)
SELECT CASE (YDFP_X%LLGP)
CASE (.TRUE.)
  YDFP_X%ISF=0
  YDFP_X%ZFK=0.0_JPRB
CASE DEFAULT
  SELECT CASE (LDSF)
  CASE (.TRUE.)
    YDFP_X%ISF=MAX(1,MIN(3,YDFP_X%ISF))
    YDFP_X%ZFK=MAX(0.0_JPRB,YDFP_X%ZFK)
  CASE DEFAULT
    YDFP_X%ISF=MIN(1,YDFP_X%ISF)
    YDFP_X%ZFK=0.0_JPRB
  END SELECT 
END SELECT

YDFP_X%CLNIL=CDNIL

IF (KCOD >= 0) THEN
  KCOD=KCOD+1
  YDFP_X%ICOD=KCOD
  IF (KCOD <= KMAX) THEN
    YDFP_DYNDS(KCOD)=YDFP_X
  ENDIF
ENDIF
IF (LHOOK) CALL DR_HOOK('FULLPOS_MIX:CTRL_TFP',1,ZHOOK_HANDLE)

END SUBROUTINE CTRL_TFP

END MODULE FULLPOS_MIX
