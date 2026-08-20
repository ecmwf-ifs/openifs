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

MODULE TYPES_FPCAT

! Purpose :
! -------
!    To define the user type "TYPES_FPCAT" for the post-processing of individual CAT/MTW indices
!    and the final combinaison of indices to compute diagnostic EDR.

!     - %CLNAME        :  the  indice name,
!     - %ZDEL[H-M-L]   :  Delta for EDR transpose
!     - %ZRAT[H-M-L]   :  Range for EDR transpose
!     - %Z[DEL-RAT]H   :  High levels
!     - %Z[DEL-RAT]M   :  Medium levels
!     - %Z[DEL-RAT]L   :  Low levels
!     - %LEDR          :  Indicate conversion into EDR unit when indice is post-processed
!     - %LCOMB[H-M-L]  :  Used for final EDR diag combinaison for each classes or layer
!     - %LCOMB         :  At least one layer needs the indice for the combinaison

!    Also define a super-type 'ALL_FPCAT_TYPES' with all TYPE_FPCAT objets

!     - %XXX           :  Individual indice of CAT
!      
!     - %LCONVEDR      : Global key to control conversion of individual indices into EDR unit when they are post-processed.
!                        When LCONVEDR =.T. if at least one individual indice LEDR is .T.
         
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
!    Olivier Jaron *METEO-FRANCE*

! Modifications :
! -------------
! Original : 2019-04-02

!-----------------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPRB, JPIM
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK

IMPLICIT NONE
SAVE

TYPE TYPE_FPCAT

LOGICAL  :: LEDR = .FALSE.

CHARACTER(LEN=16) :: CLNAME = 'XXXXXXXX'
INTEGER(KIND=JPIM):: ID = 0_JPIM
REAL(KIND=JPRB)   :: ZDELH = 0._JPRB
REAL(KIND=JPRB)   :: ZRATH = 1._JPRB
REAL(KIND=JPRB)   :: ZDELM = 0._JPRB
REAL(KIND=JPRB)   :: ZRATM = 1._JPRB
REAL(KIND=JPRB)   :: ZDELL = 0._JPRB
REAL(KIND=JPRB)   :: ZRATL = 1._JPRB

LOGICAL :: LCOMBH = .FALSE.
LOGICAL :: LCOMBM = .FALSE.
LOGICAL :: LCOMBL = .FALSE.
LOGICAL :: LCOMB  = .FALSE.

END TYPE TYPE_FPCAT

TYPE ALL_FPCAT_TYPES

LOGICAL :: LCONVEDR  = .FALSE. 
!REAL(KIND=JPRB):: PPRESH = 14483._JPRB  ! Pression de sommet de la couche HIGH  # FL450
REAL(KIND=JPRB):: PPRESH =  7170._JPRB  ! Pression de sommet de la couche HIGH  # FL600
REAL(KIND=JPRB):: PPRESM = 46569._JPRB  ! Pression de transition MEDIUM -> HIGH # FL200
REAL(KIND=JPRB):: PPRESL = 69686._JPRB  ! Pression de transition LOW -> MEDIUM  # FL100
REAL(KIND=JPRB):: PPDELTA = 5000._JPRB  ! Epaisseur (Pa) de la transition

TYPE(TYPE_FPCAT) :: BR1
TYPE(TYPE_FPCAT) :: BR1R
TYPE(TYPE_FPCAT) :: BR1PV
TYPE(TYPE_FPCAT) :: BR2
TYPE(TYPE_FPCAT) :: DUT
TYPE(TYPE_FPCAT) :: LAZ
TYPE(TYPE_FPCAT) :: TI1
TYPE(TYPE_FPCAT) :: TI1PV
TYPE(TYPE_FPCAT) :: TI2
TYPE(TYPE_FPCAT) :: TI2PV
TYPE(TYPE_FPCAT) :: CP
TYPE(TYPE_FPCAT) :: RI
TYPE(TYPE_FPCAT) :: INVRI
TYPE(TYPE_FPCAT) :: INVRIM
TYPE(TYPE_FPCAT) :: HSR
TYPE(TYPE_FPCAT) :: LRT
TYPE(TYPE_FPCAT) :: DF
TYPE(TYPE_FPCAT) :: DFR
TYPE(TYPE_FPCAT) :: DFPV
TYPE(TYPE_FPCAT) :: INVRITW
TYPE(TYPE_FPCAT) :: DVR
TYPE(TYPE_FPCAT) :: AGI
TYPE(TYPE_FPCAT) :: F2D
TYPE(TYPE_FPCAT) :: F2DR
TYPE(TYPE_FPCAT) :: F2DPV
TYPE(TYPE_FPCAT) :: NG1
TYPE(TYPE_FPCAT) :: NG2
TYPE(TYPE_FPCAT) :: NG1R
TYPE(TYPE_FPCAT) :: NG2R
TYPE(TYPE_FPCAT) :: TKER
TYPE(TYPE_FPCAT) :: TGR
TYPE(TYPE_FPCAT) :: IAW
TYPE(TYPE_FPCAT) :: IAWPV
TYPE(TYPE_FPCAT) :: IAWR
TYPE(TYPE_FPCAT) :: SV
TYPE(TYPE_FPCAT) :: BV
TYPE(TYPE_FPCAT) :: BVM
TYPE(TYPE_FPCAT) :: VSQR
TYPE(TYPE_FPCAT) :: DPV
TYPE(TYPE_FPCAT) :: VVSQRI
TYPE(TYPE_FPCAT) :: VVSQ
TYPE(TYPE_FPCAT) :: WND
TYPE(TYPE_FPCAT) :: EDR
TYPE(TYPE_FPCAT) :: TKE
TYPE(TYPE_FPCAT) :: LUN
TYPE(TYPE_FPCAT) :: VVSQMW
TYPE(TYPE_FPCAT) :: F2DMW
TYPE(TYPE_FPCAT) :: WNDMW
TYPE(TYPE_FPCAT) :: DVRMW
TYPE(TYPE_FPCAT) :: NG1MW
TYPE(TYPE_FPCAT) :: IAWMW
TYPE(TYPE_FPCAT) :: DFMW
! Final CAT/MTW/MAX combinaison
TYPE(TYPE_FPCAT) :: EDRDC
TYPE(TYPE_FPCAT) :: EDRDW
TYPE(TYPE_FPCAT) :: EDRDM


END TYPE ALL_FPCAT_TYPES
!-----------------------

CONTAINS

SUBROUTINE SUCATIND(YDCAT_X,CDNAME,LDEDR,PDELTA,PRATIO,LDCOMB,KID)


! Purpose :
! -------
!    To set default values to the type TYPE_FPCAT

!     PDELRA  : Delta for the three classes of level
!     LDCOMB  : Indicates if indece should be combined, for the three classes od level

TYPE(TYPE_FPCAT) , INTENT(INOUT) :: YDCAT_X
CHARACTER(LEN=*), INTENT(IN) :: CDNAME
INTEGER(KIND=JPIM) , INTENT(IN), OPTIONAL ::  KID
LOGICAL , INTENT(IN) :: LDEDR 
REAL(KIND=JPRB) , INTENT(IN) :: PDELTA(1:3)
REAL(KIND=JPRB) , INTENT(IN) :: PRATIO(1:3)
LOGICAL , INTENT(IN) :: LDCOMB(1:3)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('TYPES_FPCAT:SUCATIND',0,ZHOOK_HANDLE)

YDCAT_X%CLNAME=CDNAME
IF (PRESENT(KID)) YDCAT_X%ID = KID
YDCAT_X%LEDR=LDEDR
YDCAT_X%ZDELH=PDELTA(1)
YDCAT_X%ZRATH=PRATIO(1)
YDCAT_X%ZDELM=PDELTA(2)
YDCAT_X%ZRATM=PRATIO(2)
YDCAT_X%ZDELL=PDELTA(3)
YDCAT_X%ZRATL=PRATIO(3)
YDCAT_X%LCOMBH=LDCOMB(1)
YDCAT_X%LCOMBM=LDCOMB(2)
YDCAT_X%LCOMBL=LDCOMB(3)
YDCAT_X%LCOMB=LDCOMB(1).OR.LDCOMB(2).OR.LDCOMB(3)

IF (LHOOK) CALL DR_HOOK('TYPES_FPCAT:SUCATIND',1,ZHOOK_HANDLE)

END SUBROUTINE SUCATIND

SUBROUTINE PRINT_CONFIG_EDRD(YDCAT_X, CDNAME)
USE YOMLUN       , ONLY : NULOUT
! Purpose :
! -------
!   Print configuration of EDRD combinaison (coefficients, canditate for final combinaison, 
!   convert or not in EDRD unit)
CHARACTER(LEN=12) , INTENT(IN) :: CDNAME
TYPE(TYPE_FPCAT) , INTENT(IN) :: YDCAT_X
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('TYPES_FPCAT:PRINT_CONFIG_EDRD',0,ZHOOK_HANDLE)

WRITE(NULOUT, FMT ='(1X,A12,'' : '',A11,1X,I3,F9.3,1X,F9.3,1X,F9.3,1X,F9.3,1X,F9.3,1X,F9.3,1X,&
  & 1X, L6,1X, L6,1X, L6,1X, L6, 1X, L5)') &
  & CDNAME,YDCAT_X%CLNAME,YDCAT_X%ID,YDCAT_X%ZDELH,YDCAT_X%ZRATH,YDCAT_X%ZDELM,YDCAT_X%ZRATM,YDCAT_X%ZDELL,YDCAT_X%ZRATL,&
  & YDCAT_X%LCOMBH,YDCAT_X%LCOMBM,YDCAT_X%LCOMBL,YDCAT_X%LCOMB,YDCAT_X%LEDR

IF (LHOOK) CALL DR_HOOK('TYPES_FPCAT:PRINT_CONFIG_EDRD',1,ZHOOK_HANDLE)

ENDSUBROUTINE PRINT_CONFIG_EDRD


END MODULE TYPES_FPCAT
