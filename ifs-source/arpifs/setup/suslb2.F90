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

SUBROUTINE SUSLB2(YDML_DYN)

!**** *SUSLB2  - Setup relative pointers of PTRSLB1 and PTRSLB2 for 2D model.

!     Purpose.
!     --------
!           Initialize PTRSLB1 and PTRSLB2, modules that contain
!           the relative pointers for the semi-lagrangian buffers.
!           These are used in CPG2 and CPG2LAG to compute the absolute
!           pointers required. Also computes the parities needed for
!           the extension of semi-lagrangian buffer 1 over the poles.
!           Compute total number of fields in semi-lagrangian buffers.

!**   Interface.
!     ----------
!        *CALL* *SUSLB2*

!        Explicit arguments : none
!        --------------------

!        Implicit arguments :
!        --------------------
!        modules PTRSLB1, PTRSLB2

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!      None
!      Is called by SUSC2.

!     Reference.
!     ----------

!     Author.
!     -------
!      K. YESSAD (after SUSLB): move all 2D calculations of SUSLB in SUSLB2.
!      Original : 98-08-04

!     Modifications.
!     --------------
!      Modified 01-08-30 by K. YESSAD: pruning and some other cleanings.
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      K. Yessad (Nov 2009): prune lpc_old.
!      T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!      K. Yessad (July 2014): Move some variables.
!     ------------------------------------------------------------------

USE MODEL_DYNAMICS_MOD , ONLY : MODEL_DYNAMICS_TYPE
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN   , ONLY : NULOUT
USE YOMMP0   , ONLY : NPRINTLEV
USE YOMCT0   , ONLY : LALLOPR, NCONF, NUNDEFLD

!     ------------------------------------------------------------------

IMPLICIT NONE

 TYPE(MODEL_DYNAMICS_TYPE),INTENT(INOUT):: YDML_DYN
INTEGER(KIND=JPIM) ::  IFLDSLB2, IPTX, IU

LOGICAL :: LLP
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUSLB2',0,ZHOOK_HANDLE)
ASSOCIATE(YDDYN=>YDML_DYN%YRDYN,YDDYNA=>YDML_DYN%YRDYNA, &
 & YDPTRSLB1=>YDML_DYN%YRPTRSLB1, YDPTRSLB15=>YDML_DYN%YRPTRSLB15, &
 & YDPTRSLB2=>YDML_DYN%YRPTRSLB2)
ASSOCIATE(NVLAG=>YDDYN%NVLAG, NWLAG=>YDDYN%NWLAG, &
 & MSLB1C0=>YDPTRSLB1%MSLB1C0, MSLB1C9=>YDPTRSLB1%MSLB1C9, &
 & MSLB1GFL9=>YDPTRSLB1%MSLB1GFL9, MSLB1GFLSP9=>YDPTRSLB1%MSLB1GFLSP9, &
 & MSLB1PD0=>YDPTRSLB1%MSLB1PD0, MSLB1PD9=>YDPTRSLB1%MSLB1PD9, &
 & MSLB1SP0=>YDPTRSLB1%MSLB1SP0, MSLB1SP9=>YDPTRSLB1%MSLB1SP9, &
 & MSLB1T0=>YDPTRSLB1%MSLB1T0, MSLB1T9=>YDPTRSLB1%MSLB1T9, &
 & MSLB1U0=>YDPTRSLB1%MSLB1U0, MSLB1U9=>YDPTRSLB1%MSLB1U9, &
 & MSLB1UR0=>YDPTRSLB1%MSLB1UR0, MSLB1UR9=>YDPTRSLB1%MSLB1UR9, &
 & MSLB1V0=>YDPTRSLB1%MSLB1V0, MSLB1V9=>YDPTRSLB1%MSLB1V9, &
 & MSLB1VD0=>YDPTRSLB1%MSLB1VD0, MSLB1VD9=>YDPTRSLB1%MSLB1VD9, &
 & MSLB1VR0=>YDPTRSLB1%MSLB1VR0, MSLB1VR9=>YDPTRSLB1%MSLB1VR9, &
 & MSLB1WR0=>YDPTRSLB1%MSLB1WR0, MSLBUF1=>YDPTRSLB1%MSLBUF1, &
 & NFLDSLB1=>YDPTRSLB1%NFLDSLB1, MSLB1SP05=>YDPTRSLB15%MSLB1SP05, &
 & MSLB1SP95=>YDPTRSLB15%MSLB1SP95, &
 & MSLB1U05=>YDPTRSLB15%MSLB1U05, MSLB1U95=>YDPTRSLB15%MSLB1U95, &
 & MSLB1UR05=>YDPTRSLB15%MSLB1UR05, MSLB1V05=>YDPTRSLB15%MSLB1V05, &
 & MSLB1V95=>YDPTRSLB15%MSLB1V95, MSLB1VR05=>YDPTRSLB15%MSLB1VR05, &
 & MSLBUF15=>YDPTRSLB15%MSLBUF15, NFLDSLB15=>YDPTRSLB15%NFLDSLB15, &
 & MSLB2PDSI=>YDPTRSLB2%MSLB2PDSI, MSLB2SPSI=>YDPTRSLB2%MSLB2SPSI, &
 & MSLB2TSI=>YDPTRSLB2%MSLB2TSI, MSLB2U15=>YDPTRSLB2%MSLB2U15, &
 & MSLB2URL=>YDPTRSLB2%MSLB2URL, MSLB2URL5=>YDPTRSLB2%MSLB2URL5, &
 & MSLB2USI=>YDPTRSLB2%MSLB2USI, MSLB2USI5=>YDPTRSLB2%MSLB2USI5, &
 & MSLB2V15=>YDPTRSLB2%MSLB2V15, MSLB2VDSI=>YDPTRSLB2%MSLB2VDSI, &
 & MSLB2VRL=>YDPTRSLB2%MSLB2VRL, MSLB2VRL5=>YDPTRSLB2%MSLB2VRL5, &
 & MSLB2VSI=>YDPTRSLB2%MSLB2VSI, MSLB2VSI5=>YDPTRSLB2%MSLB2VSI5, &
 & MSLB2VVEL=>YDPTRSLB2%MSLB2VVEL, MSLB2WRL=>YDPTRSLB2%MSLB2WRL, &
 & MSLBUF2=>YDPTRSLB2%MSLBUF2, NFLDSLB2=>YDPTRSLB2%NFLDSLB2, &
 & RPARSL1=>YDPTRSLB1%RPARSL1, RPARSL15=>YDPTRSLB15%RPARSL15)

!     ------------------------------------------------------------------

!*    1.   COMPUTE RELATIVE POINTERS FOR SEMI LAGRANGIAN BUFFERS.
!          -----------------------------------------------------

LLP = NPRINTLEV >= 1.OR. LALLOPR
IU = NULOUT
IF (LLP) WRITE (IU,*)
IF (LLP) WRITE (IU,'(A,/)') ' SUSLB2 PRINTOUT : POINTERS,NFIELDS'

!*       1.1  COMPUTE POINTERS FOR SLBUF2

!*       Preliminary set pointers to prevent use of uninitialised variables.
!*       (All the pointers, including those used only in the 3D model).

MSLBUF2  =NUNDEFLD
MSLB2USI =NUNDEFLD
MSLB2VSI =NUNDEFLD
MSLB2TSI =NUNDEFLD
MSLB2PDSI=NUNDEFLD
MSLB2VDSI=NUNDEFLD
MSLB2SPSI=NUNDEFLD
MSLB2VVEL=NUNDEFLD
MSLB2URL =NUNDEFLD
MSLB2VRL =NUNDEFLD
MSLB2WRL =NUNDEFLD
MSLB2URL5=NUNDEFLD
MSLB2VRL5=NUNDEFLD
MSLB2USI5=NUNDEFLD
MSLB2VSI5=NUNDEFLD
MSLB2U15 =NUNDEFLD
MSLB2V15 =NUNDEFLD

!*       Definite computation of useful pointers.

IF (LLP) WRITE (IU,'(A)') ' 1) BUFFER SLBUF2:'
IPTX    = 1
MSLBUF2 = IPTX

IF (YDDYNA%LSLAG) THEN
  MSLB2USI  = IPTX
  IF (LLP) WRITE(IU,'(1X,A12,2I4)') 'MSLB2USI  = ',IPTX,1
  IPTX      = IPTX + 1
  MSLB2VSI  = IPTX
  IF (LLP) WRITE(IU,'(1X,A12,2I4)') 'MSLB2VSI  = ',IPTX,1
  IPTX      = IPTX + 1
  MSLB2URL  = IPTX
  IF (LLP) WRITE(IU,'(1X,A12,2I4)') 'MSLB2URL  = ',IPTX,1
  IPTX      = IPTX + 1
  MSLB2VRL  = IPTX
  IF (LLP) WRITE(IU,'(1X,A12,2I4)') 'MSLB2VRL  = ',IPTX,1
  IPTX      = IPTX + 1
  IF (NCONF == 421.OR.NCONF == 521) THEN
    MSLB2URL5 = IPTX
    IF (LLP) WRITE(IU,'(1X,A12,2I4)') 'MSLB2URL5 = ',IPTX,1
    IPTX      = IPTX + 1
    MSLB2VRL5 = IPTX
    IF (LLP) WRITE(IU,'(1X,A12,2I4)') 'MSLB2VRL5 = ',IPTX,1
    IPTX      = IPTX + 1
    MSLB2USI5 = IPTX
    IF (LLP) WRITE(IU,'(1X,A12,2I4)') 'MSLB2USI5 = ',IPTX,1
    IPTX      = IPTX + 1
    MSLB2VSI5 = IPTX
    IF (LLP) WRITE(IU,'(1X,A12,2I4)') 'MSLB2VSI5 = ',IPTX,1
    IPTX      = IPTX + 1
    MSLB2U15  = IPTX
    IF (LLP) WRITE(IU,'(1X,A12,2I4)') 'MSLB2U15  = ',IPTX,1
    IPTX      = IPTX + 1
    MSLB2V15  = IPTX
    IF (LLP) WRITE(IU,'(1X,A12,2I4)') 'MSLB2V15  = ',IPTX,1
    IPTX      = IPTX + 1
  ENDIF
ENDIF

IFLDSLB2 = IPTX-1
NFLDSLB2 = IFLDSLB2

!*       1.2  COMPUTE POINTERS FOR SLBUF1

!*       Preliminary set pointers to prevent use of uninitialised variables.
!*       (All the pointers, including those used only in the 3D model).

MSLBUF1  =NUNDEFLD
MSLBUF15 =NUNDEFLD
MSLB1U9  =NUNDEFLD
MSLB1V9  =NUNDEFLD
MSLB1T9  =NUNDEFLD
MSLB1GFL9=NUNDEFLD
MSLB1GFLSP9 =NUNDEFLD
MSLB1PD9 =NUNDEFLD
MSLB1VD9 =NUNDEFLD
MSLB1UR0 =NUNDEFLD
MSLB1VR0 =NUNDEFLD
MSLB1WR0 =NUNDEFLD
MSLB1UR9 =NUNDEFLD
MSLB1VR9 =NUNDEFLD
MSLB1U0  =NUNDEFLD
MSLB1V0  =NUNDEFLD
MSLB1T0  =NUNDEFLD
MSLB1PD0 =NUNDEFLD
MSLB1VD0 =NUNDEFLD
MSLB1C9  =NUNDEFLD
MSLB1SP9 =NUNDEFLD
MSLB1SP0 =NUNDEFLD
MSLB1C0  =NUNDEFLD
MSLB1UR05=NUNDEFLD
MSLB1VR05=NUNDEFLD
MSLB1U05 =NUNDEFLD
MSLB1V05 =NUNDEFLD
MSLB1SP05=NUNDEFLD
MSLB1U95 =NUNDEFLD
MSLB1V95 =NUNDEFLD
MSLB1SP95=NUNDEFLD

!*       Definite computation of useful pointers.

!        Order is:
!         - first: high-order interpolated quantities.
!         - second: linearly interpolated quantities.
!         - then repeat for trajectory variables in case of TL/AD

IF(YDDYNA%LSLAG) THEN

  IF (LLP) WRITE (IU,*)
  IF (LLP) WRITE (IU,'(A)') ' 2) BUFFER SLBUF1:'

  IPTX    = 1
  MSLBUF1 = IPTX

!        * Buffers for high order interpolations.

  MSLB1U9  = IPTX
  IF (LLP) WRITE(IU,'(1X,A14,2I4)') 'MSLB1U9     = ',IPTX,1
  IPTX     = IPTX + 1
  MSLB1V9  = IPTX
  IF (LLP) WRITE(IU,'(1X,A14,2I4)') 'MSLB1V9     = ',IPTX,1
  IPTX     = IPTX + 1
  MSLB1SP9 = IPTX
  IF (LLP) WRITE(IU,'(1X,A14,2I4)') 'MSLB1SP9    = ',IPTX,1
  IPTX     = IPTX + 1

!        * Buffers for linear interpolations.

  MSLB1UR0 = IPTX
  IF (LLP) WRITE(IU,'(1X,A14,2I4)') 'MSLB1UR0    = ',IPTX,1
  IPTX     = IPTX + 1
  MSLB1VR0 = IPTX
  IF (LLP) WRITE(IU,'(1X,A14,2I4)') 'MSLB1VR0    = ',IPTX,1
  IPTX     = IPTX + 1
  IF(NWLAG == 3) THEN
    MSLB1U0  = IPTX
    IF (LLP)WRITE(IU,'(1X,A14,2I4)') 'MSLB1U0     = ',IPTX,1
    IPTX     = IPTX + 1
    MSLB1V0  = IPTX
    IF (LLP)WRITE(IU,'(1X,A14,2I4)') 'MSLB1V0     = ',IPTX,1
    IPTX     = IPTX + 1
  ENDIF
  IF (NVLAG == 3) THEN
    MSLB1SP0 = IPTX
    IF (LLP)WRITE(IU,'(1X,A14,2I4)') 'MSLB1SP0    = ',IPTX,1
    IPTX     = IPTX + 1
  ENDIF
  NFLDSLB1 = IPTX - 1

  IF (NCONF == 421.OR.NCONF == 521) THEN

    IPTX    = 1
    MSLBUF15 = IPTX
!          * Buffers for high order interpolations.

    MSLB1U95 = IPTX
    IF (LLP) WRITE(IU,'(1X,A14,2I4)') 'MSLB1U95    = ',IPTX,1
    IPTX     = IPTX + 1
    MSLB1V95 = IPTX
    IF (LLP) WRITE(IU,'(1X,A14,2I4)') 'MSLB1V95    = ',IPTX,1
    IPTX     = IPTX + 1
    MSLB1SP95= IPTX
    IF (LLP) WRITE(IU,'(1X,A14,2I4)') 'MSLB1SP95   = ',IPTX,1
    IPTX     = IPTX + 1

!          * Buffers for linear interpolations.

    MSLB1UR05= IPTX
    IF (LLP) WRITE(IU,'(1X,A14,2I4)') 'MSLB1UR05   = ',IPTX,1
    IPTX     = IPTX + 1
    MSLB1VR05= IPTX
    IF (LLP) WRITE(IU,'(1X,A14,2I4)') 'MSLB1VR05   = ',IPTX,1
    IPTX     = IPTX + 1
    IF(NWLAG == 3) THEN
      MSLB1U05 = IPTX
      IF (LLP)WRITE(IU,'(1X,A14,2I4)') 'MSLB1U05    = ',IPTX,1
      IPTX     = IPTX + 1
      MSLB1V05 = IPTX
      IF (LLP)WRITE(IU,'(1X,A14,2I4)') 'MSLB1V05    = ',IPTX,1
      IPTX     = IPTX + 1
    ENDIF
    IF (NVLAG == 3) THEN
      MSLB1SP05= IPTX
      IF (LLP)WRITE(IU,'(1X,A14,2I4)') 'MSLB1SP05   = ',IPTX,1
      IPTX     = IPTX + 1
    ENDIF
    NFLDSLB15 = IPTX - 1
  ELSE
    NFLDSLB15 = 1
  ENDIF

ELSE

  NFLDSLB1  = 1
  NFLDSLB15 = 1

ENDIF

IF (LLP) WRITE (IU,*)

!     ------------------------------------------------------------------

!*    2.   COMPUTE PARITIES FOR EXTENSION OVER POLE.
!          -----------------------------------------

IF (YDDYNA%LSLAG) THEN
  ALLOCATE(YDPTRSLB1%RPARSL1(NFLDSLB1))
  IF(LLP)WRITE(IU,9) 'RPARSL1  ',SIZE(YDPTRSLB1%RPARSL1),SHAPE(YDPTRSLB1%RPARSL1)
  ALLOCATE(YDPTRSLB15%RPARSL15(YDPTRSLB15%NFLDSLB15))
  IF(LLP)WRITE(IU,9) 'RPARSL1  ',SIZE(YDPTRSLB15%RPARSL15),SHAPE(YDPTRSLB15%RPARSL15)
ENDIF
9 FORMAT(1X,'ARRAY ',A10,' ALLOCATED ',8I8)

IF (YDDYNA%LSLAG) THEN

  RPARSL1(MSLB1U9)=-1.0_JPRB
  RPARSL1(MSLB1V9)=-1.0_JPRB
  RPARSL1(MSLB1SP9)=+1.0_JPRB
  RPARSL1(MSLB1UR0)=-1.0_JPRB
  RPARSL1(MSLB1VR0)=-1.0_JPRB
  IF(NWLAG == 3) THEN
    RPARSL1(MSLB1U0)=-1.0_JPRB
    RPARSL1(MSLB1V0)=-1.0_JPRB
  ENDIF
  IF (NVLAG == 3) THEN
    RPARSL1(MSLB1SP0)=+1.0_JPRB
  ENDIF
  IF (NCONF == 421.OR.NCONF == 521) THEN
    RPARSL15(MSLB1U95)=-1.0_JPRB
    RPARSL15(MSLB1V95)=-1.0_JPRB
    RPARSL15(MSLB1SP95)=+1.0_JPRB
    RPARSL15(MSLB1UR05)=-1.0_JPRB
    RPARSL15(MSLB1VR05)=-1.0_JPRB
    IF(NWLAG == 3) THEN
      RPARSL15(MSLB1U05)=-1.0_JPRB
      RPARSL15(MSLB1V05)=-1.0_JPRB
    ENDIF
    IF (NVLAG == 3) THEN
      RPARSL15(MSLB1SP05)=+1.0_JPRB
    ENDIF
  ENDIF

ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUSLB2',1,ZHOOK_HANDLE)
END SUBROUTINE SUSLB2
