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

SUBROUTINE SUALDYN(YDDIMV,YDDYNA,YDDYN)

!**** *SUALDYN * - Routine to allocate space for dynamics

!     Purpose.
!     --------
!           Allocate space for the dynamics calc.
!           Not LELAM-dependent arrays.

!**   Interface.
!     ----------
!        *CALL* *SUALDYN*

!     Explicit arguments :  None
!     --------------------

!     Implicit arguments :
!     --------------------

!     Method.
!     -------

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      Mats Hamrud *ECMWF*
!      Original : 95-11-20  From SUALLO

!     Modifications.
!     --------------
!      K. Yessad (Aug 2009): remove LSITRIC option
!      K. Yessad (Jan 2010): remove useless variables.
!      M. Fisher   7-March-2012 Move SIBI out of Jb (to allow late Jb setup)
!      K. Yessad (July 2014): some reorganisation in the set-up.
!      K. Yessad (June 2017): Introduce NHQE model.
!     ------------------------------------------------------------------

USE YOMDIMV  , ONLY : TDIMV
USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
USE YOMMP0   , ONLY : NPRINTLEV
USE YOMCT0   , ONLY : LALLOPR
USE YOMLUN   , ONLY : NULOUT
USE YOMDYN   , ONLY : TDYN
USE YOMDYNA  , ONLY : TDYNA

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(TDIMV), INTENT(IN)   :: YDDIMV
TYPE(TDYNA) ,INTENT(IN)   :: YDDYNA
TYPE(TDYN)  ,INTENT(INOUT):: YDDYN

INTEGER(KIND=JPIM) :: IU

LOGICAL :: LLP
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUALDYN',0,ZHOOK_HANDLE)
ASSOCIATE(NFLEVG=>YDDIMV%NFLEVG)
!     ------------------------------------------------------------------

!*       1.    ALLOCATE SPACE FOR ARRAYS.
!              --------------------------

LLP = NPRINTLEV >= 1.OR. LALLOPR
IU = NULOUT

! * Semi-implicit scheme.
ALLOCATE(YDDYN%SIALPH(NFLEVG))
IF(LLP)WRITE(IU,9) 'SIALPH   ',SIZE(YDDYN%SIALPH),SHAPE(YDDYN%SIALPH)
ALLOCATE(YDDYN%SILNPR(NFLEVG))
IF(LLP)WRITE(IU,9) 'SILNPR   ',SIZE(YDDYN%SILNPR),SHAPE(YDDYN%SILNPR)
ALLOCATE(YDDYN%SIDELP(NFLEVG))
IF(LLP)WRITE(IU,9) 'SIDELP   ',SIZE(YDDYN%SIDELP),SHAPE(YDDYN%SIDELP)
ALLOCATE(YDDYN%SIRDEL(NFLEVG))
IF(LLP)WRITE(IU,9) 'SIRDEL   ',SIZE(YDDYN%SIRDEL),SHAPE(YDDYN%SIRDEL)
ALLOCATE(YDDYN%SITLAH(0:NFLEVG))
IF(LLP)WRITE(IU,9) 'SITLAH   ',SIZE(YDDYN%SITLAH),SHAPE(YDDYN%SITLAH)
ALLOCATE(YDDYN%SITLAF(NFLEVG))
IF(LLP)WRITE(IU,9) 'SITLAF   ',SIZE(YDDYN%SITLAF),SHAPE(YDDYN%SITLAF)
ALLOCATE(YDDYN%SIDPHI(NFLEVG))
IF(LLP)WRITE(IU,9) 'SIDPHI   ',SIZE(YDDYN%SIDPHI),SHAPE(YDDYN%SIDPHI)
ALLOCATE(YDDYN%SIB(NFLEVG,NFLEVG))
IF(LLP)WRITE(IU,9) 'SIB      ',SIZE(YDDYN%SIB),SHAPE(YDDYN%SIB)
ALLOCATE(YDDYN%SIMO(NFLEVG,NFLEVG))
IF(LLP)WRITE(IU,9) 'SIMO     ',SIZE(YDDYN%SIMO),SHAPE(YDDYN%SIMO)
ALLOCATE(YDDYN%SIMI(NFLEVG,NFLEVG))
IF(LLP)WRITE(IU,9) 'SIMI     ',SIZE(YDDYN%SIMI),SHAPE(YDDYN%SIMI)
ALLOCATE(YDDYN%SIVP(NFLEVG))
IF(LLP)WRITE(IU,9) 'SIVP     ',SIZE(YDDYN%SIVP),SHAPE(YDDYN%SIVP)
ALLOCATE(YDDYN%SIBI(NFLEVG,NFLEVG))
IF(LLP)WRITE(IU,9) 'SIBI     ',SIZE(YDDYN%SIBI),SHAPE(YDDYN%SIBI)
ALLOCATE(YDDYN%SITRAM(NFLEVG))
IF(LLP)WRITE(IU,9) 'SITRAM   ',SIZE(YDDYN%SITRAM),SHAPE(YDDYN%SITRAM)

IF (YDDYNA%LNHEE) THEN
  ALLOCATE(YDDYN%SIFAC(NFLEVG,NFLEVG))
  IF(LLP)WRITE(IU,9) 'SIFAC    ',SIZE(YDDYN%SIFAC),SHAPE(YDDYN%SIFAC)
  ALLOCATE(YDDYN%SIFACI(NFLEVG,NFLEVG))
  IF(LLP)WRITE(IU,9) 'SIFACI   ',SIZE(YDDYN%SIFACI),SHAPE(YDDYN%SIFACI)
ENDIF

IF (YDDYNA%LNHQE) THEN
  ALLOCATE(YDDYN%SI_ILAPKSSI(NFLEVG,NFLEVG,2))
  IF(LLP)WRITE(IU,8) 'SI_ILAPKSSI',SIZE(YDDYN%SI_ILAPKSSI),SHAPE(YDDYN%SI_ILAPKSSI)
ENDIF

! * Intermediate arrays to change variable for T-eqn.
ALLOCATE(YDDYN%RCORDIT(NFLEVG))
IF(LLP)WRITE(IU,9) 'RCORDIT  ',SIZE(YDDYN%RCORDIT),SHAPE(YDDYN%RCORDIT)
ALLOCATE(YDDYN%RCORDIH(0:NFLEVG))
IF(LLP)WRITE(IU,9) 'RCORDIH  ',SIZE(YDDYN%RCORDIH),SHAPE(YDDYN%RCORDIH)
ALLOCATE(YDDYN%RCORDIF(NFLEVG))
IF(LLP)WRITE(IU,9) 'RCORDIF  ',SIZE(YDDYN%RCORDIF),SHAPE(YDDYN%RCORDIF)

! * Rayleigh friction.
IF(YDDYN%LRFRIC) THEN
  ALLOCATE(YDDYN%RKRF(NFLEVG))
  IF(LLP)WRITE(IU,9) 'RKRF   ',SIZE(YDDYN%RKRF),SHAPE(YDDYN%RKRF)
ENDIF

!     ------------------------------------------------------------------
9 FORMAT(1X,'ARRAY ',A10,' ALLOCATED ',8I8)
8 FORMAT(1X,'ARRAY ',A11,' ALLOCATED ',8I8)

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUALDYN',1,ZHOOK_HANDLE)
END SUBROUTINE SUALDYN
