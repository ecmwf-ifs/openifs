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

SUBROUTINE SPNORMB(YDLAP,YDDIM,PX,KLEV,PSN,PSM,PMET,KFLEV,KFLSUR,LDNWAVE)

!**** *SPNORMB* - Compute norms in spectral space

!     Purpose.
!     --------
!        Compute the norm in spectral space for a given level and
!        a given n (and m)

!**   Interface.
!     ----------
!        *CALL* *SPNORMB(...)

!        Explicit arguments :
!        --------------------
!        PX      : Input array
!        PSN     : spectrum at constant n
!        PSM     : spectrum at constant m
!        KLEV    : number of levels of computation
!        PMET    : Metric
!        KFLEV   : first dimensioning of output arrays.
!        KFLSUR  : first dimensioning of input array.
!        LDNWAVE : .TRUE. to compute spectrum at constant n

!        Implicit arguments :  none.
!        --------------------

!     Method.
!     -------

!     Externals. None
!     ----------
!      Called by SPNORMBM.

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      Philippe Courtier  *ECMWF*
!      Original : 92-12-20

!     Modifications.
!     --------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      R. El Khatib 04-08-09 Loops on levels and local wave numbers
!     ------------------------------------------------------------------

USE YOMLAP   , ONLY : TLAP
USE YOMDIM   , ONLY : TDIM
USE PARKIND1 , ONLY : JPIM, JPRB, JPRD
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK

IMPLICIT NONE

TYPE(TLAP)        , INTENT(IN)   :: YDLAP
TYPE(TDIM)        , INTENT(IN)   :: YDDIM
INTEGER(KIND=JPIM), INTENT(IN)   :: KFLEV 
INTEGER(KIND=JPIM), INTENT(IN)   :: KFLSUR 
REAL(KIND=JPRB)   , INTENT(IN)   :: PX(:,:)
INTEGER(KIND=JPIM), INTENT(IN)   :: KLEV 
REAL(KIND=JPRB)   , INTENT(OUT)  :: PSN(:,:)
REAL(KIND=JPRB)   , INTENT(OUT)  :: PSM(:,:)
REAL(KIND=JPRB)   , INTENT(IN)   :: PMET(0:) 
LOGICAL           , INTENT(IN)   :: LDNWAVE 
INTEGER(KIND=JPIM) :: INM, ISP, ISP0, JM, JN, JLEV, JNML
REAL(KIND=JPRD) :: ZSUM, ZEPS
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SPNORMB',0,ZHOOK_HANDLE)
ASSOCIATE(NSMAX=>YDDIM%NSMAX, NUMP=>YDDIM%NUMP, MYMS=>YDLAP%MYMS, NASM0=>YDLAP%NASM0)

ZEPS=100.0_JPRB*TINY(1.0_JPRB)

!$OMP PARALLEL DO SCHEDULE(STATIC,1)  PRIVATE(JNML,JLEV,INM,JM,ISP,ISP0,JN,ZSUM)
DO JNML=1,NUMP
  INM=MYMS(JNML)
  DO JLEV=1,KLEV

!*       1.    COMPUTE for a given n (=KNM).
!              -----------------------------

    IF (LDNWAVE) THEN
      ZSUM=0.0_JPRD
      DO JM=INM,1,-1
        ISP=NASM0(JM)+(INM-JM)*2
        ZSUM = ZSUM &
         & +2.0_JPRD*(REAL(PX(JLEV,ISP),JPRD)**2+REAL(PX(JLEV,ISP+1),JPRD)**2)
      ENDDO
      ISP0=NASM0(0)+INM*2
      ZSUM =(ZSUM+REAL(PX(JLEV,ISP0),JPRD)**2)*REAL(PMET(INM),JPRD)
      PSN(JLEV,JNML)=MAX(ZEPS,ZSUM)
    ENDIF

!     ------------------------------------------------------------------

!*       2.    COMPUTE for a given m (=KNM).
!              -----------------------------

    ZSUM = 0._JPRD
    IF(INM == 0)THEN
      DO JN=0,NSMAX
        ISP=NASM0(0)+JN*2
        ZSUM=ZSUM+REAL(PMET(JN),JPRD)*REAL(PX(JLEV,ISP),JPRD)**2
      ENDDO
      PSM(JLEV,JNML)=MAX(ZEPS,ZSUM)
    ELSE
      DO JN=INM,NSMAX
        ISP=NASM0(INM)+(JN-INM)*2
        ZSUM=ZSUM+2.0_JPRD*REAL(PMET(JN),JPRD)*&
         & (REAL(PX(JLEV,ISP),JPRD)**2+REAL(PX(JLEV,ISP+1),JPRD)**2)
      ENDDO
      PSM(JLEV,JNML)=MAX(ZEPS,ZSUM)
    ENDIF

  ENDDO
ENDDO
!$OMP END PARALLEL DO

!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SPNORMB',1,ZHOOK_HANDLE)
END SUBROUTINE SPNORMB
