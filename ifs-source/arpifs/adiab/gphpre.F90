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

SUBROUTINE GPHPRE(KPROMA,KFLEV,KSTART,KPROF,YDVAB,YDCVER,PRESH,PXYB,PRESF)

!**** *GPHPRE* - Computes half and full level pressure
!                Modern version of former GPPRE.
!                Modern version of former GPPREH+GPXYB+GPPREF

!     Purpose.
!     --------
!           Computes pressures at half and full model levels.

!**   Interface.
!     ----------
!        *CALL* *GPHPRE(...)

!        Explicit arguments :
!        --------------------

!          KPROMA    : horizontal dimensioning                                (in)
!          KFLEV     : vertical dimensioning                                  (in)
!          KSTART    : start of work                                          (in)
!          KPROF     : depth of work                                          (in)
!          YDVAB     : contains information about hybrid vertical coordinate  (in)
!          PRESH     : half level pressure                                    (inout)
!          PXYB      : contains pressure depth, "delta", "alpha"              (opt out)
!          PRESF     : full level pressure                                    (opt out)

!        Implicit arguments :  NONE.
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.  None.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      K. YESSAD (Sep 2011) after GPPRE, GPPREH, GPXYB and GPPREF.

!     Modifications.
!     --------------
!   K. Yessad (Dec 2016): Prune obsolete options.
!   K. Yessad (Mar 2017): Introduce NDLNPR=2 for NHQE model.
!     ------------------------------------------------------------------

USE PARKIND1  , ONLY : JPIM, JPRB
USE YOMHOOK   , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCVER   , ONLY : TCVER
USE YOMCST    , ONLY : RD, RCVD
USE YOMVERT   , ONLY : TVAB, TOPPRES
USE INTDYN_MOD, ONLY : YYTXYB

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM)         ,INTENT(IN)    :: KPROMA 
INTEGER(KIND=JPIM)         ,INTENT(IN)    :: KFLEV 
INTEGER(KIND=JPIM)         ,INTENT(IN)    :: KSTART 
INTEGER(KIND=JPIM)         ,INTENT(IN)    :: KPROF 
TYPE(TVAB)                 ,INTENT(IN)    :: YDVAB
TYPE(TCVER)                ,INTENT(IN)    :: YDCVER
REAL(KIND=JPRB)            ,INTENT(INOUT) :: PRESH(KPROMA,0:KFLEV)
REAL(KIND=JPRB),OPTIONAL   ,INTENT(OUT)   :: PXYB(KPROMA,KFLEV,YYTXYB%NDIM)
REAL(KIND=JPRB),OPTIONAL   ,INTENT(OUT)   :: PRESF(KPROMA,KFLEV)

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: IFIRST, JLEV, JLON, JJ, JTEMP, JM
REAL(KIND=JPRB) :: ZPRESF,ZPRESFD
REAL(KIND=JPRB) :: ZPRESTOP(KPROMA)
REAL(KIND=JPRB) :: ZRPRES(KPROMA,2)
REAL(KIND=JPRB) :: ZXYB(KPROMA,KFLEV,YYTXYB%NDIM)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('GPHPRE',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

!     ------------------------------------------------------------------

!*       1.    PRELIMINARY CALCULATIONS
!              ------------------------

IF (PRESENT(PXYB).OR.PRESENT(PRESF)) THEN

  ! This is introduced to get rid of the implicit
  ! assumption that the top level input for pressure is 0 hPa.
  ! The first block if is for economy (no do loop start up) and the second for safety.
  DO JLON=KSTART,KPROF
    ZPRESTOP(JLON)=YDVAB%VAH(0)+YDVAB%VBH(0)*PRESH(JLON,KFLEV)
  ENDDO
  IF(ZPRESTOP(KSTART) <= TOPPRES)THEN
    IFIRST=2
  ELSE
    IFIRST=1
    DO JLON=KSTART,KPROF
      IF(ZPRESTOP(JLON) <= TOPPRES)THEN
        IFIRST=2
        EXIT
      ENDIF
    ENDDO
  ENDIF

ENDIF

!     ------------------------------------------------------------------

!*       2.    COMPUTES HALF LEVEL PRESSURES
!              -----------------------------

! * compute upper-air PRESH from surface PRESH.
DO JLEV=0,KFLEV-1
  DO JLON=KSTART,KPROF
    PRESH(JLON,JLEV)=YDVAB%VAH(JLEV)+YDVAB%VBH(JLEV)*PRESH(JLON,KFLEV)
  ENDDO
ENDDO

!     ------------------------------------------------------------------

!*       3.    COMPUTES PXYB
!              -------------

IF (PRESENT(PXYB).OR.PRESENT(PRESF)) THEN

  ZXYB(:,:,:)=0.0_JPRB

  IF(YDCVER%LVERTFE) THEN

    DO JLEV=1,KFLEV
      DO JLON=KSTART,KPROF
        ZXYB(JLON,JLEV,YYTXYB%M_DELP)=YDVAB%VDELA(JLEV)+YDVAB%VDELB(JLEV)*PRESH(JLON,KFLEV)
        ZXYB(JLON,JLEV,YYTXYB%M_RDELP)=1.0_JPRB/ZXYB(JLON,JLEV,YYTXYB%M_DELP)
        ZPRESF=YDVAB%VAF(JLEV)+YDVAB%VBF(JLEV)*PRESH(JLON,KFLEV)
        ZPRESFD=1.0_JPRB/ZPRESF
        ZXYB(JLON,JLEV,YYTXYB%M_LNPR)=ZXYB(JLON,JLEV,YYTXYB%M_DELP)*ZPRESFD
        ZXYB(JLON,JLEV,YYTXYB%M_RTGR)=YDVAB%VBF(JLEV)*ZPRESFD
        ZXYB(JLON,JLEV,YYTXYB%M_ALPH)=(PRESH(JLON,JLEV)-ZPRESF)*ZPRESFD
      ENDDO
    ENDDO

  ELSE

    IF(YDCVER%NDLNPR == 0) THEN

      JJ=1
      JM=2
      DO JLON=KSTART,KPROF
        ZRPRES(JLON,JM)=1.0_JPRB/PRESH(JLON,IFIRST-1)
      ENDDO
      DO JLEV=IFIRST,KFLEV
        DO JLON=KSTART,KPROF
          ZRPRES(JLON,JJ)=1.0_JPRB/PRESH(JLON,JLEV)
          ZXYB(JLON,JLEV,YYTXYB%M_DELP)=PRESH(JLON,JLEV)-PRESH(JLON,JLEV-1)
          ZXYB(JLON,JLEV,YYTXYB%M_RDELP)=1.0_JPRB/ZXYB(JLON,JLEV,YYTXYB%M_DELP)
          ZXYB(JLON,JLEV,YYTXYB%M_LNPR)=LOG(PRESH(JLON,JLEV)*ZRPRES(JLON,JM))
          ZXYB(JLON,JLEV,YYTXYB%M_RPRE)=ZRPRES(JLON,JJ)
          ZXYB(JLON,JLEV,YYTXYB%M_ALPH)=1.0_JPRB-PRESH(JLON,JLEV-1)*ZXYB(JLON,JLEV,YYTXYB%M_RDELP)&
           & *ZXYB(JLON,JLEV,YYTXYB%M_LNPR)
          ZXYB(JLON,JLEV,YYTXYB%M_RPP)=ZRPRES(JLON,JJ)*ZRPRES(JLON,JM)
          ZXYB(JLON,JLEV,YYTXYB%M_RTGR)=ZXYB(JLON,JLEV,YYTXYB%M_RDELP)&
           & *(YDVAB%VDELB(JLEV)+YDVAB%VC(JLEV)*ZXYB(JLON,JLEV,YYTXYB%M_LNPR)*ZXYB(JLON,JLEV,YYTXYB%M_RDELP))
        ENDDO
        JTEMP=JM
        JM=JJ
        JJ=JTEMP
      ENDDO
      DO JLEV=1,IFIRST-1
        DO JLON=KSTART,KPROF
          ZXYB(JLON,JLEV,YYTXYB%M_DELP)=PRESH(JLON,JLEV)-PRESH(JLON,JLEV-1)
          ZXYB(JLON,JLEV,YYTXYB%M_RDELP)=1.0_JPRB/ZXYB(JLON,JLEV,YYTXYB%M_DELP)
          ZXYB(JLON,JLEV,YYTXYB%M_LNPR)=LOG(PRESH(JLON,1)/TOPPRES)
          ZXYB(JLON,JLEV,YYTXYB%M_RPRE)=1.0_JPRB/PRESH(JLON,1)
          ZXYB(JLON,JLEV,YYTXYB%M_ALPH)=YDCVER%RHYDR0
          ZXYB(JLON,JLEV,YYTXYB%M_RPP)=1.0_JPRB/(PRESH(JLON,1)*TOPPRES)
          ZXYB(JLON,JLEV,YYTXYB%M_RTGR)=ZXYB(JLON,JLEV,YYTXYB%M_RDELP)*YDVAB%VDELB(JLEV)
        ENDDO
      ENDDO

    ELSEIF(YDCVER%NDLNPR == 1 .OR. YDCVER%NDLNPR == 2) THEN

      DO JLEV=IFIRST,KFLEV
        DO JLON=KSTART,KPROF
          ZXYB(JLON,JLEV,YYTXYB%M_DELP)=PRESH(JLON,JLEV)-PRESH(JLON,JLEV-1)
          ZXYB(JLON,JLEV,YYTXYB%M_RDELP)=1.0_JPRB/ZXYB(JLON,JLEV,YYTXYB%M_DELP)
          ZXYB(JLON,JLEV,YYTXYB%M_RPP)=1.0_JPRB/(PRESH(JLON,JLEV)*PRESH(JLON,JLEV-1))
          ZXYB(JLON,JLEV,YYTXYB%M_LNPR)=ZXYB(JLON,JLEV,YYTXYB%M_DELP)*SQRT(ZXYB(JLON,JLEV,YYTXYB%M_RPP))
          ZXYB(JLON,JLEV,YYTXYB%M_ALPH)=1.0_JPRB-PRESH(JLON,JLEV-1)*ZXYB(JLON,JLEV,YYTXYB%M_RDELP)&
           & *ZXYB(JLON,JLEV,YYTXYB%M_LNPR)
          ZXYB(JLON,JLEV,YYTXYB%M_RTGR)=ZXYB(JLON,JLEV,YYTXYB%M_RDELP)&
           & *(YDVAB%VDELB(JLEV)+YDVAB%VC(JLEV)*ZXYB(JLON,JLEV,YYTXYB%M_LNPR)*ZXYB(JLON,JLEV,YYTXYB%M_RDELP))
          ZXYB(JLON,JLEV,YYTXYB%M_RPRE)=1.0_JPRB/PRESH(JLON,JLEV)
        ENDDO
      ENDDO

      IF(YDCVER%NDLNPR == 1) THEN
        DO JLEV=1,IFIRST-1
          DO JLON=KSTART,KPROF
            ZXYB(JLON,JLEV,YYTXYB%M_DELP)=PRESH(JLON,JLEV)
            ZXYB(JLON,JLEV,YYTXYB%M_LNPR)=2.0_JPRB+RCVD/RD
          ENDDO
        ENDDO
      ELSE
        DO JLEV=1,IFIRST-1
          DO JLON=KSTART,KPROF
            ZXYB(JLON,JLEV,YYTXYB%M_DELP)=PRESH(JLON,JLEV)
            ZXYB(JLON,JLEV,YYTXYB%M_LNPR)=1.0_JPRB+ZXYB(JLON,JLEV+1,YYTXYB%M_LNPR) &
             & *(ZXYB(JLON,JLEV,YYTXYB%M_DELP)/ZXYB(JLON,JLEV+1,YYTXYB%M_DELP))*SQRT(PRESH(JLON,JLEV+1)/PRESH(JLON,JLEV))
          ENDDO
        ENDDO
      ENDIF

      DO JLEV=1,IFIRST-1
        DO JLON=KSTART,KPROF
          ZXYB(JLON,JLEV,YYTXYB%M_RDELP)=1.0_JPRB/ZXYB(JLON,JLEV,YYTXYB%M_DELP)
          ZXYB(JLON,JLEV,YYTXYB%M_ALPH)=1.0_JPRB
          ZXYB(JLON,JLEV,YYTXYB%M_RTGR)=ZXYB(JLON,JLEV,YYTXYB%M_RDELP)*YDVAB%VDELB(JLEV)
          ZXYB(JLON,JLEV,YYTXYB%M_RPRE)=1.0_JPRB/PRESH(JLON,1)
          ZXYB(JLON,JLEV,YYTXYB%M_RPP)=(ZXYB(JLON,JLEV,YYTXYB%M_LNPR)*ZXYB(JLON,JLEV,YYTXYB%M_RDELP))**2
        ENDDO
      ENDDO

    ENDIF ! NDLNPR

  ENDIF ! LVERTFE

  IF (PRESENT(PXYB)) PXYB(:,:,:)=ZXYB(:,:,:)

ENDIF

!     ------------------------------------------------------------------

!*       4.    COMPUTES FULL LEVEL PRESSURES
!              -----------------------------

IF (PRESENT(PRESF)) THEN

  IF (YDCVER%LVERTFE) THEN
    DO JLEV=1,KFLEV
      PRESF(KSTART:KPROF,JLEV)=YDVAB%VAF(JLEV)+YDVAB%VBF(JLEV)*PRESH(KSTART:KPROF,KFLEV)
    ENDDO
  ELSE
    IF (YDCVER%NDLNPR == 0) THEN
      IF (YDCVER%LAPRXPK) THEN
        DO JLEV=1,KFLEV
          DO JLON=KSTART,KPROF
            PRESF(JLON,JLEV)=(PRESH(JLON,JLEV-1)+PRESH(JLON,JLEV))*0.5_JPRB
          ENDDO
        ENDDO
      ELSE
        DO JLEV=1,KFLEV
          DO JLON=KSTART,KPROF
            PRESF(JLON,JLEV)=EXP(-ZXYB(JLON,JLEV,YYTXYB%M_ALPH))*PRESH(JLON,JLEV)
          ENDDO
        ENDDO
      ENDIF
    ELSEIF (YDCVER%NDLNPR == 1 .OR. YDCVER%NDLNPR == 2) THEN
      DO JLEV=IFIRST,KFLEV
        DO JLON=KSTART,KPROF
          PRESF(JLON,JLEV)=(1.0_JPRB-ZXYB(JLON,JLEV,YYTXYB%M_ALPH))*PRESH(JLON,JLEV)
        ENDDO
      ENDDO
      DO JLEV=1,IFIRST-1
        DO JLON=KSTART,KPROF
          PRESF(JLON,JLEV)=PRESH(JLON,JLEV)/ZXYB(JLON,JLEV,YYTXYB%M_LNPR)
        ENDDO
      ENDDO
    ENDIF ! NDLNPR
  ENDIF ! LVERTFE

ENDIF

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GPHPRE',1,ZHOOK_HANDLE)
END SUBROUTINE GPHPRE
