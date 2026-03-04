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

SUBROUTINE GPGRXYB(KPROMA,KD,KF,KFLEV,LDCOEF,YDVAB,YDCVER,PREL,PREM,PXYB,PXYBDER)

!**** *GPGRXYB* - Complement to routine "GPXYB".
!                 Computation of the horizontal gradient of quantities
!                 "alpha" and "delta" at model levels.

!     Purpose.
!     --------

!     "alpha" and "delta" are computed at model levels in routine "GPXYB",
!     but not their horizontal gradient. So this routine provides the
!     horizontal gradients at full levels. Quantity
!     "(grad(alpha)) + (grad(prehyd)/prehyd)"
!     is also provided separately (for case LVERTFE=.F.
!     its "NDLNPR=0" expression is simpler than the expressions
!     of "grad(alpha)" and "grad(prehyd)/prehyd").
!     Discretisation depends on variables "NDLNPR" and "LVERTFE".

!     For LVERTFE=.F., NDLNPR=0, discretisations are:

!      (grad(delta))[l] =
!      - (A[lbar]*B[lbar-1]-A[lbar-1]*B[lbar])/(prehyd[lbar]*prehyd[lbar-1])
!      * (grad prehyds)

!      (grad(alpha))[l] + (grad(prehyd)/prehyd)[l] =
!      B[lbar]/prehyd[lbar] * (grad prehyds)

!      Quantity "(grad(alpha))[l]" is computed by substracting
!      "(grad(prehyd)/prehyd)[l]" from
!      "(grad(alpha))[l] + (grad(prehyd)/prehyd)[l]"

!     For LVERTFE=.F., NDLNPR=1 or 2, discretisations are:

!      (grad(delta))[l] =
!      - delta[l] * (A[lbar]*B[lbar-1]-A[lbar-1]*B[lbar])
!      * (1/sqrt(prehyd[lbar]*prehyd[lbar-1])) * (1/(delta prehyd[l]))
!      * (grad prehyds)

!      (grad(alpha))[l] =
!      - alpha[l] * (A[lbar]*B[lbar-1]-A[lbar-1]*B[lbar])
!      * (1/sqrt(prehyd[lbar]*prehyd[lbar-1])) * (1/(delta prehyd[l]))
!      * (grad prehyds)

!      (grad(prehyd)/prehyd)[l] = prtgr[l] * (grad prehyds)
!      where "prtgr[l]" is computed in routine "gpxyb" as:
!      prtgr[l] = { (delta B)[l]
!      + delta[l] * (A[lbar]*B[lbar-1]-A[lbar-1]*B[lbar])/(delta prehyd[l]) }
!      * { 1/(delta prehyd[l]) }

!      In this case "(grad(alpha))[l]" is computed prior to
!      "(grad(alpha))[l] + (grad(prehyd)/prehyd)[l]"

!     For LVERTFE=.T., NDLNPR=0, discretisations are:

!      (grad(delta))[l] =
!      delta[l] * ((Delta B)[l]/(Delta prehyd)[l] - B[l]/prehyd[l])
!      * (grad prehyds)

!      grad(alpha) is useless in this case.

!     Notations:
!      - "grad" is the horizontal gradient operator.
!        (grad X = vnabla X = M vnabla' X)
!      - "prehyd" is the hydrostatic pressure.
!      - "prehyds" is the surface hydrostatic pressure.

!**   Interface.
!     ----------
!        *CALL* *GPGRXYB(...)

!        Explicit arguments :
!        --------------------
!         * INPUT:
!           KPROMA       : horizontal dimension
!           KD           : start of work
!           KF           : working length
!           KFLEV        : number of levels
!           LDCOEF       : if T, stores ZCOEFD, ZCOEFA, ZCOEFAPL in PXYBDER.
!           YDVAB        : contains information about hybrid vertical coordinate
!           PREL         : zonal component of "grad prehyds"
!           PREM         : meridian component of "grad prehyds"
!           PXYB         : contains pressure depth, "delta", "alpha".

!         * OUTPUT:
!           PXYBDER      : contains grad(delta), grad(alpha), grad(alpha + log prehyd)

!        Implicit arguments :   None.
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.    None.
!     ----------

!     Reference.
!     ----------

!     Author.
!     -------
!        K. YESSAD
!        Original : 00-08-11

!     Modifications.
!     --------------
!        K. Yessad (Dec 2008): remove dummy CDLOCK
!        K. Yessad (Jan 2011): introduce INTDYN_MOD structures.
!        K. Yessad (Dec 2011): use YDVAB.
!        K. Yessad (June 2017): introduce NDLNPR=2 (for NHQE model).
!     ------------------------------------------------------------------

USE PARKIND1  , ONLY : JPIM, JPRB
USE YOMHOOK   , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCVER   , ONLY : TCVER
USE YOMVERT   , ONLY : TVAB
USE INTDYN_MOD, ONLY : YYTXYB, YYTXYBDER

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KD 
INTEGER(KIND=JPIM),INTENT(IN)    :: KF 
LOGICAL           ,INTENT(IN)    :: LDCOEF
TYPE(TVAB)        ,INTENT(IN)    :: YDVAB
TYPE(TCVER)       ,INTENT(IN)    :: YDCVER
REAL(KIND=JPRB)   ,INTENT(IN)    :: PREL(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PREM(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PXYB(KPROMA,KFLEV,YYTXYB%NDIM) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PXYBDER(KPROMA,KFLEV,YYTXYBDER%NDIM) 

!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: ZCOEFA(KPROMA), ZCOEFAPL(KPROMA), ZCOEFD(KPROMA)

INTEGER(KIND=JPIM) :: JLEV, JROF
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GPGRXYB',0,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

!*    1/ Calculation of "grad delta" at full levels.

IF(YDCVER%LVERTFE) THEN
  DO JLEV=1,KFLEV
    DO JROF=KD,KF
      ZCOEFD(JROF)=(YDVAB%VDELB(JLEV)*PXYB(JROF,JLEV,YYTXYB%M_RDELP)-PXYB(JROF,JLEV,YYTXYB%M_RTGR))&
       & *PXYB(JROF,JLEV,YYTXYB%M_LNPR)  
      PXYBDER(JROF,JLEV,YYTXYBDER%M_LNPRL)=ZCOEFD(JROF)*PREL(JROF)
      PXYBDER(JROF,JLEV,YYTXYBDER%M_LNPRM)=ZCOEFD(JROF)*PREM(JROF)
    ENDDO
    IF (LDCOEF) PXYBDER(KD:KF,JLEV,YYTXYBDER%M_COEFD)=ZCOEFD(KD:KF)
  ENDDO
ELSE
  DO JLEV=1,KFLEV
    DO JROF=KD,KF
      ZCOEFD(JROF)=-YDVAB%VC(JLEV)*PXYB(JROF,JLEV,YYTXYB%M_RPP)
      PXYBDER(JROF,JLEV,YYTXYBDER%M_LNPRL)=ZCOEFD(JROF)*PREL(JROF)
      PXYBDER(JROF,JLEV,YYTXYBDER%M_LNPRM)=ZCOEFD(JROF)*PREM(JROF)
    ENDDO
    IF (LDCOEF) PXYBDER(KD:KF,JLEV,YYTXYBDER%M_COEFD)=ZCOEFD(KD:KF)
  ENDDO
ENDIF

!*    2/ Calculation of "grad (alpha + log prehyd)" at full levels
!        discretised as "grad alpha + (grad prehyd) / prehyd ",
!        and calculation of "grad alpha" at full levels.

IF(.NOT.YDCVER%LVERTFE) THEN
  IF(YDCVER%NDLNPR == 0) THEN
    DO JLEV=1,KFLEV
      DO JROF=KD,KF
        ZCOEFAPL(JROF)=YDVAB%VBH(JLEV)*PXYB(JROF,JLEV,YYTXYB%M_RPRE)
        ZCOEFA(JROF)=ZCOEFAPL(JROF)-PXYB(JROF,JLEV,YYTXYB%M_RTGR)
        PXYBDER(JROF,JLEV,YYTXYBDER%M_ALPHPLL)=ZCOEFAPL(JROF)*PREL(JROF)
        PXYBDER(JROF,JLEV,YYTXYBDER%M_ALPHPLM)=ZCOEFAPL(JROF)*PREM(JROF)
        PXYBDER(JROF,JLEV,YYTXYBDER%M_ALPHL)=ZCOEFA(JROF)*PREL(JROF)
        PXYBDER(JROF,JLEV,YYTXYBDER%M_ALPHM)=ZCOEFA(JROF)*PREM(JROF)
      ENDDO
      IF (LDCOEF) PXYBDER(KD:KF,JLEV,YYTXYBDER%M_COEFA)=ZCOEFA(KD:KF)
      IF (LDCOEF) PXYBDER(KD:KF,JLEV,YYTXYBDER%M_COEFAPL)=ZCOEFAPL(KD:KF)
    ENDDO
  ELSEIF(YDCVER%NDLNPR == 1 .OR. YDCVER%NDLNPR == 2) THEN
    DO JLEV=1,KFLEV
      DO JROF=KD,KF
        ZCOEFA(JROF)=-YDVAB%VC(JLEV)*PXYB(JROF,JLEV,YYTXYB%M_RPP)*PXYB(JROF,JLEV,YYTXYB%M_ALPH) &
         & /PXYB(JROF,JLEV,YYTXYB%M_LNPR)
        ZCOEFAPL(JROF)=ZCOEFA(JROF)+PXYB(JROF,JLEV,YYTXYB%M_RTGR)
        PXYBDER(JROF,JLEV,YYTXYBDER%M_ALPHPLL)=ZCOEFAPL(JROF)*PREL(JROF)
        PXYBDER(JROF,JLEV,YYTXYBDER%M_ALPHPLM)=ZCOEFAPL(JROF)*PREM(JROF)
        PXYBDER(JROF,JLEV,YYTXYBDER%M_ALPHL)=ZCOEFA(JROF)*PREL(JROF)
        PXYBDER(JROF,JLEV,YYTXYBDER%M_ALPHM)=ZCOEFA(JROF)*PREM(JROF)
      ENDDO
      IF (LDCOEF) PXYBDER(KD:KF,JLEV,YYTXYBDER%M_COEFA)=ZCOEFA(KD:KF)
      IF (LDCOEF) PXYBDER(KD:KF,JLEV,YYTXYBDER%M_COEFAPL)=ZCOEFAPL(KD:KF)
    ENDDO
  ENDIF
ENDIF

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GPGRXYB',1,ZHOOK_HANDLE)
END SUBROUTINE GPGRXYB
