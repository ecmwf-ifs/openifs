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

SUBROUTINE GPGRXYB(KPROMA,KD,KF,KFLEV,LDCOEF,YDVAB,YDCVER,PREL,PREM,PXYB,PXYBDER, &
                 & PDELP, PLNPR, PRDELP, PALPH, PRTGR, PRPRE, PRPP, &
                 & PCOEFD_DER, PLNPRL_DER, PLNPRM_DER, PCOEFA_DER, &
                 & PCOEFAPL_DER, PALPHPLL_DER, PALPHPLM_DER, PALPHL_DER, PALPHM_DER)

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
!        H. Petithomme (Dec 2020): use of pointers for optimisation
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
REAL(KIND=JPRB),TARGET,OPTIONAL,INTENT(IN) :: PXYB(KPROMA,KFLEV,YYTXYB%NDIM) 
REAL(KIND=JPRB),TARGET,OPTIONAL,INTENT(OUT):: PXYBDER(KPROMA,KFLEV,YYTXYBDER%NDIM) 
REAL(KIND=JPRB),TARGET,OPTIONAL,INTENT(IN) :: PDELP(KPROMA,KFLEV), PLNPR(KPROMA,KFLEV), PRDELP(KPROMA,KFLEV), PALPH(KPROMA,KFLEV)
REAL(KIND=JPRB),TARGET,OPTIONAL,INTENT(IN) :: PRTGR(KPROMA,KFLEV), PRPRE(KPROMA,KFLEV), PRPP(KPROMA,KFLEV)
REAL(KIND=JPRB),TARGET,OPTIONAL,INTENT(OUT) :: PCOEFD_DER(KPROMA,KFLEV),   PLNPRL_DER(KPROMA,KFLEV),   PLNPRM_DER(KPROMA,KFLEV) 
REAL(KIND=JPRB),TARGET,OPTIONAL,INTENT(OUT) :: PCOEFA_DER(KPROMA,KFLEV),   PCOEFAPL_DER(KPROMA,KFLEV), PALPHPLL_DER(KPROMA,KFLEV)
REAL(KIND=JPRB),TARGET,OPTIONAL,INTENT(OUT) :: PALPHPLM_DER(KPROMA,KFLEV), PALPHL_DER(KPROMA,KFLEV),   PALPHM_DER(KPROMA,KFLEV)

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JLEV, JROF
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
REAL(KIND=JPRB),TARGET :: ZCOEFA0(KPROMA),ZCOEFAPL0(KPROMA),ZCOEFD0(KPROMA)
REAL(KIND=JPRB),CONTIGUOUS,POINTER :: ZCOEFA(:),ZCOEFAPL(:),ZCOEFD(:)
REAL(KIND=JPRB),CONTIGUOUS,POINTER :: ZDELP (:,:), ZLNPR (:,:), ZRDELP (:,:), ZALPH (:,:), ZRTGR (:,:), ZRPRE (:,:), ZRPP (:,:)  
REAL(KIND=JPRB),CONTIGUOUS,POINTER :: ZCOEFD_DER(:,:), ZLNPRL_DER(:,:), ZLNPRM_DER(:,:), ZCOEFA_DER(:,:)
REAL(KIND=JPRB),CONTIGUOUS,POINTER :: ZCOEFAPL_DER(:,:), ZALPHPLL_DER(:,:), ZALPHPLM_DER(:,:), ZALPHL_DER(:,:), ZALPHM_DER(:,:)

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GPGRXYB',0,ZHOOK_HANDLE)

ZDELP  => NULL ()
ZLNPR  => NULL ()
ZRDELP => NULL ()
ZALPH  => NULL ()
ZRTGR  => NULL ()
ZRPRE  => NULL ()
ZRPP   => NULL ()

IF (PRESENT (PXYB)) THEN
  ZDELP  => PXYB (:,:,YYTXYB%M_DELP )
  ZLNPR  => PXYB (:,:,YYTXYB%M_LNPR )
  ZRDELP => PXYB (:,:,YYTXYB%M_RDELP)
  ZALPH  => PXYB (:,:,YYTXYB%M_ALPH )
  ZRTGR  => PXYB (:,:,YYTXYB%M_RTGR )
  ZRPRE  => PXYB (:,:,YYTXYB%M_RPRE )
  ZRPP   => PXYB (:,:,YYTXYB%M_RPP  )
ENDIF

IF (PRESENT (PDELP )) ZDELP  => PDELP 
IF (PRESENT (PLNPR )) ZLNPR  => PLNPR 
IF (PRESENT (PRDELP)) ZRDELP => PRDELP
IF (PRESENT (PALPH )) ZALPH  => PALPH 
IF (PRESENT (PRTGR )) ZRTGR  => PRTGR 
IF (PRESENT (PRPRE )) ZRPRE  => PRPRE 
IF (PRESENT (PRPP  )) ZRPP   => PRPP  

ZCOEFD_DER   => NULL ()
ZLNPRL_DER   => NULL ()
ZLNPRM_DER   => NULL ()
ZCOEFA_DER   => NULL ()
ZCOEFAPL_DER => NULL ()
ZALPHPLL_DER => NULL ()
ZALPHPLM_DER => NULL ()
ZALPHL_DER   => NULL ()
ZALPHM_DER   => NULL ()

IF (PRESENT (PXYBDER)) THEN
  ZCOEFD_DER   => PXYBDER(:,:,YYTXYBDER%M_COEFD)
  ZLNPRL_DER   => PXYBDER(:,:,YYTXYBDER%M_LNPRL)
  ZLNPRM_DER   => PXYBDER(:,:,YYTXYBDER%M_LNPRM)
  ZCOEFA_DER   => PXYBDER(:,:,YYTXYBDER%M_COEFA)
  ZCOEFAPL_DER => PXYBDER(:,:,YYTXYBDER%M_COEFAPL)
  ZALPHPLL_DER => PXYBDER(:,:,YYTXYBDER%M_ALPHPLL)
  ZALPHPLM_DER => PXYBDER(:,:,YYTXYBDER%M_ALPHPLM)
  ZALPHL_DER   => PXYBDER(:,:,YYTXYBDER%M_ALPHL)
  ZALPHM_DER   => PXYBDER(:,:,YYTXYBDER%M_ALPHM)
ENDIF

IF (PRESENT (PCOEFD_DER  ))  ZCOEFD_DER   => PCOEFD_DER   
IF (PRESENT (PLNPRL_DER  ))  ZLNPRL_DER   => PLNPRL_DER   
IF (PRESENT (PLNPRM_DER  ))  ZLNPRM_DER   => PLNPRM_DER   
IF (PRESENT (PCOEFA_DER  ))  ZCOEFA_DER   => PCOEFA_DER   
IF (PRESENT (PCOEFAPL_DER))  ZCOEFAPL_DER => PCOEFAPL_DER 
IF (PRESENT (PALPHPLL_DER))  ZALPHPLL_DER => PALPHPLL_DER 
IF (PRESENT (PALPHPLM_DER))  ZALPHPLM_DER => PALPHPLM_DER 
IF (PRESENT (PALPHL_DER  ))  ZALPHL_DER   => PALPHL_DER   
IF (PRESENT (PALPHM_DER  ))  ZALPHM_DER   => PALPHM_DER   

!     ------------------------------------------------------------------

!*    1/ Calculation of "grad delta" at full levels.

IF (.NOT.LDCOEF) ZCOEFD => ZCOEFD0(:)

! optim: compilers may report dependence between pxyb and ydvab, ignored with ivdep/nodep

IF(YDCVER%LVERTFE) THEN
    DO JLEV=1,KFLEV
    IF (LDCOEF) ZCOEFD => ZCOEFD_DER(:,JLEV)
    !DIR$ IVDEP
    !CDIR NODEP
    DO JROF=KD,KF
      ZCOEFD(JROF)=(YDVAB%VDELB(JLEV)*ZRDELP(JROF,JLEV)&
       & -ZRTGR(JROF,JLEV))*ZLNPR(JROF,JLEV)  
      ZLNPRL_DER(JROF,JLEV)=ZCOEFD(JROF)*PREL(JROF)
      ZLNPRM_DER(JROF,JLEV)=ZCOEFD(JROF)*PREM(JROF)
    ENDDO
  ENDDO
ELSE
  DO JLEV=1,KFLEV
    IF (LDCOEF) ZCOEFD => ZCOEFD_DER(:,JLEV)
    !DIR$ IVDEP
    !CDIR NODEP
    DO JROF=KD,KF
      ZCOEFD(JROF)=-YDVAB%VC(JLEV)*ZRPP(JROF,JLEV)
      ZLNPRL_DER(JROF,JLEV)=ZCOEFD(JROF)*PREL(JROF)
      ZLNPRM_DER(JROF,JLEV)=ZCOEFD(JROF)*PREM(JROF)
    ENDDO
  ENDDO

!*    2/ Calculation of "grad (alpha + log prehyd)" at full levels
!        discretised as "grad alpha + (grad prehyd) / prehyd ",
!        and calculation of "grad alpha" at full levels.

  IF (.NOT.LDCOEF) THEN
   ZCOEFA => ZCOEFA0(:)
   ZCOEFAPL => ZCOEFAPL0(:)
  ENDIF
  IF(YDCVER%NDLNPR == 0) THEN
    DO JLEV=1,KFLEV
      IF (LDCOEF) THEN
        ZCOEFA => ZCOEFA_DER(:,JLEV)
        ZCOEFAPL => ZCOEFAPL_DER(:,JLEV)
      END IF

      !DIR$ IVDEP
      !CDIR NODEP
      DO JROF=KD,KF
        ZCOEFAPL(JROF)=YDVAB%VBH(JLEV)*ZRPRE(JROF,JLEV)
        ZALPHPLL_DER(JROF,JLEV)=ZCOEFAPL(JROF)*PREL(JROF)
        ZALPHPLM_DER(JROF,JLEV)=ZCOEFAPL(JROF)*PREM(JROF)

        ZCOEFA(JROF)=ZCOEFAPL(JROF)-ZRTGR(JROF,JLEV)
        ZALPHL_DER(JROF,JLEV)=ZCOEFA(JROF)*PREL(JROF)
        ZALPHM_DER(JROF,JLEV)=ZCOEFA(JROF)*PREM(JROF)
      ENDDO
    ENDDO
  ELSEIF(YDCVER%NDLNPR == 1 .OR. YDCVER%NDLNPR == 2) THEN
    DO JLEV=1,KFLEV
      IF (LDCOEF) THEN
        ZCOEFA => ZCOEFA_DER(:,JLEV)
        ZCOEFAPL => ZCOEFAPL_DER(:,JLEV)
      END IF

      !DIR$ IVDEP
      !CDIR NODEP
      DO JROF=KD,KF
        ZCOEFA(JROF)=-YDVAB%VC(JLEV)*ZRPP(JROF,JLEV)&
          & *ZALPH(JROF,JLEV)/ZLNPR(JROF,JLEV)
        ZALPHL_DER(JROF,JLEV)=ZCOEFA(JROF)*PREL(JROF)
        ZALPHM_DER(JROF,JLEV)=ZCOEFA(JROF)*PREM(JROF)

        ZCOEFAPL(JROF)=ZCOEFA(JROF)+ZRTGR(JROF,JLEV)
        ZALPHPLL_DER(JROF,JLEV)=ZCOEFAPL(JROF)*PREL(JROF)
        ZALPHPLM_DER(JROF,JLEV)=ZCOEFAPL(JROF)*PREM(JROF)
      ENDDO
    ENDDO
  ENDIF
ENDIF

IF (LHOOK) CALL DR_HOOK('GPGRXYB',1,ZHOOK_HANDLE)
END SUBROUTINE GPGRXYB
