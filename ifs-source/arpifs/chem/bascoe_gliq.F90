! (C) Copyright 2009- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE BASCOE_GLIQ(PNDFSA,PRADIUS,PTEMP,PAIR,PH2O,PHCL,PGHCL,PGH2O)


!**   DESCRIPTION 
!-----------------------------------------------------------------------
! calculate "sticking" coefficients for
!  (Hanson & Ravishankara J. Phys. Chem. 1994, 98, 5728-5735)
!            ClONO2 +H2O -> HOCl + HNO3
!            ClONO2 +HCl -> Cl2 + HNO3
! output values ghcl and gh2o are passed to unimolecular reaction rate
! calculations by common surface
! ----------------------------------------------------------------------
!     AUTHOR.
!     -------
!        Quentin Errera     *BIRA*
!
!     MODIFICATIONS.
!     --------------
!        Vincent Huijnen           : 2014-03-03
!        implementation of BASCOE scheme            
!
!
USE PARKIND1  ,    ONLY : JPIM,   JPRB
USE YOMHOOK   ,    ONLY : LHOOK,  DR_HOOK, JPHOOK
USE BASCOE_MODULE, ONLY :  NBINS
! USE YOMLUN       , ONLY : NULOUT ,NULERR
USE YOMCST       , ONLY : RPI 
IMPLICIT NONE

!-----------------------------------------------------------------------
!*       0.1  ARGUMENTS
!             ---------
REAL   (KIND=JPRB),INTENT(IN)      :: PNDFSA(NBINS),PRADIUS(NBINS)
REAL   (KIND=JPRB),INTENT(IN)      :: PTEMP,PAIR,PH2O,PHCL
REAL   (KIND=JPRB),INTENT(OUT)     :: PGH2O,PGHCL


!-----------------------------------------------------------------------
!  Local variables 
!-----------------------------------------------------------------------

REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE
REAL(KIND=JPRB)    :: ZNT
INTEGER(KIND=JPIM) :: JL
REAL(KIND=JPRB), PARAMETER :: ZALPHA=0.3
REAL(KIND=JPRB), PARAMETER :: ZRO=2.0e3
REAL(KIND=JPRB), PARAMETER :: ZKSUR=576.
REAL(KIND=JPRB)            :: ZAH2O,ZG0,ZHSTAR,ZGS,ZP,ZGCALC,ZADIVL,ZF,ZGE,ZGHCLI,ZGH2OI
REAL(KIND=JPRB)            :: ZCONV,ZSURFACE, ZTEMP
!VH      INTEGER(KIND=JPIM)  :: PAR_TYPE(NBINS)
REAL   (KIND=JPRB)         ::ZPPHCL,ZPPH2O

REAL(KIND=JPRB), PARAMETER :: ZEXPLIMIT = LOG(HUGE(ZEXPLIMIT))

!-------------------------------------------------------------------
#ifdef WITH_COMPO_DR_HOOK
IF (LHOOK) CALL DR_HOOK('BASCOE_GLIQ',0,ZHOOK_HANDLE )
#endif

      ZNT=0.0_JPRB
      ZTEMP=MAX(185._JPRB,PTEMP)
      DO JL=1,NBINS
        ZNT=ZNT+PNDFSA(JL)
      ENDDO
      PGH2O=0.0_JPRB
      PGHCL=0.0_JPRB
      IF( ZNT < 1.E-1_JPRB) THEN
#ifdef WITH_COMPO_DR_HOOK
        IF (LHOOK) CALL DR_HOOK('BASCOE_GLIQ',1,ZHOOK_HANDLE )
#endif
        RETURN
      ENDIF

      ZCONV=1.E-2_JPRB*PAIR*28.96_JPRB/(ZTEMP*8.3E3_JPRB)

      ZPPH2O=PH2O*PAIR/101300._JPRB
      ZPPHCL=PHCL*PAIR/101300._JPRB
      ZAH2O=1013.25_JPRB*ZPPH2O/(10**(9.217_JPRB-2190._JPRB/(ZTEMP-12.7_JPRB)))
      ZG0=1.18E-4_JPRB+(9.1E-3_JPRB)*ZAH2O+0.5_JPRB*ZAH2O*ZAH2O
      ZHSTAR=EXP(6250._JPRB/ZTEMP-10.414_JPRB)*ZAH2O**3.49_JPRB
      ZGS=MAX(1.E-20_JPRB, ZAH2O*ZKSUR*ZHSTAR*ZPPHCL)
      ZP=ZRO*ZHSTAR*ZPPHCL/ZAH2O
!	write(NULOUT,*)">>>>>> ",ZP,ZAH2O,ZPPHCL,HCL
      ZGCALC=ZG0*SQRT(1._JPRB+ZP)

!VH - obsolete?!      IF(ZNT > 0.0_JPRB) THEN
      DO JL=1,NBINS
!     IF(PAR_TYPE(I) .EQ. 0 ) THEN
        ZSURFACE=4*RPI*(PRADIUS(JL))**2
        ZADIVL=100.0_JPRB*PRADIUS(JL)/(1.4E-6_JPRB*SQRT(1._JPRB/ZAH2O))
        IF(ZADIVL < ZEXPLIMIT) THEN
          ZF=(EXP(ZADIVL)+EXP(-ZADIVL))/(EXP(ZADIVL)-EXP(-ZADIVL))-1._JPRB/ZADIVL
        ELSE
          ZF=1._JPRB
        ENDIF
        ZGE=1._JPRB/(1._JPRB/(ZGS+ZF*ZGCALC)+1._JPRB/ZALPHA)
        ZGHCLI=ZGE*(ZGS+ZF*ZGCALC*ZP/(1.+ZP))/(ZGS+ZF*ZGCALC)
        ZGH2OI=ZGE-ZGHCLI
        PGHCL=PGHCL+ZGHCLI*(PNDFSA(JL))*ZSURFACE*ZCONV
        PGH2O=PGH2O+ZGH2OI*(PNDFSA(JL))*ZSURFACE*ZCONV
!     ENDIF
      ENDDO
!VH       ENDIF

      PGHCL=MAX(1.E-20_JPRB,PGHCL)
      PGH2O=MAX(1.E-20_JPRB,PGH2O)

#ifdef WITH_COMPO_DR_HOOK
IF (LHOOK) CALL DR_HOOK('BASCOE_GLIQ',1,ZHOOK_HANDLE )
#endif
END SUBROUTINE BASCOE_GLIQ
