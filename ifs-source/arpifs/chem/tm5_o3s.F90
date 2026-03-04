! (C) Copyright 2009- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE TM5_O3S(YDCHEM,KCHEM,PDT,PRR,PRJ,PY)

!**   DESCRIPTION 
!     ----------
!
!--------------------------------------------------------------------------
!   Compute O3 loss in troposphere, following CB05-type chemistry
!   To be used for O3S tracer.
!--------------------------------------------------------------------------
!
!
!
!**   INTERFACE.
!     ----------
!          *TM5_O3S* IS CALLED FROM *CHEM_TM5 / CHEM_BASCOETM5*.

! INPUTS:
! -------
! KIDIA :  Start of Array  
! KFDIA :  End  of Array 
! KCHEM :  Length of tracer array 
! KLON  :  Length of Arrays 
! PDT   :  Time step in seconds 
! PRR  (KLON,NREAC)      : reaction rates
! PRJ  (KLON,NPHOTO)     : photolysis rates
!
!
! OUTPUTS:
! -------
! PY (KLON,KCHEM)        : final   volume ratios OF TRACERS           (mol/mol)
!
! LOCAL:
! -------
!
!
!     AUTHOR.
!     -------
!        VINCENT HUIJNEN    *KNMI*
!
!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 2022-03-14



USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCHEM  , ONLY : TCHEM
USE TM5_CHEM_MODULE, ONLY : NREAC, IO3S_CB05=>IO3S, IHO2_CB05=> IHO2,IOH_CB05=>IOH, IETH_CB05=>IETH, IOLE_CB05=> IOLE, &
   & IISOP_CB05=> IISOP, IC3H6_CB05=>IC3H6 ,ITERP_CB05=>ITERP, IISPD_CB05=> IISPD, &
   & ISO2_CB05=>ISO2, ITOL_CB05=>ITOL, IXYL_CB05=>IXYL, &
!* reaction rates
   &   KO3HO2,KO3OH,KC62,KC58,KC77,KO3C3H6,KO3TERP, KO3ISPD, &
   &    KSO2O3G,KO3TOL,KO3XYL

USE BASCOETM5_MODULE, ONLY :  IO3S_CBA=>IO3S, IHO2_CBA=> IHO2,IOH_CBA=>IOH, IETH_CBA=>IETH, IOLE_CBA=> IOLE, &
   & IISOP_CBA=> IISOP, IC3H6_CBA=>IC3H6 ,ITERP_CBA=>ITERP, IISPD_CBA=> IISPD, &
   & ISO2_CBA=>ISO2, ITOL_CBA=>ITOL, IXYL_CBA=>IXYL

USE TM5_PHOTOLYSIS , ONLY : NPHOTO,&
  !* photolysis rates
   & JO3D  

IMPLICIT NONE

!-----------------------------------------------------------------------
!*       0.1  ARGUMENTS
!             ---------

TYPE(TCHEM),    INTENT(IN)    :: YDCHEM
INTEGER(KIND=JPIM),INTENT(IN) :: KCHEM 
REAL(KIND=JPRB),INTENT(IN)    :: PDT
REAL(KIND=JPRB),INTENT(IN)    :: PRR(NREAC)
REAL(KIND=JPRB),INTENT(IN)    :: PRJ(NPHOTO)   
REAL(KIND=JPRB),INTENT(INOUT) :: PY(KCHEM)  ! final concentration; only O3S updated

! * LOCAL 
REAL(KIND=JPHOOK)  :: ZHOOK_HANDLE

REAL (KIND=JPRB)   ::  ZXL3
! * counters and indices
INTEGER(KIND=JPIM) :: IO3S,IHO2, IOH,IETH,IOLE, IISOP, IC3H6,ITERP, IISPD, ISO2, IXYL, ITOL

#include "abor1.intfb.h"
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('TM5_O3S',0,ZHOOK_HANDLE)


! Find correct indices to tracer variables, depending on scheme. For now only 
! support CB05 / CB05BASCOE chemistry


SELECT CASE (TRIM(YDCHEM%CHEM_SCHEME))
  CASE ("tm5")         
    IO3S=IO3S_CB05
    IHO2=IHO2_CB05
    IOH=IOH_CB05
    IETH=IETH_CB05
    IOLE=IOLE_CB05
    IISOP=IISOP_CB05
    IC3H6=IC3H6_CB05
    ITERP=ITERP_CB05
    IISPD=IISPD_CB05
    ISO2=ISO2_CB05
    ITOL=ITOL_CB05
    IXYL=IXYL_CB05
  CASE ("bascoetm5")         
    IO3S=IO3S_CBA
    IHO2=IHO2_CBA
    IOH=IOH_CBA
    IETH=IETH_CBA
    IOLE=IOLE_CBA
    IISOP=IISOP_CBA
    IC3H6=IC3H6_CBA
    ITERP=ITERP_CBA
    IISPD=IISPD_CBA
    ISO2=ISO2_CBA
    ITOL=ITOL_CBA
    IXYL=IXYL_CBA
 CASE DEFAULT
    CALL ABOR1("Chemistry scheme not supported")
END SELECT


   ! Special treatment for O3S: Only loss in troposphere.
IF (TRIM(YDCHEM%CHEM_SCHEME)=="tm5" .AND. YDCHEM%KCHEM_SOLVE==1_JPIM) THEN
   ! Old chemistry: No xyl/tol fields defined
   ZXL3= PRR(KO3HO2)*PY(IHO2)&
      & + PRR(KO3OH)*PY(IOH)&
      ! & + PRR(KO3PO3)& VH not safe for now.
      & + PRJ(JO3D)&
      & + PRR(KC62)*PY(IETH)&
      & + PRR(KC58)*PY(IOLE)&
      & + PRR(KC77)*PY(IISOP)&
      & + PRR(KO3C3H6)*PY(IC3H6)&
      & + PRR(KO3TERP)*PY(ITERP)&
      & + PRR(KO3ISPD)*PY(IISPD)&
      & + PRR(KSO2O3G)*PY(ISO2)
ELSE
   ! Standard (kpp-based) chemistry
   ZXL3= PRR(KO3HO2)*PY(IHO2)&
      & + PRR(KO3OH)*PY(IOH)&
      ! & + PRR(KO3PO3)& VH not safe for now.
      & + PRJ(JO3D)&
      & + PRR(KC62)*PY(IETH)&
      & + PRR(KC58)*PY(IOLE)&
      & + PRR(KC77)*PY(IISOP)&
      & + PRR(KO3C3H6)*PY(IC3H6)&
      & + PRR(KO3TERP)*PY(ITERP)&
      & + PRR(KO3ISPD)*PY(IISPD)&
      & + PRR(KSO2O3G)*PY(ISO2)&
      & + PRR(KO3TOL)*PY(ITOL)&
      & + PRR(KO3XYL)*PY(IXYL)
ENDIF
PY(IO3S)=PY(IO3S)/(1.+ZXL3*PDT)

IF (LHOOK) CALL DR_HOOK('TM5_O3S',1,ZHOOK_HANDLE)
END SUBROUTINE TM5_O3S
