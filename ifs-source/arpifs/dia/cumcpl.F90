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

SUBROUTINE CUMCPL(YDDIM,YDRIP,YDPHY,KLEVCO,KSTA,KEND,KSTGLO,KTS,&
 & PSTRTU,PSTRTV,PFRTH,PFRSO,PFCLL,PFCLN,PFCS,PFEVL,PFEVN,&
 & PFPLCL,PFPLCN,PFPLSL,PFPLSN,&
 & PRUISS,PRUISP,PTSUR,PDNSHF,PALB,&
 & PUCLS, PVCLS)

!**** *CUMCPL* - CUMUL OF COUPLED FIELDS.

!     Purpose.
!     --------
!           DIAGNOSTICS OF PHYSICAL SURFACE FLUXES.

!**   Interface.
!     ----------
!        *CALL* *CUMCPL*

!        Explicit arguments :
!        --------------------

!       KLEVCO                 - dimension for PFRTH,PFRSO.
!       KSTA                   - first element of work.
!       KEND                   - last element of work.
!       KSTGLO                 - global offset.
!       KTS                    - NUMBER OF SURFACE TEMPERATURES (INPUT)
!       FLUXES COMING FROM THE PHYSICAL PARAMETERIZATIONS:

!       PSTRTU(NPROMM)              - FLUX TURBULENT DE QTE DE MVT "U" (INPUT)
!       PSTRTV(NPROMM)              - FLUX TURBULENT DE QTE DE MVT "V" (INPUT)
!       PFRTH (NPROMM,0:KLEVCO,KTS+1) - FLUX DE RAYONNEMENT THERMIQUE    (INPUT)
!       PFRSO (NPROMM,0:KLEVCO,KTS+1) - FLUX DE RAYONNEMENT SOLAIRE.     (INPUT)
!       PFCLL (NPROMM,KTS+1) - FLUX DE CHALEUR LATENTE EAU.     (INPUT)
!       PFCLN (NPROMM,KTS+1) - FLUX DE CHALEUR LATENTE NEIGE    (INPUT)
!       PFCS  (NPROMM,KTS+1) - FLUX DE CHALEUR SENSIBLE         (INPUT)
!       PFEVL (NPROMM,KTS+1) - FLUX D'EVAPORATION EAU           (INPUT)
!       PFEVN (NPROMM,KTS+1) - FLUX D'EVAPORATION NEIGE         (INPUT)
!       PFPLCL(NPROMM)              - PRECIP. CONVECTIVES LIQUIDES.    (INPUT)
!       PFPLCN(NPROMM)              - PRECIP. CONVECTIVES NEIGEUSES.   (INPUT)
!       PFPLSL(NPROMM)              - PRECIP. STRATIFORMES LIQUIDES.   (INPUT)
!       PFPLSN(NPROMM)              - PRECIP. STRATIFORMES NEIGEUSES.  (INPUT)
!       PRUISS(NPROMM)              - RUISSELLEMENT SURFACE            (INPUT)
!       PRUISP(NPROMM)              - RUISSELLEMENT PROFONDEUR         (INPUT)
!       PTSUR(NPROMM)               - TEMPERATURE SURFACE              (INPUT)
!       PDNSHF(NPROMM)              - DERIVATIVE OF NONSOLAR FLUX      (INPUT)
!       PALB(NPROMM)                - SURFACE ALBEDO                   (INPUT)
!       PUCLS(NPROMM)               - U-COMPONENT OF WIND AT 10 METERS (INPUT)
!       PVCLS(NPROMM)               - V-COMPONENT OF WIND AT 10 METERS (INPUT)
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        CNRM/GMGEC/EAC/JPh Piedelievre
!        Original : 00-11-29 from CUMCODM

!     Modifications.
!     --------------
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!     K. Yessad (Dec 2008): remove dummy CDLOCK + cleanings
!     A.Alias    13-Oct-2009 Add fields for IPCC AR5 (E.Maisonnave)
!     K. Yessad (July 2014): Move some variables.
!     ------------------------------------------------------------------

USE YOMDIM   , ONLY : TDIM
USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMCST   , ONLY : RTT      
USE YOMRIP   , ONLY : TRIP
USE YOMCPL   , ONLY : FCSTSU   ,FCSTSV    ,FCCHAS    ,&
 &                    FCRSOS   ,FCHUMS   ,FCRUIS    ,FCCHSS    ,&
 &                    FCHUML   ,FCHUMN   ,FCTSUR2   ,FCTSTS    ,FCTSFL  ,&
 &                    FCWIMO   ,FCSUBL
USE YOMPHY   , ONLY : TPHY

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(TDIM)        , INTENT(IN) :: YDDIM
TYPE(TPHY)        , INTENT(IN) :: YDPHY
TYPE(TRIP)        , INTENT(IN) :: YDRIP
INTEGER(KIND=JPIM), INTENT(IN) :: KLEVCO 
INTEGER(KIND=JPIM), INTENT(IN) :: KTS 
INTEGER(KIND=JPIM), INTENT(IN) :: KSTA 
INTEGER(KIND=JPIM), INTENT(IN) :: KEND 
INTEGER(KIND=JPIM), INTENT(IN) :: KSTGLO 
REAL(KIND=JPRB)   , INTENT(IN) :: PSTRTU(YDDIM%NPROMM) 
REAL(KIND=JPRB)   , INTENT(IN) :: PSTRTV(YDDIM%NPROMM) 
REAL(KIND=JPRB)   , INTENT(IN) :: PFRTH(YDDIM%NPROMM,0:KLEVCO,KTS+1) 
REAL(KIND=JPRB)   , INTENT(IN) :: PFRSO(YDDIM%NPROMM,0:KLEVCO,KTS+1) 
REAL(KIND=JPRB)   , INTENT(IN) :: PFCLL(YDDIM%NPROMM,KTS+1) 
REAL(KIND=JPRB)   , INTENT(IN) :: PFCLN(YDDIM%NPROMM,KTS+1) 
REAL(KIND=JPRB)   , INTENT(IN) :: PFCS(YDDIM%NPROMM,KTS+1) 
REAL(KIND=JPRB)   , INTENT(IN) :: PFEVL(YDDIM%NPROMM,KTS+1) 
REAL(KIND=JPRB)   , INTENT(IN) :: PFEVN(YDDIM%NPROMM,KTS+1) 
REAL(KIND=JPRB)   , INTENT(IN) :: PFPLCL(YDDIM%NPROMM) 
REAL(KIND=JPRB)   , INTENT(IN) :: PFPLCN(YDDIM%NPROMM) 
REAL(KIND=JPRB)   , INTENT(IN) :: PFPLSL(YDDIM%NPROMM) 
REAL(KIND=JPRB)   , INTENT(IN) :: PFPLSN(YDDIM%NPROMM) 
REAL(KIND=JPRB)   , INTENT(IN) :: PRUISS(YDDIM%NPROMM) 
REAL(KIND=JPRB)   , INTENT(IN) :: PRUISP(YDDIM%NPROMM) 
REAL(KIND=JPRB)   , INTENT(IN) :: PTSUR(YDDIM%NPROMM) 
REAL(KIND=JPRB)   , INTENT(IN) :: PDNSHF(YDDIM%NPROMM) 
REAL(KIND=JPRB)   , INTENT(IN) :: PALB(YDDIM%NPROMM) 
REAL(KIND=JPRB)   , INTENT(IN) :: PUCLS(YDDIM%NPROMM)
REAL(KIND=JPRB)   , INTENT(IN) :: PVCLS(YDDIM%NPROMM)

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: IBEG, IEND, ISTA, JROF, JTS

LOGICAL :: LL_LMCCIC
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('CUMCPL',0,ZHOOK_HANDLE)
ASSOCIATE(TSTEP=>YDRIP%TSTEP, &
 & NPROMM=>YDDIM%NPROMM, &
 & LSOLV=>YDPHY%LSOLV)
!     ------------------------------------------------------------------

LL_LMCCIC=.TRUE.

!     CUMUL

ISTA=KSTA
IEND=KEND
IBEG=KSTGLO-1
DO JROF=ISTA,IEND
  FCSTSU(IBEG+JROF)=FCSTSU(IBEG+JROF)+PSTRTU(JROF)*TSTEP
  FCSTSV(IBEG+JROF)=FCSTSV(IBEG+JROF)+PSTRTV(JROF)*TSTEP
  IF(LSOLV) THEN
    FCRUIS(IBEG+JROF)=FCRUIS(IBEG+JROF)+PRUISP(JROF)*TSTEP
  ELSE
    FCRUIS(IBEG+JROF)=FCRUIS(IBEG+JROF)+(PRUISS(JROF)+PRUISP(JROF))*TSTEP
  ENDIF
  FCHUML(IBEG+JROF)=FCHUML(IBEG+JROF)+(PFPLCL(JROF)+PFPLSL(JROF))*TSTEP
  FCHUMN(IBEG+JROF)=FCHUMN(IBEG+JROF)+(PFPLCN(JROF)+PFPLSN(JROF))*TSTEP
  FCTSUR2(IBEG+JROF)=FCTSUR2(IBEG+JROF)+(PTSUR(JROF)-RTT)*TSTEP
  IF(LL_LMCCIC) THEN
! we store here the albedo and derivative of non solar heat flux
    FCTSTS(IBEG+JROF)=FCTSTS(IBEG+JROF)+PALB(JROF)*TSTEP
    FCTSFL(IBEG+JROF)=FCTSFL(IBEG+JROF)+PDNSHF(JROF)*TSTEP
  ELSE
! we store the temperature variance and covariance with heat flux
    FCTSTS(IBEG+JROF)=FCTSTS(IBEG+JROF)+(PTSUR(JROF)-RTT)&
     & *(PTSUR(JROF)-RTT)*TSTEP  
    FCTSFL(IBEG+JROF)=FCTSFL(IBEG+JROF)+(PTSUR(JROF)-RTT)&
     & *(PFRTH(JROF,KLEVCO,1)+PFCLL(JROF,1)+PFCLN(JROF,1)&
     & +PFCS(JROF,1))*TSTEP  
  ENDIF
  FCWIMO(IBEG+JROF)=FCWIMO(IBEG+JROF)+TSTEP*SQRT(PUCLS(JROF)*PUCLS(JROF)+&
                                             &  PVCLS(JROF)*PVCLS(JROF))
  FCSUBL(IBEG+JROF)=FCSUBL(IBEG+JROF)+TSTEP*PFEVN(JROF,1)
ENDDO

DO JTS=1,KTS+1
  DO JROF=ISTA,IEND
    FCCHAS(IBEG+JROF,JTS)=FCCHAS(IBEG+JROF,JTS)&
     & +(PFRTH(JROF,KLEVCO,JTS)+PFCLL(JROF,JTS)&
     & + PFCLN(JROF,JTS)+PFCS(JROF,JTS))*TSTEP  
    FCRSOS(IBEG+JROF,JTS)=FCRSOS(IBEG+JROF,JTS)+PFRSO(JROF,KLEVCO,JTS)*TSTEP
    FCHUMS(IBEG+JROF,JTS)=FCHUMS(IBEG+JROF,JTS)&
     & +(PFPLCL(JROF)+PFPLCN(JROF)+PFPLSL(JROF)+PFPLSN(JROF)&
     & + PFEVL(JROF,JTS)+PFEVN(JROF,JTS))*TSTEP  
    IF(LL_LMCCIC) THEN
! CA2 modif Laurent Terray=====> we store here the evaporation
      FCCHSS(IBEG+JROF,JTS)=FCCHSS(IBEG+JROF,JTS)&
       & +(PFEVL(JROF,JTS)+PFEVN(JROF,JTS))*TSTEP  
    ELSE
! we store the sensible heat flux
      FCCHSS(IBEG+JROF,JTS)=FCCHSS(IBEG+JROF,JTS)+PFCS(JROF,JTS)*TSTEP
    ENDIF
  ENDDO
ENDDO

!     ------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CUMCPL',1,ZHOOK_HANDLE)
END SUBROUTINE CUMCPL
