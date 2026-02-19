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

SUBROUTINE CUMCOE(YDRIP,KPROMA,KSTA,KEND,KSTGLO,&
 & PSTRTU,PSTRTV,&
 & PFRTH ,PFRSO ,PFCLL ,PFCLN ,PFCS,&
 & PDIFTQ,&
 & PFPLCL,PFPLCN,PFPLSL,PFPLSN,&
 & PRUISS,PRUISP)  

!**** *CUMCOE* - CUMUL OF EC COUPLED FIELDS.

!     Purpose.
!     --------
!           DIAGNOSTICS OF PHYSICAL SURFACE FLUXES.

!**   Interface.
!     ----------
!        *CALL* *CUMCOE*

!        Explicit arguments :
!        --------------------

!       KPROMA         - HORIZONTAL DIMENSION                 (INPUT)
!       KSTA,KEND      - START AND END INDEX OF NPROMA CHUNK  (INPUT)
!       KSTGLO         - START ADDRESS OF CHUNK               (INPUT)
!       PSTRTU(KPROMA) - FLUX TURBULENT DE QTE DE MVT "U".    (INPUT)
!       PSTRTV(KPROMA) - FLUX TURBULENT DE QTE DE MVT "V".    (INPUT)
!       PFRTH (KPROMA,1) - FLUX DE RAYONNEMENT THERMIQUE.     (INPUT)
!       PFRSO (KPROMA,1) - FLUX DE RAYONNEMENT SOLAIRE.       (INPUT)
!       PFCLL (KPROMA,1) - FLUX DE CHALEUR LATENTE EAU.       (INPUT)
!       PFCLN (KPROMA,1) - FLUX DE CHALEUR LATENTE NEIGE      (INPUT)
!       PFCS  (KPROMA,1) - FLUX DE CHALEUR SENSIBLE           (INPUT)
!       PDIFTQ(KPROMA) - MOISTURE FLUX AT SURFACE             (INPUT)
!       PFPLCL(KPROMA) - PRECIP. CONVECTIVES LIQUIDES.        (INPUT)
!       PFPLCN(KPROMA) - PRECIP. CONVECTIVES NEIGEUSES.       (INPUT)
!       PFPLSL(KPROMA) - PRECIP. STRATIFORMES LIQUIDES.       (INPUT)
!       PFPLSN(KPROMA) - PRECIP. STRATIFORMES NEIGEUSES.      (INPUT)
!       PRUISS(KPROMA) - RUISSELLEMENT SURFACE                (INPUT)
!       PRUISP(KPROMA) - RUISSELLEMENT PROFONDEUR             (INPUT)

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
!      T. Stockdale.    Based on CUMCO.

!     Modifications.
!     --------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      K. Yessad (July 2014): Move some variables.
!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMRIP   , ONLY : TRIP
USE YOMGCO   , ONLY : STRSU, STRSV, FCHAS, FRSOS, &
 &                    FHUMS, RUIST, FCHLL, FCHLN, FCHSS, FHUML, FHUMN  

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(TRIP)        ,INTENT(IN)    :: YDRIP
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KEND 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTGLO 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSTRTU(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSTRTV(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFRTH(KPROMA,1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFRSO(KPROMA,1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFCLL(KPROMA,1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFCLN(KPROMA,1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFCS(KPROMA,1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDIFTQ(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFPLCL(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFPLCN(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFPLSL(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFPLSN(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRUISS(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRUISP(KPROMA) 

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: IBEG, IL, JROF
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CUMCOE',0,ZHOOK_HANDLE)
ASSOCIATE(TSTEP=>YDRIP%TSTEP)
!     ------------------------------------------------------------------

IBEG=KSTGLO-1
DO JROF=KSTA,KEND
  IL=IBEG+JROF
  STRSU(IL)  =STRSU(IL)  + PSTRTU(JROF)*TSTEP
  STRSV(IL)  =STRSV(IL)  + PSTRTV(JROF)*TSTEP
  RUIST(IL)  =RUIST(IL)  +(PRUISS(JROF)+PRUISP(JROF))*TSTEP
  FHUML(IL)  =FHUML(IL)  +(PFPLCL(JROF)+PFPLSL(JROF))*TSTEP
  FHUMN(IL)  =FHUMN(IL)  +(PFPLCN(JROF)+PFPLSN(JROF))*TSTEP
  FRSOS(IL,1)=FRSOS(IL,1)+ PFRSO(JROF,1)*TSTEP
  FCHLL(IL,1)=FCHLL(IL,1)+ PFCLL(JROF,1)*TSTEP
  FCHLN(IL,1)=FCHLN(IL,1)+ PFCLN(JROF,1)*TSTEP
  FCHSS(IL,1)=FCHSS(IL,1)+ PFCS (JROF,1)*TSTEP
  FCHAS(IL,1)=FCHAS(IL,1)&
   & +(PFRTH(JROF,1)+PFCLL(JROF,1)+PFCLN(JROF,1)+PFCS(JROF,1))*TSTEP
  FHUMS(IL,1)=FHUMS(IL,1)&
   & +(PFPLCL(JROF)+PFPLCN(JROF)+PFPLSL(JROF)+PFPLSN(JROF)+PDIFTQ(JROF))*TSTEP
ENDDO

!     ------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CUMCOE',1,ZHOOK_HANDLE)
END SUBROUTINE CUMCOE
