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

!----------------------------------------------------------------------
SUBROUTINE UPDSST ( YDMCC,YDPHY1,KPROMA, KSTART, KPROF, KFNUDG, KCSS, KBLK,&
 !----------------------------------------------------------------------
 ! - INPUT 1D
 & PXNUDG,&
 ! - IN/OUT 1D
 & PTS,  PTP,  PWS,  PWP,  PSNS,&
 & PALN, PDEN, PWSG, PWPG,&
 & PLSM, PALB, PEMI, PZ0F, PZ0H, PD2,  PIVEG)  
!****
!     ------------------------------------------------------------------

!     UPDATING SST AND MISCELLANEOUS IN CASE OF NUDGING

!     INDEED, UPDCLI DEFINES ONCE A DAY THE OPEN SEA/SEA ICE PROPERTIES
!     FROM THE CONTENT OF BCOND FILES. HERE WE USE AT EACH TIME STEP
!     THE CONTENT OF THE NUDGING FILES, WHICH MAY BE QUITE DIFFERENT

!     SORRY WILLIAM, THE OTHER COMMENTS USE MOLIERE LANGUAGE

!     ------------------------------------------------------------------

!     ARGUMENTS D ENTREE
!     ------------------
!       KPROMA : DIMENSION HORIZONTALE
!       KSTART : BORNE INITIALE HORIZONTALE DES CALCULS
!       KPROF : BORNE FINALE HORIZONTALE DES CALCULS
!       KFNUDG  : NBRE DE PAS DE TEMPS D'ANALYSE A INTERPOLER
!       KCSS    : NBRE DE NIVEAUX DANS LE SOL PROFOND
!       KBLK  : OFFSET (KBLK=(KSTGLO-1)/NPROMA+1)       

!     ENTREES
!     -------
!       PXNUDG(KFNUDG) : COEFFICIENTS TEMPORELS

!     ARGUMENTS IMPLICITES
!     --------------------
!       COEFF DE RAPPEL = COMMON /YOMNUD /

!     SORTIES
!     -------
!       PTS   : TEMPERATURE DE SURFACE
!       PTP   : TEMPERATURE PROFONDE
!       PWS   : RESERVOIR DE SURFACE (LIQUIDE)
!       PWP   : RESERVOIR PROFOND   (LIQUIDE)
!       PWSG  : RESERVOIR DE SURFACE (SOLIDE)
!       PWPG  : RESERVOIR PROFOND   (SOLIDE)
!       PSNS  : RESERVOIR DE NEIGE  
!       PLSM  : MASQUE TERRE-MER
!       PALB  : ALBEDO
!       PEMIS : EMISSIVITE
!       PZ0F  : LONGUEUR DE RUGOSITE
!       PZ0H  : LONGUEUR DE RUGOSITE THERMIQUE
!       PD2   : PROFONDEUR RACINAIRE
!       PIVEG : TYPE DE VEGETATION
!       PALN  : ALBEDO NEIGE
!       PDEN  : DENSITE NEIGE

!     IMPROVISE PAR
!     ------------- 
!       MICHEL DEQUE

!     AVEC LES MODIFS DE
!     ------------------
!       2001-07-04: Michel Deque: no update when coupled
!       2006-08-30: A. Alias    : RZHMER replaced by RZHZ0M*RZ0MER
!       2006-09-05: A. Alias    : RZHGLA replaced by RZHZ0G*RZ0GLA
!       2006-10-26: A. Alias    : Change in arguments list following
!                   introduction of new surface fields (CY31R2)
!       2007-10-31: A. Alias    : Test modified to proceed with nudging
!                                 Test removed (field PZ0H exists when LSOLV)
!       2009-08   : M. Deque    : set in agreement with external updcli
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMNUD   , ONLY : XNUDST, XVUST
USE YOMCST   , ONLY : RG
USE YOMPHY1  , ONLY : TPHY1
USE YOMMCC   , ONLY : TMCC

IMPLICIT NONE

TYPE(TMCC)        ,INTENT(INOUT) :: YDMCC
TYPE(TPHY1)       ,INTENT(INOUT) :: YDPHY1
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFNUDG 
INTEGER(KIND=JPIM),INTENT(IN)    :: KBLK
INTEGER(KIND=JPIM),INTENT(IN)    :: KCSS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTART 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PXNUDG(KFNUDG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTS(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTP(KPROMA,KCSS) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PWS(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PWP(KPROMA,KCSS) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PWSG(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PWPG(KPROMA,KCSS) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PALN(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDEN(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSNS(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PLSM(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PALB(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PEMI(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PZ0F(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PZ0H(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PD2(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PIVEG(KPROMA) 
INTEGER(KIND=JPIM), PARAMETER :: JPKZ = KIND(PIVEG)

REAL(KIND=JPRB) :: ZREF(KPROMA)

INTEGER(KIND=JPIM) :: JCSS, JROF, JSTEP

LOGICAL :: LLMEGL

REAL(KIND=JPRB) :: ZWXGLA, ZWXMER, ZWXSUR
REAL(KIND=JPRB) :: ZRANO , ZRDNO
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!*
!     ------------------------------------------------------------------
!     1.- MISE A JOUR SUR MER/BANQUISE

IF (LHOOK) CALL DR_HOOK('UPDSST',0,ZHOOK_HANDLE)
ASSOCIATE(NTVGLA=>YDPHY1%NTVGLA, NTVMER=>YDPHY1%NTVMER, RZ0GLA=>YDPHY1%RZ0GLA, &
 & RZHZ0M=>YDPHY1%RZHZ0M, EMMMER=>YDPHY1%EMMMER, ALBMER=>YDPHY1%ALBMER, &
 & RZHZ0G=>YDPHY1%RZHZ0G, RD2GLA=>YDPHY1%RD2GLA, RD1=>YDPHY1%RD1, &
 & RD2MER=>YDPHY1%RD2MER, TMERGL=>YDPHY1%TMERGL, ALBGLA=>YDPHY1%ALBGLA, &
 & RZ0MER=>YDPHY1%RZ0MER, EMMGLA=>YDPHY1%EMMGLA, &
 & LMCC02=>YDMCC%LMCC02, LMCC03=>YDMCC%LMCC03)
!       PAS DE MISE A JOUR DANS LE CAS DU NUDGING AVEC CORRECTION
!       PAS DE MISE A JOUR DANS LE CAS DU COUPLAGE
IF(XNUDST > 0.0_JPRB.AND.(.NOT.LMCC03)) THEN
  ZWXMER=1000._JPRB*RD2MER
  ZWXGLA=1000._JPRB*RD2GLA
  ZWXSUR=1000._JPRB*RD1
  ZRANO=0.70_JPRB
  ZRDNO=0.30_JPRB
  !      ZREF :  TEMPERATURE DE SURFACE DE RAPPEL
  DO JROF = KSTART,KPROF
    ZREF(JROF)=0.0_JPRB
  ENDDO
  DO JSTEP=1,KFNUDG
    DO JROF = KSTART,KPROF
      ZREF(JROF)=ZREF(JROF)+XVUST(JROF,JSTEP,KBLK)*PXNUDG(JSTEP)
    ENDDO
  ENDDO
  !      BOUCLE SUR LES POINTS DE MER OU BANQUISE      
  DO JROF = KSTART,KPROF
    LLMEGL=(NINT(PIVEG(JROF)) == NTVMER.OR.&
       & (NINT(10*PIVEG(JROF))) == (10*NTVGLA+1))  
    IF(LLMEGL)THEN
      IF(ZREF(JROF) <= TMERGL) THEN
!     - - - - -
!     Sea Ice 
!     - - - - -
        IF (NINT(PIVEG(JROF)) == NTVMER) THEN
!       New formed Sea Ice
          PTS(JROF)=ZREF(JROF)
          DO JCSS=1,KCSS
            PTP(JROF,JCSS)=TMERGL
          ENDDO
          PALN(JROF)=ZRANO
          PDEN(JROF)=ZRDNO
          IF(LMCC02)THEN
            PLSM(JROF)=1.0_JPRB
          ELSE
            PLSM(JROF)=0.0_JPRB
          ENDIF
          PZ0F(JROF)=RZ0GLA*RG
          PALB(JROF)=ALBGLA
          PEMI(JROF)=EMMGLA
          PWS(JROF)=0.0_JPRB
          DO JCSS=1,KCSS
            PWP(JROF,JCSS)=0.0_JPRB
            PWPG(JROF,JCSS)=ZWXGLA
          ENDDO
          PWSG(JROF)=ZWXSUR
          PSNS(JROF)=0.0_JPRB
          PIVEG(JROF)=REAL(NTVGLA,JPKZ)+.1_JPRB
          PD2(JROF)=RD2GLA
          PZ0H(JROF)=RZHZ0G*RZ0GLA*RG
        ELSE
!       Old Sea Ice
          IF(LMCC02)THEN
            PLSM(JROF)=1.0_JPRB
          ELSE
            PTS(JROF)=ZREF(JROF)
            PLSM(JROF)=0.0_JPRB
          ENDIF
          PZ0F(JROF)=RZ0GLA*RG
          PALB(JROF)=ALBGLA
          PEMI(JROF)=EMMGLA
          PWS(JROF)=0.0_JPRB
          DO JCSS=1,KCSS
            PWP(JROF,JCSS)=0.0_JPRB
            PWPG(JROF,JCSS)=ZWXGLA
          ENDDO
          PWSG(JROF)=ZWXSUR
          PIVEG(JROF)=REAL(NTVGLA,JPKZ)+.1_JPRB
          PD2(JROF)=RD2GLA
          PZ0H(JROF)=RZHZ0G*RZ0GLA*RG
        ENDIF
      ELSE
!     - - - - -
!     Open Sea
!     - - - - -
        PLSM(JROF)=0.0_JPRB
  ! don't reset Z0 because of Charnock formula
  !           PZ0F(JROF)=RZ0MER*RG
        PALB(JROF)=ALBMER
        PEMI(JROF)=EMMMER
        PTS(JROF)=ZREF(JROF)
        DO JCSS=1,KCSS
          PTP(JROF,JCSS)=ZREF(JROF)
        ENDDO
        PALN(JROF)=ZRANO
        PDEN(JROF)=ZRDNO
        PWS(JROF)=ZWXSUR
        DO JCSS=1,KCSS
          PWP(JROF,JCSS)=ZWXMER
          PWPG(JROF,JCSS)=0.0_JPRB
        ENDDO
        PWSG(JROF)=0.0_JPRB
        PSNS(JROF)=0.0_JPRB
        PIVEG(JROF)=REAL(NTVMER,JPKZ)
        PD2(JROF)=RD2MER
        PZ0H(JROF)=RZHZ0M*RZ0MER*RG
      ENDIF
    ENDIF
  ENDDO

ENDIF

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('UPDSST',1,ZHOOK_HANDLE)
END SUBROUTINE UPDSST
