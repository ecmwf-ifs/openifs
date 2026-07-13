MODULE SRFT_MOD
CONTAINS
SUBROUTINE SRFT(KIDIA  , KFDIA , KLON   , KTILES, KLEVS ,&
 & PTMST  , PTSAM1M, PWSAM1M  , KSOTY, &
 & PFRTI  , PAHFSTI, PEVAPTI ,&
 & PSLRFLTI , PSSRFLTI, PGSN   ,&
 & PCTSA  , LDLAND ,&
 & YDCST  , YDSOIL , YDFLAKE, YDURB,&
 & PTSA   , PTSDFL , PDHTTS)  
  
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOS_CST  , ONLY : TCST
USE YOS_SOIL , ONLY : TSOIL
USE YOS_FLAKE, ONLY : TFLAKE
USE YOS_URB  , ONLY : TURB

USE SRFWDIF_MOD

! (C) Copyright 1993- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

!**** *SRFT* - Computes temperature changes in soil.

!     PURPOSE.
!     --------
!**   Computes temperature changes in soil due to 
!**   surface heat flux and diffusion.  

!**   INTERFACE.
!     ----------
!          *SRFT* IS CALLED FROM *SURF*.

!     PARAMETER   DESCRIPTION                                    UNITS
!     ---------   -----------                                    -----

!     INPUT PARAMETERS (INTEGER):
!    *KIDIA*      START POINT
!    *KFDIA*      END POINT
!    *KLON*       NUMBER OF GRID POINTS PER PACKET
!    *KLEVS*      NUMBER OF SOIL LAYERS
!    *KTILES*     NUMBER OF SURFACE TILES
!    *KDHVTTS*    Number of variables for soil energy budget
!    *KDHFTTS*    Number of fluxes for soil energy budget
!    *KSOTY*      SOIL TYPE                                   (1-7)

!     INPUT PARAMETERS (REAL):
!    *PTMST*      TIME STEP                                      S

!     INPUT PARAMETERS (LOGICAL):
!    *LDLAND*     LAND/SEA MASK (TRUE/FALSE)

!     INPUT PARAMETERS AT T-1 OR CONSTANT IN TIME (REAL):
!    *PTSAM1M*    SOIL TEMPERATURE                               K
!    *PWSAM1M*    SOIL MOISTURE                                m**3/m**3
!    *PSLRFLTI*   Tiled NET LONGWAVE  RADIATION AT THE SURFACE        W/M**2
!    *PGSN*       GROUND HEAT FLUX FROM SNOW DECK TO SOIL       W/M2
!    *PCTSA*      VOLUMETRIC HEAT CAPACITY                      J/K/M**3
!    *PFRTI*      TILE FRACTIONS                              (0-1)
!            1 : WATER                  5 : SNOW ON LOW-VEG+BARE-SOIL
!            2 : ICE                    6 : DRY SNOW-FREE HIGH-VEG
!            3 : WET SKIN               7 : SNOW UNDER HIGH-VEG
!            4 : DRY SNOW-FREE LOW-VEG  8 : BARE SOIL
!            9 : LAKE                  10 : URBAN
!    *PAHFSTI*    TILE SURFACE SENSIBLE HEAT FLUX                 W/M2
!    *PEVAPTI*    TILE SURFACE MOISTURE FLUX                     KG/M2/S
!    *PSSRFLTI*   TILE NET SHORTWAVE RADIATION FLUX AT SURFACE    W/M2

!     UPDATED PARAMETERS AT T+1 (UNFILTERED,REAL):
!    *PTSA*       SOIL TEMPERATURE                               K
!    *PTSDFL*     UPWARD FLUX BETWEEN SURFACE AND DEEP LAYER     W/M**2
!    *PDHTTS*     Diagnostic array for soil T (see module yomcdh)

!     METHOD.
!     -------

!          Parameters are set and the tridiagonal solver is called.

!     EXTERNALS.
!     ----------
!     *SRFWDIF*

!     REFERENCE.
!     ----------
!          See documentation.

!     Original :
!     P.VITERBO      E.C.M.W.F.     16/03/93
!     A.Beljaars     E.C.M.W.F.     12/03/1999
!                    (implicit solution of diffusion equation)
!     P.VITERBO      E.C.M.W.F.     17/05/2000
!                    (Surface DDH for TILES)
!     J.F. Estrade *ECMWF* 03-10-01 move in surf vob
!     P. Viterbo    24-05-2004      Change surface units
!     G. Balsamo    08-01-2006      Include Van Genuchten Hydro.
!     E. Dutra      21/11/2008      change ground heat flux ! account for new lake tile
!     E. Dutra      10/10/2014  net longwave tiled
!     ------------------------------------------------------------------

IMPLICIT NONE

! Declaration of arguments

INTEGER(KIND=JPIM), INTENT(IN)   :: KIDIA
INTEGER(KIND=JPIM), INTENT(IN)   :: KFDIA
INTEGER(KIND=JPIM), INTENT(IN)   :: KLON
INTEGER(KIND=JPIM), INTENT(IN)   :: KTILES
INTEGER(KIND=JPIM), INTENT(IN)   :: KLEVS
INTEGER(KIND=JPIM), INTENT(IN)   :: KSOTY(:)

REAL(KIND=JPRB),    INTENT(IN)   :: PTMST
REAL(KIND=JPRB),    INTENT(IN)   :: PTSAM1M(:,:)
REAL(KIND=JPRB),    INTENT(IN)   :: PWSAM1M(:,:)
REAL(KIND=JPRB),    INTENT(IN)   :: PFRTI(:,:)
REAL(KIND=JPRB),    INTENT(IN)   :: PAHFSTI(:,:)
REAL(KIND=JPRB),    INTENT(IN)   :: PEVAPTI(:,:)
REAL(KIND=JPRB),    INTENT(IN)   :: PSLRFLTI(:,:)
REAL(KIND=JPRB),    INTENT(IN)   :: PSSRFLTI(:,:)
REAL(KIND=JPRB),    INTENT(IN)   :: PGSN(:)
REAL(KIND=JPRB),    INTENT(IN)   :: PCTSA(:,:)
LOGICAL,            INTENT(IN)   :: LDLAND(:)
TYPE(TCST),         INTENT(IN)   :: YDCST
TYPE(TSOIL),        INTENT(IN)   :: YDSOIL
TYPE(TFLAKE),       INTENT(IN)   :: YDFLAKE
TYPE(TURB),         INTENT(IN)   :: YDURB

REAL(KIND=JPRB),    INTENT(OUT)  :: PTSA(:,:)
REAL(KIND=JPRB),    INTENT(OUT)  :: PTSDFL(:)
REAL(KIND=JPRB),    INTENT(OUT)  :: PDHTTS(:,:,:)

!      LOCAL VARIABLES

REAL(KIND=JPRB) :: ZSURFL(KLON)
REAL(KIND=JPRB) :: ZRHS(KLON,KLEVS), ZCDZ(KLON,KLEVS),&
 & ZLST(KLON,KLEVS), ZDIF(KLON,KLEVS),&
 & ZTSA(KLON,KLEVS)  
LOGICAL :: LLDOSOIL(KLON)

INTEGER(KIND=JPIM) :: JK, JL, JS

REAL(KIND=JPRB) :: ZCONS1, ZCONS2, ZSLRFL, ZSSRFL, ZTHFL, ZTMST,&
 & ZFF, ZWU, ZLSM, ZLIC, ZLWT, ZLAMBDASAT, ZKERSTEN, ZINVWSAT
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "fcsurf.h"

!*         1. SET UP SOME CONSTANTS.
!             --- -- ---- ----------
!*    COMPUTATIONAL CONSTANTS.
!     ------------- ----------

IF (LHOOK) CALL DR_HOOK('SRFT_MOD:SRFT',0,ZHOOK_HANDLE)
ASSOCIATE(RLVTT=>YDCST%RLVTT, &
 & LEFLAKE=>YDFLAKE%LEFLAKE, &
 & LEVGEN=>YDSOIL%LEVGEN, RDAT=>YDSOIL%RDAT, RFRSMALL=>YDSOIL%RFRSMALL, &
 & RTF1=>YDSOIL%RTF1, RTF2=>YDSOIL%RTF2, RTF3=>YDSOIL%RTF3, RTF4=>YDSOIL%RTF4, &
 & RKERST1=>YDSOIL%RKERST1, RKERST2=>YDSOIL%RKERST2, RKERST3=>YDSOIL%RKERST3, &
 & RLAMBDADRY=>YDSOIL%RLAMBDADRY, RLAMBDADRYM=>YDSOIL%RLAMBDADRYM, &
 & RLAMBDAICE=>YDSOIL%RLAMBDAICE, RLAMBDAWAT=>YDSOIL%RLAMBDAWAT, &
 & RLAMSAT1=>YDSOIL%RLAMSAT1, RLAMSAT1M=>YDSOIL%RLAMSAT1M, &
 & RRCSOIL=>YDSOIL%RRCSOIL, RSIMP=>YDSOIL%RSIMP, RWSAT=>YDSOIL%RWSAT, &
 & RWSATM=>YDSOIL%RWSATM, RURBTC1=>YDURB%RURBTC1,RURBVHC=>YDURB%RURBVHC)
ZTMST=1.0_JPRB/PTMST
ZCONS1=PTMST*RSIMP
ZCONS2=1.0_JPRB-1.0_JPRB/RSIMP

!*         2. Compute net heat flux at the surface.
!             -------------------------------------

DO JL=KIDIA,KFDIA
  LLDOSOIL(JL)=LDLAND(JL)
  IF (LLDOSOIL(JL)) THEN

!         In principle this should be fractional averaging,  
!         but since the land sea mask is used fractions 3,4,5,6
!         7 and 8 add up to 1 for land. PGSN(JL) is already 
!         smeared out over the entire grid square. In future 
!         (when fractional land is used), ZSSRFL+ZSLRFL+ZTHFL 
!         should be divided by the sum of fractions 3,4,6 and 8. 

    ZSSRFL=PFRTI(JL,3)*PSSRFLTI(JL,3)&
     & +PFRTI(JL,4)*PSSRFLTI(JL,4)&
     & +PFRTI(JL,6)*PSSRFLTI(JL,6)&
     & +PFRTI(JL,8)*PSSRFLTI(JL,8)
    ZSLRFL=PFRTI(JL,3)*PSLRFLTI(JL,3)&
     &+PFRTI(JL,4)*PSLRFLTI(JL,4)&
     &+PFRTI(JL,6)*PSLRFLTI(JL,6)&
     &+PFRTI(JL,8)*PSLRFLTI(JL,8)
    ZTHFL=PFRTI(JL,3)*(PAHFSTI(JL,3)+RLVTT*PEVAPTI(JL,3))&
     & +PFRTI(JL,4)*(PAHFSTI(JL,4)+RLVTT*PEVAPTI(JL,4))&
     & +PFRTI(JL,6)*(PAHFSTI(JL,6)+RLVTT*PEVAPTI(JL,6))&
     & +PFRTI(JL,8)*(PAHFSTI(JL,8)+RLVTT*PEVAPTI(JL,8))  
     
    IF ( KTILES .GT. 9 ) THEN
     ZSSRFL=ZSSRFL+PFRTI(JL,10)*PSSRFLTI(JL,10)
     ZSLRFL=ZSLRFL+PFRTI(JL,10)*PSLRFLTI(JL,10)
     ZTHFL=ZTHFL+PFRTI(JL,10)*(PAHFSTI(JL,10)+RLVTT*PEVAPTI(JL,10))
    ENDIF

    ZSURFL(JL)=ZSSRFL+ZSLRFL+ZTHFL+PGSN(JL)

    IF ( LEFLAKE ) THEN
      IF ( PFRTI(JL,9) > RFRSMALL ) THEN
        ZSURFL(JL)=PGSN(JL)
        IF ( (PFRTI(JL,3)+PFRTI(JL,4)+PFRTI(JL,6)+PFRTI(JL,8)) > RFRSMALL ) THEN
          ZSURFL(JL)=PGSN(JL)+(ZSSRFL+ZSLRFL+ZTHFL) & 
                  & / (PFRTI(JL,3)+PFRTI(JL,4)+PFRTI(JL,6)+PFRTI(JL,8))
        ENDIF
        IF ( KTILES .GT. 9 ) THEN
          IF ( (PFRTI(JL,3)+PFRTI(JL,4)+PFRTI(JL,6)+PFRTI(JL,8)+PFRTI(JL,10)) > RFRSMALL ) THEN
            ZSURFL(JL)=PGSN(JL)+(ZSSRFL+ZSLRFL+ZTHFL) &
                    & / (PFRTI(JL,3)+PFRTI(JL,4)+PFRTI(JL,6)+PFRTI(JL,8)+PFRTI(JL,10))
          ENDIF
        ENDIF
      ENDIF
    ENDIF
    
  ELSE
    ZSURFL(JL)=0.0_JPRB
  ENDIF
ENDDO

!*         3. Compute exchange coeff. layer by layer
!             --------------------------------------

DO JK=1,KLEVS
  DO JL=KIDIA,KFDIA
    IF (LLDOSOIL(JL)) THEN
! added fix to be consistent with Peters-Lidard et al. 1998 
      IF(PTSAM1M(JL,JK) < RTF1.AND.PTSAM1M(JL,JK) > RTF2) THEN
        ZFF=0.5_JPRB*(1.0_JPRB-SIN(RTF4*(PTSAM1M(JL,JK)-RTF3)))
      ELSEIF (PTSAM1M(JL,JK) <= RTF2) THEN
        ZFF=1.0_JPRB
      ELSE
        ZFF=0.0_JPRB
      ENDIF
      !ZFF=0.0_JPRB
     ! ZWU=PWSAM1M(JL,JK)*(1.0_JPRB-ZFF)
     ! ZLWT=RLAMBDAWAT**ZWU
      IF(LEVGEN)THEN
        JS=KSOTY(JL)
        ZLSM=RLAMSAT1M(JS)

        ZWU=RWSATM(JS)*(1.0_JPRB-ZFF)
        ZLWT=RLAMBDAWAT**ZWU

        ZINVWSAT=1.0_JPRB/RWSATM(JS)
        ZLIC=RLAMBDAICE**(RWSATM(JS)-ZWU)

        ZLAMBDASAT=ZLSM*ZLIC*ZLWT

        ZKERSTEN=RKERST2*LOG10(MAX(RKERST1,PWSAM1M(JL,JK)*ZINVWSAT))+RKERST3
        ZDIF(JL,JK)=RLAMBDADRYM(JS)+ZKERSTEN*(ZLAMBDASAT-RLAMBDADRYM(JS))
        IF (YDSOIL%LESOILCOND ) THEN
          ZDIF(JL,JK)=FSOILTCOND(PWSAM1M(JL,JK),ZFF,JS)
!         print*,'soil cond',jk,PWSAM1M(JL,JK),js,ZDIF(JL,JK),FSOILTCOND(PWSAM1M(JL,JK),JS),ZDIF(JL,JK)-FSOILTCOND(PWSAM1M(JL,JK),JS)
        ENDIF
        IF ( KTILES .GT. 9 ) THEN
          IF (JK == 1) THEN
            ZDIF(JL,JK) = (1.0_JPRB-PFRTI(JL,10))*ZDIF(JL,JK) + (PFRTI(JL,10)*RURBTC1) !EXCHANGE DIFFERS FOR URBAN
          ENDIF
        ENDIF

      ELSE
        ZLSM=RLAMSAT1

        ZWU=RWSAT*(1.0_JPRB-ZFF)
        ZLWT=RLAMBDAWAT**ZWU

        ZINVWSAT=1.0_JPRB/RWSAT
        ZLIC=RLAMBDAICE**(RWSAT-ZWU)
        ZLAMBDASAT=ZLSM*ZLIC*ZLWT
        ZKERSTEN=RKERST2*LOG10(MAX(RKERST1,PWSAM1M(JL,JK)*ZINVWSAT))+RKERST3
        ZDIF(JL,JK)=RLAMBDADRY+ZKERSTEN*(ZLAMBDASAT-RLAMBDADRY)
      ENDIF
    ELSE
      ZDIF(JL,JK)=0.0_JPRB
    ENDIF
  ENDDO
ENDDO

!*         4. Set arrays
!             ----------
!     Layer 1
JK=1
DO JL=KIDIA,KFDIA
  IF (LLDOSOIL(JL)) THEN
    ZLST(JL,JK)=ZCONS1*(ZDIF(JL,JK)+ZDIF(JL,JK+1))/(RDAT(JK)+RDAT(JK+1))
    ZCDZ(JL,JK)=PCTSA(JL,JK)*RDAT(JK)
    ZRHS(JL,JK)=PTMST*ZSURFL(JL)/ZCDZ(JL,JK)
  ENDIF
ENDDO

!     Layers 2 to KLEVS-1
DO JK=2,KLEVS-1
  DO JL=KIDIA,KFDIA
    IF (LLDOSOIL(JL)) THEN
      ZLST(JL,JK)=ZCONS1*(ZDIF(JL,JK)+ZDIF(JL,JK+1))/(RDAT(JK)+RDAT(JK+1))
      ZCDZ(JL,JK)=PCTSA(JL,JK)*RDAT(JK)
      ZRHS(JL,JK)=0.0_JPRB
    ENDIF
  ENDDO
ENDDO

!     Layers KLEVS
JK=KLEVS
DO JL=KIDIA,KFDIA
  IF (LLDOSOIL(JL)) THEN
    ZLST(JL,JK)=0.0_JPRB
    ZCDZ(JL,JK)=PCTSA(JL,JK)*RDAT(JK)
    ZRHS(JL,JK)=0.0_JPRB
  ENDIF
ENDDO

!*         5. Call tridiagonal solver
!             -----------------------
CALL SRFWDIF(KIDIA,KFDIA,KLON,KLEVS,PTSAM1M,ZLST,ZRHS,ZCDZ,YDSOIL,ZTSA,LLDOSOIL)

!*         6. Flux between layer 1 and 2
!             --------------------------
DO JL=KIDIA,KFDIA
  IF (LLDOSOIL(JL)) THEN
    PTSDFL(JL)=(ZTSA(JL,2)-ZTSA(JL,1))*ZLST(JL,1)*ZTMST
  ELSE
    PTSDFL(JL)=0.0_JPRB
  ENDIF
ENDDO

!*         7. New temperatures
!             ----------------
DO JK=1,KLEVS
  DO JL=KIDIA,KFDIA
    IF (LLDOSOIL(JL)) THEN
      PTSA(JL,JK)=PTSAM1M(JL,JK)*ZCONS2+ZTSA(JL,JK)
    ELSE
      PTSA(JL,JK)=PTSAM1M(JL,JK)
    ENDIF
  ENDDO
ENDDO

!*         7a. DDH diagnostics
!              ---------------
IF (SIZE(PDHTTS) > 0) THEN
  DO JL=KIDIA,KFDIA
    IF (LLDOSOIL(JL)) THEN
! Sensible heat flux
      PDHTTS(JL,1,9)=PFRTI(JL,3)*PAHFSTI(JL,3)+&
        & PFRTI(JL,4)*PAHFSTI(JL,4)+&
        & PFRTI(JL,6)*PAHFSTI(JL,6)+&
        & PFRTI(JL,8)*PAHFSTI(JL,8)  
! Latent heat flux
      PDHTTS(JL,1,10)=RLVTT*(PFRTI(JL,3)*PEVAPTI(JL,3)+&
        & PFRTI(JL,4)*PEVAPTI(JL,4)+&
        & PFRTI(JL,6)*PEVAPTI(JL,6)+&
        & PFRTI(JL,8)*PEVAPTI(JL,8))  
! Flux snow-soil
      PDHTTS(JL,1,13)=PGSN(JL)
! Ground heat flux
      PDHTTS(JL,1,14)=-ZSURFL(JL)
! Urban
      IF ( KTILES .GT. 9 ) THEN
        PDHTTS(JL,1,9)=PDHTTS(JL,1,9)+PFRTI(JL,10)*PAHFSTI(JL,10)
        PDHTTS(JL,1,10)=PDHTTS(JL,1,10)+PFRTI(JL,10)*PEVAPTI(JL,10)
      ENDIF    
    ELSE
      PDHTTS(JL,1,9:10)=0.0_JPRB
      PDHTTS(JL,1,13)=0.0_JPRB
    ENDIF
  ENDDO
  DO JK=2,KLEVS
    DO JL=KIDIA,KFDIA
      PDHTTS(JL,JK,9)=0.0_JPRB
      PDHTTS(JL,JK,10)=0.0_JPRB
      PDHTTS(JL,JK,13)=0.0_JPRB
      PDHTTS(JL,JK,14)=0.0_JPRB 
    ENDDO
  ENDDO

  DO JK=1,KLEVS
    DO JL=KIDIA,KFDIA
      IF (LLDOSOIL(JL)) THEN
! Heat capacity per unit surface
        PDHTTS(JL,JK,1)=RRCSOIL*RDAT(JK)
        IF ( KTILES .GT. 9 ) THEN
          PDHTTS(JL,JK,1)=(RRCSOIL*(1.0_JPRB-PFRTI(JL,10)) + (PFRTI(JL,10)*RURBVHC))*RDAT(JK)
        ENDIF
! Soil temperature
        PDHTTS(JL,JK,2)=PTSAM1M(JL,JK)
! Layer energy per unit surface
        PDHTTS(JL,JK,3)=RRCSOIL*RDAT(JK)*PTSAM1M(JL,JK)
! Layer depth
        PDHTTS(JL,JK,4)=RDAT(JK)
! Soil water phase changes
        PDHTTS(JL,JK,12)= -(PCTSA(JL,JK)-RRCSOIL)* &
          & RDAT(JK)*(PTSA(JL,JK)-PTSAM1M(JL,JK))*ZTMST 
        PDHTTS(JL,JK,15)=RRCSOIL*RDAT(JK)*(ZTSA(JL,JK)-PTSAM1M(JL,JK))/PTMST
      ELSE
        PDHTTS(JL,JK,1:4)=0.0_JPRB
        PDHTTS(JL,JK,12)=0.0_JPRB
! This is needed due to uninitialized values being copied around
        PDHTTS(JL,JK,14)=0.0_JPRB
        PDHTTS(JL,JK,15)=0.0_JPRB
      ENDIF
    ENDDO
  ENDDO

  DO JK=1,KLEVS-1
    DO JL=KIDIA,KFDIA
      IF (LLDOSOIL(JL)) THEN
! (Ground) heat flux "between soil layers"
        PDHTTS(JL,JK,11)=RSIMP*(ZTSA(JL,JK)-ZTSA(JL,JK+1))*ZLST(JL,JK)*&
          & ZTMST  
      ELSE
        PDHTTS(JL,JK,11)=0.0_JPRB
      ENDIF
    ENDDO
  ENDDO
  DO JL=KIDIA,KFDIA
    PDHTTS(JL,KLEVS,11)=0.0_JPRB
  ENDDO
ENDIF

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SRFT_MOD:SRFT',1,ZHOOK_HANDLE)

END SUBROUTINE SRFT
END MODULE SRFT_MOD
