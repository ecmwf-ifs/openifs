MODULE SRFTS_MOD
CONTAINS
SUBROUTINE SRFTS(KIDIA  , KFDIA  , KLON   , KLEVS ,&
 & PTMST  , PTSAM1M, PWSAM1M, KSOTY, &
 & PFRTI  , PAHFSTI, PEVAPTI,&
 & PSLRFL ,PSSRFLTI, PGSN   ,&
 & PCTSA  , PTSA   , LDLAND ,&
 & YDCST  , YDSOIL , YDFLAKE)  
  
USE PARKIND1  , ONLY : JPIM, JPRB
USE YOMHOOK   , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOS_CST   , ONLY : TCST
USE YOS_SOIL  , ONLY : TSOIL
USE YOS_FLAKE , ONLY : TFLAKE

USE SRFWDIFS_MOD

#ifdef DOC
! (C) Copyright 2011- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

!**** *SRFTS* - Computes temperature changes in soil.

!     PURPOSE.
!     --------
!**   Computes temperature changes in soil due to 
!**   surface heat flux and diffusion.  

!**   INTERFACE.
!     ----------
!          *SRFTS* IS CALLED FROM *SURFTSTPS*.

!     PARAMETER   DESCRIPTION                                    UNITS
!     ---------   -----------                                    -----

!     INPUT PARAMETERS (INTEGER):
!    *KIDIA*      START POINT
!    *KFDIA*      END POINT
!    *KLON*       NUMBER OF GRID POINTS PER PACKET
!    *KLEVS*      NUMBER OF SOIL LAYERS
!    *KSOTY*      SOIL TYPE                                   (1-7)

!     INPUT PARAMETERS (REAL):
!    *PTMST*      TIME STEP                                      S

!     INPUT PARAMETERS (LOGICAL):
!    *LDLAND*     LAND/SEA MASK (TRUE/FALSE)

!     INPUT PARAMETERS AT T-1 OR CONSTANT IN TIME (REAL):
!    *PTSAM1M*    SOIL TEMPERATURE                               K
!    *PWSAM1M*    SOIL MOISTURE                                m**3/m**3
!    *PSLRFL*     NET LONGWAVE  RADIATION AT THE SURFACE        W/M**2
!    *PGSN*       GROUND HEAT FLUX FROM SNOW DECK TO SOIL       W/M2
!    *PCTSA*      VOLUMETRIC HEAT CAPACITY                      J/K/M**3
!    *PFRTI*      TILE FRACTIONS                              (0-1)
!            1 : WATER                  5 : SNOW ON LOW-VEG+BARE-SOIL
!            2 : ICE                    6 : DRY SNOW-FREE HIGH-VEG
!            3 : WET SKIN               7 : SNOW UNDER HIGH-VEG
!            4 : DRY SNOW-FREE LOW-VEG  8 : BARE SOIL
!    *PAHFSTI*    TILE SURFACE SENSIBLE HEAT FLUX                 W/M2
!    *PEVAPTI*    TILE SURFACE MOISTURE FLUX                     KG/M2/S
!    *PSSRFLTI*   TILE NET SHORTWAVE RADIATION FLUX AT SURFACE    W/M2

!     UPDATED PARAMETERS AT T+1 (UNFILTERED,REAL):
!    *PTSA*       SOIL TEMPERATURE                               K


!     METHOD.
!     -------

!          Parameters are set and the tridiagonal solver is called.

!     EXTERNALS.
!     ----------
!     *SRFWDIFS*

!     REFERENCE.
!     ----------
!          See documentation.

!     Original
!     --------
!          Simplified version based on SRFT
!     M. Janiskova   E.C.M.W.F.     26-07-2011  

!     Modifications
!     -------------

!     ------------------------------------------------------------------
#endif

IMPLICIT NONE

! Declaration of arguments

INTEGER(KIND=JPIM), INTENT(IN)   :: KIDIA
INTEGER(KIND=JPIM), INTENT(IN)   :: KFDIA
INTEGER(KIND=JPIM), INTENT(IN)   :: KLON
INTEGER(KIND=JPIM), INTENT(IN)   :: KLEVS
INTEGER(KIND=JPIM), INTENT(IN)   :: KSOTY(:)

REAL(KIND=JPRB),    INTENT(IN)   :: PTMST
REAL(KIND=JPRB),    INTENT(IN)   :: PTSAM1M(:,:)
REAL(KIND=JPRB),    INTENT(IN)   :: PWSAM1M(:,:)
REAL(KIND=JPRB),    INTENT(IN)   :: PFRTI(:,:)
REAL(KIND=JPRB),    INTENT(IN)   :: PAHFSTI(:,:)
REAL(KIND=JPRB),    INTENT(IN)   :: PEVAPTI(:,:)
REAL(KIND=JPRB),    INTENT(IN)   :: PSLRFL(:)
REAL(KIND=JPRB),    INTENT(IN)   :: PSSRFLTI(:,:)
REAL(KIND=JPRB),    INTENT(IN)   :: PGSN(:)
REAL(KIND=JPRB),    INTENT(IN)   :: PCTSA(:,:)

LOGICAL,            INTENT(IN)   :: LDLAND(:)

TYPE(TCST),         INTENT(IN)   :: YDCST
TYPE(TSOIL),        INTENT(IN)   :: YDSOIL
TYPE(TFLAKE),       INTENT(IN)   :: YDFLAKE

REAL(KIND=JPRB),    INTENT(OUT)  :: PTSA(:,:)

!      LOCAL VARIABLES

REAL(KIND=JPRB) :: ZSURFL(KLON)
REAL(KIND=JPRB) :: ZRHS(KLON,KLEVS), ZCDZ(KLON,KLEVS),&
 & ZLST(KLON,KLEVS), ZDIF(KLON,KLEVS),&
 & ZTSA(KLON,KLEVS)  
LOGICAL ::LLALLAYS

INTEGER(KIND=JPIM) :: JK, JL, JS, IKLEVS

REAL(KIND=JPRB) :: ZCONS1, ZCONS2, ZSLRFL, ZSSRFL, ZTHFL,&
 & ZFF, ZWU, ZLIC, ZLWT, ZLAMBDASAT, ZKERSTEN, ZINVWSAT, ZCOND

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! -------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SRFTS_MOD:SRFTS',0,ZHOOK_HANDLE)
ASSOCIATE(RLVTT=>YDCST%RLVTT, &
 & LEFLAKE=>YDFLAKE%LEFLAKE, &
 & LEVGEN=>YDSOIL%LEVGEN, RDAT=>YDSOIL%RDAT, RFRSMALL=>YDSOIL%RFRSMALL, &
 & RKERST1=>YDSOIL%RKERST1, RKERST2=>YDSOIL%RKERST2, RKERST3=>YDSOIL%RKERST3, &
 & RLAMBDADRY=>YDSOIL%RLAMBDADRY, RLAMBDADRYM=>YDSOIL%RLAMBDADRYM, &
 & RLAMBDAICE=>YDSOIL%RLAMBDAICE, RLAMBDAWAT=>YDSOIL%RLAMBDAWAT, &
 & RLAMSAT1=>YDSOIL%RLAMSAT1, RLAMSAT1M=>YDSOIL%RLAMSAT1M, RSIMP=>YDSOIL%RSIMP, &
 & RWSAT=>YDSOIL%RWSAT, RWSATM=>YDSOIL%RWSATM)

!*    0. INITIALIZATION
!     ------------------

DO JL=KIDIA,KFDIA
  ZSURFL(JL) = 0.0_JPRB
ENDDO

DO JK=1,KLEVS
  DO JL=KIDIA,KFDIA
    ZDIF(JL,JK) = 0.0_JPRB
    ZLST(JL,JK) = 0.0_JPRB
    ZCDZ(JL,JK) = 0.0_JPRB
    ZRHS(JL,JK) = 0.0_JPRB
  ENDDO
ENDDO

!* Computation done for only top or all soil layers

LLALLAYS = .TRUE.    ! done for all layers
!LLALLAYS = .FALSE.   ! done for top layer only

IF (LLALLAYS) THEN
  IKLEVS = KLEVS
ELSE
  IKLEVS = 1
ENDIF

!*         1. SET UP SOME CONSTANTS.
!             --- -- ---- ----------
!*    COMPUTATIONAL CONSTANTS.
!     ------------- ----------

ZCONS1=PTMST*RSIMP
ZCONS2=1.0_JPRB-1.0_JPRB/RSIMP

!*         2. Compute net heat flux at the surface.
!             -------------------------------------

DO JL=KIDIA,KFDIA
  IF (LDLAND(JL)) THEN

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
    ZSLRFL=(PFRTI(JL,3)+PFRTI(JL,4)+PFRTI(JL,6)+PFRTI(JL,8))*PSLRFL(JL)
    ZTHFL=PFRTI(JL,3)*(PAHFSTI(JL,3)+RLVTT*PEVAPTI(JL,3))&
     & +PFRTI(JL,4)*(PAHFSTI(JL,4)+RLVTT*PEVAPTI(JL,4))&
     & +PFRTI(JL,6)*(PAHFSTI(JL,6)+RLVTT*PEVAPTI(JL,6))&
     & +PFRTI(JL,8)*(PAHFSTI(JL,8)+RLVTT*PEVAPTI(JL,8))  
     
    ZSURFL(JL)=ZSSRFL+ZSLRFL+ZTHFL+PGSN(JL)

    IF ( LEFLAKE ) THEN
      IF ( PFRTI(JL,9) > RFRSMALL ) THEN
        ZSURFL(JL)=PGSN(JL)
        IF ( (PFRTI(JL,3)+PFRTI(JL,4)+PFRTI(JL,6)+PFRTI(JL,8)) > RFRSMALL ) THEN
          ZSURFL(JL)=PGSN(JL)+(ZSSRFL+ZSLRFL+ZTHFL) & 
           & / (PFRTI(JL,3)+PFRTI(JL,4)+PFRTI(JL,6)+PFRTI(JL,8))
        ENDIF
      ENDIF
    ENDIF
    
  ELSE
    ZSURFL(JL)=0.0_JPRB
  ENDIF
ENDDO

!*         3. Compute exchange coeff. layer by layer
!             --------------------------------------

DO JK=1,IKLEVS
  DO JL=KIDIA,KFDIA
    IF (LDLAND(JL)) THEN
      ZFF=0.0_JPRB
      ZWU=PWSAM1M(JL,JK)*(1.0_JPRB-ZFF)
      ZLWT=RLAMBDAWAT**ZWU
      IF(LEVGEN)THEN
        JS=KSOTY(JL)
        ZINVWSAT=1.0_JPRB/RWSATM(JS)
        ZLIC=RLAMBDAICE**(RWSATM(JS)-ZWU)
        ZLAMBDASAT=RLAMSAT1M(JS)*ZLIC*ZLWT
        ZCOND=MAX(RKERST1,PWSAM1M(JL,JK)*ZINVWSAT)
        ZKERSTEN=RKERST2*LOG10(ZCOND)+RKERST3
        ZDIF(JL,JK)=RLAMBDADRYM(JS)+ZKERSTEN*(ZLAMBDASAT-RLAMBDADRYM(JS))
      ELSE
        ZINVWSAT=1.0_JPRB/RWSAT
        ZLIC=RLAMBDAICE**(RWSAT-ZWU)
        ZLAMBDASAT=RLAMSAT1*ZLIC*ZLWT
        ZCOND=MAX(RKERST1,PWSAM1M(JL,JK)*ZINVWSAT)
        ZKERSTEN=RKERST2*LOG10(ZCOND)+RKERST3
        ZDIF(JL,JK)=RLAMBDADRY+ZKERSTEN*(ZLAMBDASAT-RLAMBDADRY)
      ENDIF
    ENDIF
  ENDDO
ENDDO

!*         4. Set arrays
!             ----------
!     Layer 1
JK=1
DO JL=KIDIA,KFDIA
  IF (LDLAND(JL)) THEN
    ZLST(JL,JK)=ZCONS1*(ZDIF(JL,JK)+ZDIF(JL,JK+1))/(RDAT(JK)+RDAT(JK+1))
    ZCDZ(JL,JK)=PCTSA(JL,JK)*RDAT(JK)
    ZRHS(JL,JK)=PTMST*ZSURFL(JL)/ZCDZ(JL,JK)
  ENDIF
ENDDO

IF (LLALLAYS) THEN

!     Layers 2 to KLEVS-1
  DO JK=2,KLEVS-1
    DO JL=KIDIA,KFDIA
      IF (LDLAND(JL)) THEN
        ZLST(JL,JK)=ZCONS1*(ZDIF(JL,JK)+ZDIF(JL,JK+1))/(RDAT(JK)+RDAT(JK+1))
        ZCDZ(JL,JK)=PCTSA(JL,JK)*RDAT(JK)
      ENDIF
    ENDDO
  ENDDO

!     Layers KLEVS
  JK=KLEVS
  DO JL=KIDIA,KFDIA
    IF (LDLAND(JL)) THEN
      ZLST(JL,JK)=0.0_JPRB
      ZCDZ(JL,JK)=PCTSA(JL,JK)*RDAT(JK)
    ENDIF
  ENDDO
ENDIF

!*         5. Call tridiagonal solver
!             -----------------------
CALL SRFWDIFS(KIDIA,KFDIA,KLON,KLEVS,PTSAM1M,ZLST,ZRHS,ZCDZ,ZTSA, &
 & LDLAND,LLALLAYS,YDSOIL)


!*         7. New temperatures
!             ----------------
DO JK=1,KLEVS
  DO JL=KIDIA,KFDIA
    IF (LDLAND(JL)) THEN
      PTSA(JL,JK)=PTSAM1M(JL,JK)*ZCONS2+ZTSA(JL,JK)
    ELSE
      PTSA(JL,JK)=PTSAM1M(JL,JK)
    ENDIF
  ENDDO
ENDDO

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SRFTS_MOD:SRFTS',1,ZHOOK_HANDLE)

END SUBROUTINE SRFTS
END MODULE SRFTS_MOD
