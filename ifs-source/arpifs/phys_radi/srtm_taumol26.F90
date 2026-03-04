SUBROUTINE SRTM_TAUMOL26&
 & (YDDIMV, KIDIA   , KFDIA    , KLEV,&
 & P_COLMOL  ,K_LAYTROP,&
 & P_SFLUXZEN, P_TAUG   , P_TAUR    , PRMU0&
 & )  

!     Written by Eli J. Mlawer, Atmospheric & Environmental Research.

!     BAND 26:  22650-29000 cm-1 (low - nothing; high - nothing)

!      PARAMETER (MG=16, MXLAY=203, NBANDS=14)

! Modifications
!        M.Hamrud      01-Oct-2003 CY28 Cleaning

!     JJMorcrette 2003-02-24 adapted to ECMWF environment
!        D.Salmond  31-Oct-2007 Vector version in the style of RRTM from Meteo France & NEC
!     JJMorcrette 20110610 Flexible configuration for number of g-points

USE YOMDIMV  , ONLY : TDIMV
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE PARSRTM  , ONLY : JPG
USE YOESRTM  , ONLY : NG26
USE YOESRTA26, ONLY : SFLUXREFC, RAYLC

IMPLICIT NONE

!-- Output
TYPE(TDIMV)       ,INTENT(IN)    :: YDDIMV
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA, KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_COLMOL(KIDIA:KFDIA,YDDIMV%NFLEVG) 
INTEGER(KIND=JPIM),INTENT(IN)    :: K_LAYTROP(KIDIA:KFDIA) 

REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_SFLUXZEN(KIDIA:KFDIA,JPG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_TAUG(KIDIA:KFDIA,YDDIMV%NFLEVG,JPG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_TAUR(KIDIA:KFDIA,YDDIMV%NFLEVG,JPG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRMU0(KIDIA:KFDIA)
!- from AER
!- from INTFAC      
!- from INTIND
!- from PRECISE             
!- from PROFDATA             
!- from SELF             
INTEGER(KIND=JPIM) :: IG, I_LAY, I_LAYSOLFR(KIDIA:KFDIA), I_NLAYERS, IPLON

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('SRTM_TAUMOL26',0,ZHOOK_HANDLE)
ASSOCIATE(NFLEVG=>YDDIMV%NFLEVG)

I_NLAYERS = KLEV

!     Compute the optical depth by interpolating in ln(pressure), 
!     temperature, and appropriate species.  Below LAYTROP, the water
!     vapor self-continuum is interpolated (in temperature) separately.  

I_LAYSOLFR(KIDIA:KFDIA) = K_LAYTROP(KIDIA:KFDIA)

DO I_LAY = 1, I_NLAYERS
  DO IPLON = KIDIA, KFDIA
    IF (PRMU0(IPLON) > 0.0_JPRB) THEN
      IF (I_LAY <= K_LAYTROP(IPLON)) THEN
        !  DO IG = 1, NG(26)
!CDIR UNROLL=NG26
        DO IG = 1 , NG26 
          !    TAUG(LAY,IG) = COLMOL(LAY) * RAYLC(IG)
          !    SSA(LAY,IG) = 1.0
          IF (I_LAY == I_LAYSOLFR(IPLON)) P_SFLUXZEN(IPLON,IG) = SFLUXREFC(IG) 
          P_TAUG(IPLON,I_LAY,IG) = 0.0_JPRB
          P_TAUR(IPLON,I_LAY,IG) = P_COLMOL(IPLON,I_LAY) * RAYLC(IG) 
        ENDDO
      ENDIF
    ENDIF
  ENDDO
ENDDO

DO I_LAY = 1, I_NLAYERS
  DO IPLON = KIDIA, KFDIA
    IF (PRMU0(IPLON) > 0.0_JPRB) THEN
      IF (I_LAY >= K_LAYTROP(IPLON)+1) THEN
        !  DO IG = 1, NG(26)
!CDIR UNROLL=NG26
        DO IG = 1 , NG26
          !    TAUG(LAY,IG) = COLMOL(LAY) * RAYLC(IG)
          !    SSA(LAY,IG) = 1.0
          P_TAUG(IPLON,I_LAY,IG) = 0.0_JPRB
          P_TAUR(IPLON,I_LAY,IG) = P_COLMOL(IPLON,I_LAY) * RAYLC(IG) 
        ENDDO
      ENDIF
    ENDIF
  ENDDO
ENDDO

!-----------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SRTM_TAUMOL26',1,ZHOOK_HANDLE)
END SUBROUTINE SRTM_TAUMOL26
