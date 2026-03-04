SUBROUTINE SRTM_TAUMOL17&
 & (YDDIMV, KIDIA   , KFDIA    , KLEV,&
 & P_FAC00   , P_FAC01  , P_FAC10   , P_FAC11,&
 & K_JP      , K_JT     , K_JT1     , P_ONEMINUS,&
 & P_COLH2O  , P_COLCO2 , P_COLMOL,&
 & K_LAYTROP , P_SELFFAC, P_SELFFRAC, K_INDSELF  , P_FORFAC, P_FORFRAC, K_INDFOR,&
 & P_SFLUXZEN, P_TAUG   , P_TAUR    , PRMU0&
 & )  

!     Written by Eli J. Mlawer, Atmospheric & Environmental Research.

!     BAND 17:  3250-4000 cm-1 (low - H2O,CO2; high - H2O,CO2)

! Modifications
!        M.Hamrud      01-Oct-2003 CY28 Cleaning

!     JJMorcrette 2003-02-24 adapted to ECMWF environment
!        D.Salmond  31-Oct-2007 Vector version in the style of RRTM from Meteo France & NEC
!     JJMorcrette 20010610 Flexible configuration for number of g-points

USE YOMDIMV  , ONLY : TDIMV
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE PARSRTM  , ONLY : JPG
USE YOESRTM  , ONLY : NG17
USE YOESRTA17, ONLY : ABSA, ABSB, FORREFC, SELFREFC, SFLUXREFC, RAYL, LAYREFFR, STRRAT  
USE YOESRTWN , ONLY : NSPA, NSPB

IMPLICIT NONE

!-- Output
TYPE(TDIMV)       ,INTENT(IN)    :: YDDIMV
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA, KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_FAC00(KIDIA:KFDIA,YDDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_FAC01(KIDIA:KFDIA,YDDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_FAC10(KIDIA:KFDIA,YDDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_FAC11(KIDIA:KFDIA,YDDIMV%NFLEVG) 
INTEGER(KIND=JPIM),INTENT(IN)    :: K_JP(KIDIA:KFDIA,YDDIMV%NFLEVG) 
INTEGER(KIND=JPIM),INTENT(IN)    :: K_JT(KIDIA:KFDIA,YDDIMV%NFLEVG) 
INTEGER(KIND=JPIM),INTENT(IN)    :: K_JT1(KIDIA:KFDIA,YDDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_ONEMINUS(KIDIA:KFDIA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_COLH2O(KIDIA:KFDIA,YDDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_COLCO2(KIDIA:KFDIA,YDDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_COLMOL(KIDIA:KFDIA,YDDIMV%NFLEVG) 
INTEGER(KIND=JPIM),INTENT(IN)    :: K_LAYTROP(KIDIA:KFDIA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_SELFFAC(KIDIA:KFDIA,YDDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_SELFFRAC(KIDIA:KFDIA,YDDIMV%NFLEVG) 
INTEGER(KIND=JPIM),INTENT(IN)    :: K_INDSELF(KIDIA:KFDIA,YDDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_FORFAC(KIDIA:KFDIA,YDDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_FORFRAC(KIDIA:KFDIA,YDDIMV%NFLEVG) 
INTEGER(KIND=JPIM),INTENT(IN)    :: K_INDFOR(KIDIA:KFDIA,YDDIMV%NFLEVG) 

REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_SFLUXZEN(KIDIA:KFDIA,JPG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_TAUG(KIDIA:KFDIA,YDDIMV%NFLEVG,JPG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_TAUR(KIDIA:KFDIA,YDDIMV%NFLEVG,JPG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRMU0(KIDIA:KFDIA)
!- from INTFAC      
!- from INTIND
!- from PRECISE             
!- from PROFDATA             
!- from SELF             
INTEGER(KIND=JPIM) :: IG, IND0, IND1, INDS, INDF, JS, I_LAY, I_LAYSOLFR, I_NLAYERS, IPLON

! REAL(KIND=JPRB) :: Z_FAC000, Z_FAC001, Z_FAC010, Z_FAC011, Z_FAC100, Z_FAC101,&
!  & Z_FAC110, Z_FAC111
REAL(KIND=JPRB) :: Z_FS, Z_SPECCOMB, Z_SPECMULT, Z_SPECPARM, Z_TAURAY  
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('SRTM_TAUMOL17',0,ZHOOK_HANDLE)
ASSOCIATE(NFLEVG=>YDDIMV%NFLEVG)

I_NLAYERS = KLEV

!     Compute the optical depth by interpolating in ln(pressure), 
!     temperature, and appropriate species.  Below LAYTROP, the water
!     vapor self-continuum is interpolated (in temperature) separately.  

DO I_LAY = 1, I_NLAYERS
  DO IPLON = KIDIA, KFDIA
    IF (PRMU0(IPLON) > 0.0_JPRB) THEN
      IF (I_LAY <= K_LAYTROP(IPLON)) THEN
        Z_SPECCOMB = P_COLH2O(IPLON,I_LAY) + STRRAT*P_COLCO2(IPLON,I_LAY)
        Z_SPECPARM = P_COLH2O(IPLON,I_LAY)/Z_SPECCOMB 
        IF (Z_SPECPARM >= P_ONEMINUS(IPLON)) Z_SPECPARM = P_ONEMINUS(IPLON)
        Z_SPECMULT = 8.*(Z_SPECPARM)
        JS = 1 + INT(Z_SPECMULT)
        Z_FS = MOD(Z_SPECMULT, 1.0_JPRB )
        ! Z_FAC000 = (1. - Z_FS) * P_FAC00(I_LAY)
        ! Z_FAC010 = (1. - Z_FS) * P_FAC10(I_LAY)
        ! Z_FAC100 = Z_FS * P_FAC00(I_LAY)
        ! Z_FAC110 = Z_FS * P_FAC10(I_LAY)
        ! Z_FAC001 = (1. - Z_FS) * P_FAC01(I_LAY)
        ! Z_FAC011 = (1. - Z_FS) * P_FAC11(I_LAY)
        ! Z_FAC101 = Z_FS * P_FAC01(I_LAY)
        ! Z_FAC111 = Z_FS * P_FAC11(I_LAY)
        IND0 = ((K_JP(IPLON,I_LAY)-1)*5+(K_JT(IPLON,I_LAY)-1))*NSPA(17) + JS
        IND1 = (K_JP(IPLON,I_LAY)*5+(K_JT1(IPLON,I_LAY)-1))*NSPA(17) + JS
        INDS = K_INDSELF(IPLON,I_LAY)
        INDF = K_INDFOR(IPLON,I_LAY)
        Z_TAURAY = P_COLMOL(IPLON,I_LAY) * RAYL

        !  DO IG = 1, NG(17)
!CDIR UNROLL=NG17
        DO IG = 1, NG17
          P_TAUG(IPLON,I_LAY,IG) = Z_SPECCOMB * &
           & (&
           & (1. - Z_FS) * ( ABSA(IND0,IG) * P_FAC00(IPLON,I_LAY) + &
           &                 ABSA(IND0+9,IG) * P_FAC10(IPLON,I_LAY) + &
           &                 ABSA(IND1,IG) * P_FAC01(IPLON,I_LAY) + &
           &                 ABSA(IND1+9,IG) * P_FAC11(IPLON,I_LAY) ) + &
           & Z_FS        * ( ABSA(IND0+1,IG) * P_FAC00(IPLON,I_LAY) + &
           &                 ABSA(IND0+10,IG) * P_FAC10(IPLON,I_LAY) + &
           &                 ABSA(IND1+1,IG) * P_FAC01(IPLON,I_LAY) + &
           &                 ABSA(IND1+10,IG) * P_FAC11(IPLON,I_LAY) ) &
           & )  + &
           & P_COLH2O(IPLON,I_LAY) * &
           & (P_SELFFAC(IPLON,I_LAY) * (SELFREFC(INDS,IG) + &
           & P_SELFFRAC(IPLON,I_LAY) * &
           & (SELFREFC(INDS+1,IG) - SELFREFC(INDS,IG))) + &
           & P_FORFAC(IPLON,I_LAY) * (FORREFC(INDF,IG) + &
           & P_FORFRAC(IPLON,I_LAY) * &
           & (FORREFC(INDF+1,IG) - FORREFC(INDF,IG))))   

          !     &           + TAURAY
          !    SSA(LAY,IG) = TAURAY/TAUG(LAY,IG)
          P_TAUR(IPLON,I_LAY,IG) = Z_TAURAY
        ENDDO
      ENDIF
    ENDIF
  ENDDO
ENDDO

I_LAYSOLFR = I_NLAYERS

DO I_LAY = 1, I_NLAYERS
  DO IPLON = KIDIA, KFDIA
    IF (PRMU0(IPLON) > 0.0_JPRB) THEN
      IF (I_LAY >= K_LAYTROP(IPLON)+1) THEN
        IF (K_JP(IPLON,I_LAY-1) < LAYREFFR .AND. K_JP(IPLON,I_LAY) >= LAYREFFR) &
         & I_LAYSOLFR = I_LAY  
        Z_SPECCOMB = P_COLH2O(IPLON,I_LAY) + STRRAT*P_COLCO2(IPLON,I_LAY)
        Z_SPECPARM = P_COLH2O(IPLON,I_LAY)/Z_SPECCOMB 
        IF (Z_SPECPARM >= P_ONEMINUS(IPLON)) Z_SPECPARM = P_ONEMINUS(IPLON)
        Z_SPECMULT = 4.*(Z_SPECPARM)
        JS = 1 + INT(Z_SPECMULT)
        Z_FS = MOD(Z_SPECMULT, 1.0_JPRB )
        ! Z_FAC000 = (1. - Z_FS) * P_FAC00(I_LAY)
        ! Z_FAC010 = (1. - Z_FS) * P_FAC10(I_LAY)
        ! Z_FAC100 = Z_FS * P_FAC00(I_LAY)
        ! Z_FAC110 = Z_FS * P_FAC10(I_LAY)
        ! Z_FAC001 = (1. - Z_FS) * P_FAC01(I_LAY)
        ! Z_FAC011 = (1. - Z_FS) * P_FAC11(I_LAY)
        ! Z_FAC101 = Z_FS * P_FAC01(I_LAY)
        ! Z_FAC111 = Z_FS * P_FAC11(I_LAY)
        IND0 = ((K_JP(IPLON,I_LAY)-13)*5+(K_JT(IPLON,I_LAY)-1))*NSPB(17) + JS
        IND1 = ((K_JP(IPLON,I_LAY)-12)*5+(K_JT1(IPLON,I_LAY)-1))*NSPB(17) + JS
        INDF = K_INDFOR(IPLON,I_LAY)
        Z_TAURAY = P_COLMOL(IPLON,I_LAY) * RAYL

        !  DO IG = 1, NG(17)
!CDIR UNROLL=NG17
        DO IG = 1, NG17
          P_TAUG(IPLON,I_LAY,IG) = Z_SPECCOMB * &
           !    & (Z_FAC000 * ABSB(IND0,IG) + &
           !    & Z_FAC100 * ABSB(IND0+1,IG) + &
           !    & Z_FAC010 * ABSB(IND0+5,IG) + &
           !    & Z_FAC110 * ABSB(IND0+6,IG) + &
           !    & Z_FAC001 * ABSB(IND1,IG) + &
           !    & Z_FAC101 * ABSB(IND1+1,IG) + &
           !    & Z_FAC011 * ABSB(IND1+5,IG) + &
           !    & Z_FAC111 * ABSB(IND1+6,IG)) + &
           & (&
           & (1. - Z_FS) * ( ABSB(IND0,IG) * P_FAC00(IPLON,I_LAY) + &
           &                 ABSB(IND0+5,IG) * P_FAC10(IPLON,I_LAY) + &
           &                 ABSB(IND1,IG) * P_FAC01(IPLON,I_LAY) + &
           &                 ABSB(IND1+5,IG) * P_FAC11(IPLON,I_LAY) ) + &
           & Z_FS        * ( ABSB(IND0+1,IG) * P_FAC00(IPLON,I_LAY) + &
           &                 ABSB(IND0+6,IG) * P_FAC10(IPLON,I_LAY) + &
           &                 ABSB(IND1+1,IG) * P_FAC01(IPLON,I_LAY) + &
           &                 ABSB(IND1+6,IG) * P_FAC11(IPLON,I_LAY) ) &
           & ) + &
           & P_COLH2O(IPLON,I_LAY) * &
           & P_FORFAC(IPLON,I_LAY) * (FORREFC(INDF,IG) + &
           & P_FORFRAC(IPLON,I_LAY) * &
           & (FORREFC(INDF+1,IG) - FORREFC(INDF,IG)))   
          !     &           + TAURAY
          !    SSA(LAY,IG) = TAURAY/TAUG(LAY,IG)
          IF (I_LAY == I_LAYSOLFR) P_SFLUXZEN(IPLON,IG) = SFLUXREFC(IG,JS) &
           & + Z_FS * (SFLUXREFC(IG,JS+1) - SFLUXREFC(IG,JS))  
          P_TAUR(IPLON,I_LAY,IG) = Z_TAURAY
        ENDDO
      ENDIF
    ENDIF
  ENDDO
ENDDO

!-----------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SRTM_TAUMOL17',1,ZHOOK_HANDLE)
END SUBROUTINE SRTM_TAUMOL17
