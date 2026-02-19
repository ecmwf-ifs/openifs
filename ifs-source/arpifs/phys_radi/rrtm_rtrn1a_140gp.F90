!option! -O extendreorder
!option! -pvctl no_on_adb
!option! -pvctl rsa_on_adb
SUBROUTINE RRTM_RTRN1A_140GP(YDDIMV,KIDIA,KFDIA,KLEV,K_ISTART,K_IEND,K_ICLDLYR,P_CLDFRAC,P_TAUCLD,P_ABSS1,&
 & P_OD,P_TAUSF1,P_TOTDFLUC,P_TOTDFLUX,P_TOTUFLUC,P_TOTUFLUX,&
 & P_TAVEL,P_TZ,P_TBOUND,PFRAC,P_SEMISS,P_SEMISLW)  

!-* This program calculates the upward fluxes, downward fluxes,
!   and heating rates for an arbitrary atmosphere.  The input to
!   this program is the atmospheric profile and all Planck function
!   information.  First-order "numerical" quadrature is used for the 
!   angle integration, i.e. only one exponential is computed per layer
!   per g-value per band.  Cloud overlap is treated with a generalized
!   maximum/random method in which adjacent cloud layers are treated
!   with maximum overlap, and non-adjacent cloud groups are treated
!   with random overlap.  For adjacent cloud layers, cloud information
!   is carried from the previous two layers.

!     AUTHOR.
!     -------
!      JJMorcrette, ECMWF

!     MODIFICATIONS.
!     --------------
!      NEC           25-Oct-2007 Optimisations
!      D. Salmond    11-Dec-2007 Optimizations
!      O. Riviere    15-Sep-2008 Deallocate Z_BBU1, Z_BBUTOT1, Z_ATOT1
!      NEC/FC        05-Oct-2009 Optimisations
!      P.Marguinaud  05-Jul-2011 Fix array bounds + Cleaning
!     JJMorcrette 20110613 flexible number of g-points
!      Y.Seity       23-Mar-2015 Fix array bounds problem
!      R. El Khatib  19-may-2016 Optimization for Intel compiler
! ---------------------------------------------------------------------------

USE YOMDIMV  , ONLY : TDIMV
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE PARRRTM  , ONLY : JPBAND
USE YOERRTM  , ONLY : JPGPT
USE YOERRTAB , ONLY : BPADE
USE YOERRTWN , ONLY : TOTPLNK, DELWAVE
USE YOERRTFTR, ONLY : NGB

IMPLICIT NONE

TYPE(TDIMV)       ,INTENT(IN)    :: YDDIMV
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: K_ISTART 
INTEGER(KIND=JPIM),INTENT(IN)    :: K_IEND 
INTEGER(KIND=JPIM),INTENT(IN)    :: K_ICLDLYR(KIDIA:KFDIA,YDDIMV%NFLEVG) ! Cloud indicator
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_CLDFRAC(KIDIA:KFDIA,YDDIMV%NFLEVG) ! Cloud fraction
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_TAUCLD(KIDIA:KFDIA,YDDIMV%NFLEVG,JPBAND) ! Spectral optical thickness
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_ABSS1(KIDIA:KFDIA,JPGPT,YDDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_OD(KIDIA:KFDIA,JPGPT,YDDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_TAUSF1(KIDIA:KFDIA,JPGPT,YDDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_TOTDFLUC(KIDIA:KFDIA,0:YDDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_TOTDFLUX(KIDIA:KFDIA,0:YDDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_TOTUFLUC(KIDIA:KFDIA,0:YDDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_TOTUFLUX(KIDIA:KFDIA,0:YDDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_TAVEL(KIDIA:KFDIA,YDDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_TZ(KIDIA:KFDIA,0:YDDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_TBOUND(KIDIA:KFDIA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFRAC(KIDIA:KFDIA,JPGPT,YDDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_SEMISS(KIDIA:KFDIA,JPBAND) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_SEMISLW(KIDIA:KFDIA)

! ---------------------------------------------------------------------------

INTEGER(KIND=JPIM) :: INDLAY(KIDIA:KFDIA,YDDIMV%NFLEVG),INDLEV(KIDIA:KFDIA,0:YDDIMV%NFLEVG)

REAL(KIND=JPRB) :: Z_BBU1(KIDIA:KFDIA,JPGPT*YDDIMV%NFLEVG),Z_BBUTOT1(KIDIA:KFDIA,JPGPT*YDDIMV%NFLEVG)
REAL(KIND=JPRB) :: Z_TLAYFRAC(KIDIA:KFDIA,YDDIMV%NFLEVG),Z_TLEVFRAC(KIDIA:KFDIA,0:YDDIMV%NFLEVG)
REAL(KIND=JPRB) :: Z_BGLEV(KIDIA:KFDIA,JPGPT)
!-- DS_000515
REAL(KIND=JPRB) :: Z_PLVL(KIDIA:KFDIA,JPBAND+1,0:YDDIMV%NFLEVG),Z_PLAY(KIDIA:KFDIA,JPBAND+1,0:YDDIMV%NFLEVG)
REAL(KIND=JPRB) :: Z_WTNUM = 0.5_JPRB
!-- DS_000515
!REAL(KIND=JPRB) :: Z_ODCLDNW(JPGPT,NFLEVG)
REAL(KIND=JPRB) :: Z_SEMIS(KIDIA:KFDIA,JPGPT),Z_RADUEMIT(KIDIA:KFDIA,JPGPT)

REAL(KIND=JPRB) :: Z_RADCLRU1(KIDIA:KFDIA,JPGPT) ,Z_RADCLRD1(KIDIA:KFDIA,JPGPT)
REAL(KIND=JPRB) :: Z_RADLU1(KIDIA:KFDIA,JPGPT)   ,Z_RADLD1(KIDIA:KFDIA,JPGPT)
!-- DS_000515
REAL(KIND=JPRB) :: Z_TRNCLD(KIDIA:KFDIA,YDDIMV%NFLEVG,JPBAND+1)
!-- DS_000515
!REAL(KIND=JPRB) :: Z_ABSCLDNW(JPGPT,NFLEVG)
REAL(KIND=JPRB) :: Z_ATOT1(KIDIA:KFDIA,JPGPT*YDDIMV%NFLEVG)

REAL(KIND=JPRB) :: Z_SURFEMIS(KIDIA:KFDIA,JPBAND),Z_PLNKEMIT(KIDIA:KFDIA,JPBAND)

! dimension of arrays required for cloud overlap calculations

REAL(KIND=JPRB) :: Z_CLRRADU(KIDIA:KFDIA,JPGPT),Z_CLDRADU(KIDIA:KFDIA,JPGPT),&
 & Z_OLDCLD(KIDIA:KFDIA,JPGPT)
REAL(KIND=JPRB) :: Z_OLDCLR(KIDIA:KFDIA,JPGPT),Z_RAD(KIDIA:KFDIA,JPGPT)
REAL(KIND=JPRB) :: Z_FACCLD1(KIDIA:KFDIA,YDDIMV%NFLEVG+1),Z_FACCLD2(KIDIA:KFDIA,YDDIMV%NFLEVG+1)
REAL(KIND=JPRB) :: Z_FACCLR1(KIDIA:KFDIA,YDDIMV%NFLEVG+1),Z_FACCLR2(KIDIA:KFDIA,YDDIMV%NFLEVG+1)
REAL(KIND=JPRB) :: Z_FACCMB1(KIDIA:KFDIA,YDDIMV%NFLEVG+1),Z_FACCMB2(KIDIA:KFDIA,YDDIMV%NFLEVG+1)
REAL(KIND=JPRB) :: Z_FACCLD1D(KIDIA:KFDIA,0:YDDIMV%NFLEVG),Z_FACCLD2D(KIDIA:KFDIA,0:YDDIMV%NFLEVG),&
                 & Z_FACCLR1D(KIDIA:KFDIA,0:YDDIMV%NFLEVG)
REAL(KIND=JPRB) :: Z_FACCLR2D(KIDIA:KFDIA,0:YDDIMV%NFLEVG),Z_FACCMB1D(KIDIA:KFDIA,0:YDDIMV%NFLEVG),&
                 & Z_FACCMB2D(KIDIA:KFDIA,0:YDDIMV%NFLEVG)
REAL(KIND=JPRB) :: Z_CLRRADD(KIDIA:KFDIA,JPGPT),Z_CLDRADD(KIDIA:KFDIA,JPGPT)
INTEGER(KIND=JPIM) :: ISTCLD(KIDIA:KFDIA,YDDIMV%NFLEVG+1),ISTCLDD(KIDIA:KFDIA,0:YDDIMV%NFLEVG)
!******

!REAL_B :: ZPLVL(JPGPT+1,NFLEVG)  ,ZPLAY(JPGPT+1,NFLEVG)
!REAL_B :: ZTRNCLD(JPGPT+1,NFLEVG),ZTAUCLD(JPGPT+1,NFLEVG)

INTEGER(KIND=JPIM) :: JBAND, ICLDDN_(KIDIA:KFDIA), IENT, INDBOUND(KIDIA:KFDIA), INDX, JI, JLEV, I_NBI, ILEV
INTEGER(KIND=JPIM) :: JLON

REAL(KIND=JPRB) :: Z_BBD, Z_BBDTOT, Z_BGLAY, Z_CLDSRC, Z_DBDTLAY, Z_DBDTLEV,&
 & Z_DELBGDN, Z_DELBGUP, Z_FACTOT1,&
 & Z_FMAX(KIDIA:KFDIA), Z_FMIN(KIDIA:KFDIA), Z_GASSRC, Z_ODSM, Z_PLANKBND,&
 & Z_RADCLD, Z_RADD, Z_RADMOD, Z_RAT1(KIDIA:KFDIA), Z_RAT2(KIDIA:KFDIA), Z_SUMPL(KIDIA:KFDIA),&
 & Z_SUMPLEM(KIDIA:KFDIA), Z_TBNDFRAC(KIDIA:KFDIA), Z_TRNS, Z_TTOT, ZEXTAU  
!REAL(KIND=JPRB) :: Z_URAD1_(KIDIA:KFDIA,KLEV), Z_URADCL1_(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB) :: Z_URADCL1__(KIDIA:KFDIA), Z_URAD1__(KIDIA:KFDIA)
!REAL(KIND=JPRB) :: Z_DRAD1__(KIDIA:KFDIA,KLEV), Z_DRADCL1__(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB) :: Z_DRAD1__(KIDIA:KFDIA), Z_DRADCL1__(KIDIA:KFDIA)
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!NEC

!--------------------------------------------------------------------------
! Input
!  NFLEVG                 ! Maximum number of model layers
!  JPGPT                 ! Total number of g-point subintervals
!  JPBAND                ! Number of longwave spectral bands
!  SECANG                ! Diffusivity angle
!  WTNUM                 ! Weight for radiance to flux conversion
!  KLEV                  ! Number of model layers
!  PAVEL(NFLEVG)          ! Mid-layer pressures (hPa)
!  TAVEL(NFLEVG)          ! Mid-layer temperatures (K)
!  TZ(0:NFLEVG)           ! Interface temperatures (K)
!  TBOUND                ! Surface temperature
!  CLDFRAC(NFLEVG)        ! Layer cloud fraction
!  TAUCLD(NFLEVG,JPBAND)  ! Layer cloud optical thickness
!  ITR
!  PFRAC(JPGPT,NFLEVG)    ! Planck function fractions
!  ICLDLYR(KIDIA:KFDIA,NFLEVG)        ! Flag for cloudy layers
!  ICLD                  ! Flag for cloudy column
!  SEMISS(JPBAND)        ! Surface spectral emissivity
!  BPADE                 ! Pade constant
!  OD                    ! Clear-sky optical thickness
!  TAUSF1                ! 
!  ABSS1                 !  

!  ABSS(JPGPT*NFLEVG)     !
!  ABSCLD(NFLEVG)         !
!  ATOT(JPGPT*NFLEVG)     !
!  ODCLR(JPGPT,NFLEVG)    ! 
!  ODCLD(JPBAND,NFLEVG)   !
!  EFCLFR1(JPBAND,NFLEVG) ! Effective cloud fraction
!  RADLU(JPGPT)          ! Upward radiance
!  URAD                  ! Spectrally summed upward radiance
!  RADCLRU(JPGPT)        ! Clear-sky upward radiance
!  CLRURAD               ! Spectrally summed clear-sky upward radiance
!  RADLD(JPGPT)          ! Downward radiance
!  DRAD                  ! Spectrally summed downward radiance
!  RADCLRD(JPGPT)        ! Clear-sky downward radiance
!  CLRDRAD               ! Spectrally summed clear-sky downward radiance

! Output
!  TOTUFLUX(0:NFLEVG)     ! Upward longwave flux
!  TOTDFLUX(0:NFLEVG)     ! Downward longwave flux
!  TOTUFLUC(0:NFLEVG)     ! Clear-sky upward longwave flux
!  TOTDFLUC(0:NFLEVG)     ! Clear-sky downward longwave flux

! Maximum/Random cloud overlap variables
! for upward radiaitve transfer
!  FACCLR2  fraction of clear radiance from previous layer that needs to 
!           be switched to cloudy stream
!  FACCLR1  fraction of the radiance that had been switched in the previous
!           layer from cloudy to clear that needs to be switched back to
!           cloudy in the current layer
!  FACCLD2  fraction of cloudy radiance from previous layer that needs to 
!           be switched to clear stream
!           be switched to cloudy stream
!  FACCLD1  fraction of the radiance that had been switched in the previous
!           layer from clear to cloudy that needs to be switched back to
!           clear in the current layer
! for downward radiaitve transfer
!  FACCLR2D fraction of clear radiance from previous layer that needs to 
!           be switched to cloudy stream
!  FACCLR1D fraction of the radiance that had been switched in the previous
!           layer from cloudy to clear that needs to be switched back to
!           cloudy in the current layer
!  FACCLD2D fraction of cloudy radiance from previous layer that needs to 
!           be switched to clear stream
!           be switched to cloudy stream
!  FACCLD1D fraction of the radiance that had been switched in the previous
!           layer from clear to cloudy that needs to be switched back to
!           clear in the current layer

!--------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('RRTM_RTRN1A_140GP',0,ZHOOK_HANDLE)
ASSOCIATE(NFLEVG=>YDDIMV%NFLEVG)

DO JLON = KIDIA, KFDIA
  !-start JJM_000511
  IF (P_TBOUND(JLON) < 339._JPRB .AND. P_TBOUND(JLON) >= 160._JPRB ) THEN
    INDBOUND(JLON) = P_TBOUND(JLON) - 159._JPRB
    Z_TBNDFRAC(JLON) = P_TBOUND(JLON) - INT(P_TBOUND(JLON))
  ELSEIF (P_TBOUND(JLON) >= 339._JPRB ) THEN
    INDBOUND(JLON) = 180
    Z_TBNDFRAC(JLON) = P_TBOUND(JLON) - 339._JPRB
  ELSEIF (P_TBOUND(JLON) < 160._JPRB ) THEN
    INDBOUND(JLON) = 1
    Z_TBNDFRAC(JLON) = P_TBOUND(JLON) - 160._JPRB
  ENDIF  
ENDDO
!-end JJM_000511
  
!Z_URAD1_ = 0.0_JPRB
!Z_URADCL1_ = 0.0_JPRB
P_TOTUFLUC = 0.0_JPRB
P_TOTDFLUC = 0.0_JPRB
P_TOTUFLUX = 0.0_JPRB
P_TOTDFLUX = 0.0_JPRB
Z_FACCLD1D = 0.0_JPRB
Z_FACCLD2D = 0.0_JPRB
Z_FACCLR1D = 0.0_JPRB
Z_FACCLR2D = 0.0_JPRB
Z_FACCMB1D = 0.0_JPRB
Z_FACCMB2D = 0.0_JPRB
Z_FACCLD1  = 0.0_JPRB
Z_FACCLD2  = 0.0_JPRB
Z_FACCLR1  = 0.0_JPRB
Z_FACCLR2 = 0.0_JPRB
Z_FACCMB1 = 0.0_JPRB
Z_FACCMB2 = 0.0_JPRB


DO JLEV = 0, KLEV
  DO JLON = KIDIA, KFDIA
    !-start JJM_000511
    IF (P_TZ(JLON,JLEV) < 339._JPRB .AND. P_TZ(JLON,JLEV) >= 160._JPRB ) THEN
      INDLEV(JLON,JLEV) = P_TZ(JLON,JLEV) - 159._JPRB
      Z_TLEVFRAC(JLON,JLEV) = P_TZ(JLON,JLEV) - INT(P_TZ(JLON,JLEV))
    ELSEIF (P_TZ(JLON,JLEV) >= 339._JPRB ) THEN
      INDLEV(JLON,JLEV) = 180
      Z_TLEVFRAC(JLON,JLEV) = P_TZ(JLON,JLEV) - 339._JPRB
    ELSEIF (P_TZ(JLON,JLEV) < 160._JPRB ) THEN
      INDLEV(JLON,JLEV) = 1
      Z_TLEVFRAC(JLON,JLEV) = P_TZ(JLON,JLEV) - 160._JPRB
    ENDIF    
    !-end JJM_000511
  ENDDO
ENDDO

  Z_URAD1__  (:) = 0.0_JPRB
  Z_URADCL1__(:) = 0.0_JPRB
  Z_RAT1(:) = 0.0_JPRB
  Z_RAT2(:) = 0.0_JPRB
  !_end_jjm 991209
  Z_SUMPL  (:) = 0.0_JPRB
  Z_SUMPLEM(:) = 0.0_JPRB
  ISTCLD (:,1   ) = 1
  ISTCLDD(:,KLEV) = 1

DO JLEV = 1, KLEV
  DO JLON = KIDIA, KFDIA
    !-- DS_000515
    !-start JJM_000511
    IF (P_TAVEL(JLON,JLEV) < 339._JPRB .AND. P_TAVEL(JLON,JLEV) >= 160._JPRB ) THEN
      INDLAY(JLON,JLEV) = P_TAVEL(JLON,JLEV) - 159._JPRB
      Z_TLAYFRAC(JLON,JLEV) = P_TAVEL(JLON,JLEV) - INT(P_TAVEL(JLON,JLEV))
    ELSEIF (P_TAVEL(JLON,JLEV) >= 339._JPRB ) THEN
      INDLAY(JLON,JLEV) = 180
      Z_TLAYFRAC(JLON,JLEV) = P_TAVEL(JLON,JLEV) - 339._JPRB
    ELSEIF (P_TAVEL(JLON,JLEV) < 160._JPRB ) THEN
      INDLAY(JLON,JLEV) = 1
      Z_TLAYFRAC(JLON,JLEV) = P_TAVEL(JLON,JLEV) - 160._JPRB
    ENDIF  
    !-end JJM_000511
  ENDDO
ENDDO
!-- DS_000515

DO JLEV = 1, KLEV
  DO JLON = KIDIA, KFDIA
    IF (K_ICLDLYR(JLON,JLEV) == 1) THEN

      !mji    
      ISTCLD(JLON,JLEV+1) = 0
      IF (JLEV  ==  KLEV) THEN
        Z_FACCLD1(JLON,JLEV+1) = 0.0_JPRB
        Z_FACCLD2(JLON,JLEV+1) = 0.0_JPRB
        Z_FACCLR1(JLON,JLEV+1) = 0.0_JPRB
        Z_FACCLR2(JLON,JLEV+1) = 0.0_JPRB
        !-- DS_000515      
        !      FACCMB1(LEV+1) = _ZERO_
      !      FACCMB2(LEV+1) = _ZERO_
        !mji      ISTCLD(JLON,LEV+1) = _ZERO_
      ELSEIF (P_CLDFRAC(JLON,JLEV+1)  >=  P_CLDFRAC(JLON,JLEV)) THEN
        Z_FACCLD1(JLON,JLEV+1) = 0.0_JPRB
        Z_FACCLD2(JLON,JLEV+1) = 0.0_JPRB
        IF (ISTCLD(JLON,JLEV)  ==  1) THEN
          !mji        ISTCLD(JLON,LEV+1) = 0
          Z_FACCLR1(JLON,JLEV+1) = 0.0_JPRB
          !mji        
          Z_FACCLR2(JLON,JLEV+1) = 0.0_JPRB
          IF (P_CLDFRAC(JLON,JLEV) < 1.0_JPRB) THEN
            Z_FACCLR2(JLON,JLEV+1) = (P_CLDFRAC(JLON,JLEV+1)-P_CLDFRAC(JLON,JLEV))/&
             & (1.0_JPRB-P_CLDFRAC(JLON,JLEV))  
          ENDIF   
        ELSE
          Z_FMAX(JLON) = MAX(P_CLDFRAC(JLON,JLEV),P_CLDFRAC(JLON,JLEV-1))
          !mji
          IF (P_CLDFRAC(JLON,JLEV+1)  >  Z_FMAX(JLON)) THEN
            Z_FACCLR1(JLON,JLEV+1) = Z_RAT2(JLON)
            Z_FACCLR2(JLON,JLEV+1) = (P_CLDFRAC(JLON,JLEV+1)-Z_FMAX(JLON))/(1.0_JPRB-Z_FMAX(JLON))
            !mji          
          ELSEIF (P_CLDFRAC(JLON,JLEV+1) < Z_FMAX(JLON)) THEN
            Z_FACCLR1(JLON,JLEV+1) = (P_CLDFRAC(JLON,JLEV+1)-P_CLDFRAC(JLON,JLEV))/&
             & (P_CLDFRAC(JLON,JLEV-1)-P_CLDFRAC(JLON,JLEV))  
            Z_FACCLR2(JLON,JLEV+1) = 0.0_JPRB
            !mji
          ELSE
            Z_FACCLR1(JLON,JLEV+1) = Z_RAT2(JLON)  
            Z_FACCLR2(JLON,JLEV+1) = 0.0_JPRB
          ENDIF
        ENDIF
        IF (Z_FACCLR1(JLON,JLEV+1) > 0.0_JPRB .OR. Z_FACCLR2(JLON,JLEV+1) > 0.0_JPRB) THEN
          Z_RAT1(JLON) = 1.0_JPRB
          Z_RAT2(JLON) = 0.0_JPRB
        ENDIF
      ELSE
        Z_FACCLR1(JLON,JLEV+1) = 0.0_JPRB
        Z_FACCLR2(JLON,JLEV+1) = 0.0_JPRB
        IF (ISTCLD(JLON,JLEV)  ==  1) THEN
          !mji        ISTCLD(JLON,LEV+1) = 0
          Z_FACCLD1(JLON,JLEV+1) = 0.0_JPRB
          Z_FACCLD2(JLON,JLEV+1) = (P_CLDFRAC(JLON,JLEV)-P_CLDFRAC(JLON,JLEV+1))/P_CLDFRAC(JLON,JLEV)
        ELSE
          Z_FMIN(JLON) = MIN(P_CLDFRAC(JLON,JLEV),P_CLDFRAC(JLON,JLEV-1))
          IF (P_CLDFRAC(JLON,JLEV+1)  <=  Z_FMIN(JLON)) THEN
            Z_FACCLD1(JLON,JLEV+1) = Z_RAT1(JLON)
            Z_FACCLD2(JLON,JLEV+1) = (Z_FMIN(JLON)-P_CLDFRAC(JLON,JLEV+1))/Z_FMIN(JLON)
          ELSE
            Z_FACCLD1(JLON,JLEV+1) = (P_CLDFRAC(JLON,JLEV)-P_CLDFRAC(JLON,JLEV+1))/&
             & (P_CLDFRAC(JLON,JLEV)-Z_FMIN(JLON))  
            Z_FACCLD2(JLON,JLEV+1) = 0.0_JPRB
          ENDIF
        ENDIF
        IF (Z_FACCLD1(JLON,JLEV+1) > 0.0_JPRB .OR. Z_FACCLD2(JLON,JLEV+1) > 0.0_JPRB) THEN
          Z_RAT1(JLON) = 0.0_JPRB
          Z_RAT2(JLON) = 1.0_JPRB
        ENDIF
      ENDIF
      !fcc
      IF (JLEV == 1) THEN
        Z_FACCMB1(JLON,JLEV+1) = 0.
        Z_FACCMB2(JLON,JLEV+1) = Z_FACCLD1(JLON,JLEV+1) * Z_FACCLR2(JLON,JLEV)
      ELSE
        Z_FACCMB1(JLON,JLEV+1) = Z_FACCLR1(JLON,JLEV+1) * Z_FACCLD2(JLON,JLEV) *P_CLDFRAC(JLON,JLEV-1)
        Z_FACCMB2(JLON,JLEV+1) = Z_FACCLD1(JLON,JLEV+1) * Z_FACCLR2(JLON,JLEV) *&
         & (1.0_JPRB - P_CLDFRAC(JLON,JLEV-1))   
      ENDIF
      !end fcc
    ELSE
      !-- DS_000515
      ISTCLD(JLON,JLEV+1) = 1
    ENDIF
  ENDDO
ENDDO

!_start_jjm 991209
  Z_RAT1(:) = 0.0_JPRB
  Z_RAT2(:) = 0.0_JPRB
!_end_jjm 991209


!-- DS_000515

DO JLEV = KLEV, 2, -1
  DO JLON = KIDIA, KFDIA
    IF (K_ICLDLYR(JLON,JLEV) == 1) THEN
      !mji
      ISTCLDD(JLON,JLEV-1) = 0  
        !mji      ISTCLDD(JLON,LEV-1) = _ZERO_
      IF (P_CLDFRAC(JLON,JLEV-1)  >=  P_CLDFRAC(JLON,JLEV)) THEN
        Z_FACCLD1D(JLON,JLEV-1) = 0.0_JPRB
        Z_FACCLD2D(JLON,JLEV-1) = 0.0_JPRB
        IF (ISTCLDD(JLON,JLEV)  ==  1) THEN
          !mji        ISTCLDD(JLON,LEV-1) = 0
          Z_FACCLR1D(JLON,JLEV-1) = 0.0_JPRB
          Z_FACCLR2D(JLON,JLEV-1) = 0.0_JPRB
          IF (P_CLDFRAC(JLON,JLEV) < 1.0_JPRB) THEN
            Z_FACCLR2D(JLON,JLEV-1) = (P_CLDFRAC(JLON,JLEV-1)-P_CLDFRAC(JLON,JLEV))/&
             & (1.0_JPRB-P_CLDFRAC(JLON,JLEV))  
          ENDIF
        ELSE
          Z_FMAX(JLON) = MAX(P_CLDFRAC(JLON,JLEV),P_CLDFRAC(JLON,JLEV+1))
          !mji
          IF (P_CLDFRAC(JLON,JLEV-1)  >  Z_FMAX(JLON)) THEN
            Z_FACCLR1D(JLON,JLEV-1) = Z_RAT2(JLON)
            Z_FACCLR2D(JLON,JLEV-1) = (P_CLDFRAC(JLON,JLEV-1)-Z_FMAX(JLON))/(1.0_JPRB-Z_FMAX(JLON))
            !mji
          ELSEIF (P_CLDFRAC(JLON,JLEV-1) < Z_FMAX(JLON)) THEN
            Z_FACCLR1D(JLON,JLEV-1) = (P_CLDFRAC(JLON,JLEV-1)-P_CLDFRAC(JLON,JLEV))/&
             & (P_CLDFRAC(JLON,JLEV+1)-P_CLDFRAC(JLON,JLEV))  
            Z_FACCLR2D(JLON,JLEV-1) = 0.0_JPRB
            !mji
          ELSE          
            Z_FACCLR1D(JLON,JLEV-1) = Z_RAT2(JLON)
            Z_FACCLR2D(JLON,JLEV-1) = 0.0_JPRB
          ENDIF
        ENDIF
        IF (Z_FACCLR1D(JLON,JLEV-1) > 0.0_JPRB .OR. Z_FACCLR2D(JLON,JLEV-1) > 0.0_JPRB)THEN
          Z_RAT1(JLON) = 1.0_JPRB
          Z_RAT2(JLON) = 0.0_JPRB
        ENDIF
      ELSE
        Z_FACCLR1D(JLON,JLEV-1) = 0.0_JPRB
        Z_FACCLR2D(JLON,JLEV-1) = 0.0_JPRB
        IF (ISTCLDD(JLON,JLEV)  ==  1) THEN
          !mji        ISTCLDD(JLON,LEV-1) = 0
          Z_FACCLD1D(JLON,JLEV-1) = 0.0_JPRB
          Z_FACCLD2D(JLON,JLEV-1) = (P_CLDFRAC(JLON,JLEV)-P_CLDFRAC(JLON,JLEV-1))/P_CLDFRAC(JLON,JLEV)
        ELSE
          Z_FMIN(JLON) = MIN(P_CLDFRAC(JLON,JLEV),P_CLDFRAC(JLON,JLEV+1))
          IF (P_CLDFRAC(JLON,JLEV-1)  <=  Z_FMIN(JLON)) THEN
            Z_FACCLD1D(JLON,JLEV-1) = Z_RAT1(JLON)
            Z_FACCLD2D(JLON,JLEV-1) = (Z_FMIN(JLON)-P_CLDFRAC(JLON,JLEV-1))/Z_FMIN(JLON)
          ELSE
            Z_FACCLD1D(JLON,JLEV-1) = (P_CLDFRAC(JLON,JLEV)-P_CLDFRAC(JLON,JLEV-1))/&
             & (P_CLDFRAC(JLON,JLEV)-Z_FMIN(JLON))  
            Z_FACCLD2D(JLON,JLEV-1) = 0.0_JPRB
          ENDIF
        ENDIF
        IF (Z_FACCLD1D(JLON,JLEV-1) > 0.0_JPRB .OR. Z_FACCLD2D(JLON,JLEV-1) > 0.0_JPRB)THEN
          Z_RAT1(JLON) = 0.0_JPRB
          Z_RAT2(JLON) = 1.0_JPRB
        ENDIF
      ENDIF
      IF(JLEV/=KLEV) THEN
        Z_FACCMB1D(JLON,JLEV-1) = Z_FACCLR1D(JLON,JLEV-1) * Z_FACCLD2D(JLON,JLEV) *P_CLDFRAC(JLON,JLEV+1)
        Z_FACCMB2D(JLON,JLEV-1) = Z_FACCLD1D(JLON,JLEV-1) * Z_FACCLR2D(JLON,JLEV) *&
       & (1.0_JPRB - P_CLDFRAC(JLON,JLEV+1))  
      ENDIF
    ELSE
      ISTCLDD(JLON,JLEV-1) = 1
    ENDIF
  ENDDO
ENDDO
!-----------
ILEV = 1
DO JLON = KIDIA, KFDIA
  IF (K_ICLDLYR(JLON,ILEV) == 1) THEN
    !mji
    ISTCLDD(JLON,ILEV-1) = 0  
    Z_FACCLD1D(JLON,ILEV-1) = 0.0_JPRB
    Z_FACCLD2D(JLON,ILEV-1) = 0.0_JPRB
    Z_FACCLR1D(JLON,ILEV-1) = 0.0_JPRB
    Z_FACCLR2D(JLON,ILEV-1) = 0.0_JPRB
    Z_FACCMB1D(JLON,ILEV-1) = 0.0_JPRB
    Z_FACCMB2D(JLON,ILEV-1) = 0.0_JPRB
    !mji      ISTCLDD(JLON,LEV-1) = _ZERO_
    Z_FACCMB1D(JLON,ILEV-1) = Z_FACCLR1D(JLON,ILEV-1) * Z_FACCLD2D(JLON,ILEV) *P_CLDFRAC(JLON,ILEV+1)
    Z_FACCMB2D(JLON,ILEV-1) = Z_FACCLD1D(JLON,ILEV-1) * Z_FACCLR2D(JLON,ILEV) *&
     & (1.0_JPRB - P_CLDFRAC(JLON,ILEV+1))  
  ELSE
    ISTCLDD(JLON,ILEV-1) = 1
  ENDIF
ENDDO
!-----------


!- Loop over frequency bands.

DO JBAND = K_ISTART, K_IEND
  DO JLON = KIDIA, KFDIA
    Z_DBDTLEV = TOTPLNK(INDBOUND(JLON)+1,JBAND)-TOTPLNK(INDBOUND(JLON),JBAND)
    Z_PLANKBND = DELWAVE(JBAND) * (TOTPLNK(INDBOUND(JLON),JBAND) + Z_TBNDFRAC(JLON) * Z_DBDTLEV)
    Z_DBDTLEV = TOTPLNK(INDLEV(JLON,0)+1,JBAND) -TOTPLNK(INDLEV(JLON,0),JBAND)
    !-- DS_000515
    Z_PLVL(JLON,JBAND,0) = DELWAVE(JBAND)&
     & * (TOTPLNK(INDLEV(JLON,0),JBAND) + Z_TLEVFRAC(JLON,0)*Z_DBDTLEV)  

    Z_SURFEMIS(JLON,JBAND) = P_SEMISS(JLON,JBAND)
    Z_PLNKEMIT(JLON,JBAND) = Z_SURFEMIS(JLON,JBAND) * Z_PLANKBND
    Z_SUMPLEM(JLON)  = Z_SUMPLEM(JLON) + Z_PLNKEMIT(JLON,JBAND)
    Z_SUMPL(JLON)    = Z_SUMPL(JLON)   + Z_PLANKBND
    !--DS
  ENDDO
ENDDO
!---

!-- DS_000515
DO JLEV = 1, KLEV

!cdir outerunroll=8
DO JBAND = K_ISTART, K_IEND
!cdir gthreorder
!cdir on_adb(INDLEV)
!cdir on_adb(INDLAY)
!cdir on_adb(Z_TLAYFRAC)
!cdir on_adb(Z_TLEVFRAC)
!cdir       nodep
    DO JLON = KIDIA, KFDIA
      !----              
      !- Calculate the integrated Planck functions for at the
      !  level and layer temperatures.
      !  Compute cloud transmittance for cloudy layers.
      Z_DBDTLEV = TOTPLNK(INDLEV(JLON,JLEV)+1,JBAND) - TOTPLNK(INDLEV(JLON,JLEV),JBAND)
      Z_DBDTLAY = TOTPLNK(INDLAY(JLON,JLEV)+1,JBAND) - TOTPLNK(INDLAY(JLON,JLEV),JBAND)
      !-- DS_000515
      Z_PLAY(JLON,JBAND,JLEV) = DELWAVE(JBAND)&
       & *(TOTPLNK(INDLAY(JLON,JLEV),JBAND)+Z_TLAYFRAC(JLON,JLEV)*Z_DBDTLAY)  
      Z_PLVL(JLON,JBAND,JLEV) = DELWAVE(JBAND)&
       & *(TOTPLNK(INDLEV(JLON,JLEV),JBAND)+Z_TLEVFRAC(JLON,JLEV)*Z_DBDTLEV)  
      !-- DS_000515
    ENDDO

  ENDDO

!cdir outerunroll=8
DO JBAND = K_ISTART, K_IEND
!cdir on_adb(K_ICLDLYR)
     DO JLON = KIDIA, KFDIA
      IF (K_ICLDLYR(JLON,JLEV) > 0) THEN
        ZEXTAU = MIN( P_TAUCLD(JLON,JLEV,JBAND), 200._JPRB)
        Z_TRNCLD(JLON,JLEV,JBAND) = EXP( -ZEXTAU )
      ENDIF
     ENDDO

  ENDDO


ENDDO



DO JLON = KIDIA, KFDIA
  P_SEMISLW(JLON) = Z_SUMPLEM(JLON) / Z_SUMPL(JLON)
ENDDO

!--DS
!O JI = 1, JPGPT
! NBI = NGB(JI)
! DO LEV =  1 , KLEV
!-- DS_000515
!   ZPLAY(JI,LEV) = PLAY(LEV,NGB(JI))
!   ZPLVL(JI,LEV) = PLVL(LEV-1,NGB(JI))
!   ZTAUCLD(JI,LEV) = TAUCLD(LEV,NGB(JI))
!   ZTRNCLD(JI,LEV) = TRNCLD(LEV,NGB(JI))
!-- DS_000515
! ENDDO
!NDDO
!----    


!- For cloudy layers, set cloud parameters for radiative transfer.
!DO JLON = KIDIA, KFDIA
!  DO JLEV = 1, KLEV
!    IF (K_ICLDLYR(JLON,JLEV) > 0) THEN
!      DO JI = 1, JPGPT
!        !--DS          
!        !            NBI = NGB(JLON,JI)
!        Z_ODCLDNW (JI,JLEV) =            P_TAUCLD(JLON,JLEV,NGB(JI))
!        Z_ABSCLDNW(JI,JLEV) = 1.0_JPRB - Z_TRNCLD(JLON,JLEV,NGB(JI))
!        !----            
!        !            EFCLFRNW(JI,LEV) = ABSCLDNW(JI,LEV) * CLDFRAC(LEV)
!      ENDDO
!    ENDIF
!  ENDDO

  !- Initialize for radiative transfer.
    Z_RADCLRD1(:,:) = 0.0_JPRB
    Z_RADLD1(:,:)   = 0.0_JPRB
  DO JI = 1, JPGPT
  DO JLON = KIDIA, KFDIA
    I_NBI = NGB(JI)
    Z_SEMIS(JLON,JI) = Z_SURFEMIS(JLON,I_NBI)
    Z_RADUEMIT(JLON,JI) = PFRAC(JLON,JI,1) * Z_PLNKEMIT(JLON,I_NBI)
    Z_BGLEV(JLON,JI) = PFRAC(JLON,JI,KLEV) * Z_PLVL(JLON,I_NBI,KLEV)
  ENDDO
  ENDDO

  !- Downward radiative transfer.
  !  *** DRAD1 holds summed radiance for total sky stream
  !  *** DRADCL1 holds summed radiance for clear sky stream

!cdir on_adb(ICLDDN_)
  ICLDDN_ = 0
  DO JLEV = KLEV, 1, -1
!cdir on_adb(Z_DRAD1__)
!cdir on_adb(Z_DRADCL1__)
    Z_DRAD1__(:)  = 0.0_JPRB
    Z_DRADCL1__(:) = 0.0_JPRB
    IENT = JPGPT * (JLEV-1)
!cdir on_adb(K_ICLDLYR)
!cdir on_adb(ICLDDN_)
  DO JLON = KIDIA, KFDIA
    IF (K_ICLDLYR(JLON,JLEV) == 1) THEN
      ICLDDN_(JLON) = 1
    ENDIF
  ENDDO

!cdir outerunroll=4
  DO JI = 1, JPGPT
    INDX = IENT + JI
!cdir on_adb(K_ICLDLYR)
!cdir on_adb(ISTCLD)
!cdir on_adb(ICLDDN_)
!cdir on_adb(Z_DRAD1__)
!cdir on_adb(Z_DRADCL1__)
!cdir on_adb(P_CLDFRAC)
!cdir on_adb(Z_FACCLR1D)
!cdir on_adb(Z_FACCLD1D)
!cdir on_adb(Z_FACCMB1D)
!cdir on_adb(Z_FACCMB2D)
!cdir on_adb(Z_FACCLR2D)
!cdir on_adb(Z_FACCLD2D)
!DEC$ VECTOR ALWAYS
      DO JLON = KIDIA, KFDIA

        !--DS            
        !            NBI = NGB(JI)
        Z_BGLAY = PFRAC(JLON,JI,JLEV) * Z_PLAY(JLON,NGB(JI),JLEV)
        !----            
        Z_DELBGUP     = Z_BGLEV(JLON,JI) - Z_BGLAY
        Z_BBU1(JLON,INDX) = Z_BGLAY + P_TAUSF1(JLON,JI,JLEV) * Z_DELBGUP
        !--DS            
        Z_BGLEV(JLON,JI) = PFRAC(JLON,JI,JLEV) * Z_PLVL(JLON,NGB(JI),JLEV-1)
        !----            
        Z_DELBGDN = Z_BGLEV(JLON,JI) - Z_BGLAY
        Z_BBD = Z_BGLAY + P_TAUSF1(JLON,JI,JLEV) * Z_DELBGDN

        IF (K_ICLDLYR(JLON,JLEV) == 1) THEN

      !  *** Cloudy layer
        IF (ISTCLDD(JLON,JLEV)  ==  1) THEN
        !***
          Z_CLDRADD(JLON,JI) = P_CLDFRAC(JLON,JLEV) * Z_RADLD1(JLON,JI)
          Z_CLRRADD(JLON,JI) = Z_RADLD1(JLON,JI) - Z_CLDRADD(JLON,JI)
          Z_OLDCLD(JLON,JI) = Z_CLDRADD(JLON,JI)
          Z_OLDCLR(JLON,JI) = Z_CLRRADD(JLON,JI)
          Z_RAD(JLON,JI) = 0.0_JPRB
        ENDIF

        !- total-sky downward flux          
        Z_ODSM = P_OD(JLON,JI,JLEV) + P_TAUCLD(JLON,JLEV,NGB(JI))
        Z_FACTOT1 = Z_ODSM / (BPADE + Z_ODSM)
        Z_BBUTOT1(JLON,INDX) = Z_BGLAY + Z_FACTOT1 * Z_DELBGUP
        Z_ATOT1(JLON,INDX) = P_ABSS1(JLON,JI,JLEV) + (1.0_JPRB - Z_TRNCLD(JLON,JLEV,NGB(JI)))&
         & - P_ABSS1(JLON,JI,JLEV) * (1.0_JPRB - Z_TRNCLD(JLON,JLEV,NGB(JI)))  
        Z_BBDTOT = Z_BGLAY + Z_FACTOT1 * Z_DELBGDN
        Z_GASSRC = Z_BBD * P_ABSS1(JLON,JI,JLEV)
        Z_TTOT = 1.0_JPRB - Z_ATOT1(JLON,INDX)
        Z_CLDSRC = Z_BBDTOT * Z_ATOT1(JLON,INDX)
      
        ! Separate RT equations for clear and cloudy streams      
        Z_CLDRADD(JLON,JI) = Z_CLDRADD(JLON,JI) * Z_TTOT + P_CLDFRAC(JLON,JLEV) * Z_CLDSRC
        Z_CLRRADD(JLON,JI) = Z_CLRRADD(JLON,JI) * (1.0_JPRB-P_ABSS1(JLON,JI,JLEV)) +&
         & (1.0_JPRB - P_CLDFRAC(JLON,JLEV)) * Z_GASSRC  

        !  Total sky downward radiance
        Z_RADLD1(JLON,JI) = Z_CLDRADD(JLON,JI) + Z_CLRRADD(JLON,JI)
        Z_DRAD1__(JLON) = Z_DRAD1__(JLON) + Z_RADLD1(JLON,JI)
      
        !  Clear-sky downward radiance          
        Z_RADCLRD1(JLON,JI) = Z_RADCLRD1(JLON,JI)+(Z_BBD-Z_RADCLRD1(JLON,JI))*P_ABSS1(JLON,JI,JLEV)
        Z_DRADCL1__(JLON) = Z_DRADCL1__(JLON) + Z_RADCLRD1(JLON,JI)

        !* Code to account for maximum/random overlap:
        !   Performs RT on the radiance most recently switched between clear and
        !   cloudy streams
        Z_RADMOD = Z_RAD(JLON,JI) * (Z_FACCLR1D(JLON,JLEV-1) * (1.0_JPRB-P_ABSS1(JLON,JI,JLEV)) +&
         & Z_FACCLD1D(JLON,JLEV-1) *  Z_TTOT) - &
         & Z_FACCMB1D(JLON,JLEV-1) * Z_GASSRC + &
         & Z_FACCMB2D(JLON,JLEV-1) * Z_CLDSRC  
       
        !   Computes what the clear and cloudy streams would have been had no
        !   radiance been switched       
        Z_OLDCLD(JLON,JI) = Z_CLDRADD(JLON,JI) - Z_RADMOD
        Z_OLDCLR(JLON,JI) = Z_CLRRADD(JLON,JI) + Z_RADMOD
      
        !   Computes the radiance to be switched between clear and cloudy.      
        Z_RAD(JLON,JI) = -Z_RADMOD + Z_FACCLR2D(JLON,JLEV-1)*Z_OLDCLR(JLON,JI) -&
         & Z_FACCLD2D(JLON,JLEV-1)*Z_OLDCLD(JLON,JI)  
        Z_CLDRADD(JLON,JI) = Z_CLDRADD(JLON,JI) + Z_RAD(JLON,JI)
        Z_CLRRADD(JLON,JI) = Z_CLRRADD(JLON,JI) - Z_RAD(JLON,JI)
        !***

    ELSE

      !  *** Clear layer
      !  *** DRAD1 holds summed radiance for total sky stream
      !  *** DRADCL1 holds summed radiance for clear sky stream

          !- total-sky downward radiance
          Z_RADLD1(JLON,JI) = Z_RADLD1(JLON,JI)+(Z_BBD-Z_RADLD1(JLON,JI))*P_ABSS1(JLON,JI,JLEV)
          Z_DRAD1__(JLON) = Z_DRAD1__(JLON) + Z_RADLD1(JLON,JI)

          IF (ICLDDN_(JLON) == 1) THEN
        
          !- clear-sky downward radiance
          !-  Set clear sky stream to total sky stream as long as layers
          !-  remain clear.  Streams diverge when a cloud is reached.
          Z_RADCLRD1(JLON,JI) = Z_RADCLRD1(JLON,JI)+(Z_BBD-Z_RADCLRD1(JLON,JI))*P_ABSS1(JLON,JI,JLEV)
          Z_DRADCL1__(JLON) = Z_DRADCL1__(JLON) + Z_RADCLRD1(JLON,JI)
            
          ELSE
        
          !- clear-sky downward flux          
          !-  Set clear sky stream to total sky stream as long as layers
          !-  remain clear.  Streams diverge when a cloud is reached.
          Z_RADCLRD1(JLON,JI) = Z_RADLD1(JLON,JI)
          ENDIF
    ENDIF
   ENDDO
  ENDDO

!cdir on_adb(K_ICLDLYR)
!cdir on_adb(Z_DRAD1__)
!cdir on_adb(Z_DRADCL1__)
  DO JLON = KIDIA, KFDIA
   IF (K_ICLDLYR(JLON,JLEV) /= 1 .AND. ICLDDN_(JLON) /= 1) THEN
        Z_DRADCL1__(JLON) = Z_DRAD1__(JLON)
   ENDIF
  ENDDO
    
!cdir on_adb(Z_DRAD1__)
!cdir on_adb(Z_DRADCL1__)
  DO JLON = KIDIA, KFDIA

    P_TOTDFLUC(JLON,JLEV-1) = Z_DRADCL1__(JLON) * Z_WTNUM
    P_TOTDFLUX(JLON,JLEV-1) = Z_DRAD1__(JLON)   * Z_WTNUM

  ENDDO
ENDDO


! l'inversion de ces 2 boucles amene de legeres diff numeriques (sommation)

  Z_URAD1__(:) = 0.0_JPRB
  Z_URADCL1__(:) = 0.0_JPRB
!cdir novector
DO JLON = KIDIA, KFDIA
  DO JI = 1, JPGPT
    ! Spectral reflectivity and reflectance
    ! Includes the contribution of spectrally varying longwave emissivity 
    ! and reflection from the surface to the upward radiative transfer.
    ! Note: Spectral and Lambertian reflections are identical for the one
    ! angle flux integration used here.

    !- Lambertian reflection.
    ! Clear-sky radiance
    !    RADCLD = _TWO_ * (RADCLRD1(JI) * WTNUM(1) )
    Z_RADCLD = Z_RADCLRD1(JLON,JI)
    Z_RADCLRU1(JLON,JI) = Z_RADUEMIT(JLON,JI) + (1.0_JPRB - Z_SEMIS(JLON,JI)) * Z_RADCLD
    Z_URADCL1__(JLON) = Z_URADCL1__(JLON) + Z_RADCLRU1(JLON,JI)

    ! Total sky radiance
    !    RADD = _TWO_ * (RADLD1(JI) * WTNUM(1) )
    Z_RADD = Z_RADLD1(JLON,JI)
    Z_RADLU1(JLON,JI) = Z_RADUEMIT(JLON,JI) + (1.0_JPRB - Z_SEMIS(JLON,JI)) * Z_RADD
    Z_URAD1__(JLON) = Z_URAD1__(JLON)+ Z_RADLU1(JLON,JI)
  ENDDO
ENDDO
DO JLON = KIDIA, KFDIA
  P_TOTUFLUC(JLON,0) = Z_WTNUM * Z_URADCL1__(JLON)
  P_TOTUFLUX(JLON,0) = Z_WTNUM * Z_URAD1__(JLON)
ENDDO



  !ELSE
  !!- Specular reflection.
  !  DO JI = 1, JPGPT
  !    RADCLU = RADUEMIT(JI)
  !    RADCLRU1(JI) = RADCLU + (_ONE_ - SEMIS(JI)) * RADCLRD1(JI)
  !    URADCL1 = URADCL1 + RADCLRU1(JI)

  !    RADU = RADUEMIT(JI)
  !    RADLU1(JI) = RADU + (_ONE_ - SEMIS(JI)) * RADLD1(JI)
  !    URAD1 = URAD1 + RADLU1(JI)
  !  ENDDO
  !  TOTUFLUC(0) = URADCL1 * WTNUM(1)
  !  TOTUFLUX(0) = URAD1   * WTNUM(1)
  !ENDIF

  !- Upward radiative transfer.
  !- *** URAD1 holds the summed radiance for total sky stream
  !- *** URADCL1 holds the summed radiance for clear sky stream
  
  
  !- Upward radiative transfer.
  !- *** URAD1 holds the summed radiance for total sky stream
  !- *** URADCL1 holds the summed radiance for clear sky stream
  DO JLEV = 1, KLEV
!cdir on_adb(Z_URAD1__)
!cdir on_adb(Z_URADCL1__)
     Z_URAD1__(:) = 0.0_JPRB
     Z_URADCL1__(:) = 0.0_JPRB
    IENT = JPGPT * (JLEV-1)
!cdir outerunroll=2
    DO JI = 1, JPGPT
      INDX = IENT + JI
!cdir on_adb(K_ICLDLYR)
!cdir on_adb(ISTCLD)
!cdir on_adb(P_CLDFRAC)
!cdir on_adb(Z_URAD1__)
!cdir on_adb(Z_URADCL1__)
!cdir on_adb(Z_FACCLR1D)
!cdir on_adb(Z_FACCLD1D)
!cdir on_adb(Z_FACCMB1D)
!cdir on_adb(Z_FACCMB2D)
!cdir on_adb(Z_FACCLR2D)
!cdir on_adb(Z_FACCLD2D)
!DEC$ VECTOR ALWAYS
     DO JLON = KIDIA, KFDIA
      ! Check flag for cloud in current layer
      IF (K_ICLDLYR(JLON,JLEV) == 1) THEN

        IF (ISTCLD(JLON,JLEV)  ==  1) THEN
          Z_CLDRADU(JLON,JI) = P_CLDFRAC(JLON,JLEV) * Z_RADLU1(JLON,JI)
          Z_CLRRADU(JLON,JI) = Z_RADLU1(JLON,JI) - Z_CLDRADU(JLON,JI)
          Z_OLDCLD(JLON,JI) = Z_CLDRADU(JLON,JI)
          Z_OLDCLR(JLON,JI) = Z_CLRRADU(JLON,JI)
          Z_RAD(JLON,JI) = 0.0_JPRB
        ENDIF

      !- *** Cloudy layer
        !- total-sky upward flux          
        Z_GASSRC = Z_BBU1(JLON,INDX) * P_ABSS1(JLON,JI,JLEV)

        !- If first cloudy layer in sequence, split up radiance into clear and
        !    cloudy streams depending on cloud fraction
        Z_TTOT = 1.0_JPRB - Z_ATOT1(JLON,INDX)
        Z_TRNS = 1.0_JPRB - P_ABSS1(JLON,JI,JLEV)
        Z_CLDSRC = Z_BBUTOT1(JLON,INDX) * Z_ATOT1(JLON,INDX)

        !- Separate RT equations for clear and cloudy streams      
        Z_CLDRADU(JLON,JI) = Z_CLDRADU(JLON,JI) * Z_TTOT + P_CLDFRAC(JLON,JLEV) * Z_CLDSRC
        Z_CLRRADU(JLON,JI) = Z_CLRRADU(JLON,JI) * Z_TRNS +(1.0_JPRB - P_CLDFRAC(JLON,JLEV)) * Z_GASSRC
        !***

        !- total sky upward flux
        Z_RADLU1(JLON,JI) = Z_CLDRADU(JLON,JI) + Z_CLRRADU(JLON,JI)
        Z_URAD1__(JLON) = Z_URAD1__(JLON) + Z_RADLU1(JLON,JI)
      
        !- clear-sky upward flux
        Z_RADCLRU1(JLON,JI) = Z_RADCLRU1(JLON,JI) + (Z_BBU1(JLON,INDX)-Z_RADCLRU1(JLON,JI))&
         & *P_ABSS1(JLON,JI,JLEV)  
        Z_URADCL1__(JLON) = Z_URADCL1__(JLON) + Z_RADCLRU1(JLON,JI)

        !* Code to account for maximum/random overlap:
        !   Performs RT on the radiance most recently switched between clear and
        !   cloudy streams
        Z_RADMOD = Z_RAD(JLON,JI) * (Z_FACCLR1(JLON,JLEV+1) * Z_TRNS +&
         & Z_FACCLD1(JLON,JLEV+1) *  Z_TTOT) - &
         & Z_FACCMB1(JLON,JLEV+1) * Z_GASSRC + &
         & Z_FACCMB2(JLON,JLEV+1) * Z_CLDSRC  
       
        !   Computes what the clear and cloudy streams would have been had no
        !   radiance been switched       
        Z_OLDCLD(JLON,JI) = Z_CLDRADU(JLON,JI) - Z_RADMOD
        Z_OLDCLR(JLON,JI) = Z_CLRRADU(JLON,JI) + Z_RADMOD
      
        !   Computes the radiance to be switched between clear and cloudy.      
        Z_RAD(JLON,JI) = -Z_RADMOD + Z_FACCLR2(JLON,JLEV+1)*Z_OLDCLR(JLON,JI) -&
         & Z_FACCLD2(JLON,JLEV+1)*Z_OLDCLD(JLON,JI)  
        Z_CLDRADU(JLON,JI) = Z_CLDRADU(JLON,JI) + Z_RAD(JLON,JI)
        Z_CLRRADU(JLON,JI) = Z_CLRRADU(JLON,JI) - Z_RAD(JLON,JI)
        !***

      ELSE

        !- *** Clear layer
        !- total-sky upward flux          
        Z_RADLU1(JLON,JI) = &
       &Z_RADLU1(JLON,JI)+(Z_BBU1(JLON,INDX)-Z_RADLU1(JLON,JI))*P_ABSS1(JLON,JI,JLEV)
        Z_URAD1__(JLON) = &
       &Z_URAD1__(JLON) + Z_RADLU1(JLON,JI)
        !- clear-sky upward flux
        !   Upward clear and total sky streams must be separate because surface
        !   reflectance is different for each.
        Z_RADCLRU1(JLON,JI) = &
       &Z_RADCLRU1(JLON,JI)+(Z_BBU1(JLON,INDX)-Z_RADCLRU1(JLON,JI))*P_ABSS1(JLON,JI,JLEV)
        Z_URADCL1__(JLON) = &
       &Z_URADCL1__(JLON) + Z_RADCLRU1(JLON,JI)
      ENDIF

    ENDDO

  ENDDO
!cdir on_adb(Z_URAD1__)
!cdir on_adb(Z_URADCL1__)
  DO JLON = KIDIA, KFDIA
    P_TOTUFLUC(JLON,JLEV) = Z_WTNUM * Z_URADCL1__(JLON)
    P_TOTUFLUX(JLON,JLEV) = Z_WTNUM * Z_URAD1__(JLON)
  ENDDO

ENDDO

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('RRTM_RTRN1A_140GP',1,ZHOOK_HANDLE)
END SUBROUTINE RRTM_RTRN1A_140GP
