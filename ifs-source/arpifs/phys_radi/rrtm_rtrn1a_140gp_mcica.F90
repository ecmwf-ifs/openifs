 SUBROUTINE RRTM_RTRN1A_140GP_MCICA&
 &(YDDIMV, KIDIA, KFDIA, KLEV, KSTART, KEND, KCOLS ,&
 &PCLDFRAC, PTAUCLD, PABSS1 ,&
 &POD , PTAUSF1 ,&
 &PTOTDFLUC,PTOTDFLUX,PTOTUFLUC,PTOTUFLUX,&
 &PTAVEL,PTZ,PTBOUND,PFRAC,PSEMISS,PSEMISLW,&
 &  PLWDERIVATIVE) 


!-* This program calculates the upward fluxes, downward fluxes,
!   and heating rates for an arbitrary atmosphere.  The input to
!   this program is the atmospheric profile and all Planck function
!   information.  First-order "numerical" quadrature is used for the 
!   angle integration, i.e. only one exponential is computed per layer
!   per g-value per band. 
!   Consistent with McICA, cloud fraction is only 0 or 1 and the routine has 
!   been cleaned of all calculations related to partial cloud cover.

!     AUTHOR.
!     -------
!      JJMorcrette, ECMWF

!     MODIFICATIONS.
!     --------------
!      McICA JJMorcrette       ECMWF, 20050110
!      NEC           25-Oct-2007 Optimisations
!      JJMorcrette 20110613  flexible number of g-points
!      RJHogan     20140519  Added profile of partial derivatives to output
!      RJHogan     20160914  Dr Hook profiling around Planck-function part
! ---------------------------------------------------------------------------

USE YOMDIMV  , ONLY : TDIMV
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE PARRRTM  , ONLY : JPBAND
USE YOERRTM  , ONLY : JPGPT
USE YOERRTAB , ONLY : BPADE
USE YOERRTWN , ONLY : TOTPLNK, DELWAVE
USE YOERRTFTR, ONLY : NGB
USE YOMJFH   , ONLY : N_VMASS
USE YOMCST   , ONLY : RSIGMA

IMPLICIT NONE

TYPE(TDIMV)       ,INTENT(IN)    :: YDDIMV
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA, KFDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTART 
INTEGER(KIND=JPIM),INTENT(IN)    :: KEND 
INTEGER(KIND=JPIM),INTENT(IN)    :: KCOLS
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCLDFRAC(KIDIA:KFDIA,KCOLS,YDDIMV%NFLEVG) ! Cloud fraction
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTAUCLD(KIDIA:KFDIA,YDDIMV%NFLEVG,KCOLS) ! Spectral optical thickness
REAL(KIND=JPRB)   ,INTENT(IN)    :: PABSS1(KIDIA:KFDIA,JPGPT*YDDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: POD(KIDIA:KFDIA,JPGPT,YDDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTAUSF1(KIDIA:KFDIA,JPGPT*YDDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTOTDFLUC(KIDIA:KFDIA,0:YDDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTOTDFLUX(KIDIA:KFDIA,0:YDDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTOTUFLUC(KIDIA:KFDIA,0:YDDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTOTUFLUX(KIDIA:KFDIA,0:YDDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTAVEL(KIDIA:KFDIA,YDDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTZ(KIDIA:KFDIA,0:YDDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTBOUND(KIDIA:KFDIA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFRAC(KIDIA:KFDIA,JPGPT,YDDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSEMISS(KIDIA:KFDIA,JPBAND) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSEMISLW(KIDIA:KFDIA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PLWDERIVATIVE(KIDIA:KFDIA,0:YDDIMV%NFLEVG) ! See below for definition
! ---------------------------------------------------------------------------

INTEGER(KIND=JPIM) :: INDLAY(KIDIA:KFDIA,YDDIMV%NFLEVG),INDLEV(KIDIA:KFDIA,0:YDDIMV%NFLEVG)

REAL(KIND=JPRB) :: ZBBU1(KIDIA:KFDIA,JPGPT*YDDIMV%NFLEVG),ZBBUTOT1(KIDIA:KFDIA,JPGPT*YDDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZTLAYFRAC(KIDIA:KFDIA,YDDIMV%NFLEVG),ZTLEVFRAC(KIDIA:KFDIA,0:YDDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZBGLEV(KIDIA:KFDIA,JPGPT)
!-- DS_000515
REAL(KIND=JPRB) :: ZZPLVL(KIDIA:KFDIA,JPBAND+1,0:YDDIMV%NFLEVG),ZZPLAY(KIDIA:KFDIA,JPBAND+1,0:YDDIMV%NFLEVG), &
 & ZWTNUM(KIDIA:KFDIA,3)
REAL(KIND=JPRB) :: ZPLVL(KIDIA:KFDIA,JPGPT+1,0:YDDIMV%NFLEVG)  ,ZPLAY(KIDIA:KFDIA,JPGPT+1,YDDIMV%NFLEVG)
!-- DS_000515
REAL(KIND=JPRB) :: ZODCLDNW(KIDIA:KFDIA,JPGPT,YDDIMV%NFLEVG), ZEFCLFRNW(KIDIA:KFDIA,JPGPT,YDDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZSEMIS(KIDIA:KFDIA,JPGPT),ZRADUEMIT(KIDIA:KFDIA,JPGPT)

REAL(KIND=JPRB) :: ZRADCLRU1(KIDIA:KFDIA,JPGPT) ,ZRADCLRD1(KIDIA:KFDIA,JPGPT)
REAL(KIND=JPRB) :: ZRADLU1(KIDIA:KFDIA,JPGPT)   ,ZRADLD1(KIDIA:KFDIA,JPGPT)

REAL(KIND=JPRB) :: ZABSCLDNW(KIDIA:KFDIA,JPGPT,YDDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZATOT1(KIDIA:KFDIA,JPGPT*YDDIMV%NFLEVG)

REAL(KIND=JPRB) :: ZSURFEMIS(KIDIA:KFDIA,JPBAND),ZPLNKEMIT(KIDIA:KFDIA,JPBAND)

INTEGER(KIND=JPIM) :: ISTCLD(KIDIA:KFDIA,YDDIMV%NFLEVG+1),ISTCLDD(KIDIA:KFDIA,0:YDDIMV%NFLEVG)
!******

!REAL_B :: ZPLVL(JPGPT+1,NFLEVG)  ,ZPLAY(JPGPT+1,NFLEVG)
!REAL_B :: ZTRNCLD(JPGPT+1,NFLEVG),ZTAUCLD(JPGPT+1,NFLEVG)

INTEGER(KIND=JPIM) :: JBAND, IENT, INDBOUND(KIDIA:KFDIA), INDEX, JI, JLAY, JLEV, INBI, JL, JLEN

REAL(KIND=JPRB) :: ZBBD, ZBBDTOT, ZBGLAY, ZDBDTLAY, ZDBDTLEV,&
 & ZDELBGDN, ZDELBGUP, ZDRAD1(KIDIA:KFDIA), ZDRADCL1(KIDIA:KFDIA), ZFACTOT1,&
 & ZGASSRC, ZODSM, ZPLANKBND, ZRADCLD(KIDIA:KFDIA), ZRADD, ZSUMPL(KIDIA:KFDIA),&
 & ZSUMPLEM(KIDIA:KFDIA), ZTBNDFRAC(KIDIA:KFDIA), ZURAD1(KIDIA:KFDIA), ZURADCL1(KIDIA:KFDIA), ZEXTAU  
! Partial derivative per g point
REAL(KIND=JPRB) :: ZLWDERIVATIVEGPOINT(KIDIA:KFDIA,KCOLS)
! Partial derivative of surface emitted flux in a band w.r.t. total flux
REAL(KIND=JPRB) :: ZDBANDFLUXUP_DTOTALFLUXUP(KIDIA:KFDIA,JPBAND)
REAL(KIND=JPRB) :: ZDTSKIN_DTOTALPLANCK(KIDIA:KFDIA) ! 4 * sigma * Tskin**3
REAL(KIND=JPRB) :: ZTMP1(JPGPT*(KFDIA-KIDIA+1)), ZTMP2(JPGPT*(KFDIA-KIDIA+1))
REAL(KIND=JPRB) :: ZTRANSMISSIVITY ! Combined transmissivity of gas and cloud
REAL(KIND=JPRB) :: ZFLUXFACTOR     ! factor to convert to flux
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE, ZHOOK_HANDLE_PLANCK

!--------------------------------------------------------------------------
! Input
!  NFLEVG                ! Maximum number of model layers
!  JPGPT                 ! Total number of g-point subintervals
!  JPBAND                ! Number of longwave spectral bands
!  SECANG                ! Diffusivity angle
!  WTNUM                 ! Weight for radiance to flux conversion
!  KLEV                  ! Number of model layers
!  PAVEL(NFLEVG)         ! Mid-layer pressures (hPa)
!  TAVEL(NFLEVG)         ! Mid-layer temperatures (K)
!  TZ(0:NFLEVG)          ! Interface temperatures (K)
!  TBOUND                ! Surface temperature
!  CLDFRAC(NFLEVG)       ! Layer cloud fraction
!  TAUCLD(NFLEVG,JPBAND) ! Layer cloud optical thickness
!  ITR
!  PFRAC(JPGPT,NFLEVG)   ! Planck function fractions
!  ICLDLYR(NFLEVG)       ! Flag for cloudy layers
!  ICLD                  ! Flag for cloudy column
!  IREFLECT              ! Flag for specular reflection
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
!  PLwDerivative(0:NFLEVG)! Partial derivative of total upward flux at
!                           each level with respect to upward flux at
!                           surface, used to correct heating rates at
!                           gridpoints/timesteps between calls to the
!                           full radiation scheme

!--------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('RRTM_RTRN1A_140GP_MCICA',0,ZHOOK_HANDLE)
ASSOCIATE(NFLEVG=>YDDIMV%NFLEVG)

IF (LHOOK) CALL DR_HOOK('RRTM_RTRN1A_140GP_MCICA:PLANCK',0,ZHOOK_HANDLE_PLANCK)
DO JL = KIDIA, KFDIA
  ZWTNUM(JL,1)=0.5_JPRB
  ZWTNUM(JL,2)=0.0_JPRB
  ZWTNUM(JL,3)=0.0_JPRB

!-start JJM_000511
  IF (PTBOUND(JL) < 339._JPRB .AND. PTBOUND(JL) >= 160._JPRB ) THEN
    INDBOUND(JL) = PTBOUND(JL) - 159._JPRB
    ZTBNDFRAC(JL) = PTBOUND(JL) - INT(PTBOUND(JL))
  ELSEIF (PTBOUND(JL) >= 339._JPRB ) THEN
    INDBOUND(JL) = 180
    ZTBNDFRAC(JL) = PTBOUND(JL) - 339._JPRB
  ELSEIF (PTBOUND(JL) < 160._JPRB ) THEN
    INDBOUND(JL) = 1
    ZTBNDFRAC(JL) = PTBOUND(JL) - 160._JPRB
  ENDIF  
ENDDO
!-end JJM_000511
  
DO JLAY = 0, KLEV
  DO JL = KIDIA, KFDIA
    PTOTUFLUC(JL,JLAY) = 0.0_JPRB
    PTOTDFLUC(JL,JLAY) = 0.0_JPRB
    PTOTUFLUX(JL,JLAY) = 0.0_JPRB
    PTOTDFLUX(JL,JLAY) = 0.0_JPRB
!-start JJM_000511
    IF (PTZ(JL,JLAY) < 339._JPRB .AND. PTZ(JL,JLAY) >= 160._JPRB ) THEN
      INDLEV(JL,JLAY) = PTZ(JL,JLAY) - 159._JPRB
      ZTLEVFRAC(JL,JLAY) = PTZ(JL,JLAY) - INT(PTZ(JL,JLAY))
    ELSEIF (PTZ(JL,JLAY) >= 339._JPRB ) THEN
      INDLEV(JL,JLAY) = 180
      ZTLEVFRAC(JL,JLAY) = PTZ(JL,JLAY) - 339._JPRB
    ELSEIF (PTZ(JL,JLAY) < 160._JPRB ) THEN
      INDLEV(JL,JLAY) = 1
      ZTLEVFRAC(JL,JLAY) = PTZ(JL,JLAY) - 160._JPRB
    ENDIF    
!-end JJM_000511
  ENDDO
ENDDO
!print 9901,ZTLEVFRAC(0),ZTLEVFRAC(KLEV)
9901 FORMAT(1X,'rrtm: ',10E12.5)

DO JL = KIDIA, KFDIA
  ZSUMPL(JL)   = 0.0_JPRB
  ZSUMPLEM(JL) = 0.0_JPRB

  ISTCLD(JL,1) = 1
  ISTCLDD(JL,KLEV) = 1
ENDDO

DO JLEV = 1, KLEV
  DO JL = KIDIA, KFDIA
!  -- DS_000515
!  -start JJM_000511
    IF (PTAVEL(JL,JLEV) < 339._JPRB .AND. PTAVEL(JL,JLEV) >= 160._JPRB ) THEN
      INDLAY(JL,JLEV) = PTAVEL(JL,JLEV) - 159._JPRB
      ZTLAYFRAC(JL,JLEV) = PTAVEL(JL,JLEV) - INT(PTAVEL(JL,JLEV))
    ELSEIF (PTAVEL(JL,JLEV) >= 339._JPRB ) THEN
      INDLAY(JL,JLEV) = 180
      ZTLAYFRAC(JL,JLEV) = PTAVEL(JL,JLEV) - 339._JPRB
    ELSEIF (PTAVEL(JL,JLEV) < 160._JPRB ) THEN
      INDLAY(JL,JLEV) = 1
      ZTLAYFRAC(JL,JLEV) = PTAVEL(JL,JLEV) - 160._JPRB
    ENDIF  
!  -end JJM_000511
  ENDDO
ENDDO
!-- DS_000515 

! Compute reciprocal of partial derivative of total Planck function
DO JL = KIDIA, KFDIA
   ZDTSKIN_DTOTALPLANCK(JL) = 1.0_JPRB / (4.0_JPRB * RSIGMA * PTBOUND(JL)**3)
ENDDO

! Factor to convert to fluxes
ZFLUXFACTOR = 0.5_JPRB * (2.0_JPRB * ASIN(1.0_JPRB)) * 2.0E4_JPRB

!- Loop over frequency bands.

DO JBAND = KSTART, KEND
  DO JL = KIDIA, KFDIA
    ZDBDTLEV = TOTPLNK(INDBOUND(JL)+1,JBAND)-TOTPLNK(INDBOUND(JL),JBAND)
    ZPLANKBND = DELWAVE(JBAND) * (TOTPLNK(INDBOUND(JL),JBAND) + ZTBNDFRAC(JL) * ZDBDTLEV)

    ZDBANDFLUXUP_DTOTALFLUXUP(JL,JBAND) = ZFLUXFACTOR &
         & * ZDTSKIN_DTOTALPLANCK(JL) * ZDBDTLEV * DELWAVE(JBAND)

    ZDBDTLEV = TOTPLNK(INDLEV(JL,0)+1,JBAND) -TOTPLNK(INDLEV(JL,0),JBAND)
!  -- DS_000515
    ZZPLVL(JL,JBAND,0) = DELWAVE(JBAND)&
     & * (TOTPLNK(INDLEV(JL,0),JBAND) + ZTLEVFRAC(JL,0)*ZDBDTLEV)  

    ZSURFEMIS(JL,JBAND) = PSEMISS(JL,JBAND)
    ZPLNKEMIT(JL,JBAND) = ZSURFEMIS(JL,JBAND) * ZPLANKBND
    ZSUMPLEM(JL)  = ZSUMPLEM(JL) + ZPLNKEMIT(JL,JBAND)
    ZSUMPL(JL)    = ZSUMPL(JL)   + ZPLANKBND
!  --DS

!    print 9880,JBAND,ZZPLVL(JBAND,0)
9880 FORMAT(1X,'ZZPLVL(JL,',I3,')=',E13.6)
  ENDDO
ENDDO
!---

!-- DS_000515
DO JLEV = 1, KLEV
  DO JBAND = KSTART, KEND
    DO JL = KIDIA, KFDIA
!  ----              
!  - Calculate the integrated Planck functions for at the
!    level and layer temperatures.
!    Compute cloud transmittance for cloudy layers.
      ZDBDTLEV = TOTPLNK(INDLEV(JL,JLEV)+1,JBAND) - TOTPLNK(INDLEV(JL,JLEV),JBAND)
      ZDBDTLAY = TOTPLNK(INDLAY(JL,JLEV)+1,JBAND) - TOTPLNK(INDLAY(JL,JLEV),JBAND)
!  -- DS_000515
      ZZPLAY(JL,JBAND,JLEV) = DELWAVE(JBAND)&
       & *(TOTPLNK(INDLAY(JL,JLEV),JBAND)+ZTLAYFRAC(JL,JLEV)*ZDBDTLAY)  
      ZZPLVL(JL,JBAND,JLEV) = DELWAVE(JBAND)&
       & *(TOTPLNK(INDLEV(JL,JLEV),JBAND)+ZTLEVFRAC(JL,JLEV)*ZDBDTLEV)  
!  -- DS_000515
    ENDDO
  ENDDO
ENDDO

IF (LHOOK) CALL DR_HOOK('RRTM_RTRN1A_140GP_MCICA:PLANCK',1,ZHOOK_HANDLE_PLANCK)

DO JL = KIDIA, KFDIA
  PSEMISLW(JL) = ZSUMPLEM(JL) / ZSUMPL(JL)
ENDDO

DO JLEV =  1 , KLEV
  DO JI = 1, JPGPT
    INBI=NGB(JI)
    DO JL = KIDIA, KFDIA
      ZPLAY(JL,JI,JLEV) = ZZPLAY(JL,INBI,JLEV)
!!!  ZPLVL(JI,JLEV) = ZZPLVL(INBI,JLEV-1)
      ZPLVL(JL,JI,JLEV) = ZZPLVL(JL,INBI,JLEV)
    ENDDO
  ENDDO
ENDDO

DO JI = 1, JPGPT
  INBI=NGB(JI)
  DO JL = KIDIA, KFDIA
    ZPLVL(JL,JI, 0  ) = ZZPLVL(JL,INBI, 0  )
  ENDDO
ENDDO

!----      

!- For cloudy layers, set cloud parameters for radiative transfer.
JLEN=KFDIA-KIDIA+1
DO JLEV = 1, KLEV
  IF( N_VMASS > 0 )THEN
    DO JI = 1, JPGPT
    DO JL = KIDIA, KFDIA
      ZTMP1(JL+(JI-1)*JLEN) = -MIN( PTAUCLD(JL,JLEV,JI), 200._JPRB )
    ENDDO
    ENDDO
    CALL VEXP(ZTMP2,ZTMP1,JPGPT*JLEN)
    DO JI = 1, JPGPT
    DO JL = KIDIA, KFDIA
      ZODCLDNW(JL,JI,JLEV) = PTAUCLD(JL,JLEV,JI)
      ZABSCLDNW(JL,JI,JLEV) = 1.0_JPRB - ZTMP2(JL+(JI-1)*JLEN)
      ZEFCLFRNW(JL,JI,JLEV) = ZABSCLDNW(JL,JI,JLEV) * PCLDFRAC(JL,JI,JLEV)
    ENDDO
    ENDDO
  ELSE
    DO JI = 1, JPGPT
    DO JL = KIDIA, KFDIA
      ZEXTAU = MIN( PTAUCLD(JL,JLEV,JI), 200._JPRB )          
      ZODCLDNW(JL,JI,JLEV) = PTAUCLD(JL,JLEV,JI)
      ZABSCLDNW(JL,JI,JLEV) = 1.0_JPRB - EXP( -ZEXTAU)    
      ZEFCLFRNW(JL,JI,JLEV) = ZABSCLDNW(JL,JI,JLEV) * PCLDFRAC(JL,JI,JLEV)
    ENDDO
    ENDDO
  ENDIF
ENDDO

!- Initialize for radiative transfer.
DO JI = 1, JPGPT
  INBI = NGB(JI)
  DO JL = KIDIA, KFDIA
    ZRADCLRD1(JL,JI) = 0.0_JPRB
    ZRADLD1(JL,JI)   = 0.0_JPRB
    ZSEMIS(JL,JI) = ZSURFEMIS(JL,INBI)
    ZRADUEMIT(JL,JI) = PFRAC(JL,JI,1) * ZPLNKEMIT(JL,INBI)
!!!  ZBGLEV(JI) = PFRAC(JI,KLEV) * ZPLVL(INBI,KLEV)
!!!  ZBGLEV(JI) = PFRAC(JI,KLEV) * ZZPLVL(INBI,KLEV)
    ZBGLEV(JL,JI) = PFRAC(JL,JI,KLEV) * ZPLVL(JL,JI,KLEV)

!  print 9890,JI,ZBGLEV(JI),PFRAC(JI,KLEV),ZPLVL(JI,KLEV)
9890 FORMAT(1X,'BGLEV',I4,3E14.7)
  ENDDO
ENDDO

!- Downward radiative transfer.
!  *** DRAD1 holds summed radiance for total sky stream
!  *** DRADCL1 holds summed radiance for clear sky stream

DO JLEV = KLEV, 1, -1
  IENT = JPGPT * (JLEV-1)
  DO JL = KIDIA, KFDIA
  ZDRAD1(JL)   = 0.0_JPRB
  ZDRADCL1(JL) = 0.0_JPRB

!  *** for a layer, either fully clear or overcast
  ENDDO
  DO JI = 1, JPGPT
    INDEX = IENT + JI
    DO JL = KIDIA, KFDIA
      ZBGLAY = PFRAC(JL,JI,JLEV) * ZPLAY(JL,JI,JLEV)
      ZDELBGUP     = ZBGLEV(JL,JI) - ZBGLAY
      ZBBU1(JL,INDEX) = ZBGLAY + PTAUSF1(JL,INDEX) * ZDELBGUP
      ZBGLEV(JL,JI) = PFRAC(JL,JI,JLEV) * ZPLVL(JL,JI,JLEV-1)
      ZDELBGDN = ZBGLEV(JL,JI) - ZBGLAY
      ZBBD = ZBGLAY + PTAUSF1(JL,INDEX) * ZDELBGDN

!    Clear-sky downward radiance          
      ZRADCLRD1(JL,JI) = ZRADCLRD1(JL,JI)+(ZBBD-ZRADCLRD1(JL,JI))*PABSS1(JL,INDEX)
      ZDRADCL1(JL) = ZDRADCL1(JL) + ZRADCLRD1(JL,JI)

!  - total-sky downward flux          
      ZODSM = POD(JL,JI,JLEV) + ZODCLDNW(JL,JI,JLEV)
      ZFACTOT1 = ZODSM / (BPADE + ZODSM)
      ZBBUTOT1(JL,INDEX) = ZBGLAY + ZFACTOT1 * ZDELBGUP
      ZATOT1(JL,INDEX) = PABSS1(JL,INDEX) + ZABSCLDNW(JL,JI,JLEV)&
        &           - PABSS1(JL,INDEX) * ZABSCLDNW(JL,JI,JLEV)  
      ZBBDTOT = ZBGLAY + ZFACTOT1 * ZDELBGDN
      ZGASSRC = ZBBD * PABSS1(JL,INDEX)

      ZRADLD1(JL,JI) = ZRADLD1(JL,JI) &
        &          - ZRADLD1(JL,JI)*(PABSS1(JL,INDEX)+ZEFCLFRNW(JL,JI,JLEV) &
        &          * (1._JPRB - PABSS1(JL,INDEX))) &
        &          + ZGASSRC &
        &          + PCLDFRAC(JL,JI,JLEV)*(ZBBDTOT*ZATOT1(JL,INDEX) - ZGASSRC)
      ZDRAD1(JL) = ZDRAD1(JL) + ZRADLD1(JL,JI)

!    if (JI == 1 .AND. Jlev <= 2) then
!      print 9891,jlev,JI, ZBGLAY,ZDELBGUP,ZBBU1(INDEX),ZBGLEV(JI),&
!      & PFRAC(JI,JLEV),ZPLVL(JI,JLEV-1),ZDELBGDN,ZBBD,PABSS1(INDEX),PTAUSF1(INDEX),ZRADCLRD1(JI)
9891  FORMAT(1X,2I4,11E12.5)
!    endif

    ENDDO
  ENDDO

!  print 9902,JLEV,ZDRADCL1
9902 FORMAT(1X,'ZDRADCL1(JL)=',I3,E13.6)

  DO JL = KIDIA, KFDIA
    PTOTDFLUC(JL,JLEV-1) = ZDRADCL1(JL) * ZWTNUM(JL,1)
    PTOTDFLUX(JL,JLEV-1) = ZDRAD1(JL)   * ZWTNUM(JL,1)
  ENDDO
ENDDO


! Spectral reflectivity and reflectance
! Includes the contribution of spectrally varying longwave emissivity 
! and reflection from the surface to the upward radiative transfer.
! Note: Spectral and Lambertian reflections are identical for the one
! angle flux integration used here.

DO JL = KIDIA, KFDIA
  ZURAD1(JL)   = 0.0_JPRB
  ZURADCL1(JL) = 0.0_JPRB
ENDDO

!IF (IREFLECT  ==  0) THEN
!- NB: Only Lambertian reflection is considered
DO JI = 1, JPGPT
  DO JL = KIDIA, KFDIA
!   Clear-sky radiance
    ZRADCLD(JL) = ZRADCLRD1(JL,JI)
    ZRADCLRU1(JL,JI) = ZRADUEMIT(JL,JI) + (1.0_JPRB - ZSEMIS(JL,JI)) * ZRADCLD(JL)
    ZURADCL1(JL) = ZURADCL1(JL) + ZRADCLRU1(JL,JI)

!   Total sky radiance
    ZRADD = ZRADLD1(JL,JI)
    ZRADLU1(JL,JI) = ZRADUEMIT(JL,JI) + (1.0_JPRB - ZSEMIS(JL,JI)) * ZRADD
    ZURAD1(JL) = ZURAD1(JL) + ZRADLU1(JL,JI)
  ENDDO
ENDDO
DO JL = KIDIA, KFDIA
  PTOTUFLUC(JL,0) = ZURADCL1(JL) * 0.5_JPRB
  PTOTUFLUX(JL,0) = ZURAD1(JL) * 0.5_JPRB
ENDDO

! Initialize the partial derivatives for each g point
DO JI = 1, JPGPT
  INBI = NGB(JI)
  DO JL = KIDIA, KFDIA
    ! In the simplest approximation at the surface these are the
    ! fraction of the total upwelling flux contributed by this g
    ! point
    !     ZLwDerivativeGpoint(JL,JI) &
    !          = ZRADLU1(JL,JI) / ZURAD1(JL)
    ! More accurately we consider the reweighting of the bands as the
    ! shape of the Planck function changes
     ZLWDERIVATIVEGPOINT(JL,JI) &
          & = PFRAC(JL,JI,1) * ZDBANDFLUXUP_DTOTALFLUXUP(JL,INBI)
  ENDDO
ENDDO

! The partial derivative of the surface flux with respect to the
! surface flux is one
DO JL = KIDIA, KFDIA
   PLWDERIVATIVE(JL,0) = 1.0_JPRB
ENDDO

!- Upward radiative transfer.
!- *** URAD1 holds the summed radiance for total sky stream
!- *** URADCL1 holds the summed radiance for clear sky stream
DO JLEV = 1, KLEV
  IENT=JPGPT * (JLEV-1)
  DO JL = KIDIA, KFDIA
    ZURAD1(JL)   = 0.0_JPRB
    ZURADCL1(JL) = 0.0_JPRB
    ! Initialize the partial derivatives to zero
    PLWDERIVATIVE(JL,JLEV) = 0.0_JPRB
  ENDDO

  DO JI = 1, JPGPT
    INDEX = IENT + JI
    DO JL = KIDIA, KFDIA

!  - clear-sky upward flux
      ZRADCLRU1(JL,JI) = ZRADCLRU1(JL,JI) + (ZBBU1(JL,INDEX)-ZRADCLRU1(JL,JI))&
        & *PABSS1(JL,INDEX)  
      ZURADCL1(JL) = ZURADCL1(JL) + ZRADCLRU1(JL,JI)

!  - total-sky upward flux          
      ZGASSRC = ZBBU1(JL,INDEX) * PABSS1(JL,INDEX)
      ! The combined transmissivity of gas and cloud is
      ! (1-emissivity_cloud)*(1-emissivity_gas):
      ZTRANSMISSIVITY = (1._JPRB - ZEFCLFRNW(JL,JI,JLEV)) &
           &           * (1._JPRB - PABSS1(JL,INDEX))
      ZRADLU1(JL,JI) = ZRADLU1(JL,JI) * ZTRANSMISSIVITY &
        &          + ZGASSRC &
        &          + PCLDFRAC(JL,JI,JLEV)*(ZBBUTOT1(JL,INDEX)*ZATOT1(JL,INDEX) - ZGASSRC)
      ZURAD1(JL) = ZURAD1(JL) + ZRADLU1(JL,JI)

! Accumulate the partial derivatives
      ZLWDERIVATIVEGPOINT(JL,JI) &
        & = ZLWDERIVATIVEGPOINT(JL,JI) * ZTRANSMISSIVITY
      PLWDERIVATIVE(JL,JLEV) &
        & = PLWDERIVATIVE(JL,JLEV) + ZLWDERIVATIVEGPOINT(JL,JI)
    ENDDO
  ENDDO

  DO JL = KIDIA, KFDIA
    PTOTUFLUC(JL,JLEV) = ZURADCL1(JL) * ZWTNUM(JL,1)
    PTOTUFLUX(JL,JLEV) = ZURAD1(JL)   * ZWTNUM(JL,1)
  ENDDO

ENDDO

!-----------------------------------------------------------------------


END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('RRTM_RTRN1A_140GP_MCICA',1,ZHOOK_HANDLE)
END SUBROUTINE RRTM_RTRN1A_140GP_MCICA
