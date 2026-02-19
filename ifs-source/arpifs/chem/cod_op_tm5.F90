! (C) Copyright 2009- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE COD_OP_TM5&
 & (YDDIMV, YDERAD,KIDIA, KFDIA, KLON , KLEV, KMXCOUNT, &
 & PQF, PTF, PAP, PPRESH, PPRESF,&
 & PLSF, PWND,  PLWCF, PIWCF,&
 &   PCLD_REFF,PCOD, PTAUC,PTAUA, PPMCLD )

!*** *COD_OP_TM5* Operator routine for cloud optical depth observations at 550 nm, 
!*** copied from COD_OP - Vincent Huijnen

!**   INTERFACE.
!     ----------
!          *COD_OP_TM5* IS CALLED FROM *CHEM_TM5*.
!     CALL COD_OP(KLON,ILEN,IMXCOUNT,ZQF5, ZTF5, ZPRESH5, ZPRESF5, ZLS5, ZCLWF5,ZCIWF5,ZXPP)
!     WHERE KIDIA, KFDIA, = START and End  (INPUT)
!           KLON       = Dimension (INPUT)
!           KMXCOUNT   = Number of channels (INPUT)
!           PQF        =  Specific humidity at observation points, model levels (INPUT)
!           PTF        =  Temperature at observation points, model levels (INPUT)
!           PPRESH,PPRESF =   Half/Full level pressure values at obs points, model levels (INPUT)
!           PLSF       = land-sea mask 
!           PLWCF      = cloud liquid water (INPUT)
!           PIWCF      = cloud ice water (INPUT)
!           ZXPP       = TOTAL Cloud optical depth at 550 nm (OUTPUT)
!
!           PCLD_REFF  = cloud effective radius
!           PCOD       = TOTAL CLOUD OPT DEPTH 
!           PTAUC      = CLOUD OPT DEPTH per layer
!           PPMCLD     = CLOUD asymmetry factor


!**   INTERFACE.
!     ----------
!          *COD_OP_TM5* IS CALLED FROM *CHEM_TM5*.

!     AUTHOR.
!     -------
!        ANGELA BENEDETTI  *ECMWF*

!     MODIFICATIONS.
!     --------------
!        20060303 - adaptation of Jean-Jacques Morcrette's cloprop550.F90
!        20071114 - output of PTAUC for each layer, naj
!        20090428 - J. Flemming: Add default formulation NRADIP=3
!        20120830 - V. Huijnen: Change standard cloud droplet size, 
!                               change output.

!     PURPOSE.
!     --------
!     Process cloud optical thickness at 0.55 um from mass mixing ratio
!     of liquid and ice water

!-----------------------------------------------------------------------
USE YOMDIMV  , ONLY : TDIMV
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMCST    ,ONLY : RD, RG, RPI, RTT
USE YOECLOP550,ONLY : RSA55, RSB55, RSC55,RSD55,RFA55
USE YOERAD    ,ONLY : TERAD
USE YOERDU    ,ONLY : REPLOG, REPSCW
USE YOETHF   , ONLY : RTICE
USE YOMLUN    ,ONLY : NULERR

IMPLICIT NONE
!-----------------------------------------------------------------------

!*       0.1   ARGUMENTS
!              ---------

TYPE(TDIMV)       ,INTENT(IN)    :: YDDIMV
TYPE(TERAD)       ,INTENT(INOUT) :: YDERAD
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA, KFDIA , KLON, KLEV , KMXCOUNT

REAL(KIND=JPRB)   ,INTENT(IN)    :: PPRESH(KLON,0:KLEV) ! PRESSURE ON HALF LEVELS FROM TRAJ
REAL(KIND=JPRB)   ,INTENT(IN)    :: PPRESF(KLON,KLEV)   ! PRESSURE ON FULL LEVELS FROM TRAJ
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTF(KLON,KLEV)    ! temperature FROM TRAJ
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQF(KLON,KLEV)    ! specific humidity FROM TRAJ
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAP(KLON,KLEV)      ! Cloud fraction
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLSF(KLON)             ! Land-sea mask 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PWND(KLON)             ! surface wind
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLWCF(KLON,KLEV)  ! liquid water content FROM TRAJ (kg/kg)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PIWCF(KLON,KLEV)  ! ice water content FROM TRAJ (kg/kg)


REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCLD_REFF(KLON,KLEV)   ! CLOUD effective radius 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCOD(KLON,KMXCOUNT)         ! TOTAL CLOUD OPT DEPTH 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTAUC(KLON,KLEV,KMXCOUNT) ! CLOUD OPT DEPTH per layer
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTAUA(KLON,KLEV,KMXCOUNT) ! CLOUD OPT DEPTH per layer - Absorption
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PPMCLD(KLON,KLEV)   ! CLOUD asymmetry factor

!*       0.5   LOCAL VARIABLES
!              ---------------

INTEGER(KIND=JPIM) :: JK, JL, JJ

REAL(KIND=JPRB) :: ZDESR(KLON), ZIWC, ZLWC
REAL(KIND=JPRB) :: ZRADIP(KLON), ZRADLP(KLON)
REAL(KIND=JPRB) :: ZNTOT, ZDEN(KLON)
!REAL(KIND=JPRB) :: ZTAUC(KLON,NFLEVG)

REAL(KIND=JPRB) :: Z1RADI, ZBETAI, ZDEFRE,&
                  & ZDPOG, ZNUM,&
                  & ZPODT, ZQIP, ZQLP, ZREFDE, ZTEMPC, ZTOI, ZTOLT,ZTOLA,ZTOLS
REAL(KIND=JPRB) :: ZRGI, ZRDI
REAL(KIND=JPRB) :: ZTCELS, ZAIWC, ZBIWC, ZFSR 

REAL(KIND=JPRB) :: ZCCNL(KLON),ZCCNO(KLON)
REAL(KIND=JPRB) :: ZK, ZQ, ZEXPO, ZEXPL, ZTEMP, ZREFAC, ZD, ZFAC
INTEGER(KIND=JPIM) :: IRADLP

! Local, best estimate of global average CCN over ocean (N/cm3)
REAL(KIND=JPRB),PARAMETER :: ZCCNSEA=40.

! Slingo E and F coefficients for asymmetry parameter for first wavelength band
REAL(KIND=JPRB) :: ZGC
REAL(KIND=JPRB),PARAMETER :: RSE1=0.829
REAL(KIND=JPRB),PARAMETER :: RSF1=2.482E-03

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"
!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('COD_OP_TM5',0,ZHOOK_HANDLE)
ASSOCIATE(NFLEVG=>YDDIMV%NFLEVG, &
 & LCCNL=>YDERAD%LCCNL, LCCNO=>YDERAD%LCCNO, RCCNLND=>YDERAD%RCCNLND, &
 & RCCNSEA=>YDERAD%RCCNSEA, RMINICE=>YDERAD%RMINICE, RRE2DE=>YDERAD%RRE2DE)
!-- Select method for cloud radius calculation
IRADLP = 3_JPIM

!VH I noted that option IRADLP=2 leads to rather low effective cloud 
!VH    droplet radius, on average 4 um, rather than ~8um.

!-- units:
!   ------
! PLWCF  is the liquid water mass mixing ratio      kg kg-1
! PIWCF  is the ice water ,ass mixing ratio         kg kg-1

ZREFDE = RRE2DE
ZDEFRE = 1.0_JPRB / ZREFDE
ZRGI = 1.0_JPRB/RG
ZRDI = 1.0_JPRB/RD

PCOD(KIDIA:KFDIA,1:KMXCOUNT) = 0._JPRB

IF(KMXCOUNT /= 1) THEN
 WRITE(NULERR,*) 'COD_OP_TM5: ROUTINE CANNOT HANDLE MORE THAN WAVELENGHT'
 WRITE(NULERR,*) 'PROGRAM WILL STOP'
 CALL ABOR1('')
ENDIF

IF (IRADLP==2) THEN

IF (LCCNO) THEN
  DO JL=KIDIA,KFDIA
    IF (PWND(JL) > 30._JPRB) THEN
      ZQ=327._JPRB
    ELSEIF (PWND(JL) > 15._JPRB) THEN
      ZQ=EXP(0.13_JPRB*PWND(JL)+1.89_JPRB)
    ELSE
      ZQ=EXP(0.16_JPRB*PWND(JL)+1.44_JPRB)
    ENDIF
    ZEXPO=1.2_JPRB+0.5_JPRB*LOG10(ZQ)
    ZCCNO(JL)=10._JPRB**ZEXPO
    ZCCNO(JL)= MIN(  287._JPRB,MAX( 32._JPRB,ZCCNO(JL)))
  ENDDO
ELSE
  DO JL=KIDIA,KFDIA
    ZCCNO(JL)=RCCNSEA
  ENDDO
ENDIF

IF (LCCNL) THEN
  DO JL=KIDIA,KFDIA
    IF (PWND(JL) <= 15._JPRB) THEN
      ZQ=EXP(0.16_JPRB*PWND(JL)+1.45_JPRB)
    ELSE
      ZQ=EXP(0.13_JPRB*PWND(JL)+1.89_JPRB)
    ENDIF
    ZEXPL=2.21_JPRB+0.3_JPRB*LOG10(ZQ)
    ZCCNL(JL)=10._JPRB**ZEXPL
    ZCCNL(JL) = MIN( 1743._JPRB,MAX(292._JPRB,ZCCNL(JL)))
  ENDDO
ELSE
  DO JL=KIDIA,KFDIA
    ZCCNL(JL)=RCCNLND
  ENDDO
ENDIF 

DO JL=KIDIA,KFDIA
  IF (PLSF(JL) < 0.5_JPRB) THEN
    ! Spectral dispersion of the droplet size distribution over sea 
    ! ZD = 0.33_JPRB
    ! Constant (k) relating volume radius to effective radius rv^3=k*re^3
    ! k = (1.0+ZD*ZD)^3)/((1.0+3.0*ZD*ZD)^2)
    ! ZK=(1.0_JPRB+ZD*ZD)**3/(1.0_JPRB+3.0_JPRB*ZD*ZD)**2
    ZK = 0.77_JPRB

    ! Cloud droplet concentration in cm-3 (activated CCN) over ocean
    ZNTOT=-1.15E-03_JPRB*ZCCNO(JL)*ZCCNO(JL)+0.963_JPRB*ZCCNO(JL)+5.30_JPRB
  ELSE
    ! Spectral dispersion (d) of the droplet size distribution over land 
    ! ZD=0.43_JPRB
    ! Constant (k) relating volume radius to effective radius rv^3=k*re^3
    ! k = (1.0+ZD*ZD)^3)/((1.0+3.0*ZD*ZD)^2)
    ! ZK=(1.0_JPRB+ZD*ZD)**3/(1.0_JPRB+3.0_JPRB*ZD*ZD)**2
    ZK = 0.69_JPRB

    ! Cloud droplet concentration in cm-3 (activated CCN) over land
    ZNTOT=-2.10E-04_JPRB*ZCCNL(JL)*ZCCNL(JL)+0.568_JPRB*ZCCNL(JL)-27.9_JPRB

  ENDIF
  ZDEN(JL)=4.0_JPRB*RPI*ZNTOT*ZK ! *density of water 1000 kg m-3
ENDDO

ENDIF ! preparation for IRADLP == 2

DO JJ=1,KMXCOUNT
  DO JK=1,KLEV
!DEC$ IVDEP
    DO JL=KIDIA,KFDIA
      IF ( PAP(JL,JK) >=0.001_JPRB ) THEN
        ZTEMP=1.0_JPRB/PAP(JL,JK)
        ZDPOG=ZRGI*(PPRESH(JL,JK)-PPRESH(JL,JK-1))

!-- cloud and ice water path in g m-2
        ZQIP=MAX(0._JPRB,1.E+03_JPRB*ZDPOG*PIWCF(JL,JK)*ZTEMP)
        ZQLP=MAX(0._JPRB,1.E+03_JPRB*ZDPOG*PLWCF(JL,JK)*ZTEMP)
        ! ZQIP=MAX(0._JPRB,1.E+03_JPRB*ZDPOG*PIWCF(JL,JK))
        ! ZQLP=MAX(0._JPRB,1.E+03_JPRB*ZDPOG*PLWCF(JL,JK))
!-- cloud and ice water content in g m-3
        ZPODT=ZRDI*PPRESF(JL,JK)/PTF(JL,JK)
        ZIWC=1.E+03_JPRB*PIWCF(JL,JK)*ZPODT*ZTEMP 
        ZLWC=1.E+03_JPRB*PLWCF(JL,JK)*ZPODT*ZTEMP

      ELSE
        ZQIP = 0._JPRB
        ZQLP = 0._JPRB
        ZLWC = 0._JPRB
        ZIWC = 0._JPRB
      ENDIF

!-- cloud and ice water content in g m-3
      !VH ZPODT = PPRESF(JL,JK)/(RD*PTF(JL,JK))
      !VH ZIWC= ZQLP*ZPODT
      !VH ZLWC= ZQIP*ZPODT


      
      IF ( IRADLP == 1 ) THEN 
        ! Very simple approximation...
 
        ! simple distinction between land (10) and ocean (13) Zhang and Rossow
        ! VH these cloud droplet sizes are globally too high.  
        ! VH Suggestion to use something lower
     
        IF (PLSF(JL) < 0.5_JPRB) THEN
          ZRADLP(JL)=10.0_JPRB
        ELSE
          ZRADLP(JL)=7.0_JPRB
        ENDIF 
      
      ELSEIF ( IRADLP == 2 ) THEN
        ZNUM  = 3.0_JPRB*ZLWC
        !VH ZNUM  = 3.0_JPRB*(ZLWC+ZRWC)
        ZTEMP = 1.0_JPRB/ZDEN(JL)
        ! So for ZTEMP in SI units
        ! g m-3 and cm-3 units cancel out with density of water 10^6/(1000*1000)
        ! Need a factor of 10^6 to convert to microns 
        ! and cubed root is factor of 100 which appears in ZRADLP equation below

        ! Adjustment to Martin_et_al(1994) effective radius 
        ! from Wood(2000) param (Eq. 19)
        !VH - switch this off... I don't have ZRWC 
        ! IF (ZLWC > REPSCW) THEN
        !  !   IF (ZLWC(JL,JK) > REPSCW .AND. PLSM(JL) < 0.5_JPRB) THEN
        !  ZRLRATIO  = ZRWC(JL,JK)/ZLWC(JL,JK)
        !  ZKL       = 0.222_JPRB
        !  ZKRATIO   = (ZKL/ZK)**0.333_JPRB
        !  ZREFACDEN = 1.0_JPRB+0.2_JPRB*ZKRATIO*ZRLRATIO
        !  ZREFAC    = ((1.0_JPRB+ZRLRATIO)**0.666_JPRB)/ZREFACDEN
        ! ELSE
         ZREFAC    = 1.0_JPRB
        ! ENDIF
 
        IF((ZNUM*ZTEMP) > REPLOG)THEN

          ! R_eff=100*(ZNUM/ZDEN)^0.333, factor of 100 to convert from m to microns
          ZRADLP(JL)=100._JPRB*ZREFAC*EXP(0.333_JPRB*LOG(ZNUM*ZTEMP))

          ! Limit effective radius to within defined range
          ZRADLP(JL)=MAX(ZRADLP(JL), 4.0_JPRB)
          ! Limit maximum to ~16. 30 seems too high.
          ZRADLP(JL)=MIN(ZRADLP(JL),16.0_JPRB)
        ELSE
          ZRADLP(JL)=4.0_JPRB
        ENDIF
!
      ELSEIF (IRADLP ==3) THEN
      
!-- droplet effective radius, Martin et al., 1994, using fixed CCN for sea and land
        IF (PLSF(JL) < 0.5_JPRB) THEN
          ! Sea
          ZD = 0.333_JPRB
          ! ZNTOT=-1.15E-03_JPRB*RCCNSEA*RCCNSEA+0.963_JPRB*RCCNSEA+5.30_JPRB
          ZNTOT=-1.15E-03_JPRB*ZCCNSEA*ZCCNSEA+0.963_JPRB*ZCCNSEA+5.30_JPRB
          ZDEN(JL)=4.0_JPRB*RPI*ZNTOT*(1.0_JPRB+ZD*ZD)**3
          ZNUM=3.0_JPRB*ZLWC*(1.0_JPRB+3.0_JPRB*ZD*ZD)**2
        ELSE
          !Land 
          ZD = 0.43_JPRB
          ZNTOT=-2.10E-04_JPRB*RCCNLND*RCCNLND+0.568_JPRB*RCCNLND-27.9_JPRB
          ZDEN(JL)=4.0_JPRB*RPI*ZNTOT*(1.0_JPRB+ZD*ZD)**3
          ZNUM=3.0_JPRB*ZLWC*(1.0_JPRB+3.0_JPRB*ZD*ZD)**2
        ENDIF
        IF((ZNUM/ZDEN(JL)) > REPLOG)THEN
          ZRADLP(JL)=100.0_JPRB*EXP(0.333_JPRB*LOG(ZNUM/ZDEN(JL)))
          ZRADLP(JL)=MAX(ZRADLP(JL), 4.0_JPRB)
          ZRADLP(JL)=MIN(ZRADLP(JL),16.0_JPRB)
        ELSE
          ZRADLP(JL)=4.0_JPRB
        ENDIF
      
      ELSEIF (IRADLP == 4) THEN
         ! Simple scaling based on liquid water content...  
         ZRADLP(JL)=4.0_JPRB + 11.0_JPRB*ZLWC
         ! Simple scaling based on liquid water path...  
         !VH - 2014-04-15 -that's wrong - it should be LWC!, see, e.g., Martin et al., AMS 1994)
         ! ZRADLP(JL)=4.0_JPRB + 11.0_JPRB*ZQLP
         ZRADLP(JL)=MIN(ZRADLP(JL),12.0_JPRB)
      ENDIF   ! ( IRADLP == 4 )

      PCLD_REFF(JL,JK)=ZRADLP(JL)

!-- ice particle effective size

      IF (PTF(JL,JK) < RTICE) THEN
        ZTEMPC=PTF(JL,JK)-RTT
      ELSE
        ZTEMPC=RTICE-RTT
      ENDIF
      ZRADIP(JL)=326.3_JPRB+ZTEMPC*(12.42_JPRB + ZTEMPC*(0.197_JPRB + ZTEMPC*&
       & 0.0012_JPRB))
      ZRADIP(JL)=MAX(ZRADIP(JL),30.0_JPRB)
      ZRADIP(JL)=MIN(ZRADIP(JL),60.0_JPRB)
      ZDESR(JL)= ZDEFRE * ZRADIP(JL)

!- default formulation NRADIP=3
      ZTEMPC=PTF(JL,JK)-83.15_JPRB
      ZTCELS=PTF(JL,JK)-RTT
      ZFSR=1.2351_JPRB+0.0105_JPRB*ZTCELS
      ZAIWC=45.8966_JPRB*ZIWC**0.2214_JPRB
      ZBIWC=0.7957_JPRB *ZIWC**0.2535_JPRB
      ZDESR(JL)=ZFSR*(ZAIWC+ZBIWC*ZTEMPC)
      ZDESR(JL)=MIN(MAX(ZDESR(JL),RMINICE), 155.0_JPRB)

!-- cloud optical properties

      ZTOLS=0.0_JPRB
      ZTOLA=0.0_JPRB
      ZTOI =0.0_JPRB
      ZGC  =0.0_JPRB
      
!-- for liquid water clouds (Slingo, 1989)

      IF (ZQLP > REPSCW) THEN
        ZTOLT = ZQLP *(RSA55 + RSB55 /ZRADLP(JL))
        ZGC   = RSE1 + RSF1*ZRADLP(JL)
!-- distribute over scattering and absorbing components
 
       ZFAC = 1._JPRB - RSC55 - RSD55*ZRADLP(JL)  
       ZFAC = MIN(ZFAC,.999999_JPRB)        

       ZTOLS = ZTOLT * ZFAC
       ZTOLA = ZTOLT * (1._JPRB-ZFAC)

      ENDIF



!-- for ice water clouds (Fu, 1996)

      IF (ZQIP > REPSCW) THEN
        Z1RADI = 1.0_JPRB / ZDESR(JL)
        ZBETAI = RFA55(1) + Z1RADI * RFA55(2)
        ZTOI = ZQIP * ZBETAI
      ENDIF

      !VH ZTAUC(JL,JK)=ZTOI+ZTOL
      PTAUC(JL,JK,JJ)=ZTOI+ZTOLS
      PTAUA(JL,JK,JJ)=ZTOLA
      PCOD(JL,JJ)=PCOD(JL,JJ)+ PTAUC(JL,JK,JJ)

      !VH TM5 requires optical depth per layer rather than integrated from above.
      !VH PTAUC(JL,JK,JJ)=PCOD(JL,JJ)

!-- Asymmetry factor

      IF (PTAUC(JL,JK,JJ) > 0.) THEN 
        PPMCLD(JL,JK)=ZGC 
      ELSE 
        PPMCLD(JL,JK)=0._JPRB
      ENDIF

   ENDDO
 ENDDO
ENDDO

! BEN - add check for optical depth above 100.
DO JJ=1,KMXCOUNT
  DO JL=KIDIA,KFDIA
   IF(PCOD(JL,JJ) > 100._JPRB ) THEN
    PCOD(JL,JJ) = 100.0_JPRB
   ENDIF
  ENDDO
ENDDO



!-----------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('COD_OP_TM5',1,ZHOOK_HANDLE)
END SUBROUTINE COD_OP_TM5

