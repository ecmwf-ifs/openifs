! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE AER_CLCLD ( & 
!  input
 & YDMODEL,KIDIA,    KFDIA,    KLON,   KLEV,&
 & PTSPHY,&
 & PT,       PQ,       PQSAT, &
 & PAPH,    PAP,&
 & PGELAM,   PGEMU,    PCLON,   PSLON,&
 & PL,       PI,       PA,&
!  output
 & PLCRIT_AER,&
 & PICRIT_AER,&
 & PRE_LIQ,&
 & PRE_ICE,&
 & PCCN,&
 & PNICE )
!  diagnostics
! & PEXTRA,   KFLDX )

!     AER_CLCLD 
!     ---------
!          Sets up all information for aerosols 
!          effects on clouds and convection

!     AUTHOR
!          A. Tompkins  E.C.M.W.F.

!     PURPOSE.
!     --------

!     INTERFACE.
!     ----------

!          *AER_CLCLD* IS CALLED FROM *CALLPAR*

!     PARAMETER     DESCRIPTION                                   UNITS
!     ---------     -----------                                   -----

! -   INPUT ARGUMENTS.
!     -------------------

! KIDIA   : START OF HORIZONTAL LOOP
! KFDIA   : END   OF HORIZONTAL LOOP
! KLON    : HORIZONTAL DIMENSION
! KLEV    : END OF VERTICAL LOOP AND VERTICAL DIMENSION
! PGELAM     : LONGITUDE
! PCLON      : COSINE OF LONGITUDE
! PSLON      : SINE   OF LONGITUDE
! PGEMU      : SINE OF LATITUDE

! -   OUTPUT ARGUMENTS.
!     -------------------
! PLCRIT_AER : critical liquid mmr for autoconversion process
! PICRIT_AER : critical liquid mmr for autoconversion process
! PRE_LIQ : liq Re
! PRE_ICE : ice Re
! PCCN    : liquid cloud condensation nuclei
! PNICE   : ice number concentration (cf. CCN)

!     Modifications:
!     --------------
!      K. Yessad (July 2014): Move some variables.
! -----------------------------------------------------------------------------

USE TYPE_MODEL , ONLY : MODEL
USE PARKIND1  ,ONLY : JPIM     ,JPRB, JPRD
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMCST   , ONLY : RG       ,RD       ,RETV     ,&
 & RLVTT    ,RLSTT    ,RTT     ,RPI
USE YOMCT3    ,ONLY : NSTEP
USE YOETHF   , ONLY : R2ES     ,R3LES    ,R3IES    ,R4LES    ,&
 & R4IES    ,R5LES    ,R5IES    ,R5ALVCP  ,R5ALSCP  ,&
 & RALVDCP  ,RALSDCP  ,RTWAT    ,&
 & RTICE    ,RTICECU  ,&
 & RTWAT_RTICE_R      ,RTWAT_RTICECU_R,&
 & RKOOP1   ,RKOOP2
USE YOEAEROP , ONLY : ALF_BC, ALF_DD, ALF_SS, ALF_SU

! -----------------------------------------------------------------------------

IMPLICIT NONE

! input variables
TYPE(MODEL)       ,INTENT(INOUT) :: YDMODEL
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSPHY 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPH(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAP(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PT(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQ(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PL(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PI(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PA(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQSAT(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGELAM(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCLON(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSLON(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEMU(KLON) 

! output
REAL(KIND=JPRB)   ,INTENT(INOUT)    :: PLCRIT_AER(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT)    :: PICRIT_AER(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT)    :: PRE_LIQ(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT)    :: PRE_ICE(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT)    :: PCCN(KLON,KLEV)     ! liquid cloud condensation nuclei
REAL(KIND=JPRB)   ,INTENT(INOUT)    :: PNICE(KLON,KLEV)    ! ice number concentration (cf. CCN)

! diagnostics
!INTEGER(KIND=JPIM),INTENT(IN)    :: KFLDX
!REAL(KIND=JPRB)   ,INTENT(INOUT) :: PEXTRA(KLON,KLEV,KFLDX) 
!! LOGICAL, PARAMETER :: LDUMPEXTRA=.FALSE.

!----------------------------------------------------------------------

! local arrays
REAL(KIND=JPRD) :: ZGEMU(KLON)

! Aerosol arrays
REAL(KIND=JPRB) :: ZQAER(KLON,6,KLEV)  ! aerosol mass
REAL(KIND=JPRB) :: ZMAER(KLON,KLEV,6)  ! mass of aerosol
REAL(KIND=JPRB) :: ZMAERMN(6) ! annual column mean mass of aerosol
REAL(KIND=JPRB) :: ZQOZ(KLON,KLEV)     ! dummy ozone array
REAL(KIND=JPRB) :: ZECPO3(KLON,KLEV)   ! dummy prognostic ozone array
!REAL(KIND=JPRB) :: ZRE_LIQ(KLON,KLEV)  ! Effective radius
!REAL(KIND=JPRB) :: ZRE_ICE(KLON,KLEV)  ! ice effect radius
REAL(KIND=JPRB) :: ZICENUCLEI(KLON,KLEV) ! number concentration of ice nuclei
REAL(KIND=JPRB) :: ZQS(KLON,KLEV)      ! saturation
REAL(KIND=JPRB) :: ZDUMAER(KLON,KLEV,12)

REAL(KIND=JPRB) :: ZS0, ZSCRITHOMO, ZSVP, ZTEMPC
REAL(KIND=JPRB) :: ZNCRIT_GIERENS, ZNCRIT_REN
REAL(KIND=JPRB) :: ZNICEHOMO
!REAL(KIND=JPRB) :: ZLIQCLD, ZICERE

REAL(KIND=JPRB) :: ZCLD
REAL(KIND=JPRB) :: ZRHO_ICE, ZRHO_LIQ ! density of pristine ice crystals and cloud
REAL(KIND=JPRB) :: ZRLIQ_CRIT, ZRICE_CRIT ! critical radii for autoconversion process

! for RH look up tables
INTEGER(KIND=JPIM) :: IRH(KLON,KLEV)
INTEGER(KIND=JPIM) :: JTYP, JTAB, IBIN
INTEGER(KIND=JPIM) :: JAERSS, JAERDD, JAEROM, JAERSU, JAERBC

! these are reduced compared to USE YOEAEROP, since 1 band only
!REAL(KIND=JPRB) :: ZALF_BC(1)     
!REAL(KIND=JPRB) :: ZALF_DD(3)   
!REAL(KIND=JPRB) :: ZALF_OM(12)  
!REAL(KIND=JPRB) :: ZALF_SS(12,3)
!REAL(KIND=JPRB) :: ZALF_SU(12)  
!REAL(KIND=JPRB) :: ZRHTAB(12)  

REAL(KIND=JPRB) :: ZALF, ZRH, ZWTOT

! general arrays
REAL(KIND=JPRB) :: ZTHF(KLON,KLEV+1)   ! T on half levels
REAL(KIND=JPRB) :: ZRHO(KLON,KLEV)     ! density

REAL(KIND=JPRB) :: ZDPR

! misc variables
REAL(KIND=JPRB) :: ZEPSEC
INTEGER(KIND=JPIM) :: IWAVL
INTEGER(KIND=JPIM) :: JK, JL

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!--------------------------------------------------------------------------

!------------------------
! interface include files
!------------------------
#include "radact.intfb.h"

#include "fcttre.func.h"
#include "fccld.func.h"

!--------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('AER_CLCLD',0,ZHOOK_HANDLE)
ASSOCIATE(YDEAERSNK=>YDMODEL%YRML_PHY_AER%YREAERSNK,YDECLDP=>YDMODEL%YRML_PHY_EC%YRECLDP, &
 & YDRIP=>YDMODEL%YRML_GCONF%YRRIP,YDERAD=>YDMODEL%YRML_PHY_RAD%YRERAD)
ASSOCIATE(RRHTAB=>YDEAERSNK%RRHTAB, &
 & LAERICEAUTO=>YDECLDP%LAERICEAUTO, LAERICESED=>YDECLDP%LAERICESED, &
 & LAERLIQAUTOLSP=>YDECLDP%LAERLIQAUTOLSP, LAERLIQCOLL=>YDECLDP%LAERLIQCOLL, &
 & RCCNSS=>YDECLDP%RCCNSS, RCCNSU=>YDECLDP%RCCNSU, RCLCRIT=>YDECLDP%RCLCRIT, &
 & RCLDMAX=>YDECLDP%RCLDMAX, RLCRITSNOW=>YDECLDP%RLCRITSNOW, &
 & RNICE=>YDECLDP%RNICE, &
 & LEPO3RA=>YDERAD%LEPO3RA, &
 & NSTART=>YDRIP%NSTART)
!--------------------------------------------------------------------------


!######################################################################
!                       1.0 Basic variables
!######################################################################
! IF (KFLDX>0 .AND. LDUMPEXTRA) PEXTRA(:,:,:)=0.0_JPRB

ZEPSEC=1.E-10_JPRB
IWAVL=8               ! reference to 550 nm

! move to cldp module
ZRHO_ICE=900.0_JPRB
ZRHO_LIQ=1000.0_JPRB
ZRICE_CRIT=60.0E-6_JPRB ! ice to snow critical radius
ZRLIQ_CRIT=9.3E-6_JPRB  ! cloud to rain critical radius
ZWTOT=0.1_JPRB ! governs critical N

! change to vector vmass 
DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    ZQS(JL,JK)=FOEEWM(PT(JL,JK))/PAP(JL,JK)
    ZQS(JL,JK)=MIN(0.5_JPRB,ZQS(JL,JK))
    ZQS(JL,JK)=ZQS(JL,JK)/(1.0_JPRB-RETV*ZQS(JL,JK))
  ENDDO
ENDDO
DO JK=2,KLEV
  DO JL=KIDIA,KFDIA
    ZTHF(JL,JK)=(PT(JL,JK-1)*PAP(JL,JK-1)&
   & *(PAP(JL,JK)-PAPH(JL,JK))&
   & +PT(JL,JK)*PAP(JL,JK)*(PAPH(JL,JK)-PAP(JL,JK-1)))&
   & *(1.0_JPRB/(PAPH(JL,JK)*(PAP(JL,JK)-PAP(JL,JK-1))))  
  ENDDO
ENDDO

DO JL=KIDIA,KFDIA
  ZTHF(JL,KLEV+1)=PT(JL,KLEV) ! should be surface temperature
  ZTHF(JL,1)=PT(JL,1)
ENDDO

!------------------------
! density
!------------------------
DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    ZRHO(JL,JK)=PAP(JL,JK)/(RD*PT(JL,JK))
  ENDDO
ENDDO
IF (LEPO3RA) THEN
  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      ZECPO3(JL,JK)=0.0_JPRB
    ENDDO
  ENDDO
ENDIF

!######################################################################
!               2. Retrieve aerosols optical depths 
!######################################################################
! line 2: KRINT=1, KDLON=KLON , P2=KLON, KSHIFT=0, 
! line 4: ozone set to dummy variable.
ZGEMU=PGEMU
CALL RADACT ( YDMODEL%YRML_PHY_RAD%YREAERD,YDERAD,YDEAERSNK,YDRIP, &
     & YDMODEL%YRML_GCONF%YRSPP_CONFIG,  & 
     & KIDIA , KFDIA, KLON , KLEV,&
     & 1    , KLON  , KLON , 0    , 1   ,&
     & PAPH , &
     & PGELAM, ZGEMU, PCLON, PSLON, ZTHF,&
     & PQ   , PQSAT , ZECPO3,&
     & ZQAER, ZDUMAER, ZQOZ  )  

! 1=sulphate+organic
! 2=sea salt
! 3=sand dust
! 4=black carbon
! 5=Volcanic
! 6=Background

!DO JK=2,KLEV
!  DO JL=KIDIA,KFDIA
!    IF (KFLDX>0 .AND. LDUMPEXTRA) THEN
!    PEXTRA(JL,JK,1)=ZQAER(JL,1,JK)
!    PEXTRA(JL,JK,2)=ZQAER(JL,2,JK)
!    PEXTRA(JL,JK,3)=ZQAER(JL,3,JK)
!    PEXTRA(JL,JK,4)=ZQAER(JL,4,JK)
!    !PEXTRA(JL,JK,5)=ZQAER(JL,5,JK)
!    !PEXTRA(JL,JK,6)=ZQAER(JL,6,JK)
!    ENDIF
!  ENDDO
!ENDDO

!######################################################################
!        3. Retrieve Aerosol mass from Tau
!######################################################################

!     Tegen order:
!     1=sulphate+organic
!     2=sea salt
!     3=sand dust
!     4=black carbon
!     5=Volcanic
!     6=Background
! set up indexes
JAERSU=1
JAERSS=2
JAERDD=3
JAERBC=4
JAEROM=1

! Table of column/annual mean mass of aerosol
! converted to microgram per m**3
ZMAERMN(JAERSU)=1.02E-09_JPRB*1.E9_JPRB 
ZMAERMN(JAERSS)=2.12E-10_JPRB*1.E9_JPRB 
ZMAERMN(JAERDD)=1.01E-09_JPRB*1.E9_JPRB 
ZMAERMN(JAERBC)=3.05E-11_JPRB*1.E9_JPRB 

!========================================
! conversion of aerosols from Tau to Mass 
!========================================
!-- define RH index from "clear-sky" (not yet!) relative humidity

! for now fix to using Bin 1, meaning the smallest SS particles
IBIN=1

! taken from YOEAEROP: 17- band data for IWAVL=8 corresponding to 550nm

DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    ZRH=100.0_JPRB*PQ(JL,JK)/PQSAT(JL,JK)
    ZRH=MIN(MAX(ZRH,1.0_JPRB),100.0_JPRB)
    DO JTAB=1,12
      IF (ZRH > RRHTAB(JTAB)) THEN
        IRH(JL,JK)=JTAB
      ENDIF
    ENDDO
  ENDDO
ENDDO

ZMAER(:,:,:) = 0.0_JPRB

DO JTYP=1,4
  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      IF      (JTYP == JAERSS) THEN
        ZALF=ALF_SS(IRH(JL,JK),IWAVL,IBIN)
      ELSEIF (JTYP == JAERDD) THEN
        ZALF=ALF_DD(IBIN,IWAVL)
      ELSEIF (JTYP == JAERSU) THEN
        ZALF=ALF_SU(IRH(JL,JK),IWAVL)
      ELSEIF (JTYP == JAERBC) THEN
        ZALF=ALF_BC(IWAVL)
      ELSEIF (JTYP == 5 .OR. JTYP == 6) THEN
        ZALF=ALF_SU(IRH(JL,JK),IWAVL)
      ENDIF

!      IF (KFLDX>0 .AND. LDUMPEXTRA) PEXTRA(JL,JK,4+JTYP)=0.0_JPRB

      ZDPR=PAPH(JL,JK+1)-PAPH(JL,JK)
      IF (ZALF /= 0.0_JPRB .AND. ZDPR /=0.0_JPRB ) THEN 
!       kg/kg
        ZMAER(JL,JK,JTYP) = ZQAER(JL,JTYP,JK)*RG/(ZDPR*ZALF*1000._JPRB)
! temporary - store in extra field
!       micrograms per m**3
        ZMAER(JL,JK,JTYP) = ZMAER(JL,JK,JTYP)*ZRHO(JL,JK)*1.E9_JPRB
      ENDIF
    ENDDO
  ENDDO
ENDDO
DO JK=1,KLEV
  IF (NSTEP == NSTART+10) THEN
    JL=(KIDIA+KFDIA)/2
  ENDIF
ENDDO

! store after adjustment to dust
!DO JTYP=1,4
!  DO JK=1,KLEV
!    DO JL=KIDIA,KFDIA
!      IF (KFLDX>0 .AND. LDUMPEXTRA) PEXTRA(JL,JK,4+JTYP)=ZMAER(JL,JK,JTYP)
!    ENDDO
!  ENDDO
!ENDDO

!######################################################################
!                4. WARM PHASE MICROPHYSICS
!######################################################################

IF (LAERLIQAUTOLSP.OR.LAERLIQCOLL) THEN

!---------------------------------------------------------------------
! Turn aerosol mass into a CCN Number concentration for warm rain
! From Menon et al JAS 2002
!---------------------------------------------------------------------
  DO JK=1,KLEV
!DEC$ IVDEP
    DO JL=KIDIA,KFDIA

! note - for the moment organic is combined with sulphate...

!     N cm**-3
      PCCN(JL,JK)=10.0_JPRB**(2.41+&
! CCNSU factor average for sulphate 0.50 and organic 0.15
                    & RCCNSU*LOG10(MAX(ZMAER(JL,JK,JAERSU),ZEPSEC)) +&
                    & RCCNSS*LOG10(MAX(ZMAER(JL,JK,JAERSS),ZEPSEC)) )
      PCCN(JL,JK)=MAX(PCCN(JL,JK),50._JPRB)
!      IF (KFLDX>0 .AND. LDUMPEXTRA) PEXTRA(JL,JK,9)=PCCN(JL,JK)

! this is effective radius calculation - not used for now

      ZCLD=PL(JL,JK)/MAX(PA(JL,JK),ZEPSEC)
      ZCLD=MIN(MAX(ZCLD,0.0_JPRB),RCLDMAX)

! number is 3/(4*pi*rho_liq*10^6)  [10^6 for N in right units]
      PRE_LIQ(JL,JK)=(2.387E-10_JPRB*ZRHO(JL,JK)*ZCLD/PCCN(JL,JK))**0.333_JPRB
!      IF (KFLDX>0 .AND. LDUMPEXTRA) PEXTRA(JL,JK,11)=1.e6_JPRB*PRE_LIQ(JL,JK) ! microns
    ENDDO
  ENDDO
!---------------------------------------------------------------------
! Turn CCN Number concentration for warm rain into critical 
! mixing ratio for autoconversion process
! Obvious, but also from Rotstayn and penner JAS 2001
!---------------------------------------------------------------------

! Number is 1.e6. 4/3 .pi.rho_l. r_crit**3
! 1e6 since CCN in units of  N cm**-3
! r_crit=9.3 microns
  DO JK=1,KLEV
!DEC$ IVDEP
    DO JL=KIDIA,KFDIA
      PLCRIT_AER(JL,JK)=1.333E6_JPRB*RPI*ZRHO_LIQ*PCCN(JL,JK)*ZRLIQ_CRIT**3.0_JPRB/ZRHO(JL,JK)
!      PLCRIT_AER(JL,JK)=3.36e-6_JPRB*PCCN(JL,JK)/ZRHO(JL,JK)
!      IF (KFLDX>0 .AND. LDUMPEXTRA) PEXTRA(JL,JK,10)=PLCRIT_AER(JL,JK)
! limit the effect to ratio of the "background" value
      PLCRIT_AER(JL,JK)=MAX(PLCRIT_AER(JL,JK),0.1_JPRB*RCLCRIT)
      PLCRIT_AER(JL,JK)=MIN(PLCRIT_AER(JL,JK),10.0_JPRB*RCLCRIT)
    ENDDO
    IF (NSTEP == NSTART+10) THEN
      JL=(KIDIA+KFDIA)/2
    ENDIF
  ENDDO
ENDIF

!######################################################################
!           5. ICE PHASE MICROPHYSICS
!######################################################################

IF (LAERICESED.OR.LAERICEAUTO) THEN

!---------------------------------------------------------------------
! Turn aerosol mass into a Ice Number concentration for ice processes
!---------------------------------------------------------------------
  DO JK=1,KLEV
!DEC$ IVDEP
    DO JL=KIDIA,KFDIA
!                       0.01_JPRB is "default" value from
! Demott et al. Ice SS=55% or Meyers et al. 1992 JAS, ISS=25%
! In a prognostic scheme this will be function of clear sky humidity

! By relating IN to Aerosol mass we are assuming that the mode of the 
! Aerosol size distribution lies in the accumulation or coarse mode

! The relationship will implicitly introduce the exponential height
! dependence that Sassen (1992) and K and Curry (1998) explicitly 
! introduced to their parametrizations.

      ZICENUCLEI(JL,JK)=0.01_JPRB*&
    &  (ZMAER(JL,JK,JAERSU)+ZMAER(JL,JK,JAERBC)+ZMAER(JL,JK,JAERDD)) &
    & /(ZMAERMN(JAERSU)    +ZMAERMN(JAERBC)    +ZMAERMN(JAERDD))
      ZICENUCLEI(JL,JK)=MAX(ZICENUCLEI(JL,JK),0.0_JPRB)
!      IF (KFLDX>0 .AND. LDUMPEXTRA) PEXTRA(JL,JK,12)=ZICENUCLEI(JL,JK)

! T in oC
      ZTEMPC=PT(JL,JK)-RTT

! Re form for ice crystals from Liou and Oort 1994
! used to derive Re(ice) as in Lohmann JC 2002 
      ZNICEHOMO=0.0_JPRB
      PRE_ICE(JL,JK)=0.5_JPRB*(326.3_JPRB+ZTEMPC* &
        & (12.42_JPRB + ZTEMPC*(0.197_JPRB + ZTEMPC*0.0012_JPRB)))
      PRE_ICE(JL,JK)=MAX(PRE_ICE(JL,JK),0.0_JPRB)

! effect Re to volume mean from S Moss or Lohmann and Kaercher papers 200?
      PRE_ICE(JL,JK)=(MAX(SQRT(5.113E6_JPRB+2.809E3_JPRB*PRE_ICE(JL,JK)**3.0_JPRB)-2.261E3_JPRB,0.0_JPRB))**0.333_JPRB
      PRE_ICE(JL,JK)=MAX(PRE_ICE(JL,JK),1.0_JPRB)  ! diameter minimum 1.0 microns

! more default values if not applying
      PNICE(JL,JK)=RNICE ! place as default

      IF (PT(JL,JK)<238._JPRB .AND. PI(JL,JK)>ZEPSEC) THEN
        ZS0=1.3_JPRB
        ZSCRITHOMO=2.349_JPRB-PT(JL,JK)/259.0_JPRB !ren form of Koop 2000 
        ZSVP=MAX(ZEPSEC,ZQS(JL,JK)*PAP(JL,JK)/0.622_JPRB)

! Klaus Gierens critical ice nuclei: Gierens (2003)
        ZNCRIT_GIERENS=2.81E11_JPRB*(10.0_JPRB**(4.0_JPRB-0.02_JPRB*PT(JL,JK)))**0.75_JPRB&
       &*(ZWTOT**1.5_JPRB)*PAP(JL,JK)**1.5_JPRB/&
       &(PT(JL,JK)**5.415_JPRB*(1.5_JPRB*ZSVP)**0.5_JPRB*(ZSCRITHOMO-ZS0)**0.75_JPRB)
        ZNCRIT_GIERENS=ZNCRIT_GIERENS/1.E6_JPRB ! cm**-3

! Ren and Mackensie QJRMS 2005 critical ice nuclei
        ZNCRIT_REN=5.4E10_JPRB*(ZWTOT**1.5_JPRB)*PAP(JL,JK)**1.5_JPRB*&
       & (ZSCRITHOMO/(ZSCRITHOMO-1.0_JPRB))**1.5_JPRB/ &
       & (PT(JL,JK)**5.415_JPRB*(1.5_JPRB*ZSVP)**0.5_JPRB)
        ZNCRIT_REN=ZNCRIT_REN/1.E6_JPRB ! cm**-3

!        IF (KFLDX>0 .AND. LDUMPEXTRA) PEXTRA(JL,JK,13)=ZNCRIT_GIERENS
!        IF (KFLDX>0 .AND. LDUMPEXTRA) PEXTRA(JL,JK,14)=ZNCRIT_REN

! from Re derive the number concentration - here ice density is 900 kg/m**3 
! Re is in microns, 1e18 factor
        ZCLD=PI(JL,JK)/MAX(PA(JL,JK),ZEPSEC)
        ZCLD=MIN(MAX(ZCLD,0.0_JPRB),RCLDMAX)
        IF (ZCLD>ZEPSEC) THEN
          ZNICEHOMO=0.75_JPRB*ZRHO(JL,JK)*ZCLD/(RPI*ZRHO_ICE*1.0E-18_JPRB*PRE_ICE(JL,JK)**3.0_JPRB)
        ENDIF
        ZNICEHOMO = ZNICEHOMO/1.E6_JPRB ! cm**-3
!        IF (KFLDX>0 .AND. LDUMPEXTRA) PEXTRA(JL,JK,16)=ZNICEHOMO

! following Ren and Mackensie, 2005, linearly interpolate to get Ice number
        IF (ZICENUCLEI(JL,JK)<ZNCRIT_REN) THEN
          PNICE(JL,JK)=ZICENUCLEI(JL,JK)+(1.0_JPRB-ZICENUCLEI(JL,JK)/ZNCRIT_REN)*ZNICEHOMO
        ELSE
          PNICE(JL,JK)=ZICENUCLEI(JL,JK) 
        ENDIF

! number is 3/(4*pi*rho_liq*10^6)  [10^6 for N in cm**-3]
        PRE_ICE(JL,JK)=(0.75_JPRB*ZRHO(JL,JK)*ZCLD/(RPI*ZRHO_ICE*1.E6_JPRB*PNICE(JL,JK)))**0.333_JPRB
        PRE_ICE(JL,JK)=PRE_ICE(JL,JK)*1.E6_JPRB
!        IF (KFLDX>0 .AND. LDUMPEXTRA) PEXTRA(JL,JK,15)=PRE_ICE(JL,JK)

!---------------------------------------------------------------------
! Turn IN Number concentration for ice into critical 
! mixing ratio for autoconversion process
! Obvious, but also from Rotstayn and penner JAS 2001
!---------------------------------------------------------------------

! Number is 1.e6. 4/3 .pi.rho_l. r_crit**3
! 1e6 since IN in units of  N cm**-3
! r_crit=9.3 microns
        PICRIT_AER(JL,JK)=1.333E6_JPRB*RPI*ZRHO_ICE*PNICE(JL,JK)*ZRICE_CRIT**3.0_JPRB/ZRHO(JL,JK)
! limit the effect to ratio of the "background" value
        PICRIT_AER(JL,JK)=MAX(PICRIT_AER(JL,JK),0.1_JPRB*RLCRITSNOW)
        PICRIT_AER(JL,JK)=MIN(PICRIT_AER(JL,JK),10.0_JPRB*RLCRITSNOW)
      ELSE ! T>233K
        PICRIT_AER(JL,JK)=RLCRITSNOW
      ENDIF
!      IF (KFLDX>0 .AND. LDUMPEXTRA) PEXTRA(JL,JK,17)=PICRIT_AER(JL,JK)
    ENDDO
  ENDDO
ENDIF

! stop crash with excessive values?
!IF (KFLDX>0 .AND. LDUMPEXTRA) THEN 
!DO JF=1,KFLDX
!  DO JK=1,KLEV
!    DO JL=KIDIA,KFDIA
!      PEXTRA(JL,JK,JF)=MAX(PEXTRA(JL,JK,JF),-1.E32_JPRB)
!      PEXTRA(JL,JK,JF)=MIN(PEXTRA(JL,JK,JF),1.E32_JPRB)
!    ENDDO
!  ENDDO
!ENDDO
!ENDIF

!----------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('AER_CLCLD',1,ZHOOK_HANDLE)
END SUBROUTINE AER_CLCLD
