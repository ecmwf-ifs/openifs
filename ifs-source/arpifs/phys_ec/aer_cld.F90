! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE AER_CLD ( & 
!  input
 & YDECLDP,KIDIA , KFDIA, KLON,  KLEV, &
 & KACTAERO, &
 & PAERO , PT , &
 & PAP  , &
 & PL    , PI   , PA   , &
!  output
 & PLCRIT_AER, PICRIT_AER, &
 & PRE_LIQ   , PRE_ICE, &
 & PCCN      , PNICE &
 & )    

!     CLOUDAER 
!     --------
!          Sets up all information for subsequent calculation 
!          of the effect of prognostic aerosols on clouds and convection
!          (2nd indirect effect)

!     AUTHOR
!          A. Tompkins/JJMorcrette  E.C.M.W.F.

!     PURPOSE.
!     --------

!     INTERFACE.
!     ----------

!          *AER_CLD* IS CALLED FROM *CALLPAR*

!     PARAMETER     DESCRIPTION                                   UNITS
!     ---------     -----------                                   -----

! -   INPUT ARGUMENTS.
!     -------------------

! KIDIA    : START OF HORIZONTAL LOOP
! KFDIA    : END   OF HORIZONTAL LOOP
! KLON     : HORIZONTAL DIMENSION
! KLEV     : END OF VERTICAL LOOP AND VERTICAL DIMENSION
! KACTAERO : NUMBER OF ACTIVE PROGNOSTIC AEROSOLS

! -   OUTPUT ARGUMENTS.
!     -------------------
! PLCRIT_AER : critical liquid mmr for autoconversion process
! PICRIT_AER : critical liquid mmr for autoconversion process
! PRE_LIQ : liq Re
! PRE_ICE : ice Re
! PCCN    : liquid cloud condensation nuclei
! PNICE   : ice number concentration (cf. CCN)


USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK    ,DR_HOOK, JPHOOK

USE YOMCST   , ONLY : RD       ,RETV     ,&
 & RLVTT    ,RLSTT    ,RTT     ,RPI
USE YOETHF   , ONLY : R2ES     ,R3LES    ,R3IES    ,R4LES    ,&
 & R4IES    ,R5LES    ,R5IES    ,R5ALVCP  ,R5ALSCP  ,&
 & RALVDCP  ,RALSDCP  ,RTWAT    ,&
 & RTICE    ,RTICECU  ,&
 & RTWAT_RTICE_R      ,RTWAT_RTICECU_R,&
 & RKOOP1   ,RKOOP2
USE YOECLDP  , ONLY : TECLDP

IMPLICIT NONE

! input variables
TYPE(TECLDP)      ,INTENT(INOUT) :: YDECLDP
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KACTAERO
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAERO(KLON,KLEV,KACTAERO)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAP(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PT(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PL(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PI(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PA(KLON,KLEV) 

! output
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PLCRIT_AER(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PICRIT_AER(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PRE_LIQ(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PRE_ICE(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PCCN(KLON,KLEV)     ! liquid cloud condensation nuclei
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PNICE(KLON,KLEV)    ! ice number concentration (cf. CCN)

! diagnostics
!!INTEGER(KIND=JPIM),INTENT(IN)    :: KFLDX
!!REAL(KIND=JPRB)   ,INTENT(INOUT) :: PEXTRA(KLON,KLEV,KFLDX) 
!LOGICAL, PARAMETER :: LDUMPEXTRA=.FALSE.

!----------------------------------------------------------------------

! local arrays

! Aerosol arrays
REAL(KIND=JPRB) :: ZMAER(KLON,KLEV,6)    ! mass of aerosol   [kg kg-1]
REAL(KIND=JPRB) :: ZMAERMN(6)            ! annual column mean mass of aerosol
REAL(KIND=JPRB) :: ZICENUCLEI(KLON,KLEV) ! number concentration of ice nuclei
REAL(KIND=JPRB) :: ZQS(KLON,KLEV)        ! saturation
REAL(KIND=JPRB) :: ZCCN0(KLON)

REAL(KIND=JPRB) :: ZEXPN, ZRELVNT, ZS0, ZSCRITHOMO, ZSVP, ZTEMP, ZTEMPC
REAL(KIND=JPRB) :: ZNCRIT_GIERENS, ZNCRIT_REN
REAL(KIND=JPRB) :: ZNICEHOMO

REAL(KIND=JPRB) :: ZCLD
REAL(KIND=JPRB) :: ZRHO_ICE, ZRHO_LIQ ! density of pristine ice crystals and cloud
REAL(KIND=JPRB) :: ZRLIQ_CRIT, ZRICE_CRIT ! critical radii for autoconversion process

INTEGER(KIND=JPIM) :: JAE, JAERSS, JAERDU, JAEROM, JAERSU, JAERBC

REAL(KIND=JPRB) :: ZWTOT

! general arrays
REAL(KIND=JPRB) :: ZRHO(KLON,KLEV)            ! density
REAL(KIND=JPRB) :: ZAERO(KLON,KLEV,KACTAERO)  ! protected aerosol amounts

! misc variables
REAL(KIND=JPRB) :: ZEPSEC, ZEPSTAU
INTEGER(KIND=JPIM) :: ISRCCN
INTEGER(KIND=JPIM) :: JK, JL

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!------------------------

#include "fcttre.func.h"
#include "fccld.func.h"

!--------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('AER_CLD',0,ZHOOK_HANDLE)
ASSOCIATE(LAERICEAUTO=>YDECLDP%LAERICEAUTO, LAERICESED=>YDECLDP%LAERICESED, &
 & LAERLIQAUTOLSP=>YDECLDP%LAERLIQAUTOLSP, LAERLIQCOLL=>YDECLDP%LAERLIQCOLL, &
 & NAECLBC=>YDECLDP%NAECLBC, NAECLDU=>YDECLDP%NAECLDU, NAECLOM=>YDECLDP%NAECLOM, &
 & NAECLSS=>YDECLDP%NAECLSS, NAECLSU=>YDECLDP%NAECLSU, RCCNOM=>YDECLDP%RCCNOM, &
 & RCCNSS=>YDECLDP%RCCNSS, RCCNSU=>YDECLDP%RCCNSU, RCLCRIT=>YDECLDP%RCLCRIT, &
 & RCLDMAX=>YDECLDP%RCLDMAX, RLCRITSNOW=>YDECLDP%RLCRITSNOW, &
 & RNICE=>YDECLDP%RNICE)
!######################################################################
!                       1.0 Basic variables
!######################################################################
!!! IF (KFLDX>0 .AND. LDUMPEXTRA) PEXTRA(:,:,:)=0.0_JPRB

ZEPSEC =1.E-10_JPRB
ZEPSTAU=1.E-20_JPRB
DO JAE=1,KACTAERO
  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      ZAERO(JL,JK,JAE)=MAX(ZEPSTAU,PAERO(JL,JK,JAE))
    ENDDO
  ENDDO
ENDDO

! move to cldp module
ZRHO_ICE=900.0_JPRB
ZRHO_LIQ=1000.0_JPRB
ZRICE_CRIT=60.0E-6_JPRB ! ice to snow critical radius
ZRLIQ_CRIT=9.3E-6_JPRB  ! cloud to rain critical radius
ZWTOT=0.1_JPRB ! governs critical N

!-- saturation mixing ratio [kg kg-1] and and density [kg m-3]
DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    ZQS(JL,JK)=FOEEWM(PT(JL,JK))/PAP(JL,JK)
    ZQS(JL,JK)=MIN(0.5_JPRB,ZQS(JL,JK))
    ZQS(JL,JK)=ZQS(JL,JK)/(1.0_JPRB-RETV*ZQS(JL,JK))
    ZRHO(JL,JK)=PAP(JL,JK)/(RD*PT(JL,JK))
  ENDDO
ENDDO

!######################################################################

JAERSS=1
JAEROM=2
JAERBC=3
JAERSU=4
JAERDU=5

ZMAERMN(JAERSS)=2.12E-10_JPRB*1.E9_JPRB 
ZMAERMN(JAERDU)=1.01E-09_JPRB*1.E9_JPRB
ZMAERMN(JAEROM)=3.05E-11_JPRB*1.E9_JPRB
ZMAERMN(JAERBC)=3.05E-11_JPRB*1.E9_JPRB
ZMAERMN(JAERSU)=1.02E-09_JPRB*1.E9_JPRB

!-- Prognostic aerosols are entered as:
!   PAERO
!     1       1- 3  sea-salt  0.03 - 0.5 -  5  - 20 microns
!     2       4- 6  dust      0.03 - 0.5 - 0.9 - 20 microns
!     3       7- 8  POM     hydrophilic, hydrophobic
!     4       9-10  BC      hydrophilic, hydrophobic
!     5      11-12  SO4/SO2 including sulfate prognostic stratospheric aerosols (SO4 is 11)

!  What is used is: 
! JAERSS=1  SS1 (0.03-0.5um)     i.e. NAECLSS=1
! JAEROM=2  POM1 (hydrophilic)   i.e. NAECLOM=7
! JAERBC=3  BC1 (hydrophilic)    i.e. NAECLBC=9
! JAERSU=4  SO4                  i.e. NAECLSU=11 
! JAERDU=5  DU1 (0.03-0.55um)    i.e. NAECLDU=4
 
! ... in principle it should be the coarse dust, i.e., DU3 (0.9-20um), i.e., ZAERO(..,.., 6)

ZMAER(:,:,:) = 0.0_JPRB

DO JK=1,KLEV
  DO JL=KIDIA,KFDIA

!  PAERO are entered in kg/kg
!  ZMAER are in ug / m3              ( micrograms per m**3 )
    ZMAER(JL,JK,1) = ZAERO(JL,JK,NAECLSS)*ZRHO(JL,JK)*1.E9_JPRB    ! sea-salt
    ZMAER(JL,JK,2) = ZAERO(JL,JK,NAECLOM)*ZRHO(JL,JK)*1.E9_JPRB    ! org.matter
    ZMAER(JL,JK,3) = ZAERO(JL,JK,NAECLBC)*ZRHO(JL,JK)*1.E9_JPRB    ! black carbon
    ZMAER(JL,JK,4) = ZAERO(JL,JK,NAECLSU)*ZRHO(JL,JK)*1.E9_JPRB    ! sulphate
    ZMAER(JL,JK,5) = ZAERO(JL,JK,NAECLDU)*ZRHO(JL,JK)*1.E9_JPRB    ! dust
  ENDDO
ENDDO

!######################################################################
!                4. WARM PHASE MICROPHYSICS
!######################################################################

IF (LAERLIQAUTOLSP.OR.LAERLIQCOLL) THEN

!---------------------------------------------------------------------
! Turn aerosol mass into a CCN Number concentration for warm rain
! From Menon et al, 2002: JAS, 59, 692-713  Eqns 1a, 1b
!---------------------------------------------------------------------
  PCCN(KIDIA:KFDIA,1:KLEV)=32.0_JPRB
  DO JK=1,KLEV
!DEC$ IVDEP
    DO JL=KIDIA,KFDIA

!     PCCN = N in cm**-3

      ZRELVNT=ZMAER(JL,JK,JAERSS)+ZMAER(JL,JK,JAEROM)+ZMAER(JL,JK,JAERSU)

! CCNxx factor average for org.matter 0.13, sea salt 0.05, and sulphate 0.50
      ZEXPN = RCCNOM*LOG10(MAX(ZMAER(JL,JK,JAEROM),ZEPSEC)) +&
            & RCCNSS*LOG10(MAX(ZMAER(JL,JK,JAERSS),ZEPSEC)) +&
            & RCCNSU*LOG10(MAX(ZMAER(JL,JK,JAERSU),ZEPSEC))
      ZCCN0(JL)=10.0_JPRB**(2.41_JPRB + ZEXPN)

      IF (ZRELVNT > 0._JPRB) THEN
!-- NB: no bounds: if either SS, OM and/or SU aerosols are present, just use the diagnosed CCN
        PCCN(JL,JK) = ZCCN0(JL)
      ENDIF

! this is effective radius calculation - not used for now
      IF (PA(JL,JK) >= 0.001_JPRB) THEN
        ZTEMP=1.0_JPRB/PA(JL,JK)
        ZCLD=MAX(0.0_JPRB, PL(JL,JK)*ZTEMP)
      ELSE
        ZCLD=0.0_JPRB
      ENDIF
! number is 3/(4*pi*rho_liq*10^6)  [10^6 for N in right units]
      PRE_LIQ(JL,JK)=(2.387E-10_JPRB*ZRHO(JL,JK)*ZCLD/PCCN(JL,JK))**0.333_JPRB
    ENDDO
  ENDDO

!---------------------------------------------------------------------
! Turn CCN Number concentration for warm rain into critical 
! mixing ratio for autoconversion process
! Obvious, but also from Rotstayn, L. and J.E. Penner, 2001: JAS 2001
!---------------------------------------------------------------------

! Number is 1.e6. 4/3 .pi.rho_l. r_crit**3
! 1e6 since CCN in units of  N cm**-3
! r_crit=9.3 microns

  ISRCCN=2

  DO JK=1,KLEV
!DEC$ IVDEP
    DO JL=KIDIA,KFDIA
      PLCRIT_AER(JL,JK)=1.333E6_JPRB*RPI*ZRHO_LIQ*PCCN(JL,JK)*ZRLIQ_CRIT**3.0_JPRB/ZRHO(JL,JK)

      IF (ISRCCN == 0) THEN
!-- PLCRIT_AER is to be used as computed above
      ELSEIF (ISRCCN == 1) THEN         
      ELSEIF (ISRCCN == 2) THEN
! limit the effect to ratio of the "background" value
        PLCRIT_AER(JL,JK)=MAX(PLCRIT_AER(JL,JK),0.1_JPRB*RCLCRIT)
        PLCRIT_AER(JL,JK)=MIN(PLCRIT_AER(JL,JK),10.0_JPRB*RCLCRIT)
      ENDIF
    ENDDO
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
! Demott et al., 1994: JAS, 51, 77-90   -->  ISS = 25 %
! Demott et al., 1997: JGR, 102, 19575-19584
! Meyers et al., 1992: JAM, 31, 708-721  --> ISS = 25 %

! Demott et al. Ice SS=55% or Meyers et al. 1992 JAM, ISS=25%
! In a prognostic scheme this will be function of clear sky humidity

! By relating IN to Aerosol mass we are assuming that the mode of the 
! aerosol size distribution lies in the accumulation or coarse mode

! The relationship will implicitly introduce the exponential height
! dependence that Sassen (1992) and K and Curry (1998) explicitly 
! introduced to their parametrizations.
! Sassen, 1992: 
! Khvorostyanov and Curry, 2000: GRL 27, 4081-4084.

      ZICENUCLEI(JL,JK)=0.01_JPRB*&
    &  (ZMAER(JL,JK,JAERSU)+ZMAER(JL,JK,JAERBC)+ZMAER(JL,JK,JAERDU))&
    & /(ZMAERMN(JAERSU)    +ZMAERMN(JAERBC)    +ZMAERMN(JAERDU))
      ZICENUCLEI(JL,JK)=MAX(ZICENUCLEI(JL,JK),0.0_JPRB)
!      IF (KFLDX>0 .AND. LDUMPEXTRA) PEXTRA(JL,JK,12)=ZICENUCLEI(JL,JK)

! T in oC
      ZTEMPC=PT(JL,JK)-RTT

! Re form for ice crystals from Liou and Oort 1994    (Ou and Liou, 1995, Atmos. Res., 35, 127-138)???
! used to derive Re(ice) as in Lohmann, 2002: JAS, 59, 647-656. Eqn 5 
      ZNICEHOMO=0.0_JPRB
      PRE_ICE(JL,JK)=0.5_JPRB*(326.3_JPRB+ZTEMPC*&
        & (12.42_JPRB + ZTEMPC*(0.197_JPRB + ZTEMPC*0.0012_JPRB)))
      PRE_ICE(JL,JK)=MAX(PRE_ICE(JL,JK),0.0_JPRB)

! effect Re to volume mean from S Moss or Lohmann and Kaercher papers 200? 2002, JGR, 107B, doi:10.1029/2001JD000767.
      PRE_ICE(JL,JK)=(MAX(SQRT(5.113E6_JPRB+2.809E3_JPRB*PRE_ICE(JL,JK)**3.0_JPRB)-2.261E3_JPRB,0.0_JPRB))**0.333_JPRB
      PRE_ICE(JL,JK)=MAX(PRE_ICE(JL,JK),1.0_JPRB)  ! diameter minimum 1.0 microns

! more default values if not applying
      PNICE(JL,JK)=RNICE ! place as default

      IF (PT(JL,JK)<238._JPRB .AND. PI(JL,JK)>ZEPSEC) THEN
        ZS0=1.3_JPRB
        ZSCRITHOMO=2.349_JPRB-PT(JL,JK)/259.0_JPRB !ren form of Koop 2000 
        ZSVP=MAX(ZEPSEC,ZQS(JL,JK)*PAP(JL,JK)/0.622_JPRB)

! Klaus Gierens critical ice nuclei: Gierens, 2003: ACP, 3, 437-446.
        ZNCRIT_GIERENS=2.81E11_JPRB*(10.0_JPRB**(4.0_JPRB-0.02_JPRB*PT(JL,JK)))**0.75_JPRB&
       &*(ZWTOT**1.5_JPRB)*PAP(JL,JK)**1.5_JPRB/&
       &(PT(JL,JK)**5.415_JPRB*(1.5_JPRB*ZSVP)**0.5_JPRB*(ZSCRITHOMO-ZS0)**0.75_JPRB)
        ZNCRIT_GIERENS=ZNCRIT_GIERENS/1.E6_JPRB ! cm**-3

! Ren and Mackensie, 2005: QJRMS, 131B, 1585-1605: critical ice nuclei
        ZNCRIT_REN=5.4E10_JPRB*(ZWTOT**1.5_JPRB)*PAP(JL,JK)**1.5_JPRB*&
       & (ZSCRITHOMO/(ZSCRITHOMO-1.0_JPRB))**1.5_JPRB/&
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

! following Ren and Mackensie, 2005, QJ 131, linearly interpolate to get Ice number
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
! Obvious, but also from Rotstayn and Penner, 2001: JClimate, 14, 2960-2975.
!---------------------------------------------------------------------

! Number is 1.e6. 4/3 .pi.rho_l. r_crit**3
! 1e6 since IN in units of  N cm**-3
! r_crit=60.0 microns
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
IF (LHOOK) CALL DR_HOOK('AER_CLD',1,ZHOOK_HANDLE)
END SUBROUTINE AER_CLD
