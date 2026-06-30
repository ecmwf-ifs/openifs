MODULE SUSURB_MOD
CONTAINS
!SUBROUTINE SUSURB(LD_LURBAN,LD_LURBUI,YDURB)
SUBROUTINE SUSURB(YDURB)
    ! (C) Copyright 2021- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

!**   *SUSURB* IS THE SET-UP ROUTINE FOR COMMON BLOCK *YOS_URB*

!     PURPOSE
!     -------
!          THIS ROUTINE INITIALIZES THE CONSTANTS IN COMMON BLOCK
!     *YOS_URB*

!     INTERFACE.
!     ----------
!     CALLLED FROM *SUSURF*

!     METHOD.
!     This routine sets up variables and, calculates canyon emissivity,
!     canyon roughness length/ heat coefficient/ resistance and
!     roof roughness length/ heat coefficient/ resistance
!     canyon albedo is calculated at each timestep to account for S.Z.angle

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     Urban scheme based on Macdonald 1998, Harman 2004, Oleson 2008 and Porson 2010

!     MODIFICATIONS
!     -------------
!==============================================================================

! Modules used : 
USE PARKIND1  , ONLY : JPIM, JPRB
USE YOMHOOK   , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOS_URB   , ONLY : TURB

IMPLICIT NONE

TYPE(TURB),     INTENT(INOUT)   :: YDURB

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!REAL (KIND = JPRB):: RHWR           ! Height-to-roadwidth (aspect) ratio
!REAL (KIND = JPRB):: RHGT           ! Average building height (m)
!REAL (KIND = JPRB):: RWRR           ! Road-to-building+road ratio (width)
!REAL (KIND = JPRB):: RCANE          ! Canyon Emissivity
!REAL (KIND = JPRB):: RCANZTM        ! Canyon Roughness Length for Momentum
!REAL (KIND = JPRB):: RCANZTH        ! Canyon Roughness Length for Heat
!REAL (KIND = JPRB):: RCANHC         ! Canyon Heat Coefficient
!REAL (KIND = JPRB):: RROOZTM        ! Roof Roughness Length for Momentum
!REAL (KIND = JPRB):: RROOZTH        ! Roof Roughness Length for Heat
!REAL (KIND = JPRB):: RROOHC         ! Roof Heat Coefficient

REAL (KIND = JPRB)::aa
REAL (KIND = JPRB)::bb
REAL (KIND = JPRB)::cc
REAL (KIND = JPRB)::dd
REAL (KIND = JPRB)::RMODW
INTEGER           ::i
INTEGER           ::j

INTEGER, PARAMETER:: matdim = 3
REAL (KIND = JPRB):: shf(matdim,matdim)         ! shape factor
REAL (KIND = JPRB):: delta(matdim,matdim)       ! diagonal
REAL (KIND = JPRB):: B(matdim,matdim)       ! diagonal

REAL (KIND = JPRB):: shfinv(matdim,matdim)   !! Inverse matrix
REAL (KIND = JPRB):: detinv

REAL (KIND = JPRB):: sc_hwr    ! Scaled HWR for all directions
REAL (KIND = JPRB):: d_h
REAL (KIND = JPRB):: disp      ! Displacement height
REAL (KIND = JPRB):: zref      ! Reference height 
REAL (KIND = JPRB):: pi
REAL (KIND = JPRB):: l_rec     ! Recirc. region length 
REAL (KIND = JPRB):: width     ! Canyon Width

REAL (KIND = JPRB):: lrw       ! rescaling resistance over edge recirc.
REAL (KIND = JPRB):: ltw       ! rescaling resistance over edge recirc.
REAL (KIND = JPRB):: hrh       ! rescaling resistance over edge recirc.
REAL (KIND = JPRB):: lttw      ! rescaling resistance over edge recirc.
REAL (KIND = JPRB):: swi
REAL (KIND = JPRB):: uct       ! wind with random orientation
REAL (KIND = JPRB):: can_road  ! wind travelling along road
REAL (KIND = JPRB):: can_wall  ! wind travelling along wall
REAL (KIND = JPRB):: l_slo     ! Length of sloping edge 
REAL (KIND = JPRB):: l_dow     ! Length of downstream wall

REAL (KIND = JPRB):: fdr       ! Flow downstream road
REAL (KIND = JPRB):: fur       ! Flow upstream road
REAL (KIND = JPRB):: fuw       ! Flow upstream wall
REAL (KIND = JPRB):: fdw       ! Flow downstream wall
REAL (KIND = JPRB):: low       ! Lower limit of flow downstream wall
REAL (KIND = JPRB):: z0x       ! Scaled circulation property Nakamura and Oke, 1988

REAL (KIND = JPRB):: r_3       ! Resistance recirc. wall
REAL (KIND = JPRB):: r_4       ! Resistance recirc. road
REAL (KIND = JPRB):: r_5       ! Resistance recirc. internal
REAL (KIND = JPRB):: r_6       ! Resistance vent. road
REAL (KIND = JPRB):: r_7       ! Resistance vent. internal
REAL (KIND = JPRB):: r_8       ! Resistance vent. wall

REAL (KIND = JPRB):: r_recirc  ! Resistance vent. road
REAL (KIND = JPRB):: r_venti   ! Resistance vent. internal
REAL (KIND = JPRB):: r_bulk    ! Resistance vent. wall

IF (LHOOK) CALL DR_HOOK('SUSURB_MOD:SUSURB',0,ZHOOK_HANDLE)
ASSOCIATE( &
 & RBUIZ0M=>YDURB%RBUIZ0M,RWALTHK=>YDURB%RWALTHK,RROOTHK=>YDURB%RROOTHK, &
 & RROATHK=>YDURB%RROATHK,RWALALB=>YDURB%RWALALB,RROOALB=>YDURB%RROOALB, &
 & RROAALB=>YDURB%RROAALB,RWALEMIS=>YDURB%RWALEMIS,RROOEMIS=>YDURB%RROOEMIS, &
 & RROAEMIS=>YDURB%RROAEMIS,RSOIEMIS=>YDURB%RSOIEMIS,RWALVHC=>YDURB%RWALVHC, &
 & RROOVHC=>YDURB%RROOVHC,RROAVHC=>YDURB%RROAVHC,RAIRVHC=>YDURB%RAIRVHC, &
 & RSOITMP=>YDURB%RSOITMP,RSOITC=>YDURB%RSOITC,RWALTC=>YDURB%RWALTC, &
 & RROOTC=>YDURB%RROOTC,RROATC=>YDURB%RROATC,RSATSH=>YDURB%RSATSH, &
 & RCDG=>YDURB%RCDG,RCDA=>YDURB%RCDA,RCHIS=>YDURB%RCHIS, &
 & RSBCONS=>YDURB%RSBCONS,RGACC=>YDURB%RGACC,RAIRRHO=>YDURB%RAIRRHO, &
 & RMODZ=>YDURB%RMODZ,RVKSQ=>YDURB%RVKSQ,REXPDR=>YDURB%REXPDR, &
 & RHWR=>YDURB%RHWR,RHGT=>YDURB%RHGT,RWRR=>YDURB%RWRR,RCANEMIS=>YDURB%RCANEMIS, &
 & RURBEMIS=>YDURB%RURBEMIS,RCANZTM=>YDURB%RCANZTM,RCANZTH=>YDURB%RCANZTH, &
 & RCANHC=>YDURB%RCANHC,RCANRES=>YDURB%RCANRES,RROORES=>YDURB%RROORES,&
 & RROOZTM=>YDURB%RROOZTM,RROOZTH=>YDURB%RROOZTH,RROOHC=>YDURB%RROOHC,&
 & RURBRES=>YDURB%RURBRES,RURBZTM=>YDURB%RURBZTM,RURBZTH=>YDURB%RURBZTH,&
 & RURBHC=>YDURB%RURBHC,RCANTC=>YDURB%RCANTC,RCANVHC=>YDURB%RCANVHC,&
 & RURBTC=>YDURB%RURBTC,RURBVHC=>YDURB%RURBVHC, RURBALP=>YDURB%RURBALP,&
 & RURBCON=>YDURB%RURBCON,RURBLAM=>YDURB%RURBLAM,RURBSAT=>YDURB%RURBSAT,&
 & RURBSRES=>YDURB%RURBSRES, RURBTC1=>YDURB%RURBTC1)

!==============================================================================
 
! Defined constants
 
RCDG    = 1.2_JPRB              ! Drag Coefficient (Macdonald 1998)
RCDA    = 4.43_JPRB             ! Coefficient For Canyon (alpha) (Macdonald 1998)
RCHIS   = 0.3_JPRB              ! Fraction of SW radiation scatter by sky
RSBCONS = 5.67E+08_JPRB         ! Stefan-Boltzmann Constant (W m-2 K-4)
RGACC   = 9.81_JPRB             ! Acceleration due to gravity (m s-2)
RVKSQ   = 0.16_JPRB             ! Von Karman constant squared
REXPDR  = 0.15_JPRB             ! Exponential decay of recirculation

! Material properties

RBUIZ0M = 0.005_JPRB            ! Rough. len. for mom. of building materials
RWALTHK = 0.15_JPRB             ! Thickness of wall (m) (not currently used)
RROOTHK = 0.15_JPRB             ! Thickness of roof (m) (not currently used)
RROATHK = 0.15_JPRB             ! Thickness of road (m) (not currently used)
RWALALB = 0.6_JPRB              ! Albedo of wall
RROOALB = 0.16_JPRB             ! Albedo of roof
RROAALB = 0.05_JPRB             ! Albedo of road
RWALEMIS= 0.96_JPRB             ! Emissivity of wall
RROOEMIS= 0.96_JPRB             ! Emissivity of roof
RROAEMIS= 0.99_JPRB             ! Emissivity of road
RWALVHC = 3.0E+06_JPRB          ! Volumetric heat capcity of wall (J m-3 K-1)
RROOVHC = 2.0E+06_JPRB          ! Volumetric heat capcity of roof (J m-3 K-1)
RROAVHC = 3.0E+06_JPRB          ! Volumetric heat capcity of road (J m-3 K-1)
RWALTC  = 20.0_JPRB             ! Thermal conductivity of wall (W m-1 K-1)(skin thick inc.*0.07)
RROOTC  = 10.0_JPRB             ! Thermal conductivity of roof (W m-1 K-1)(replace 0.035)
RROATC  = 20.0_JPRB             ! Thermal conductivity of road (W m-1 K-1)(with dz-1)

! Mapped variables which are currently set as a single value
RHGT    = 8.0_JPRB             ! Average building height (m) (estimated. from london see UK https://buildingheights.emu-analytics.net)
RHWR    = 0.5_JPRB              ! Height-to-roadwidth (aspect) ratio
RWRR    = 0.5_JPRB              ! Road-to-building+road ratio (width)

RSOIEMIS= 0.95_JPRB             ! Emissivity of soil (not currently used)
RSOITMP = 280.0_JPRB            ! Constant sub-soil temperature (K) (not currently used)
RSOITC  = 0.22_JPRB             ! Thermal conductivity of soil (W m-1 K-1) (not currently used)
RSATSH  = 0.00_JPRB             ! Saturated specific humidity (kg kg-1) (not currently used)

!RFLOTHK = 1.00_JPRB             ! Thickness of floor (m) (not currently used)
!RFLOEMIS= 0.95_JPRB             ! Floor emissivity (not currently used)
!RFLOVHC = 1.37E+06_JPRB         ! Volumetric heat capactiy of floor (K m-3 K-1) (not currently used)

! Atmospheric properties (can be taken from model)
RMODZ   = 10.0_JPRB               ! Lowest model height (m)
RMODW   = 10.0_JPRB              ! Lowest model height - wind (m)
RAIRRHO = 1.225_JPRB             ! Density of air (Kg m-3)
RAIRVHC = RAIRRHO*1.005E+03_JPRB ! Heat capactiy of (moist) air (K m-3K-1) - Check units/values

! Switches
!LURBUI = LD_LURBUI              ! Building heat ON/OFF (not currently used)

! MOISTURE VARIABLES

RURBALP  = 7.65E-01_JPRB 
RURBCON  = 6.67E-13_JPRB
RURBLAM  = 35.2_JPRB
RURBSAT  = 0.15_JPRB
RURBSRES = 0.001_JPRB
!==============================================================================

! CALCULATE CANYON EMISSIVITY BASED ON HWR - Based on Harman 2004 and Porson 2010

! Exchange matrix

 aa = (1.0_JPRB + RHWR**2.0_JPRB)**0.5_JPRB - RHWR                      ! sky-road
 bb = (1.0_JPRB + (1.0_JPRB/RHWR)**2.0_JPRB)**0.5_JPRB - 1.0_JPRB/RHWR  ! wall-wall
 cc = (1.0_JPRB - aa)                                                   ! road-wall
 dd = (1.0_JPRB - bb)/2.0_JPRB                                          ! wall-sky

 shf(1,1) = 0.0_JPRB                      ! 1 - Sky, 2 - Wall, 3 - Road
 shf(1,2) = cc
 shf(1,3) = aa 
 shf(2,1) = dd
 shf(2,2) = bb
 shf(2,3) = dd
 shf(3,1) = aa
 shf(3,2) = cc
 shf(3,3) = 0.0_JPRB

  ! Exchange
  
 DO i = 1,matdim
  shf(1,i) = shf(1,i)*(1.0_JPRB - RROAEMIS)
  shf(2,i) = shf(2,i)*(1.0_JPRB - RWALEMIS)
  shf(3,i) = 0.0_JPRB
 ENDDO
  
 DO i = 1,matdim
  DO j = 1,matdim
   IF (i == j) THEN
    delta(i,j) = 1.0_JPRB
   ELSE
    delta(i,j) = 0.0_JPRB
   END IF
  END DO
 END DO

 shf = delta-shf
 
!Invert shf

!Calculate the inverse determinant of the matrix
 detinv = 1/(shf(1,1)*shf(2,2)*shf(3,3) - shf(1,1)*shf(2,3)*shf(3,2)&
          - shf(1,2)*shf(2,1)*shf(3,3) + shf(1,2)*shf(2,3)*shf(3,1)&
          + shf(1,3)*shf(2,1)*shf(3,2) - shf(1,3)*shf(2,2)*shf(3,1))

! Calculate the inverse of the matrix
 B(1,1) = +detinv * (shf(2,2)*shf(3,3) - shf(2,3)*shf(3,2))
 B(2,1) = -detinv * (shf(2,1)*shf(3,3) - shf(2,3)*shf(3,1))
 B(3,1) = +detinv * (shf(2,1)*shf(3,2) - shf(2,2)*shf(3,1))
 B(1,2) = -detinv * (shf(1,2)*shf(3,3) - shf(1,3)*shf(3,2))
 B(2,2) = +detinv * (shf(1,1)*shf(3,3) - shf(1,3)*shf(3,1))
 B(3,2) = -detinv * (shf(1,1)*shf(3,2) - shf(1,2)*shf(3,1))
 B(1,3) = +detinv * (shf(1,2)*shf(2,3) - shf(1,3)*shf(2,2))
 B(2,3) = -detinv * (shf(1,1)*shf(2,3) - shf(1,3)*shf(2,1))
 B(3,3) = +detinv * (shf(1,1)*shf(2,2) - shf(1,2)*shf(2,1))
 
 RCANEMIS = ((RROAEMIS/(1.0_JPRB-RROAEMIS))*B(1,3)) + (2.0_JPRB*RHWR*(RWALEMIS/(1.0_JPRB-RWALEMIS))*B(2,3))

 RURBEMIS = RROOEMIS*(1.0_JPRB/(1.0_JPRB+RWRR)) + RCANEMIS*(RWRR/(1.0_JPRB+RWRR))

!==============================================================================
 
 
 
! CALCULATE ROUGHNESS LENGTH AND HEAT COEFFICIENT BASED ON HGT, HWR and WRR

 ! CANYON
 !
 ! The canyon roughness is dependant on building height, road width and canyon width
 ! The formulation here is mainly based on Harman 2004 and Porson 2010
 !


 sc_hwr  = 0.5_JPRB*(RHWR/2.0_JPRB*ATAN(1.0_JPRB))
 d_h     = max(1.0_JPRB-RWRR*(RCDA**(RWRR-1.0_JPRB)),0.0_JPRB)
 disp    = 1.4_JPRB
 width   = RHGT/RHWR
 z0x     = 0.1_JPRB * RBUIZ0M


 RCANZTM = (RCDG*(1.0_JPRB-d_h)*sc_hwr*RWRR/RVKSQ)**(-0.5_JPRB)
 RCANZTM = (1.0_JPRB-d_h)*EXP(-RCANZTM)
 RCANZTM = RCANZTM*RHGT
 RCANZTM = max(RCANZTM, RBUIZ0M)

! ----- Bulk Resistance Calculations -----
 pi      = 4.0_JPRB*ATAN(1.0_JPRB)
 zref    = 0.1_JPRB*RHGT
 l_rec   = 3.0_JPRB*RHGT
 
! Factors relating to recirculation region edge 
! from Harman 2004 (equation 13), Grimmond and Oke 1999a and MacDonald 1998
 lrw = 3.0_JPRB*(2.0_JPRB*rhwr/pi)
 ltw = 1.5_JPRB*(2.0_JPRB*rhwr/pi)
 hrh = (lrw-1.0_JPRB)/(lrw-ltw)
 lrw = MIN(lrw, 1.0_JPRB) ! Limits equal to that of Porson 2010.
 ltw = MIN(ltw, 1.0_JPRB)
 hrh = MAX(hrh, 0.0_JPRB)
 hrh = MIN(hrh, 1.0_JPRB)
 lttw = lrw+(1.0_JPRB-hrh)*SQRT((lrw-ltw)**2.0_JPRB+RHWR**2.0_JPRB) ! Friction velocity

! Wind components
!w_ori = SQRT( ruwind**2.0 + rvwind**2.0)*2.0/pi
!w_hyp = SQRT( ruwind**2.0 + rvwind**2.0)

 swi = MAX(RHGT-disp, RCANZTM)
 uct = MAX((2.0_JPRB/pi)*LOG(swi/RCANZTM)/LOG(RMODW/RCANZTM+1.0_JPRB),0.1_JPRB)
!This contains assumptions about normalised wind, 2.0/pi assumes random orientation

 
! Harman 2004 (Oke 1987) flow regime with length of sloing edge (l_slo) and downstream wall (l_dow)

! Three regimes 
! 1 - HWR LT 1/3 = Isolated
! 2 - HWR GT 1/3 and LT 2/3 = Interference
! 3 - HWR GT 2/3 = Skimming. 
! Min. calculations taken from Porson 2010.
 
 IF (RHWR <= 1.0_JPRB/3.0_JPRB) THEN 
! Isolated flow see Figure 1a Harman 2004
  l_slo = RHGT*SQRT(3.25_JPRB)
  l_dow = RHGT
 ELSE IF (1.0_JPRB/3.0_JPRB < RHWR .AND. RHWR <= 2.0_JPRB/3.0_JPRB) THEN 
! Wake interference see Figure 1b Harman 2004
  l_slo = SQRT(((width-l_rec/2.0_JPRB)**2.0_JPRB)*((2.0_JPRB*RHGT/l_rec)**2.0_JPRB+1.0_JPRB))
  l_dow = 2.0_JPRB*(width-l_rec/2.0_JPRB)*RHGT/l_rec
 ELSE 
! Skimming flow see Figure 1c Harman 2004
  l_slo = 0.0_JPRB
  l_dow = 0.0_JPRB
 END IF
  l_slo = 1.25_JPRB*l_slo

!!!! Flow regime = Isolated

 IF (RHWR <= 1.0_JPRB/3.0_JPRB) THEN
  IF (l_rec == width) THEN
   l_rec = l_rec - 1.0E-7_JPRB ! For computational reasons to prevent 0.0 difference
  END IF

 ! Recirculation region near wake of building
  fur = uct*RHGT*EXP(-RCDA*l_slo/RHGT)*(1.0_JPRB-EXP(-RCDA*l_rec/RHGT))
  fur = fur/(RCDA*l_rec)
   
  fuw = uct*EXP(-RCDA*(l_slo+l_rec)/RHGT)*(1.0_JPRB-EXP(-RCDA))
  fuw = fuw / RCDA
 
 ! Ventilated region downstream of recirculation region
  fdr = RHGT*EXP(-RCDA*l_slo/RHGT)*(1.0_JPRB-EXP(-RCDA*(width-l_rec)/RHGT))
  fdr = fdr/(RCDA*(width-l_rec)) ! width =/= l_rec

 ! Lower flow limit
  low = uct*LOG(zref/RBUIZ0M)/LOG(RHGT/RBUIZ0M) ! Harman eq 21
  fdr = MAX(fdr,low)

  fdw = uct*EXP(-RCDA*(l_slo+width-l_rec)/RHGT)*(1.0_JPRB-EXP(-RCDA))/RCDA
 ! Lower flow limit
  low = uct*((1.0_JPRB/(1.0_JPRB-(RBUIZ0M/RHGT)))-(1.0_JPRB/LOG(RHGT/RBUIZ0M))) ! Harman eq 22
  fdw = MAX(fdw, low)

!!!! Flow regime = Wake Interference 

 ELSE IF ( 1.0_JPRB/3.0_JPRB < RHWR .AND. RHWR <= 2.0_JPRB/3.0_JPRB) THEN
  IF (l_dow == RHGT) THEN
   l_dow = l_dow - 1.0E-7_JPRB ! For computational reasons to prevent 0.0 difference
  END IF

 ! Upstream Recirculation region near wake of building
  fur = uct*RHGT*EXP(-RCDA*(l_slo+RHGT-l_dow)/RHGT)*(1.0_JPRB-EXP(-RCDA*width/RHWR))
  fur = fur/(RCDA*width)

  fuw = uct * EXP(-RCDA*(l_slo+RHGT-l_dow+width)/RHGT)*(1.0_JPRB-EXP(-RCDA))
  fuw = fuw/RCDA
  
 ! Downstream Ventilated region of recirculation region
  fdw = uct*RHGT*EXP(-RCDA * l_slo/RHGT)*(1.0_JPRB-EXP(-RCDA*l_dow/RHGT))
  fdw = fdw/(RCDA*l_dow)

  low = MAX((RHGT - l_dow), RBUIZ0M)
  low = uct*(RHGT/l_dow - 1.0_JPRB / LOG(RHGT/RBUIZ0M)-(RHGT - l_dow)*LOG(low/RBUIZ0M)/(l_dow*LOG(RHGT/RBUIZ0M)))
  fdw = MAX(fdw, low)

 ! Downstream Recirculation region near wake of building impacts wall only
  fdr = RHGT*EXP(-RCDA*l_slo/RHGT)*(1.0_JPRB-EXP(-RCDA*(RHGT-l_dow)/RHGT))
  fdr = fdr*uct/(RCDA*(RHGT - l_dow))
  fdw = (l_dow*fdw+(RHGT-l_dow)*fdr)/RHGT ! Both recirc. and vent.
 

 !!!! Flow regime = Skimming
 
 ELSE IF (RHWR >= 2.0_JPRB/3.0_JPRB) THEN

  low = l_rec/(2.0_JPRB*width) ! Adjusts RCDA (flow resistance) - low used as spare

  fdw = uct*(1.0_JPRB-EXP(-RCDA))/RCDA
  fdr = fdw ! Not needed - provides Harman R8 value

  fur = uct*RHGT*EXP(-low*RCDA)*(1.0_JPRB-EXP(-low*RCDA*width/RHGT))
  fur = fur/(low*RCDA*width)

  fuw = uct*EXP(-low*RCDA*(RHGT+width)/RHGT)*(1.0_JPRB-EXP(-low*RCDA))
  fuw = fuw/(low*RCDA)

END IF

! Using material roughness to account for roughness along road and walls 
 can_road = uct*LOG(zref/RBUIZ0M)/LOG(RHGT/RBUIZ0M)
 can_wall = uct*(1.0_JPRB/(1.0_JPRB-(RBUIZ0M/RHGT))-1.0_JPRB/LOG(RHGT/RBUIZ0M))

! Combine roughness "tangential" and alongside canyon elements. 
fur = SQRT(fur**2.0_JPRB + can_road**2.0_JPRB)
fdr = SQRT(fdr**2.0_JPRB + can_road**2.0_JPRB)
fuw = SQRT(fuw**2.0_JPRB + can_wall**2.0_JPRB)
fdw = SQRT(fdw**2.0_JPRB + can_wall**2.0_JPRB)

! Harman resistances calculated from representative 
! wind speeds of the turbulent transport (f**) 
! See Harman 2004 figure 3 for detailed plan

 low = LOG(zref/RBUIZ0M) * LOG(zref/z0x)
 r_3 = low/(RVKSQ*fuw) ! Harman eq 15    
 r_4 = low/(RVKSQ*fur) ! Harman eq 16
!R5 represents the transfer across top of recirculation region
!Road used because it is more representative of canyon (Porson 2010)
 r_5 = ((1.0_JPRB-fur)*((LOG(RMODW/RCANZTM+1.0_JPRB))/SQRT(RVKSQ))**2.0_JPRB) / (1.0_JPRB+lttw-ltw) ! Adjusted Harman eq17
 r_6 = low/(RVKSQ*fdr) ! Harman eq 23
 r_7 = (1.0_JPRB-fdw)*((LOG(RMODW/RCANZTM+1.0_JPRB))/SQRT(RVKSQ))**2.0_JPRB ! Harman eq 24 adjusted to Porson 2010
 r_8 = low/(RVKSQ*fdw) ! Harman eq 25
 
 r_3 = r_3 / RHWR
 r_4 = r_4 / (MIN(l_rec, width) / width )
 r_5 = r_5 / (MIN(l_rec/2.0_JPRB, width) / width )

 low = MAX(1E-3_JPRB, (RHGT-l_dow)/width) ! Recirc. using Porson estimation
 swi = r_8/low ! reuse swi, spare variable name
 
 low = MAX(1E-3_JPRB, (width-MIN(l_rec, width))/width)
 r_6 = r_6/low
 
 low = MAX(1E-3_JPRB, (width-MIN(l_rec/2.0_JPRB, width))/width)
 r_7 = r_7/low
 
 low = MAX(1E-3_JPRB,l_dow/width) ! Venti.
 r_8 = r_8/low

 r_recirc = ((r_3*r_4*swi)/((r_3*r_4)+(r_3*swi)+(r_4*swi))) + r_5
 r_venti  = ((r_6*r_8)/(r_6+r_8)) + r_7
 r_bulk   = 1.0_JPRB/((1.0_JPRB/r_recirc) + (1.0_JPRB/r_venti)) ! Porson eq A11
 
 RCANRES = r_bulk
 RCANZTH = RHGT*(RMODZ+RCANZTM)/(RHGT*EXP( RVKSQ*r_bulk/LOG(RMODW/RCANZTM + 1 ))) ! Porson eq 18
! RCANZTH = MAX( RCANZTH, 1.0E-10_JPRB ) ! Porson limit
 IF ( RHWR < (1.0_JPRB/3.0_JPRB) ) THEN ! Porson limit
  RCANZTH = MIN(RCANZTH, 0.1_JPRB*RBUIZ0M)
 END IF

 RCANHC  = 1.0_JPRB/( (1.0_JPRB/RVKSQ)*LOG((RMODW+RCANZTM)/RCANZTM)*LOG((RMODW+RCANZTM)/RCANZTH)  ) ! Porson eq 19
 
 ! ROOF


 RROOZTM = RCANZTM

! Porson R_bulk,f (APPENDIX)

 low = MAX((1.1_JPRB*RHGT - disp), RROOZTM)
 swi = MAX(LOG(low/RROOZTM)/LOG(RMODW/RROOZTM+1.0_JPRB),1.0E-7_JPRB)

! r_int = r_recirc for computational efficiency
 r_recirc = (LOG(zref/RBUIZ0M+1.0_JPRB)*LOG((zref+RBUIZ0M)/(0.1_JPRB*RBUIZ0M)))/(RVKSQ*swi) ! Porson eq A8

! r_inert = r_venti for computational efficiency
 r_venti = (1.0_JPRB-swi)*LOG(RMODW/RROOZTM + 1.0_JPRB)/RVKSQ ! Porson eq A10

 r_bulk = r_recirc + r_venti ! Porson eq A12

 RROORES = r_bulk
 RROOZTH = RHGT*(RMODZ+RROOZTM)/( RHGT*EXP(RVKSQ*r_bulk/LOG(RMODW/RROOZTM+1.0_JPRB)) ) ! Porson eq 18

 !RROOZTH = MAX(RROOZTH, 1.0E-10_JPRB) ! Porson limit
 IF ( RHWR < (1.0_JPRB/3.0_JPRB) ) THEN   ! Porson limit
  RROOZTH = MIN(RROOZTH, 0.1_JPRB*RBUIZ0M)
 END IF

RROOHC  = 1.0_JPRB/( (1.0_JPRB/RVKSQ)*LOG((RMODW+RROOZTM)/RCANZTM)*LOG((RMODW+RROOZTM)/RCANZTH)  ) ! Porson eq 19

!URBAN PROPERTIES
RURBRES = (1.0_JPRB/(1.0_JPRB+RWRR))*RROORES + (RWRR/(1.0_JPRB+RWRR))*RCANRES
RURBZTM = (1.0_JPRB/(1.0_JPRB+RWRR))*RROOZTM + (RWRR/(1.0_JPRB+RWRR))*RCANZTM
!RURBZTM = 1.9_JPRB
RURBZTH = (1.0_JPRB/(1.0_JPRB+RWRR))*RROOZTH + (RWRR/(1.0_JPRB+RWRR))*RCANZTH
!RURBZTH = 1.0E-2_JPRB
RURBHC  = (1.0_JPRB/(1.0_JPRB+RWRR))*RROOHC  + (RWRR/(1.0_JPRB+RWRR))*RCANHC

RCANTC  = 2.0_JPRB*RHWR*RWALTC    + RROATC
!RCANTC  = (RHWR/(RHWR+1.0_JPRB))*RWALTC + (1.0_JPRB/(RHWR+1.0_JPRB))*RROATC

RCANVHC = 2.0_JPRB*RHWR*RWALVHC   + RROAVHC
!RCANVHC  = (RHWR/(RHWR+1.0_JPRB))*RWALVHC + (1.0_JPRB/(RHWR+1.0_JPRB))*RROAVHC

!RURBTC  = (1.0_JPRB-RWRR)*RROOTC  + RWRR*RCANTC
RURBTC  = (RWRR/(1.0_JPRB+RWRR))*RCANTC + (1.0_JPRB/(1.0_JPRB+RWRR))*RROOTC
RURBTC1 = RURBTC*0.00001_JPRB


!RURBVHC = (1-RWRR)*RROOVHC + RWRR*RCANVHC
RURBVHC  = (RWRR/(1.0_JPRB+RWRR))*RCANVHC + (1.0_JPRB/(1.0_JPRB+RWRR))*RROOVHC

!==============================================================================

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUSURB_MOD:SUSURB',1,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

END SUBROUTINE SUSURB
END MODULE SUSURB_MOD
