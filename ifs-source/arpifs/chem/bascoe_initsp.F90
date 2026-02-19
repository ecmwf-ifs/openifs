! (C) Copyright 2009- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE BASCOE_INITSP(PNTOTAL,PMEDIAN,PGSTDEV,PTAIR,PAIR,PNDSA)
!
!*****************************************************************************
!
!     This subroutine calculates the initial values of particle number
!     density of a stratospheric sulfuric acid particle ensemble with an
!     initial log-normal size distribution (parameters given below).
!     If the input parameter PNTOTAL holds a positive real value, this
!     value will be used as the total number density of sulfate aerosols;
!     if PNTOTAL is equal to zero on input, a typical observed vertical
!     profile will be used (cf. Hofmann et al. GRL 13,1252,1986).
!     The parameters PMEDIAN and PGSTDEV hold the values of lognormal
!     median radius and geometric standard deviation of the initial
!     sulfate particle size distribution.
!     Subroutine SETBIN must be called once prior to the use of this
!     subroutine.
!
!
!     Input/output variables:
!     REAL(KIND=8)    PNTOTAL,PTAIR,PMEDIAN,PGSTDEV,PAIR
!     It uses:        NBINS, PTSIZE(NBINS,8),PNDSA(NBINS)
!
!     Input: 
!         NBINS:                  Number of particle radii bins
!         PNTOTAL:  (par/cm**3)    Total number density of particles
!         PMEDIAN:  (m)            Lognormal median particle radius
!         PGSTDEV:                 Lognormal geometric standard deviation
!         PTAIR:    (K)            Ambient air temperature
!         PAIR:    (Pa)           Ambient air pressure
!         PTSIZE:  (m,m**2,m**3)  Particle radii, surface, volume
!
!     Output:
!         PNDSA:    (par/kg air)   Number density of sulfate particles
!
USE BASCOE_MODULE, ONLY : NBINS, PTSIZE, D1
USE PARKIND1  ,    ONLY : JPIM,   JPRB
USE YOMHOOK   ,    ONLY : LHOOK,  DR_HOOK, JPHOOK
USE YOMCST    ,    ONLY : R , RMD     

IMPLICIT NONE
REAL(KIND=JPRB), INTENT(IN)  :: PNTOTAL,PMEDIAN,PGSTDEV,PTAIR,PAIR
REAL(KIND=JPRB), INTENT(OUT) :: PNDSA(NBINS)

!     Physical constants:
!REAL, PARAMETER( &
!       Molar weight of dry air (kg/mole)
!     &          MAIR=28.9644E-3, &         ! VH - use RMD * 1e-3
!       Universal gas constant (J/(mole K))
!     &          RGAS=8.31441)              ! VH -use R
!
!     Auxiliary local variables:
REAL(KIND=JPRB)     :: ZP,ZNDTOT,ZD2,ZD3,ZD4
INTEGER(KIND=JPIM)  :: JL

REAL(KIND=JPHOOK)     :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('BASCOE_INITSP',0,ZHOOK_HANDLE )
!----------------------------------------------------------------------------
IF(PNTOTAL==0.0_JPRB) THEN
!         Total number density of sulfate particles a function of pressure:
  ZP=MAX(10.0E2_JPRB,MIN(PAIR,400.0E2_JPRB))
  IF(ZP <= 100.0E2_JPRB) THEN
    ZNDTOT=3.3448E-5_JPRB*ZP+6.6552E-1_JPRB
  ELSEIF(ZP <= 300.0E2_JPRB) THEN
    ZNDTOT=3.0103E-5_JPRB*ZP+6.9897E-1_JPRB
  ELSE
    ZNDTOT=1.3522E-4_JPRB*ZP-2.4545_JPRB
  ENDIF
  ZNDTOT=10.0_JPRB**ZNDTOT
  ZNDTOT=1.0E6_JPRB*ZNDTOT
ELSEIF(PNTOTAL > 0.0_JPRB) THEN
  ZNDTOT=PNTOTAL*1.0E6_JPRB
ELSE
  ZNDTOT=0.0_JPRB
ENDIF
!----------------------------------------------------------------------------
!*  Conversion to number of particles per kg air:

ZNDTOT=ZNDTOT * R * PTAIR/(RMD*1e-3_JPRB*PAIR)

!*  Calculate log-normal distribution:

!VH  - Explicit comutation below:
! CALL BASCOE_LGNDST(NBINS,ZNDTOT,PMEDIAN,PGSTDEV,PTSIZE,PNDSA)
!
!     This calculates a log-normal distribution with
!     total number density of particles ZNDTOT, median radius PMEDIAN,
!     and geometric standard deveation PGSTDEV. 
!     NBINS:		      Number of particle radii bins
!     NTOTAL:  (par/kg air)   Total number density of particles
!     MEDIAN:  (m)	      Median particle radius
!     GSTDEV:		      Geometric standard deveation
!     PTSIZE:  (m,m**2,m**3)  Particle radii, surface, volume
!
!     Output:
!     PNDSA:      (par/kg air)   Number density of particles

ZD3=LOG(PGSTDEV)
ZD2=ZNDTOT*D1/ZD3
ZD3=-0.5_JPRB/(ZD3*ZD3)
ZD4=LOG(PMEDIAN)
!
DO JL=1,NBINS
   PNDSA(JL)=ZD2*EXP(ZD3*(PTSIZE(JL,7)-ZD4)**2)
ENDDO


IF (LHOOK) CALL DR_HOOK('BASCOE_INITSP',1,ZHOOK_HANDLE )
END SUBROUTINE BASCOE_INITSP
