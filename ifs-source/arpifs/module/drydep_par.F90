! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction
! 
! (C) Copyright 1989- Meteo-France.
! 

MODULE DRYDEP_PAR

!!****  * MODULE DRYDEP_PAR
!!
!!     declare wesely parameterization constants
!!
!!    AUTHOR
!!    ------
!!    J. FLEMMING based on SUMO code J. Arteta
!!
!!
!!    MODIFICATIONS
!!    -------------
!!    Original   June 2014
!!
!!-------------------------------------------------------------------------------

USE  PARKIND1 ,ONLY : JPIM, JPRB

IMPLICIT NONE
SAVE

  INTEGER (KIND=JPIM), PARAMETER :: JPVEGNB_GC = 11_JPIM  ! number of GeosChem vegetation types
  INTEGER (KIND=JPIM), PARAMETER :: JPVEGNB_WE = 11_JPIM  ! number of Wesely vegetation types
  INTEGER (KIND=JPIM), PARAMETER :: JPSEANB_WE = 5_JPIM   ! number of Wesely season types

  REAL (KIND=JPRB), PARAMETER :: ZLUSW = 50.0_JPRB
  REAL (KIND=JPRB), PARAMETER :: ZGSS_W = 50.0_JPRB

  REAL (KIND=JPRB), PARAMETER :: ZEPS = 1.E-5_JPRB
  REAL (KIND=JPRB), PARAMETER :: ZEPS1 = 1.E-4_JPRB
  REAL (KIND=JPRB), PARAMETER :: ZEPS2 = 10.0_JPRB

! DRY DEPOSITION TYPES:
!      WESELY ORIGINAL         
!  1)   URBAN                  
!  2)   AGRICULTURAL           
!  3)   RANGE LAND             
!  4)   DECIDUOUS FOREST       
!  5)   CONIFEROUS FOREST      
!  6)   MIXED FOREST & WETLAND 
!  7)   WATER                  
!  8)   BARREN/DESERT          
!  9)   NON-FORESTED WETLAND   
! 10)   MIXED AGRI & RANGE     
! 11)   ROCKY OPEN & SHRUBS    

REAL (KIND=JPRB), PARAMETER, DIMENSION(JPVEGNB_WE, JPSEANB_WE) :: ZRLU =  & ! Wesely rlu
&  RESHAPE (SOURCE = (/9999, 2000, 2000, 2000, 2000, 2000, 9999, 9999, 2500, 2000, 4000,   &
&                      9999, 9000, 9000, 9000, 4000, 8000, 9999, 9999, 9000, 9000, 9000,   &
&                      9999, 9999, 9000, 9000, 4000, 8000, 9999, 9999, 9000, 9000, 9000,   &
&                      9999, 9999, 9999, 9999, 6000, 9000, 9999, 9999, 9000, 9000, 9000,   &
&                      9999, 4000, 4000, 4000, 2000, 3000, 9999, 9999, 4000, 4000, 8000/), &
&           SHAPE = (/JPVEGNB_WE, JPSEANB_WE/))

REAL (KIND=JPRB), PARAMETER, DIMENSION(JPVEGNB_WE, JPSEANB_WE) :: ZRAC = & ! Wesely rac
&  RESHAPE (SOURCE = (/100, 200, 100, 2000, 2000, 2000, 0, 0, 300, 150, 200,   &
&                      100, 150, 100, 1500, 2000, 1700, 0, 0, 200, 120, 140,   &
&                      100,  10, 100, 1000, 2000, 1500, 0, 0, 100,  50, 120,   &
&                      100,  10,  10, 1000, 2000, 1500, 0, 0,  50,  10,  50,   &
&                      100,  50,  80, 1200, 2000, 1500, 0, 0, 200,  60, 120/), &
&           SHAPE = (/JPVEGNB_WE, JPSEANB_WE/))

REAL (KIND=JPRB), PARAMETER, DIMENSION(JPVEGNB_WE, JPSEANB_WE) :: ZRGSS = & ! Wesely rgsS
&  RESHAPE (SOURCE = (/400, 150, 350, 500, 500, 100, 0, 1000,   0, 220, 400,   &
&                      400, 200, 350, 500, 500, 100, 0, 1000,   0, 300, 400,   &
&                      400, 150, 350, 500, 500, 200, 0, 1000,   0, 200, 400,   &
&                      100, 100, 100, 100, 100, 100, 0, 1000, 100, 100,  50,   &
&                      500, 150, 350, 500, 500, 200, 0, 1000,   0, 250, 400/), &
&           SHAPE = (/JPVEGNB_WE, JPSEANB_WE/))

REAL (KIND=JPRB), PARAMETER, DIMENSION(JPVEGNB_WE, JPSEANB_WE) :: ZRGSO = & ! Wesely rgso
&  RESHAPE (SOURCE = (/300,  150,  200,  200,  200,  300, 2000, 400, 1000,  180,  200,   &
&                      300,  150,  200,  200,  200,  300, 2000, 400,  800,  180,  200,   &
&                      300,  150,  200,  200,  200,  300, 2000, 400, 1000,  180,  200,   &
&                      600, 3500, 3500, 3500, 3500, 3500, 2000, 400, 3500, 3500, 3500,   &
&                      300,  150,  200,  200,  200,  300, 2000, 400, 1000,  180,  200/), &
&           SHAPE = (/JPVEGNB_WE, JPSEANB_WE/))

REAL (KIND=JPRB), PARAMETER, DIMENSION(JPVEGNB_WE, JPSEANB_WE) :: ZRCLS = & ! Wesely rclS
&  RESHAPE (SOURCE = (/9999, 2000, 2000, 2000, 2000, 2000, 9999, 9999, 2500, 2000, 4000,   &
&                      9999, 9000, 9000, 9000, 2000, 4000, 9999, 9999, 9000, 9000, 9000,   &
&                      9999, 9999, 9000, 9000, 3000, 6000, 9999, 9999, 9000, 9000, 9000,   &
&                      9999, 9999, 9999, 9000,  200,  400, 9999, 9999, 9000, 9999, 9000,   &
&                      9999, 4000, 4000, 4000, 2000, 3000, 9999, 9999, 4000, 4000, 8000/), &
&           SHAPE = (/JPVEGNB_WE, JPSEANB_WE/))

REAL (KIND=JPRB), PARAMETER, DIMENSION(JPVEGNB_WE, JPSEANB_WE) :: ZRCLO = & ! Wesely rclO
&  RESHAPE (SOURCE = (/9999, 1000, 1000, 1000, 1000, 1000, 9999, 9999, 1000, 1000, 1000,   &
&                      9999,  400,  400,  400, 1000,  600, 9999, 9999,  400,  400,  400,   &
&                      9999, 1000,  400,  400, 1000,  600, 9999, 9999,  800,  600,  600,   &
&                      9999, 1000, 1000,  400, 1500,  600, 9999, 9999,  800, 1000,  800,   &
&                      9999, 1000,  500,  500, 1500,  700, 9999, 9999,  600,  800,  800/), &
&           SHAPE = (/JPVEGNB_WE, JPSEANB_WE/))

INTEGER (KIND=JPIM), PARAMETER :: JPVEGNB_AR=4   ! Number of vegetation types:
                                      ! 1: sea; 2: ice; 3: low vegetation;
                                      ! 4: forests for the vegetation


!REAL(KIND=JPRB) , PARAMETER :: GF1=0.55
REAL(KIND=JPRB) , PARAMETER :: RSMAX=5000.
REAL(KIND=JPRB) , PARAMETER :: GCONV=1.E3
                   ! FACTEUR DE CHANGEMENT D'UNITE (MASSE VOL. DE L'EAU LIQUIDE)
REAL(KIND=JPRB) , PARAMETER :: EWFC=0.35
REAL(KIND=JPRB) , PARAMETER :: GWWILT=37.1342E-3
REAL(KIND=JPRB) , PARAMETER :: EWWILT=0.5
REAL(KIND=JPRB) , PARAMETER :: GWLMX=0.2_JPRB

REAL(KIND=JPRB) , PARAMETER :: GWFC=89.0467E-3

REAL(KIND=JPRB) , PARAMETER :: EDB=5.       !determines PCD in near the neutral regime
REAL(KIND=JPRB) , PARAMETER :: EDC=5.       !determines PCD in an unstable asymptotic regime
REAL(KIND=JPRB) , PARAMETER :: EDD=5.       !determines PCD in a stable asymptotic regime
REAL(KIND=JPRB) , PARAMETER :: Z2B=2.*EDB
REAL(KIND=JPRB) , PARAMETER :: Z3BC=3.*EDB*EDC
REAL(KIND=JPRB) , PARAMETER :: Z3B=3._JPRB*EDB


!REAL (KIND=JPRB), DIMENSION(0:3,4) :: ZGZ ! ISBA - roughness length
!REAL (KIND=JPRB) :: GWLEX


REAL (KIND=JPRB) , PARAMETER  :: PPDIMOH2O=0.15/0.595  ! molecular diffusivity of H2O (cm2/s)


REAL(KIND=JPRB) , PARAMETER               :: PPRAND = 0.72 ! Van Karman constant
!
REAL(KIND=JPRB) , PARAMETER               :: PPKARMAN = 0.4 ! Van Karman constant, in RKAP, set in ifs/phys_ec/suvdf.F90

! Symbolic constants extracted from suphy1.F
REAL(KIND=JPRB) , PARAMETER :: WCRIN=10.  !used in the computation of proportion
                              ! of grid cell covered by snow
REAL(KIND=JPRB) , PARAMETER               :: PPRMAX = 1E+5 ! resistance maximum value
                              ! used for surface, aerodynamic and laminar
                              ! resistances
REAL(KIND=JPRB), PARAMETER               :: PPINVI = 0.15
                              ! kinematic viscosity of air (cm2/s)

!REAL(KIND=JPRB),  PARAMETER ::      RTT=273.16
!REAL(KIND=JPRB),  PARAMETER ::      RLVTT=2.5008E+6
!REAL(KIND=JPRB),  PARAMETER ::      RLSTT=2.8345E+6

!  SATURATED VAPOUR PRESSURE

!REAL(KIND=JPRB), PARAMETER ::  RESTT=611.14
!REAL(KIND=JPRB), PARAMETER ::  RGAMW=(RCW-RCPV)/RV
!REAL(KIND=JPRB), PARAMETER ::  RBETW=RLVTT/RV+RGAMW*RTT
!REAL(KIND=JPRB),           ::      RALPW
!REAL(KIND=JPRB), PARAMETER ::  RGAMS=(RCS-RCPV)/RV
!REAL(KIND=JPRB), PARAMETER ::  RBETS=RLSTT/RV+RGAMS*RTT
!REAL(KIND=JPRB),           ::      RALPS
!REAL(KIND=JPRB), PARAMETER ::  RGAMD=RGAMS-RGAMW
!REAL(KIND=JPRB), PARAMETER ::  RBETD=RBETS-RBETW
!REAL(KIND=JPRB),           ::  RALPD


! Species - specific paramaters
! set in chem_ini

TYPE TDRYDEP
REAL(KIND=JPRB),   ALLOCATABLE :: RDIMO(:)   ! Molecular Diffusivity

REAL(KIND=JPRB) ,  ALLOCATABLE :: RCHEN(:)   ! efficicent Henry coefficient  1

REAL(KIND=JPRB) ,  ALLOCATABLE :: RCHENXP(:) ! efficicent Henry coefficient 2

REAL(KIND=JPRB) ,  ALLOCATABLE :: RCF0(:)    ! Reactivity
END TYPE TDRYDEP

END MODULE DRYDEP_PAR

