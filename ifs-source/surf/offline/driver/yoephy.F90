MODULE YOEPHY
! (C) Copyright 2005- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

USE PARKIND1  ,ONLY : JPRB, JPIM

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOEPHY* - SWITCHES RELATED TO DIABATIC PROCESSES
!     -----------------------------------------------------------------

!        * E.C.M.W.F. PHYSICS PACKAGE *

LOGICAL :: LEPHYS
LOGICAL :: LECOND
LOGICAL :: LECUMF
LOGICAL :: LEDCLD
LOGICAL :: LEEVAP
LOGICAL :: LEGWDG
LOGICAL :: LEOZOC
LOGICAL :: LEQNGT
LOGICAL :: LERADI
LOGICAL :: LERADS
LOGICAL :: LESHCV
LOGICAL :: LESICE
LOGICAL :: LESURF
LOGICAL :: LEVDIF
LOGICAL :: LAGPHY
LOGICAL :: LEPCLD
LOGICAL :: LEO3CH
LOGICAL :: LNEEONLINE
LOGICAL :: LBUD23
LOGICAL :: LEMETHOX
LOGICAL :: LERA40
LOGICAL :: LECURR
LOGICAL :: LVDFTRAC
LOGICAL :: LMFTRAC
LOGICAL :: LERAIN
LOGICAL :: LEVGEN
LOGICAL :: LESSRO
INTEGER(KIND=JPIM) :: NALBEDOSCHEME
INTEGER(KIND=JPIM) :: NEMISSSCHEME
LOGICAL :: LEMWAVE
LOGICAL :: LEOCWA
LOGICAL :: LEOCCO
LOGICAL :: LEOCSA
LOGICAL :: LEOCLA
REAL(KIND=JPRB) :: RTHRFRTI
INTEGER(KIND=JPIM) :: NPHYINT
INTEGER(KIND=JPIM) :: NPHPROMA

LOGICAL :: LEFLAKE  ! FLake model
LOGICAL :: LWCOU
LOGICAL :: LWCOU2W
LOGICAL :: LWCOUHMF
LOGICAL :: LEOCML   ! KPP
LOGICAL :: LELAIV   ! LAI Climatology
LOGICAL :: LESN09  ! snow 2009 
LOGICAL :: LECTESSEL ! CTESSEL
LOGICAL :: LEAGS    ! AGS for CO2&Evap
LOGICAL :: LEFARQUHAR ! Farquhar photosynthesis model
LOGICAL :: LEOPTSURF ! Read optimized parameters from namelist
LOGICAL :: LEC4MAP  ! MAP FOR C3/C4 PHOTOSYNTHESIS TYPE
LOGICAL :: LEAIRCO2COUP    ! Variable air CO2 in photosynthesis
REAL(KIND=JPRB) :: RLAIINT
LOGICAL :: LECLIM10D ! 10-day clim interpolation
LOGICAL :: LESNML  ! Multi-layer snow activated 
LOGICAL :: LEURBAN ! Urban tile active
LOGICAL :: LEINTWIND  ! Interpolate wind to match T/Q level 
LOGICAL :: LEWARMSTART ! Apply warm start to surf prognostics (only snow)
LOGICAL :: LECOLDSTART ! Apply cold start to surf prognostics (only snow)
INTEGER(KIND=JPIM) :: NSNMLWS ! Type of warm start to use (1,2,3)
REAL(KIND=JPRB) :: RALFMINPSN ! Albedo of permanent snow
LOGICAL :: LECMF1WAY  ! 1 way coupling with cama-flood
INTEGER(KIND=JPIM) :: LECMF2LAKEC ! 2 way coupling with cama-flood updated lake cover 
                                  ! 0 == OFF
                                  ! 1 == lake cover with flood fraction
                                  ! 2 == add flood fraction to lake cover 

!
!     REFERENCE.
!     ----------

!     J.-J. MORCRETTE       E.C.M.W.F.      91/07/14

!     MODIFICATIONS
!     -------------

!     P. Viterbo   ECMWF   03-12-2004  Include user-defined RTHRFRTI
!     G. Balsamo   ECMWF   08-01-2006  Include Van Genuchten Hydro LEVGEN
!     G. Balsamo   ECMWF   11-01-2006  Include sub-grid sf. runoff LESSRO
!     S. Boussetta/G.Balsamo May 2010  Include CTESSEL switch LECTESSEL 
!     G. Balsamo/S. Boussetta June 2011 Include LEAGS switch (modularity of CO2&Evap)
!     S.Boussetta  Nov 2013  Include 10-day clim interpolation switch LECLIM10D
!     R. Hogan     ECMWF   14-01-2019  Replace LE4ALB with NALBEDOSCHEME
!     A. Agusti-Panareda ECMWF 18-11-2020 Include LEAIRCO2COUP (use variable air CO2 in photosynthesis)
!     A. Agusti-Panareda ECMWF 02-06-2021 Include photosynthesis parameters that are optimized with observations
!     ------------------------------------------------------------------

!  NAME     TYPE     PURPOSE
!  ----  :  ----   : ---------------------------------------------------
! LEPHYS : LOGICAL : SWITCH THE FULL E.C.M.W.F. PHYSICS PACKAGE ON
! LAGPHY : LOGICAL : IF TRUE, PHYSICS PACKAGE CALLED IN LAGGED MODE
! LECOND : LOGICAL : TURN THE LARGE-SCALE CONDENSATION ON
! LECUMF : LOGICAL : TURN THE MASS-FLUX CUMULUS CONVECTION SCHEME ON
! LEDCLD : LOGICAL : TURN THE DIAGNOSTIC CLOUD SCHEME ON
! LEPCLD : LOGICAL : TURN THE PROGNOSTIC CLOUD SCHEME ON
! LEEVAP : LOGICAL : TURN THE EVAPORATION OF PRECIPITATION ON
! LEGWDG : LOGICAL : TURN THE GRAVITY WAVE DRAG ON
! LEOZOC : LOGICAL : TURN THE CLIMATOLOGICAL OZONE ON
! LEQNGT : LOGICAL : TURN THE NEGATIVE HUMIDITY FIXER ON
! LERADI : LOGICAL : TURN THE RADIATION SCHEME ON
! LERADS : LOGICAL : TURN THE INTERACTIVE SURFACE RADIATIVE PROPERTIESON
! LESHCV : LOGICAL : TURN THE SHALLOW CONV. IN THE MASS-FLUX SCHEME ON
! LESICE : LOGICAL : TURN THE INTERACTIVE SEA ICE PROCESSES ON
! LESURF : LOGICAL : TURN THE INTERACTIVE SURFACE PROCESSES ON
! LEVGEN : LOGICAL : TURN THE VAN GENUCHTEN HYDROLOGY ON
! LESSRO : LOGICAL : TURN THE SUB-GRID SURFACE RUNOFF ON
! LEVDIF : LOGICAL : TURN THE VERTICAL DIFFUSION ON
! LEO3CH : LOGICAL : TURN THE O3 CHEMISTRY ON (for EC prog. ozone)
! LNEEONLINE: LOGICAL: USE ON-LINE CTESSEL IF TRUE 
! LBUD23 : LOGICAL : SWITCH FOR 3 AND 2 DIMENSIONAL BUDGETS 
! LEMETHOX: LOGICAL: TURN THE METHANE OXIDATION ON
! LERA40 : LOGICAL : EXTRA PHYSICS DIAGNOSTICS FOR ERA40
! LECURR : LOGICAL : IF TRUE, OCEAN CURRENT BOUNDARY CONDITION IS USED
! LVDFTRAC: LOGICAL: TURN TRACER TRANSPORT BY VERTICAL DIFFUSION ON 
! LMFTRAC: LOGICAL : TURN TRACER TRANSPORT BY MASS FLUX CONVECTION ON
! LERAIN : LOGICAL : RAIN ASSIMILATION
! LEOCWA : LOGICAL : WARM OCEAN LAYER PARAMETRIZATION
! LEOCCO : LOGICAL : COOL OCEAN SKIN PARAMETRIZATION
! LEOCSA : LOGICAL : SALINTY EFFECT ON SATURATION AT OCEAN SURFACE
! LEOCLA : LOGICAL : LANGMUIR CIURCULATION EFFECT IN VOSKIN 
! RTHRFRTI : INTEGER : MINIMUM FRACTION FOR ALL SURFACE TILES
! NALBEDOSCHEME : INTEGER : Surface albedo, (0) ERBE,
!       (1) MODIS 4 component (UV-Vis+NIR)x(direct+diffuse), (2) MODIS 6 component
! NEMISSSCHEME : INTEGER : (0) 2-band emissivity (window, non-window), 6-band emissivity
! LEFLAKE: LOGICAL : IF TRUE USE FLAKE OVER LAKES 
! LEOCML : LOGICAL : IF TRUE USE OCEAN MIXED LAYER MODEL
! LEOCML : LOGICAL : IF TRUE USE LAI MONTHLY CLIMATOLOGY
! LESN09 : LOGICAL  : IF TRUE use snow 2009 
! LECTESSEL : LOGICAL  : IF TRUE USE CTESSEL surface scheme CO2 components
! LEAGS  : LOGICAL  : IF TRUE USE CTESSEL for CO2 and Evap components
! LEFARQUHAR  : LOGICAL  : IF TRUE USE FARQUHAR MODEL FOR PHOTOSYNTHESIS
! LEOPTSURF : LOGICAL : IF TRUE READ OPTIMIZED PARAMETERS IN FARQUHAR MODEL FROM NAMELIST, OTHERWISE USE DEFAULT VALUES
! LEC4MAP  : LOGICAL  : IF TRUE USE C3/C4 PHOTOSYNTHESIS MAP FROM CLIMATE FIELD
! LEAIRCO2COUP  : LOGICAL  : IF TRUE USE CTESSEL for CO2 and Evap components
! RLAIINT : REAL: Relaxation factor between interactive LAI and climatological LAI (1:fully interactive, 0:climatological LAI is used, )
! LECLIM10D: Logical: IF TRUE interpolate between 10-day climate values (for albedo and LAI) 

!     -----------------------------------------------------------------
END MODULE YOEPHY
