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

MODULE YOE_CUCONVCA

! Modified : 
!    R. El Khatib 09-Mar-2012 Allocate RCUCONVCA/RNLCONVCA for safe bound checkings later
!    L. Bengtsson 05-Aug-2014 Correction if using with GOL
!    L. Gerard    31-Mar-2016 Add RCADELX
USE PARKIND1 , ONLY : JPIM, JPRB
USE RANDOM_NUMBERS_MIX, ONLY: RANDOMNUMBERSTREAM

IMPLICIT NONE

SAVE

!     ------------------------------------------------------

!*    control parameters for cellular automaton for convection

!      yd_random_stream_CA: random number stream for CA update
!      NIJH         : Multiplicative for number of points in high resolution grid
!      NBmax        : maximum possible number of neighbours for every gridcell
!      RCADELX      : Length of individual CA cell (meters) to compute NIJH
!      RCUCONVCA    : Array for interaction with the physics
!      RNLCONVCA    : Array for interaction with the physics
!      RCAPECONVCA  : Array to store "good" CAPE between its calculation [MOD(KSTEP*TSPHY,3600._JPRB)]
!      NCELLCU      : Lat-lon grid for cellular automaton
!      NFERTCU      : Lat-lon grid for cellular automaton (fertile cell indicator)
!      NBLS         : neighbours for large scale grid (= model red gaussian grid)
!      NBSS         : neighbours for small scale grid (= NIJH*NIJH cells in every LS cell)
!      NRNBLS       : number of neighbours for large scale cells
!      NRNBSS       : number of neighbours for small scale cells
!      NBNS         : neighbour north/south or on same lat (1/-1/0)
!      NBEW         : neighbour east/west or on same long (1/-1/0)
!      RWASALIVE    : Array indicating grid-points where CA was alive before
!      RWGHTCU      : Weighted and smoothed CA pattern
!      LCUCONV_CA   : TRUE IF CELLULAR AUTOMATON USED IN CONVECTION
!      LCA_ADVECT   : switch for "kind of" semi-Lagrangian advection
!      LCA_GLOBAL   : switch for global CA instead ov "convective" CA
!      LCA_SMOOTH   : swith for smoothing CA on large scale grid
!      LCA_RANTROP  : switch for random coupling of CA to deep convection (concerns initialization)
!      LCA_TEST     : switch for initialize CA at single point
!      LCA_ADVTEST  : switch for advection test (CA only evolved first 5 steps)
!      NTESTPROC    : set on which processor the single point should lie
!      NTESTGP      : set which gridpoint the single point is on
!      NSPINUP      : set number of spin-up cycles for global pattern
!      CA_FORC      : switch to choose convective forcing
!      CA_WIND      : switch to choose CA-wind (real/idealized)
!      CA_PROB      : switch to choose probabilities
!      NLIVES       : switch to choose number of lives (scaled by CAPE in callpar)
!      NFERTYRS     : switch to choose max number of steps a cell can be fertile
!      NFRCASEED    : Frequency of seeding CA where physics diagnoses deep convection
!      RCA_SEEDPROB : Probability of random seeding for global CA

!      RPROB_SURVIVE: probabilities for cell survival
!      RPROB_BIRTH: probabilities for cell birth
!      RPROB_FUN_FERT_S: function for specifying probabilities for survival as function of upwind(uw)+downwind(dw) neighbours
!      RPROB_FUN_UW_DW_S: function for specifying probabilities for survival as function of upwind(uw)-downwind(dw) neighbours
!      RPROB_FUN_FERT_B: function for specifying probabilities for birth as function of upwind(uw)+downwind(dw) neighbours
!      RPROB_FUN_UW_DW_B: function for specifying probabilities for birth as function of upwind(uw)-downwind(dw) neighbours

!      SLLONG(YRSL%NASLB1)       : Array containing the longitude for each gridpoint in a semi/lagrangian buffer array
!      SLLAT(YRSL%NASLB1)        : Array containing the latitude for each gridpoint in a semi/lagrangian buffer array
!      SLDLONG(YRSL%NASLB1)      : Array containing the width of each gridbox in longitudinal direction 
!                             in a semi/lagrangian buffer array
!      SLDLAT(YRSL%NASLB1)       : Array containing the width of each gridbox in latitudinal direction 
!                             in a semi/lagrangian buffer array
!      SLDDLAT(YRSL%NASLB1)      : Array containing the difference of the gridbox's northern boundary latitude 
!                             and the latitude of the gp (necessary because latitudes are not uniformly distributed)

INTEGER(KIND=JPIM)            :: NIJH=4
INTEGER(KIND=JPIM), PARAMETER :: NBMAX=8

TYPE :: TECUCONVCA
TYPE (RANDOMNUMBERSTREAM)     :: YD_RANDOM_STREAM_CA
INTEGER(KIND=JPIM)            :: NLIVES
INTEGER(KIND=JPIM)            :: NFERTYRS
INTEGER(KIND=JPIM)            :: NSPINUP

INTEGER(KIND=JPIM),ALLOCATABLE:: NCELLCU(:,:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NFERTCU(:,:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NBLS(:,:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NBSS(:,:,:,:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NRNBLS(:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NRNBSS(:,:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NBEW(:,:,:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NBNS(:,:,:)
REAL(KIND=JPRB)   ,ALLOCATABLE:: RCUCONVCA(:)
REAL(KIND=JPRB)   ,ALLOCATABLE:: RNLCONVCA(:)
REAL(KIND=JPRB)   ,ALLOCATABLE:: RCAPECONVCA(:)
REAL(KIND=JPRB)   ,ALLOCATABLE:: RWASALIVE(:)
REAL(KIND=JPRB)   ,ALLOCATABLE:: RWGHTCU(:)
REAL(KIND=JPRB)   ,ALLOCATABLE:: RLATLONNBLS(:,:,:)
REAL(KIND=JPRB)   ,ALLOCATABLE:: RLATLONNBSS(:,:,:,:)

REAL(KIND=JPRB)               :: RPROB_SURVIVE(0:8,0:8)
REAL(KIND=JPRB)               :: RPROB_BIRTH(0:8,0:8)
REAL(KIND=JPRB)               :: RPROB_FUN_FERT_S(0:16)
REAL(KIND=JPRB)               :: RPROB_FUN_UW_DW_S(0:16)
REAL(KIND=JPRB)               :: RPROB_FUN_FERT_B(0:16)
REAL(KIND=JPRB)               :: RPROB_FUN_UW_DW_B(0:16)

REAL(KIND=JPRB)   ,ALLOCATABLE:: RLONDEP(:,:,:)
REAL(KIND=JPRB)   ,ALLOCATABLE:: RLATDEP(:,:,:)

REAL(KIND=JPRB)               :: RCA_SEEDPROB
REAL(KIND=JPRB)               :: RCADELX

LOGICAL                       :: LCUCONV_CA
LOGICAL                       :: LCA_ADVECT
LOGICAL                       :: LCA_GLOBAL
LOGICAL                       :: LCA_SMOOTH
LOGICAL                       :: LCA_TEST,LCA_ADVTEST
LOGICAL                       :: LCA_RANTROP
LOGICAL                       :: LCA_NBDEBUG
LOGICAL                       :: LCA_EXTRACT
CHARACTER(LEN=10)             :: CA_FORC
CHARACTER(LEN=10)             :: CA_WIND
CHARACTER(LEN=10)             :: CA_PROB
INTEGER(KIND=JPIM)            :: NFRCASEED
INTEGER(KIND=JPIM)            :: NTESTPROC, NTESTGP
INTEGER(KIND=JPIM)            :: NDXUNREAL, NDYUNREAL

REAL(KIND=JPRB)   ,ALLOCATABLE:: SLLONG(:), SLLAT(:)
REAL(KIND=JPRB)   ,ALLOCATABLE:: SLDLONG(:), SLDLAT(:), SLDDLAT(:)
!----------------------------------------------------------------------------
CONTAINS
  PROCEDURE, PASS :: PRINT => PRINT_CONFIGURATION 
END TYPE TECUCONVCA


!     ------------------------------------------------------

CONTAINS

SUBROUTINE PRINT_CONFIGURATION(SELF, KDEPTH, KOUTNO)
  IMPLICIT NONE
  CLASS(TECUCONVCA), INTENT(IN) :: SELF
  INTEGER          , INTENT(IN) :: KDEPTH
  INTEGER          , INTENT(IN) :: KOUTNO

  INTEGER :: IDEPTHLOC

  IDEPTHLOC = KDEPTH+2
  
  WRITE(KOUTNO,*) REPEAT(' ',KDEPTH   ) // 'model%yrml_phy_ec%yrecuconvca : '
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NLIVES = ', SELF%NLIVES
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NFERTYRS = ', SELF%NFERTYRS
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NSPINUP = ', SELF%NSPINUP
  IF (ALLOCATED(SELF%NCELLCU)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NCELLCU ALLOCATED OF SHAPE ',SHAPE(SELF%NCELLCU), &
    &        ' SUM ',SUM(SELF%NCELLCU)
  IF (ALLOCATED(SELF%NFERTCU)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NFERTCU ALLOCATED OF SHAPE ',SHAPE(SELF%NFERTCU), &
    &        ' SUM ',SUM(SELF%NFERTCU)
  IF (ALLOCATED(SELF%NBLS)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NBLS ALLOCATED OF SHAPE ',SHAPE(SELF%NBLS), &
    &        ' SUM ',SUM(SELF%NBLS)
  IF (ALLOCATED(SELF%NBSS)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NBSS ALLOCATED OF SHAPE ',SHAPE(SELF%NBSS), &
    &        ' SUM ',SUM(SELF%NBSS)
  IF (ALLOCATED(SELF%NRNBLS)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NRNBLS ALLOCATED OF SHAPE ',SHAPE(SELF%NRNBLS), &
    &        ' SUM ',SUM(SELF%NRNBLS)
  IF (ALLOCATED(SELF%NRNBSS)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NRNBSS ALLOCATED OF SHAPE ',SHAPE(SELF%NRNBSS), &
    &        ' SUM ',SUM(SELF%NRNBSS)
  IF (ALLOCATED(SELF%NBEW)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NBEW ALLOCATED OF SHAPE ',SHAPE(SELF%NBEW), &
    &        ' SUM ',SUM(SELF%NBEW)
  IF (ALLOCATED(SELF%NBNS)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NBNS ALLOCATED OF SHAPE ',SHAPE(SELF%NBNS), &
    &        ' SUM ',SUM(SELF%NBNS)
  IF (ALLOCATED(SELF%RCUCONVCA)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RCUCONVCA ALLOCATED OF SHAPE ',SHAPE(SELF%RCUCONVCA), &
    &        ' SUM ',SUM(SELF%RCUCONVCA)
  IF (ALLOCATED(SELF%RNLCONVCA)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RNLCONVCA ALLOCATED OF SHAPE ',SHAPE(SELF%RNLCONVCA), &
    &        ' SUM ',SUM(SELF%RNLCONVCA)
  IF (ALLOCATED(SELF%RCAPECONVCA)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RCAPECONVCA ALLOCATED OF SHAPE ',SHAPE(SELF%RCAPECONVCA), &
    &        ' SUM ',SUM(SELF%RCAPECONVCA)
  IF (ALLOCATED(SELF%RWASALIVE)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RWASALIVE ALLOCATED OF SHAPE ',SHAPE(SELF%RWASALIVE), &
    &        ' SUM ',SUM(SELF%RWASALIVE)
  IF (ALLOCATED(SELF%RWGHTCU)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RWGHTCU ALLOCATED OF SHAPE ',SHAPE(SELF%RWGHTCU), &
    &        ' SUM ',SUM(SELF%RWGHTCU)
  IF (ALLOCATED(SELF%RLATLONNBLS)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RLATLONNBLS ALLOCATED OF SHAPE ',&
    &        SHAPE(SELF%RLATLONNBLS), ' SUM ',SUM(SELF%RLATLONNBLS)
  IF (ALLOCATED(SELF%RLATLONNBSS)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RLATLONNBSS ALLOCATED OF SHAPE ',&
    &        SHAPE(SELF%RLATLONNBSS), ' SUM ',SUM(SELF%RLATLONNBSS)
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RPROB_SURVIVE sum ', SUM(SELF%RPROB_SURVIVE)
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RPROB_BIRTH sum ', SUM(SELF%RPROB_BIRTH)
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RPROB_FUN_FERT_S sum ', SUM(SELF%RPROB_FUN_FERT_S)
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RPROB_FUN_UW_DW_S sum ', SUM(SELF%RPROB_FUN_UW_DW_S)
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RPROB_FUN_FERT_B sum ', SUM(SELF%RPROB_FUN_FERT_B)
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RPROB_FUN_UW_DW_B sum ', SUM(SELF%RPROB_FUN_UW_DW_B)
  IF (ALLOCATED(SELF%RLONDEP)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RLONDEP ALLOCATED OF SHAPE ',SHAPE(SELF%RLONDEP), &
    &        ' SUM ',SUM(SELF%RLONDEP)
  IF (ALLOCATED(SELF%RLATDEP)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RLATDEP ALLOCATED OF SHAPE ',SHAPE(SELF%RLATDEP), &
    &        ' SUM ',SUM(SELF%RLATDEP)
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RCA_SEEDPROB = ', SELF%RCA_SEEDPROB
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RCADELX = ', SELF%RCADELX
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LCUCONV_CA = ', SELF%LCUCONV_CA
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LCA_ADVECT = ', SELF%LCA_ADVECT
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LCA_GLOBAL = ', SELF%LCA_GLOBAL
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LCA_SMOOTH = ', SELF%LCA_SMOOTH
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LCA_TEST = ', SELF%LCA_TEST
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LCA_ADVTEST = ', SELF%LCA_ADVTEST
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LCA_RANTROP = ', SELF%LCA_RANTROP
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LCA_NBDEBUG = ', SELF%LCA_NBDEBUG
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LCA_EXTRACT = ', SELF%LCA_EXTRACT
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'CA_FORC = ', SELF%CA_FORC
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'CA_WIND = ', SELF%CA_WIND
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'CA_PROB = ', SELF%CA_PROB
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NFRCASEED = ', SELF%NFRCASEED
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NTESTPROC = ', SELF%NTESTPROC
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NTESTGP = ', SELF%NTESTGP
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NDXUNREAL = ', SELF%NDXUNREAL
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NDYUNREAL = ', SELF%NDYUNREAL
  IF (ALLOCATED(SELF%SLLONG)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SLLONG ALLOCATED OF SHAPE ',SHAPE(SELF%SLLONG), &
    &        ' SUM ',SUM(SELF%SLLONG)
  IF (ALLOCATED(SELF%SLLAT)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SLLAT ALLOCATED OF SHAPE ',SHAPE(SELF%SLLAT), &
    &        ' SUM ',SUM(SELF%SLLAT)
  IF (ALLOCATED(SELF%SLDLONG)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SLDLONG ALLOCATED OF SHAPE ',SHAPE(SELF%SLDLONG), &
    &        ' SUM ',SUM(SELF%SLDLONG)
  IF (ALLOCATED(SELF%SLDLAT)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SLDLAT ALLOCATED OF SHAPE ',SHAPE(SELF%SLDLAT), &
    &        ' SUM ',SUM(SELF%SLDLAT)
  IF (ALLOCATED(SELF%SLDDLAT)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SLDDLAT ALLOCATED OF SHAPE ',SHAPE(SELF%SLDDLAT), &
    &        ' SUM ',SUM(SELF%SLDDLAT)

END SUBROUTINE PRINT_CONFIGURATION

!     ------------------------------------------------------
!     CA-initialization
!       - read NAMELIST
!       - if CA active initialize constants and allocate arrays
!     ------------------------------------------------------
SUBROUTINE INI_CUCONVCA(YDGEOMETRY,YDDYNA,YDECUCONVCA,YDSL)

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMDYNA      , ONLY : TDYNA
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCT0       , ONLY : LECMWF
USE YOMCST       , ONLY : RPI
USE YOMLUN       , ONLY : NULOUT, NULNAM
USE YOMGRIB      , ONLY : NENSFNB
USE EINT_MOD     , ONLY : SL_STRUCT

IMPLICIT NONE
TYPE(GEOMETRY)  , INTENT(IN)   :: YDGEOMETRY
TYPE(TDYNA)     , INTENT(IN)   :: YDDYNA
TYPE(TECUCONVCA), INTENT(INOUT),TARGET :: YDECUCONVCA
TYPE(SL_STRUCT) , INTENT(INOUT):: YDSL
INTEGER(KIND=JPIM) :: JI,JJ, JIJ
REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

CHARACTER(LEN=10) :: CA_FORC, CA_WIND, CA_PROB
LOGICAL, POINTER ::  LCUCONV_CA, LCA_ADVECT, LCA_GLOBAL, LCA_SMOOTH, LCA_TEST,&
 & LCA_RANTROP, LCA_ADVTEST, LCA_EXTRACT, LCA_NBDEBUG
INTEGER(KIND=JPIM), POINTER ::  NFRCASEED, NLIVES, NFERTYRS, NTESTPROC,&
 & NTESTGP, NDXUNREAL, NDYUNREAL, NSPINUP
REAL(KIND=JPRB), POINTER ::  RCA_SEEDPROB, RCADELX

#include "namca.nam.h"
#include "posnam.intfb.h"
#include "setran.intfb.h"

IF (LHOOK) CALL DR_HOOK('YOE_CUCONVCA:INI_CUCONVCA',0,ZHOOK_HANDLE)

! Associate pointers for variables in namelist
LCUCONV_CA   => YDECUCONVCA%LCUCONV_CA
LCA_ADVECT   => YDECUCONVCA%LCA_ADVECT
LCA_GLOBAL   => YDECUCONVCA%LCA_GLOBAL
LCA_SMOOTH   => YDECUCONVCA%LCA_SMOOTH
LCA_TEST     => YDECUCONVCA%LCA_TEST
NFRCASEED    => YDECUCONVCA%NFRCASEED
NLIVES       => YDECUCONVCA%NLIVES
NFERTYRS     => YDECUCONVCA%NFERTYRS
RCA_SEEDPROB => YDECUCONVCA%RCA_SEEDPROB
RCADELX      => YDECUCONVCA%RCADELX
NTESTPROC    => YDECUCONVCA%NTESTPROC
NTESTGP      => YDECUCONVCA%NTESTGP
NDXUNREAL    => YDECUCONVCA%NDXUNREAL
NDYUNREAL    => YDECUCONVCA%NDYUNREAL
NSPINUP      => YDECUCONVCA%NSPINUP
LCA_RANTROP  => YDECUCONVCA%LCA_RANTROP
LCA_ADVTEST  => YDECUCONVCA%LCA_ADVTEST
LCA_EXTRACT  => YDECUCONVCA%LCA_EXTRACT
LCA_NBDEBUG  => YDECUCONVCA%LCA_NBDEBUG

!     ----------------
!     Setup namelist variables
!     ----------------

LCUCONV_CA=.FALSE.! use Cellular Automaton 
IF(LECMWF)THEN
  LCA_ADVECT=.TRUE.
ELSE
  LCA_ADVECT=.FALSE. ! use pseudo semi-lag advection of CA
ENDIF
IF(LECMWF)THEN
  NIJH=4
ELSE
  NIJH=2
ENDIF
RCADELX=0._JPRB
LCA_GLOBAL=.FALSE.
LCA_SMOOTH=.TRUE.
LCA_RANTROP=.FALSE.
LCA_EXTRACT =.FALSE.
CA_FORC="CONV_TQ"
CA_WIND="REAL"
CA_PROB="WIND"
IF(LECMWF)THEN
  RCA_SEEDPROB=0.995
ELSE
  RCA_SEEDPROB=0.75
ENDIF
NFRCASEED=1
NLIVES=10
NFERTYRS=3
LCA_TEST=.FALSE.
LCA_ADVTEST=.FALSE.
NTESTPROC=8
NTESTGP=1
NDXUNREAL=-1
NDYUNREAL=0
NSPINUP=1000

!     ----------------
!     Read namelist
!     ----------------
CALL POSNAM(NULNAM,'NAMCA')
READ(NULNAM,NAMCA)

! These are not a pointers into YDECUCONVCA to workaround a PGI compiler bug,
! so we have to update YDECUCONVCA components after reading the namelist
YDECUCONVCA%CA_FORC = CA_FORC 
YDECUCONVCA%CA_WIND = CA_WIND
YDECUCONVCA%CA_PROB = CA_PROB
IF (RCADELX>1._JPRB) THEN
   NIJH=MAX(1,NINT(YDGEOMETRY%YRGEM%RDELXN/RCADELX))
ENDIF

IF (.NOT. YDDYNA%LSLAG) LCUCONV_CA=.FALSE.

IF (.NOT. LCUCONV_CA) THEN
  WRITE(NULOUT,*) 'convective CA not active!'
ELSE
  WRITE(NULOUT,*) 'Initializing convective CA:'
  WRITE(NULOUT,*) 'CA settings are:'
  WRITE(NULOUT,*) 'Number of CA cells in grid cell: ',NIJH
  WRITE(NULOUT,*) 'Number of lives: ',NLIVES
  WRITE(NULOUT,*) 'Number of fertile Years: ',NFERTYRS
  WRITE(NULOUT,*) 'Frequency of seeding: ',NFRCASEED
  WRITE(NULOUT,*) '"semi-lagrangian" advection of CA: ',LCA_ADVECT
  WRITE(NULOUT,*) 'Wind used: ', CA_WIND
  WRITE(NULOUT,*) 'Forcing applied to convection: ', CA_FORC
  WRITE(NULOUT,*) 'Probabilities used: ', CA_PROB
  WRITE(NULOUT,*) 'Ensemble member number: ', NENSFNB

  !   ----------------
  !   Setup CA grid and allocate arrays
  !   ----------------

 
  ALLOCATE(YDECUCONVCA%RCUCONVCA(YDGEOMETRY%YRGEM%NGPTOT))
  ALLOCATE(YDECUCONVCA%RNLCONVCA(YDGEOMETRY%YRGEM%NGPTOT))
  ALLOCATE(YDECUCONVCA%RCAPECONVCA(YDGEOMETRY%YRGEM%NGPTOT))
  ALLOCATE(YDECUCONVCA%RWASALIVE(YDGEOMETRY%YRGEM%NGPTOT))
  ALLOCATE(YDECUCONVCA%NCELLCU(YDGEOMETRY%YRGEM%NGPTOT,NIJH*NIJH))
  ALLOCATE(YDECUCONVCA%NFERTCU(YDGEOMETRY%YRGEM%NGPTOT,NIJH*NIJH))
  ALLOCATE(YDECUCONVCA%NBLS(YDSL%NASLB1,NBMAX))
  ALLOCATE(YDECUCONVCA%NBSS(YDSL%NASLB1,2,NBMAX,NIJH*NIJH))
  ALLOCATE(YDECUCONVCA%NBNS(YDSL%NASLB1,NBMAX,NIJH*NIJH))
  ALLOCATE(YDECUCONVCA%NBEW(YDSL%NASLB1,NBMAX,NIJH*NIJH))
  ALLOCATE(YDECUCONVCA%NRNBLS(YDSL%NASLB1))
  ALLOCATE(YDECUCONVCA%NRNBSS(YDSL%NASLB1,NIJH*NIJH))
  ALLOCATE(YDECUCONVCA%RWGHTCU(YDGEOMETRY%YRGEM%NGPTOT))

  ASSOCIATE(NGPBLKS=>YDGEOMETRY%YRDIM%NGPBLKS, NPROMA=>YDGEOMETRY%YRDIM%NPROMA, &
   & NFLEVG=>YDGEOMETRY%YRDIMV%NFLEVG, &
   & RCUCONVCA=>YDECUCONVCA%RCUCONVCA,  RNLCONVCA=>YDECUCONVCA%RNLCONVCA, &
   & RCAPECONVCA=>YDECUCONVCA%RCAPECONVCA,  RWASALIVE=>YDECUCONVCA%RWASALIVE, &
   & NCELLCU=>YDECUCONVCA%NCELLCU, NFERTCU=>YDECUCONVCA%NFERTCU,  &
   & NBLS=>YDECUCONVCA%NBLS, NBSS=>YDECUCONVCA%NBSS, NBNS=>YDECUCONVCA%NBNS, &
   & NBEW=>YDECUCONVCA%NBEW, NRNBLS=>YDECUCONVCA%NRNBLS, &
   & NRNBSS=>YDECUCONVCA%NRNBSS, RWGHTCU=>YDECUCONVCA%RWGHTCU, &
   & RPROB_FUN_FERT_S=>YDECUCONVCA%RPROB_FUN_FERT_S, &
   & RPROB_FUN_UW_DW_S=>YDECUCONVCA%RPROB_FUN_UW_DW_S, &
   & RPROB_FUN_FERT_B=>YDECUCONVCA%RPROB_FUN_FERT_B, &
   & RPROB_FUN_UW_DW_B=>YDECUCONVCA%RPROB_FUN_UW_DW_B, &
   & RPROB_SURVIVE=>YDECUCONVCA%RPROB_SURVIVE, &
   & RPROB_BIRTH=>YDECUCONVCA%RPROB_BIRTH, &
   & YD_RANDOM_STREAM_CA=>YDECUCONVCA%YD_RANDOM_STREAM_CA)
  RCUCONVCA=0.0
  RNLCONVCA=0.0 !NLIVES

  RCAPECONVCA=0.0
  RWASALIVE=0.0

  NCELLCU=0
  NFERTCU=0

  NBLS=0
  NBSS=0
  NRNBLS=0
  NRNBSS=0
  NBNS=0
  NBEW=0

  !initialize random number stream (as function of ensemble member number and the date)
  CALL SETRAN (NENSFNB,YD_RANDOM_STREAM_CA)

  !initialize probabilities for CA
  IF ( CA_PROB == 'WIND' ) THEN
    ! wind dominated probabilities
    !                        0      1      2      3      4      5      6      7      8      9
    RPROB_FUN_FERT_S  = (/ 0.000, 0.500, 0.950, 0.950, 0.950, 0.700, 0.550, 0.000, 0.000, 0.000,&
    !                       10     11     12     13     14     15     16
                         & 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /)
    !                        0      1      2      3      4      5      6      7      8      9
    RPROB_FUN_UW_DW_S = (/ 0.000, 0.000, 0.000, 0.000, 0.000, 0.005, 0.100, 0.250, 0.700, 0.950,&
    !                       10     11     12     13     14     15     16
                         & 0.950, 0.900, 0.850, 0.750, 0.700, 0.650, 0.650 /)
    !                        0      1      2      3      4      5      6      7      8      9
    RPROB_FUN_FERT_B  = (/ 0.000, 0.950, 0.950, 0.700, 0.600, 0.300, 0.150 , 0.000, 0.000, 0.000,&
    !                       10     11     12     13     14     15     16
                         & 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /)
    !                        0      1      2      3      4      5      6      7      8      9
    RPROB_FUN_UW_DW_B = (/ 0.000, 0.000, 0.000, 0.000, 0.000, 0.005, 0.100, 0.250, 0.700, 0.950,&
    !                       10     11     12     13     14     15     16
                         & 0.950, 0.900, 0.850, 0.750, 0.700, 0.650, 0.650 /)
  ELSEIF ( CA_PROB == 'NOWIND' ) THEN
    ! equal probabilities for UW/DW
    !                        0      1      2      3      4      5      6      7      8      9
    RPROB_FUN_FERT_S  = (/ 0.000, 0.250, 0.600, 0.950, 0.950, 0.950, 0.800, 0.650, 0.500, 0.000,&
    !                       10     11     12     13     14     15     16
                         & 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /)
    !                        0      1      2      3      4      5      6      7      8      9
    RPROB_FUN_UW_DW_S = (/ 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,&
    !                       10     11     12     13     14     15     16
                         & 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000 /)
    !                        0      1      2      3      4      5      6      7      8      9
    RPROB_FUN_FERT_B  = (/ 0.000, 0.250, 0.950, 0.950, 0.700, 0.600, 0.300, 0.150 , 0.000, 0.000,&
    !                       10     11     12     13     14     15     16
                         & 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /)
    !                        0      1      2      3      4      5      6      7      8      9
    RPROB_FUN_UW_DW_B = (/ 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,&
    !                       10     11     12     13     14     15     16
                         & 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000 /)
  ELSEIF ( CA_PROB == 'GOL' ) THEN
    ! Deterministic rules, following 'the game of life'
    !                        0      1      2      3      4      5      6      7      8      9
    RPROB_FUN_FERT_S  = (/ 0.000, 0.000, 1.000, 1.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,&
    !                       10     11     12     13     14     15     16
                         & 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /)
    !                        0      1      2      3      4      5      6      7      8      9
    RPROB_FUN_UW_DW_S = (/ 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,&
    !                       10     11     12     13     14     15     16
                         & 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000 /)
    !                        0      1      2      3      4      5      6      7      8      9
    RPROB_FUN_FERT_B  = (/ 0.000, 0.000, 0.000, 1.000, 0.000, 0.000, 0.000, 0.000 , 0.000, 0.000,&
    !                       10     11     12     13     14     15     16
                         & 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /)
    !                        0      1      2      3      4      5      6      7      8      9
    RPROB_FUN_UW_DW_B = (/ 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,&
    !                       10     11     12     13     14     15     16
                         & 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000 /)
  ENDIF

  ! This loop has been manually collapsed to gain a bit longer vector
  DO JIJ=0,80
    JI=INT(JIJ/9,JPIM)  !DO JI=0,8
    JJ=MOD(JIJ,9)       !DO JJ=0,8
    RPROB_SURVIVE(JI,JJ)=RPROB_FUN_FERT_S(JI+JJ)*RPROB_FUN_UW_DW_S(JJ-JI+8)
    RPROB_BIRTH(JI,JJ)=RPROB_FUN_FERT_B(JI+JJ)*RPROB_FUN_UW_DW_B(JJ-JI+8)
  ENDDO

  !setup LAT/LONG, DLAT, DLONG arrays
  CALL SETUP_LATLONGHELP(YDGEOMETRY,YDECUCONVCA,YDSL)
  !caluclate neighbours for CA cells on reduced gaussian grid
  CALL CALCULATE_NEIGHBOURS(YDECUCONVCA,YDGEOMETRY,YDSL)

  IF (LCA_ADVECT) THEN
    ALLOCATE(YDECUCONVCA%RLONDEP(NPROMA,NFLEVG,NGPBLKS))
    ALLOCATE(YDECUCONVCA%RLATDEP(NPROMA,NFLEVG,NGPBLKS))
  ENDIF

  END ASSOCIATE
ENDIF

IF (LHOOK) CALL DR_HOOK('YOE_CUCONVCA:INI_CUCONVCA',1,ZHOOK_HANDLE)
END SUBROUTINE INI_CUCONVCA


!     ------------------------------------------------------
!     prepare arrays containing latitude, longitude and 
!     gridbox size for every gridbox in SL-buffer
!     ------------------------------------------------------
SUBROUTINE SETUP_LATLONGHELP(YDGEOMETRY,YDECUCONVCA,YDSL)

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK 
USE YOMCST       , ONLY : RPI

!stuff for SL-halo
USE YOMMP0   , ONLY : MY_REGION_NS, MY_REGION_EW
USE EINT_MOD , ONLY : SL_STRUCT

IMPLICIT NONE
TYPE(GEOMETRY)  , INTENT(IN)    :: YDGEOMETRY
TYPE(TECUCONVCA), INTENT(INOUT) :: YDECUCONVCA
TYPE(SL_STRUCT) , INTENT(INOUT) :: YDSL
REAL(KIND=JPHOOK)              :: ZHOOK_HANDLE
INTEGER(KIND=JPIM)           :: J,JLS, JLAT, JK
INTEGER(KIND=JPIM)           :: IGP, ILG
INTEGER(KIND=JPIM)           :: INSTANLON, INRLON
REAL(KIND=JPRB), ALLOCATABLE :: ZDLONG(:),ZDLAT(:),ZDDLAT(:)
REAL(KIND=JPRB)              :: ZDLO, ZDLA, ZDDLA

!stuff for SL-halo
LOGICAL                      :: LLINC
INTEGER(KIND=JPIM)           :: IDUMARR(2)
INTEGER(KIND=JPIM)           :: IFIXSFLD(2)
INTEGER(KIND=JPIM)           :: IFLDSLB1
REAL(KIND=JPRB), ALLOCATABLE :: ZPB1(:,:)

#include "slcomm.intfb.h"

IF (LHOOK) CALL DR_HOOK('YOE_CUCONVCA:SETUP_LATLONGHELP',0,ZHOOK_HANDLE)

ALLOCATE(YDECUCONVCA%SLLONG(YDSL%NASLB1))
ALLOCATE(YDECUCONVCA%SLLAT(YDSL%NASLB1))
ALLOCATE(YDECUCONVCA%SLDLONG(YDSL%NASLB1))
ALLOCATE(YDECUCONVCA%SLDLAT(YDSL%NASLB1))
ALLOCATE(YDECUCONVCA%SLDDLAT(YDSL%NASLB1))

ASSOCIATE(YDVAB=>YDGEOMETRY%YRVAB,YDVETA=>YDGEOMETRY%YRVETA,YDVFE=>YDGEOMETRY%YRVFE,YDSTA=>YDGEOMETRY%YRSTA, &
 & YDLAP=>YDGEOMETRY%YRLAP,YDCSGLEG=>YDGEOMETRY%YRCSGLEG, &
 & YDGSGEOM_NB=>YDGEOMETRY%YRGSGEOM_NB, YDSPGEOM=>YDGEOMETRY%YSPGEOM)
ASSOCIATE(NLOENG=>YDGEOMETRY%YRGEM%NLOENG, NGPTOT=>YDGEOMETRY%YRGEM%NGPTOT, &
 & NDGENL=>YDGEOMETRY%YRDIM%NDGENL, NDGLG=>YDGEOMETRY%YRDIM%NDGLG, NPROMA=>YDGEOMETRY%YRDIM%NPROMA, &
 & MYLATS=>YDGEOMETRY%YRMP%MYLATS, NONL=>YDGEOMETRY%YRMP%NONL, NPTRFRSTLAT=>YDGEOMETRY%YRMP%NPTRFRSTLAT, &
 & SLLONG=>YDECUCONVCA%SLLONG, SLLAT=>YDECUCONVCA%SLLAT, &
 & SLDLONG=>YDECUCONVCA%SLDLONG, SLDLAT=>YDECUCONVCA%SLDLAT, &
 & SLDDLAT=>YDECUCONVCA%SLDDLAT)
ALLOCATE(ZDLONG(NGPTOT),ZDLAT(NGPTOT),ZDDLAT(NGPTOT))

!calculate dlong, dlat
IGP=1
DO JLAT=1,NDGENL
  ILG=MYLATS(JLAT)
  ZDLO=2*RPI/REAL(YDGEOMETRY%YRGEM%NLOENG(ILG),JPRB)         !calculate delta longitude between gp on lat
  IF (ILG == 1) THEN
    ZDLA=0.5_JPRB*RPI-0.5_JPRB*(YDCSGLEG%RLATIG(ILG)+YDCSGLEG%RLATIG(ILG+1))
    ZDDLA=0.5_JPRB*RPI-YDCSGLEG%RLATIG(ILG)
  ELSEIF (ILG == NDGLG) THEN
    ZDLA= 0.5_JPRB*(YDCSGLEG%RLATIG(ILG)+YDCSGLEG%RLATIG(ILG-1))+0.5_JPRB*RPI
    ZDDLA=0.5_JPRB*(YDCSGLEG%RLATIG(ILG-1)-YDCSGLEG%RLATIG(ILG))
  ELSE
    ZDLA=0.5_JPRB*(YDCSGLEG%RLATIG(ILG-1)-YDCSGLEG%RLATIG(ILG+1))
    ZDDLA=0.5_JPRB*(YDCSGLEG%RLATIG(ILG-1)-YDCSGLEG%RLATIG(ILG))
  ENDIF
  INSTANLON=NPTRFRSTLAT(MY_REGION_NS)+JLAT-1      !index for NSTA/NONL arrays
  INRLON=NONL(INSTANLON,MY_REGION_EW)             !number of gridpoints on this processor for current latitude
  IF (INRLON < 1 .OR. INRLON > NLOENG(ILG)) THEN
    CYCLE                    !go to next lat, if there are no points on this lat
  ENDIF
  DO J=1,INRLON                             !loop over gp on this latitude
    ZDLONG(IGP)=ZDLO
    ZDLAT(IGP)=ZDLA
    ZDDLAT(IGP)=ZDDLA
    IGP=IGP+1
  ENDDO
ENDDO

!SL-stuff
LLINC=.FALSE.
IFIXSFLD(:)=0
IFLDSLB1=5
ALLOCATE(ZPB1(YDSL%NASLB1,IFLDSLB1))
ZPB1=0

!fill buffer
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JK,JLS)
DO JK=1,NGPTOT,NPROMA
  DO JLS=JK,JK+MIN(NPROMA,NGPTOT-JK+1)-1
    ZPB1(YDSL%NSLCORE(JLS),1)=YDGSGEOM_NB%GELAM(JLS)
    ZPB1(YDSL%NSLCORE(JLS),2)=YDGSGEOM_NB%GELAT(JLS)
    ZPB1(YDSL%NSLCORE(JLS),3)=ZDLONG(JLS)
    ZPB1(YDSL%NSLCORE(JLS),4)=ZDLAT(JLS)
    ZPB1(YDSL%NSLCORE(JLS),5)=ZDDLAT(JLS)
  ENDDO
ENDDO
!$OMP END PARALLEL DO

!get halo
CALL SLCOMM(YDSL,IDUMARR,IFLDSLB1,LLINC,0,ZPB1)

!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JLS,JK)
DO JK=1,YDSL%NASLB1,NPROMA
  DO JLS=JK,JK+MIN(NPROMA,YDSL%NASLB1-JK+1)-1
    SLLONG(JLS)=ZPB1(JLS,1)
    SLLAT(JLS)=ZPB1(JLS,2)
    SLDLONG(JLS)=ZPB1(JLS,3)
    SLDLAT(JLS)=ZPB1(JLS,4)
    SLDDLAT(JLS)=ZPB1(JLS,5)
  ENDDO
ENDDO
!$OMP END PARALLEL DO

DEALLOCATE(ZPB1,ZDLONG,ZDLAT,ZDDLAT)

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('YOE_CUCONVCA:SETUP_LATLONGHELP',1,ZHOOK_HANDLE)

END SUBROUTINE SETUP_LATLONGHELP


!     ------------------------------------------------------
!     calculate neighbours on reduced gaussian grid
!     for each gridcell on this processor
!       - neighbour given in terms of index in the SL-buffer
!     ------------------------------------------------------
SUBROUTINE CALCULATE_NEIGHBOURS(YDECUCONVCA,YDGEOMETRY,YDSL)

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK 
USE YOMMP0  , ONLY : MY_REGION_NS, MY_REGION_EW
USE YOMCST  , ONLY : RPI
USE EINT_MOD, ONLY : SL_STRUCT
USE YOMCT0  , ONLY : LRPLANE

IMPLICIT NONE

TYPE(TECUCONVCA),INTENT(INOUT) :: YDECUCONVCA
TYPE(GEOMETRY)  , INTENT(IN)    :: YDGEOMETRY
TYPE(SL_STRUCT) , INTENT(INOUT) :: YDSL
INTEGER(KIND=JPIM)         :: JLAT, J, JLL,JI,JXSS,JYSS, JNB
INTEGER(KIND=JPIM)         :: ILG, ILATNB, ILGNB
INTEGER(KIND=JPIM)         :: INSTANLON
INTEGER(KIND=JPIM)         :: INRLON
INTEGER(KIND=JPIM)         :: INRNB
INTEGER(KIND=JPIM)         :: IGP, IGPNB, IGPSS
INTEGER(KIND=JPIM)         :: IINCR
INTEGER(KIND=JPIM)         :: ILON,ILONNB
INTEGER(KIND=JPIM)         :: IYNB, IXNB
INTEGER(KIND=JPIM)         :: IEW

REAL(KIND=JPRB)            :: ZDLONG, ZDLONGNB, ZDLONGMAX,ZDELTA
REAL(KIND=JPRB)            :: ZLONGSSGP, ZLONGSSGPNB
INTEGER(KIND=JPIM)         :: ILATGP(YDSL%NASLB1),ILONGP(YDSL%NASLB1)
REAL(KIND=JPRB)            :: ZTWOPI,ZEPSILON

LOGICAL                    :: LLCOMPLI

REAL(KIND=JPHOOK)            :: ZHOOK_HANDLE


IF (LHOOK) CALL DR_HOOK('YOE_CUCONVCA:CALCULATE_NEIGHBOURS',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP)
ASSOCIATE(NLOENG=>YDGEM%NLOENG, &
 & NDGENL=>YDDIM%NDGENL, NDGLG=>YDDIM%NDGLG, &
 & NSTA=>YDMP%NSTA, MYLATS=>YDMP%MYLATS, NONL=>YDMP%NONL, &
 & NPTRFRSTLAT=>YDMP%NPTRFRSTLAT, &
 & SLLONG=>YDECUCONVCA%SLLONG, NBLS=>YDECUCONVCA%NBLS, NBSS=>YDECUCONVCA%NBSS, &
 & NBNS=>YDECUCONVCA%NBNS, NBEW=>YDECUCONVCA%NBEW, NRNBLS=>YDECUCONVCA%NRNBLS, &
 & NRNBSS=>YDECUCONVCA%NRNBSS)

ZTWOPI=2*RPI
!work out "accuracy" of (x-ZTWOPI)
ZEPSILON=NEAREST(ZTWOPI,-1._JPRB)
ZEPSILON=ZTWOPI-ZEPSILON
ILATGP=-99
ILONGP=-99

!calculate longitude for every point in SL-array
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JLAT,ILG,INRLON,ZDLONG,J,IGP,ILON)
DO JLAT=0,NDGENL+1
  ILG=MYLATS(1)+JLAT-1           !global latitude number
  IF (ILG < 1) THEN
    ILG=1
  ELSEIF (ILG > NDGLG) THEN
    ILG=NDGLG
  ENDIF

  INRLON=YDSL%NSLONL(JLAT)             !number of gridpoints on this processor for current latitude
  DO J=1,INRLON
    IGP=YDSL%NSLOFF(JLAT)+J
    ILON=YDSL%NSLSTA(JLAT)+J-1-1       !second -1 because zero meridian = index one
    ILONGP(IGP)=ILON
    ILATGP(IGP)=JLAT
  ENDDO
ENDDO
!$OMP END PARALLEL DO

! Following loop remains to be vectorized. Ideally the code is rewritten
!  in a way to move the 1:INRLON loop to become the innermost one.

!$OMP PARALLEL DO SCHEDULE(STATIC) &
!$OMP& PRIVATE(JLAT,ILG,INSTANLON,INRLON,J,INRNB,IGP,JLL,ILATNB,IINCR,ILGNB) &
!$OMP& PRIVATE(LLCOMPLI,JI,IGPNB,ZDLONG,ZDLONGNB,ZDLONGMAX,ILON,ILONNB,ZDELTA) &
!$OMP& PRIVATE(JYSS,JXSS,IGPSS,IYNB,IXNB,JNB,ZLONGSSGP,ZLONGSSGPNB,IEW)
DO JLAT=1,NDGENL                            !loop over latitudes on this processor
  ILG=MYLATS(JLAT)                          !global latitude number
  INSTANLON=NPTRFRSTLAT(MY_REGION_NS)+JLAT-1      !index for NSTA/NONL arrays
  INRLON=NONL(INSTANLON,MY_REGION_EW)             !number of gridpoints on this processor for current latitude
  IF (INRLON < 1 .OR. INRLON > NLOENG(ILG)) THEN
    CYCLE                    !go to next lat, if there are no points on this lat (why <= 1 and not <1 ????)
  ENDIF
  DO J=1,INRLON                             !loop over gp on this latitude
    INRNB=0                                !reset number of neighbours found
    IGP=YDSL%NSLOFF(JLAT)+NSTA(INSTANLON,MY_REGION_EW)-YDSL%NSLSTA(JLAT)+J     !index of current gp in SL-buffer !
                                                                 !(!!!! potentially false-> inconsistent with IFS DOCU)
    DO JLL=-1,1                            !loop over lat above, own, below
      ILATNB=JLAT+JLL                     !local lat index of neighbour
      IINCR=1                             !east-west gridpoint increment
      IF (JLL == 0) IINCR=2               !on own latitude JI=0 is the point itself

      ILGNB=ILG+JLL                        !global lat index of neighbour

      ! LAM security
      IF (LRPLANE .AND. (ILGNB < 1)) ILGNB=1 
      IF (LRPLANE .AND. (ILGNB > NDGLG)) ILGNB=NDGLG

      !because there is no SL buffer in EW direction on lats containing all gp's there has to be done
      !this complicated logical
      LLCOMPLI=.FALSE.
      IF (NLOENG(ILGNB) == YDSL%NSLONL(ILATNB) .AND. (J == 1 .OR. J == INRLON)) THEN
        LLCOMPLI=.TRUE.
      ENDIF

      IF (NLOENG(ILG) == NLOENG(ILGNB)) THEN ! easy case: same number of longitudes on both lats 
                                             !(problematic at poles becuase NLOENG outside of array range??)
        DO JI=-1,1,IINCR
          IGPNB=YDSL%NSLOFF(ILATNB)+NSTA(INSTANLON,MY_REGION_EW)+J-YDSL%NSLSTA(ILATNB)+JI    !index of neighbour in SL-buffer
          IF (LLCOMPLI) THEN
            IF (IGPNB == YDSL%NSLOFF(ILATNB)) IGPNB=YDSL%NSLOFF(ILATNB)+YDSL%NSLONL(ILATNB)
            IF (IGPNB == YDSL%NSLOFF(ILATNB)+YDSL%NSLONL(ILATNB)+1) IGPNB=YDSL%NSLOFF(ILATNB)+1
          ENDIF
          INRNB=INRNB+1
          NBLS(IGP,INRNB)=IGPNB
        ENDDO
      ELSE                  !complicated case, different number of longs on lats
        ZDLONG=ZTWOPI/REAL(NLOENG(ILG),JPRB)       !calculate delta longitude between gp on lat
        ZDLONGNB=ZTWOPI/REAL(NLOENG(ILGNB),JPRB)
        ZDLONGMAX=0.5_JPRB*ZDLONG+0.5_JPRB*ZDLONGNB+ZEPSILON   !max delta in long for a gp to be a neighbour

        DO JI=1,YDSL%NSLONL(ILATNB)
          IGPNB=YDSL%NSLOFF(ILATNB)+JI                     !index of neighbour in SL-buffer

          ZDELTA=ABS(SLLONG(IGP)-SLLONG(IGPNB))
          DO WHILE (ZDELTA >= ZTWOPI)
            ZDELTA=ZDELTA-ZTWOPI
          ENDDO
          IF (ZDELTA > RPI) ZDELTA=ZTWOPI-ZDELTA
          IF (ZDELTA <= ZDLONGMAX .AND. IGPNB /= IGP) THEN
            INRNB=INRNB+1
            NBLS(IGP,INRNB)=IGPNB
          ENDIF
        ENDDO
      ENDIF
    ENDDO ! JLL
    NRNBLS(IGP)=INRNB

    !find neighbours for small scale grid
    DO JYSS=1,NIJH
      DO JXSS=1,NIJH
        IGPSS=(JYSS-1)*NIJH+JXSS             !index of small scale cell in large scale cell
        INRNB=0                              !reset number of neighbours found
        DO JLL=-1,1                          !loop over lat above, own, below

          IYNB=JYSS+JLL                      !Y ob neibouring row
          ILATNB=JLAT                          !large scale latitude of neighbour
          ILGNB=ILG
          IF (IYNB < 1) THEN                 !if SS nb-row is outside of LS cell (to the north)
            ILATNB=JLAT-1
            ILGNB=ILG-1
            IYNB=NIJH
          ENDIF
          IF (IYNB > NIJH) THEN                 !if SS nb-row is outside of LS cell (to the south)
            ILATNB=JLAT+1
            ILGNB=ILG+1
            IYNB=1
          ENDIF

          ! LAM security
          IF (LRPLANE .AND. (ILGNB < 1)) ILGNB=1
          IF (LRPLANE .AND. (ILGNB > NDGLG)) ILGNB=NDGLG

          !because there is no SL buffer in EW direction on lats containing all gp's there has to be done
          !this complicated logical
          LLCOMPLI=.FALSE.
          IF (NLOENG(ILGNB) == YDSL%NSLONL(ILATNB) .AND. (J == 1 .OR. J == INRLON)) THEN
            LLCOMPLI=.TRUE.
          ENDIF

          IINCR=1                             !east-west gridpoint increment
          IF (JLL == 0) IINCR=2               !on own latitude JI=0 is the point itself

          IF (ILATNB == JLAT) THEN            !easy case: NB and GP on same LS Lat
            DO JI=-1,1,IINCR
              IXNB=JXSS+JI
              IGPNB=IGP
              IF (IXNB < 1) THEN                 !if SS nb-column is outside of LS cell (to the west)
                IGPNB=IGP-1
                IXNB=NIJH
                IF (LLCOMPLI .AND. IGPNB == YDSL%NSLOFF(ILATNB)) IGPNB=YDSL%NSLOFF(ILATNB)+YDSL%NSLONL(ILATNB)
              ENDIF
              IF (IXNB > NIJH) THEN              !if SS nb-column is outside of LS cell (to the east)
                IGPNB=IGP+1
                IXNB=1
                IF (LLCOMPLI .AND. IGPNB == YDSL%NSLOFF(ILATNB)+YDSL%NSLONL(ILATNB)+1) IGPNB=YDSL%NSLOFF(ILATNB)+1
              ENDIF
              INRNB=INRNB+1
              NBSS(IGP,1,INRNB,IGPSS)=IGPNB
              NBSS(IGP,2,INRNB,IGPSS)=(IYNB-1)*NIJH+IXNB
              NBNS(IGP,INRNB,IGPSS)=-JLL
              NBEW(IGP,INRNB,IGPSS)=JI
            ENDDO
          ELSEIF (NLOENG(ILG) == NLOENG(ILGNB)) THEN      !NB and GP on different LS Lat, but same number of LS gp on lat
                                                          !necessarry because below's complicated case doesn't reliable detect
                                                          !neighbours for NLOENG(ILG) == NLOENG(ILGNB) (roundoff errors)
            ILON=NSTA(INSTANLON,MY_REGION_EW)+J-1-1
            DO JI=-1,1,IINCR
              IXNB=JXSS+JI
              ILONNB=ILON
              IF (IXNB < 1) THEN                 !if SS nb-column is outside of LS cell (to the west)
                ILONNB=ILON-1
                IXNB=NIJH
                IF (LLCOMPLI .AND. ILONNB == YDSL%NSLSTA(ILATNB)-1-1) ILONNB=ILONNB+YDSL%NSLONL(ILATNB)
              ENDIF
              IF (IXNB > NIJH) THEN              !if SS nb-column is outside of LS cell (to the east)
                ILONNB=ILON+1
                IXNB=1
                IF (LLCOMPLI .AND. ILONNB == YDSL%NSLSTA(ILATNB)-1-1+YDSL%NSLONL(ILATNB)+1) ILONNB=YDSL%NSLSTA(ILATNB)-1-1+1
              ENDIF
              DO JNB=1,NRNBLS(IGP)                 !loop over LS neighbours
                IGPNB=NBLS(IGP,JNB)
                IF (ILATNB == ILATGP(IGPNB) .AND. ILONNB == ILONGP(IGPNB)) THEN
                  INRNB=INRNB+1
                  NBSS(IGP,1,INRNB,IGPSS)=IGPNB
                  NBSS(IGP,2,INRNB,IGPSS)=(IYNB-1)*NIJH+IXNB
                  NBNS(IGP,INRNB,IGPSS)=-JLL
                  NBEW(IGP,INRNB,IGPSS)=JI
                ENDIF
              ENDDO
            ENDDO
          ELSE                                !complicated case: NB and GP on different LS Lat
            ZDLONG=ZTWOPI/REAL(NLOENG(ILG),JPRB)            !calculate delta longitude between LS gp on lat
            ZDLONGNB=ZTWOPI/REAL(NLOENG(ILGNB),JPRB)
            ZDLONGMAX=0.5_JPRB*ZDLONG/REAL(NIJH,JPRB)+0.5_JPRB*ZDLONGNB/REAL(NIJH,JPRB)+ZEPSILON !max delta in long 
                                                                                                 !for a gp to be a nbr
            DO JNB=1,NRNBLS(IGP)                 !loop over LS neighbours
              IGPNB=NBLS(IGP,JNB)
              IF (ILATNB == ILATGP(IGPNB)) THEN
                DO JI=1,NIJH                       !E-W loop over SS cells in LS cell
                  ZLONGSSGP=SLLONG(IGP)-(0.5_JPRB-0.5_JPRB/REAL(NIJH,JPRB)-REAL(JXSS-1,JPRB)/REAL(NIJH,JPRB))*ZDLONG
                  ZLONGSSGPNB=SLLONG(IGPNB)-(0.5_JPRB-0.5_JPRB/REAL(NIJH,JPRB)-REAL(JI-1,JPRB)/REAL(NIJH,JPRB))*ZDLONGNB
                  ZDELTA=ABS(ZLONGSSGP-ZLONGSSGPNB)
                  DO WHILE (ZDELTA >= ZTWOPI)
                    ZDELTA=ZDELTA-ZTWOPI
                  ENDDO
                  IF (ZDELTA > RPI) ZDELTA=ZTWOPI-ZDELTA
                  IF (ZDELTA <= ZDLONGMAX) THEN
                    INRNB=INRNB+1
                    NBSS(IGP,1,INRNB,IGPSS)=IGPNB
                    NBSS(IGP,2,INRNB,IGPSS)=(IYNB-1)*NIJH+JI

                    IF (ZLONGSSGP >= ZTWOPI) ZLONGSSGP=ZLONGSSGP-ZTWOPI
                    IF (ZLONGSSGP <    0._JPRB) ZLONGSSGP=ZLONGSSGP+ZTWOPI
                    IF (ZLONGSSGPNB >= ZTWOPI) ZLONGSSGPNB=ZLONGSSGPNB-ZTWOPI
                    IF (ZLONGSSGPNB <    0._JPRB) ZLONGSSGPNB=ZLONGSSGPNB+ZTWOPI
                    ZDELTA=ZLONGSSGPNB-ZLONGSSGP
                    IEW=NINT(SIGN(1._JPRB,ZDELTA))
                    IF (ABS(ZDELTA) > RPI) THEN
                      IEW=-1*IEW
                    ENDIF

                    NBNS(IGP,INRNB,IGPSS)=-JLL
                    NBEW(IGP,INRNB,IGPSS)=IEW
                  ENDIF
                ENDDO
              ENDIF
            ENDDO
          ENDIF

        ENDDO !JLL
        NRNBSS(IGP,IGPSS)=INRNB
      ENDDO   !JXSS
    ENDDO     !JYSS
    !end of finding small scale grid

  ENDDO ! J
ENDDO
!$OMP END PARALLEL DO

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('YOE_CUCONVCA:CALCULATE_NEIGHBOURS',1,ZHOOK_HANDLE)
END SUBROUTINE CALCULATE_NEIGHBOURS


!     ------------------------------------------------------
!     initialize CA state
!       - set KCELL=KLIVES and KFERT=NFERTYRS where KINI == 1
!       - used for initialization at first timestep and for
!         seeding the CA with "new" cells later on
!     ------------------------------------------------------
SUBROUTINE INITIALIZE_CELLS(YDECUCONVCA,YDGEOMETRY,KINI,KLIVES,KCELL,KFERT)

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK 
USE YOMLUN   , ONLY : NULOUT

IMPLICIT NONE

TYPE(TECUCONVCA), INTENT(INOUT):: YDECUCONVCA
TYPE(GEOMETRY)    , INTENT(IN)    :: YDGEOMETRY
INTEGER(KIND=JPIM), INTENT(IN)    :: KLIVES(YDGEOMETRY%YRGEM%NGPTOT)
INTEGER(KIND=JPIM), INTENT(IN)    :: KINI(YDGEOMETRY%YRGEM%NGPTOT)
INTEGER(KIND=JPIM), INTENT(INOUT) :: KCELL(YDGEOMETRY%YRGEM%NGPTOT,NIJH*NIJH)
INTEGER(KIND=JPIM), INTENT(INOUT) :: KFERT(YDGEOMETRY%YRGEM%NGPTOT,NIJH*NIJH)

INTEGER(KIND=JPIM)               :: JLS,JSS,JK
INTEGER(KIND=JPIM)               :: INIJH2
REAL(KIND=JPHOOK)                  :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('YOE_CUCONVCA:INITIALIZE_CELLS',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP)
ASSOCIATE(NGPTOT=>YDGEM%NGPTOT, &
 & NPROMA=>YDDIM%NPROMA, &
 & NFERTYRS=>YDECUCONVCA%NFERTYRS)

INIJH2=NIJH*NIJH

!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JLS,JSS,JK)
DO JK=1,NGPTOT,NPROMA
  DO JSS=1,INIJH2
    DO JLS=JK,JK+MIN(NPROMA,NGPTOT-JK+1)-1
      IF(KINI(JLS)==1) THEN
        KCELL(JLS,JSS)=KLIVES(JLS)
        KFERT(JLS,JSS)=NFERTYRS
      ENDIF
    ENDDO
  ENDDO
ENDDO
!$OMP END PARALLEL DO

CALL FLUSH(NULOUT)

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('YOE_CUCONVCA:INITIALIZE_CELLS',1,ZHOOK_HANDLE)
END SUBROUTINE INITIALIZE_CELLS


!     ------------------------------------------------------
!     update CA
!       - evolve CA according to the probabilities set in the
!         initial setup
!       - advect CA if LD_ADVECT is set
!       - average to coarse CA grid
!     ------------------------------------------------------
SUBROUTINE UPDCELAUT_RGG(YDGEOMETRY,YDECUCONVCA,YDECUMF,YDSL, & 
 &                       KLIVE,KDX,KDY,KFERT,KCELL,PFERTIN,PCELLIN,PWGHT,PRAND1D,LD_ADVECT)

USE YOECUMF      , ONLY : TECUMF
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK 
USE YOMCT3       , ONLY : NSTEP
USE EINT_MOD     , ONLY : SL_STRUCT

IMPLICIT NONE

TYPE(GEOMETRY)    , INTENT(IN)    :: YDGEOMETRY
TYPE(TECUCONVCA)  , INTENT(INOUT) :: YDECUCONVCA
TYPE(TECUMF)      , INTENT(INOUT) :: YDECUMF
TYPE(SL_STRUCT)   , INTENT(INOUT) :: YDSL
INTEGER(KIND=JPIM), INTENT(INOUT) :: KLIVE(YDGEOMETRY%YRGEM%NGPTOT)                !max number of lives for GP
INTEGER(KIND=JPIM), INTENT(IN)    :: KDX(YDGEOMETRY%YRGEM%NGPTOT),KDY(YDGEOMETRY%YRGEM%NGPTOT)   !wind increments
REAL(KIND=JPRB)   , INTENT(IN)    :: PRAND1D(NIJH*NIJH*YDGEOMETRY%YRGEM%NGPTOT)    !random numbers
REAL(KIND=JPRB)   , INTENT(INOUT) :: PFERTIN(YDSL%NASLB1,NIJH*NIJH)    !fertile status of CA including SL halo
REAL(KIND=JPRB)   , INTENT(INOUT) :: PCELLIN(YDSL%NASLB1,NIJH*NIJH)    !status of CA including SL halo

INTEGER(KIND=JPIM),INTENT(OUT) :: KFERT(YDGEOMETRY%YRGEM%NGPTOT,NIJH*NIJH)      !fertile status of CA
INTEGER(KIND=JPIM),INTENT(OUT) :: KCELL(YDGEOMETRY%YRGEM%NGPTOT,NIJH*NIJH)      !status of CA
REAL(KIND=JPRB)   ,INTENT(OUT) :: PWGHT(YDGEOMETRY%YRGEM%NGPTOT)                !CA smoothed to LS

LOGICAL,INTENT(IN)               :: LD_ADVECT

REAL(KIND=JPHOOK)                  :: ZHOOK_HANDLE

INTEGER(KIND=JPIM)               :: INIJH2
INTEGER(KIND=JPIM)               :: JLS,JSS,JNB
INTEGER(KIND=JPIM)               :: IFERT_UW(YDGEOMETRY%YRDIM%NPROMA)
INTEGER(KIND=JPIM)               :: IFERT_DW(YDGEOMETRY%YRDIM%NPROMA)
INTEGER(KIND=JPIM)               :: ISUM(YDGEOMETRY%YRDIM%NPROMA)

INTEGER(KIND=JPIM)               :: JNBX,JNBMAX,JK, JIND
INTEGER(KIND=JPIM)               :: INBLS, INBSS, IMEAN, INBMAX
INTEGER(KIND=JPIM)               :: IDX, IDY
INTEGER(KIND=JPIM)               :: IEW, INS
INTEGER(KIND=JPIM)               :: ILSSL, IDET

IF (LHOOK) CALL DR_HOOK('YOE_CUCONVCA:UPDCELAUT_RGG',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP)
ASSOCIATE(NGPTOT=>YDGEM%NGPTOT, &
 & NPROMA=>YDDIM%NPROMA, &
 & CA_PROB=>YDECUCONVCA%CA_PROB, LCA_ADVTEST=>YDECUCONVCA%LCA_ADVTEST, &
 & RPROB_SURVIVE=>YDECUCONVCA%RPROB_SURVIVE, RPROB_BIRTH=>YDECUCONVCA%RPROB_BIRTH, &
 & NBSS=>YDECUCONVCA%NBSS, NBNS=>YDECUCONVCA%NBNS, NBEW=>YDECUCONVCA%NBEW, &
 & NRNBSS=>YDECUCONVCA%NRNBSS, NFERTYRS=>YDECUCONVCA%NFERTYRS)
INIJH2=NIJH*NIJH

IF ( CA_PROB == 'GOL' ) THEN
  IDET=0
ELSE
  IDET=1
ENDIF

! advect cells
IF (LD_ADVECT .AND. (.NOT.LCA_ADVTEST)) CALL ADVECT_CA_SOPH(YDGEOMETRY,YDECUCONVCA,YDECUMF,YDSL,KCELL,KFERT,PFERTIN,PCELLIN)

!following lines for testing advection
IF (LCA_ADVTEST .AND. LD_ADVECT .AND. NSTEP > 5) THEN
  CALL ADVECT_CA_SOPH(YDGEOMETRY,YDECUCONVCA,YDECUMF,YDSL,KCELL,KFERT,PFERTIN,PCELLIN)
  IF (NSTEP == 6) THEN
    RPROB_SURVIVE(:,:)=1
    RPROB_BIRTH(:,:)=-1
  ENDIF
ENDIF

! get maximum dimension for easier vectorization
JNBMAX=0
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) &
!$OMP& PRIVATE(JLS,JSS,JK,ILSSL,INBMAX) REDUCTION(MAX:JNBMAX)
DO JK=1,NGPTOT,NPROMA
  INBMAX=0
  DO JSS= 1,INIJH2
    DO JLS= JK,JK+MIN(NPROMA,NGPTOT-JK+1)-1
      ILSSL=YDSL%NSLCORE(JLS)
      INBMAX= MAX(INBMAX,NRNBSS(ILSSL,JSS))
    ENDDO
  ENDDO
  JNBMAX=MAX(INBMAX,JNBMAX)
ENDDO
!$OMP END PARALLEL DO

! evolve CA
!$OMP  PARALLEL DO SCHEDULE(DYNAMIC,1) &
!$OMP& PRIVATE (JK,JLS,JSS,JNB,JNBX,JIND) &
!$OMP& PRIVATE (IMEAN,ILSSL,INBLS,INBSS,IEW,INS,IDX,IDY) &
!$OMP& PRIVATE (IFERT_UW,IFERT_DW,ISUM) FIRSTPRIVATE(IDET)
DO JK= 1,NGPTOT,NPROMA
  DO JSS= 1,INIJH2
    IFERT_UW(:)= 0    !IFERT_UW(1:NPROMA)= 0
    IFERT_DW(:)= 0    !IFERT_DW(1:NPROMA)= 0
    ISUM(:)=0         !ISUM(1:NPROMA)=0

    DO JNBX= 1,JNBMAX
!cdir nodep
      DO JLS= JK,JK+MIN(NPROMA,NGPTOT-JK+1)-1

        JIND=JLS-JK+1

        !count fertile neighbours (upwind and downwind)
        !     ------------------------------------------------------
        !     count those neighbours that haven't lost a life and 
        !     are younger than NFERTYRS
        !     (separate counts for upwind and downwind -> not yet
        !      implemented for redGG CA)
        !     ------------------------------------------------------

        ILSSL=YDSL%NSLCORE(JLS)
        IMEAN=MIN(1,MAX(0,SIGN(JNBX,NRNBSS(ILSSL,JSS)-JNBX))) ! 0 when out of bounds, 1 elsewhere
        JNB=MIN(NRNBSS(ILSSL,JSS),JNBX)
        INBLS=NBSS(ILSSL,1,JNB,JSS)
        INBSS=NBSS(ILSSL,2,JNB,JSS)
        IEW=NBEW(ILSSL,JNB,JSS)
        INS=NBNS(ILSSL,JNB,JSS)
 
        IF ((PFERTIN(INBLS,INBSS) > 0._JPRB) .AND. (IMEAN == 1)) THEN
          IDX=KDX(JLS)
          IDY=KDY(JLS)
          IF (&
            & (IEW == 0 .AND. IDY /= 0 .AND. INS == SIGN(1,IDY))&
            & .OR. (INS == 0 .AND. IDX /= 0 .AND. IEW /= SIGN(1,IDX))&
            & .OR. (INS /= 0 .AND. IEW /= 0 .AND.&
            & ((IDX /= 0 .AND. IDY /= 0 .AND. IEW /= SIGN(1,IDX) .AND.&
            & INS == SIGN(1,IDY))&
            & .OR. (IDX /= 0 .AND. IDY == 0 .AND. IEW /= SIGN(1,IDX) )&
            & .OR. (IDX == 0 .AND. IDY /= 0 .AND. INS == SIGN(1,IDY))))&
            & ) THEN                       !check if neighbour is upwind
            IFERT_UW(JIND)= IFERT_UW(JIND)+1
          ELSE
            IFERT_DW(JIND)= IFERT_DW(JIND)+1
          ENDIF
        ISUM(JIND)=ISUM(JIND)+PCELLIN(INBLS,INBSS)
        ENDIF
      ENDDO ! JLS
    ENDDO ! JNBX

!cdir nodep
    DO JLS= JK,JK+MIN(NPROMA,NGPTOT-JK+1)-1
      JIND=JLS-JK+1
      ! check for living cell
      IF (KCELL(JLS,JSS)/=0) THEN
        IF (PRAND1D((JLS-1)*INIJH2+JSS)<=RPROB_SURVIVE(IFERT_DW(JIND),IFERT_UW(JIND))) THEN ! survival test
          !KCELL(JLS,JSS)=KCELL(JLS,JSS) ! deterministic GOL rules
          KFERT(JLS,JSS)= MAX(KFERT(JLS,JSS)-IDET,0)
        ELSE  !  cells failing the survival test -> lose a life
          KCELL(JLS,JSS)=  KCELL(JLS,JSS) - 1
          KFERT(JLS,JSS)=  (1-IDET)*(KFERT(JLS,JSS)-1)
        ENDIF
      ELSE        ! dead cell
        IF (PRAND1D((JLS-1)*INIJH2+JSS)<=RPROB_BIRTH(IFERT_DW(JIND),IFERT_UW(JIND))) THEN  ! check for birth condition
          KCELL(JLS,JSS)= KLIVE(JLS)
          KFERT(JLS,JSS)= (1-IDET)* KLIVE(JLS)&
           & -IDET*MIN(0,SIGN(NFERTYRS,-KLIVE(JLS)))
        ELSE                    !  cells failing the birth condition
          KCELL(JLS,JSS)=  0
          KFERT(JLS,JSS)=  0
        ENDIF
      ENDIF
    ENDDO ! JLS

  ENDDO ! JSS
ENDDO ! JK

!$OMP END PARALLEL DO

! average blocks of cells to a coarser grid
CALL WEIGHTING_FIELD(YDGEOMETRY,YDECUCONVCA,YDSL,KCELL,PWGHT)
    
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('YOE_CUCONVCA:UPDCELAUT_RGG',1,ZHOOK_HANDLE)
END SUBROUTINE UPDCELAUT_RGG


!     ------------------------------------------------------
!     - average over NIJHxNIJH cell blocks, smooth and weight
!     ------------------------------------------------------
SUBROUTINE WEIGHTING_FIELD(YDGEOMETRY,YDECUCONVCA,YDSL,KCELL,PWGHT)
!     ------------------------------------------------------

!     square the KCELL values (on range 0->31); average over
!     4x4 KCELL blocks, and then normalize  
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK 
USE EINT_MOD     , ONLY : SL_STRUCT
USE ORDER_INDEPENDENT_SUMMATION_MOD, ONLY : ORDER_INDEP_GLOBAL_SUM

IMPLICIT NONE

TYPE(GEOMETRY)    , INTENT(IN)    :: YDGEOMETRY
TYPE(TECUCONVCA)  , INTENT(INOUT) :: YDECUCONVCA
TYPE(SL_STRUCT)   , INTENT(INOUT) :: YDSL
INTEGER(KIND=JPIM), INTENT(IN)    :: KCELL(YDGEOMETRY%YRGEM%NGPTOT,NIJH*NIJH)
REAL(KIND=JPRB)   , INTENT(INOUT) :: PWGHT(YDGEOMETRY%YRGEM%NGPTOT)
INTEGER(KIND=JPIM)               :: JLS,JSS,JK
INTEGER(KIND=JPIM)               :: INIJH2
REAL(KIND=JPRB)                  :: ZDIV, ZSUM
REAL(KIND=JPHOOK)                  :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('YOE_CUCONVCA:WEIGHTING_FIELD',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP)
ASSOCIATE(NGPTOTG=>YDGEM%NGPTOTG, NGPTOT=>YDGEM%NGPTOT, &
 & NPROMA=>YDDIM%NPROMA, &
 & LCA_GLOBAL=>YDECUCONVCA%LCA_GLOBAL, LCA_SMOOTH=>YDECUCONVCA%LCA_SMOOTH)
INIJH2=NIJH*NIJH
ZDIV=1._JPRB/REAL(INIJH2,JPRB)

!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JLS,JK)
DO JK=1,NGPTOT,NPROMA
  DO JLS=JK,JK+MIN(NPROMA,NGPTOT-JK+1)-1
    PWGHT(JLS)=0.
  ENDDO
ENDDO
!$OMP END PARALLEL DO

!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JLS,JSS,JK)
DO JK=1,NGPTOT,NPROMA
  DO JSS=1,INIJH2
    DO JLS=JK,JK+MIN(NPROMA,NGPTOT-JK+1)-1
      PWGHT(JLS)= PWGHT(JLS)+ZDIV*KCELL(JLS,JSS)
    ENDDO
  ENDDO
ENDDO
!$OMP END PARALLEL DO

IF (LCA_SMOOTH) CALL SMOOTH_CA(YDGEOMETRY,YDECUCONVCA,YDSL,PWGHT)

IF (LCA_GLOBAL) THEN    !scale global field to zero mean

  ZSUM=ORDER_INDEP_GLOBAL_SUM(PWGHT)
  ZSUM=ZSUM/REAL(NGPTOTG,JPRB)

  IF (ZSUM /= 0._JPRB ) THEN
    !$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JK,JLS)
    DO JK=1,NGPTOT,NPROMA
      DO JLS=JK,JK+MIN(NPROMA,NGPTOT-JK+1)-1
        PWGHT(JLS)=PWGHT(JLS)/ZSUM - 1._JPRB
      ENDDO
    ENDDO
    !$OMP END PARALLEL DO
  ENDIF

ENDIF

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('YOE_CUCONVCA:WEIGHTING_FIELD',1,ZHOOK_HANDLE)
END SUBROUTINE WEIGHTING_FIELD


!     ------------------------------------------------------
!       apply smoothing
!     ------------------------------------------------------
SUBROUTINE SMOOTH_CA(YDGEOMETRY,YDECUCONVCA,YDSL,PFLD)

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK 

!stuff for SL-halo
USE EINT_MOD , ONLY : SL_STRUCT

IMPLICIT NONE
TYPE(GEOMETRY)   , INTENT(IN)    :: YDGEOMETRY
TYPE(TECUCONVCA) , INTENT(INOUT) :: YDECUCONVCA
TYPE(SL_STRUCT)  , INTENT(INOUT) :: YDSL
REAL(KIND=JPRB)  , INTENT(INOUT) :: PFLD(YDGEOMETRY%YRGEM%NGPTOT)
INTEGER(KIND=JPIM)               :: JLS,JNB,INB,JK
REAL(KIND=JPRB)                  :: ZFLD_SMOOTH(YDGEOMETRY%YRGEM%NGPTOT)
REAL(KIND=JPHOOK)                  :: ZHOOK_HANDLE

!stuff for SL-halo
LOGICAL                          :: LLINC
INTEGER(KIND=JPIM)               :: IDUMARR(2)
INTEGER(KIND=JPIM)               :: IFIXSFLD(2)
INTEGER(KIND=JPIM)               :: IFLDSLB1
REAL(KIND=JPRB), ALLOCATABLE     :: ZPB1(:,:)

#include "slcomm.intfb.h"

IF (LHOOK) CALL DR_HOOK('YOE_CUCONVCA:SMOOTH_CA',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP)
ASSOCIATE(NGPTOT=>YDGEM%NGPTOT, &
 & NPROMA=>YDDIM%NPROMA, &
 & NRNBLS=>YDECUCONVCA%NRNBLS, NBLS=>YDECUCONVCA%NBLS)
!*       3.2  get halo for CA
!     ------------------------------------------------------

LLINC=.FALSE.
IFIXSFLD(:)=0
IFLDSLB1=1
ALLOCATE(ZPB1(YDSL%NASLB1,IFLDSLB1))
ZPB1=0

!fill buffer
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JLS,JK)
DO JK=1,NGPTOT,NPROMA
  DO JLS=JK,JK+MIN(NPROMA,NGPTOT-JK+1)-1
    ZPB1(YDSL%NSLCORE(JLS),1)=PFLD(JLS)
  ENDDO
ENDDO
!$OMP END PARALLEL DO

!get halo
CALL SLCOMM(YDSL,IDUMARR,IFLDSLB1,LLINC,0,ZPB1)

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JLS,INB,JNB)
DO JLS=1,NGPTOT
  ZFLD_SMOOTH(JLS)=4._JPRB*PFLD(JLS)
  INB=4
  DO JNB=1,NRNBLS(YDSL%NSLCORE(JLS))
    ZFLD_SMOOTH(JLS)=ZFLD_SMOOTH(JLS)+2._JPRB*ZPB1(NBLS(YDSL%NSLCORE(JLS),JNB),1)
    INB=INB+2
  ENDDO
  ZFLD_SMOOTH(JLS)=ZFLD_SMOOTH(JLS)/(REAL(INB,JPRB))
  PFLD(JLS)= ZFLD_SMOOTH(JLS)
ENDDO
!$OMP END PARALLEL DO

DEALLOCATE(ZPB1)

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('YOE_CUCONVCA:SMOOTH_CA',1,ZHOOK_HANDLE)
END SUBROUTINE SMOOTH_CA


!     ------------------------------------------------------
!     perform advection (using departure points from SL-advection)
!     ------------------------------------------------------
SUBROUTINE ADVECT_CA_SOPH(YDGEOMETRY,YDECUCONVCA,YDECUMF,YDSL,KCELL,KFERT,PFERTIN,PCELLIN)

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK 
USE YOMLUN   , ONLY : NULERR
USE YOECUMF  , ONLY : TECUMF
USE YOMCST   , ONLY : RPI
USE EINT_MOD , ONLY : SL_STRUCT

IMPLICIT NONE

TYPE(GEOMETRY)    , INTENT(IN)    :: YDGEOMETRY
TYPE(TECUCONVCA)  , INTENT(INOUT) :: YDECUCONVCA
TYPE(TECUMF)      , INTENT(INOUT) :: YDECUMF
TYPE(SL_STRUCT)   , INTENT(INOUT) :: YDSL
REAL(KIND=JPRB)   , INTENT(INOUT) :: PFERTIN(YDSL%NASLB1,NIJH*NIJH)    !fertile status of CA including SL halo
REAL(KIND=JPRB)   , INTENT(INOUT) :: PCELLIN(YDSL%NASLB1,NIJH*NIJH)    !status of CA including SL halo
INTEGER(KIND=JPIM), INTENT(INOUT) :: KFERT(YDGEOMETRY%YRGEM%NGPTOT,NIJH*NIJH)      !fertile status of CA
INTEGER(KIND=JPIM), INTENT(INOUT) :: KCELL(YDGEOMETRY%YRGEM%NGPTOT,NIJH*NIJH)

INTEGER(KIND=JPIM)                :: IDEPLS(YDGEOMETRY%YRGEM%NGPTOT,NIJH*NIJH)
INTEGER(KIND=JPIM)                :: IDEPSS(YDGEOMETRY%YRGEM%NGPTOT,NIJH*NIJH)

INTEGER(KIND=JPIM)                :: ICEND,IBL, IOFF
INTEGER(KIND=JPIM)                :: IGP, IJKT, IGPDEP
INTEGER(KIND=JPIM)                :: IXSS,IYSS
INTEGER(KIND=JPIM)                :: JK, JI, JLAT, JLS, JSS
INTEGER(KIND=JPIM)                :: JXSS,JYSS
REAL(KIND=JPRB)                   :: ZLATDEP, ZLONDEP
INTEGER(KIND=JPIM)                :: ILATDEP, ILATDEPG

INTEGER(KIND=JPIM)                :: ILONG, ILON
REAL(KIND=JPRB)                   :: ZLATBOUND(1:YDGEOMETRY%YRDIM%NDGLG-1) !,ZDLATG(1:YDGEOMETRY%YRDIM%NDGLG)
REAL(KIND=JPRB)                   :: ZDLONG, ZRADDEG
REAL(KIND=JPRB)                   :: ZDLONGDEP, ZDLATDEP

REAL(KIND=JPRB)   ,ALLOCATABLE    :: ZLONDEPC(:,:,:),ZLDLAT(:),ZLDLONG(:)
REAL(KIND=JPRB)   ,ALLOCATABLE    :: ZLATDEPC(:,:,:)

REAL(KIND=JPHOOK)                   :: ZHOOK_HANDLE

!stuff for SL-halo
LOGICAL                           :: LLINC
INTEGER(KIND=JPIM)                :: IDUMARR(2)
INTEGER(KIND=JPIM)                :: IFIXSFLD(2)
INTEGER(KIND=JPIM)                :: IFLDSLB1
REAL(KIND=JPRB), ALLOCATABLE      :: ZPB1(:,:)

#include "slcomm.intfb.h"

IF (LHOOK) CALL DR_HOOK('YOE_CUCONVCA:ADVECT_CA_SOPH',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP,  &
 & YDGSGEOM_NB=>YDGEOMETRY%YRGSGEOM_NB, YDCSGLEG=>YDGEOMETRY%YRCSGLEG,  &
 & YDVAB=>YDGEOMETRY%YRVAB, YDVETA=>YDGEOMETRY%YRVETA, YDVFE=>YDGEOMETRY%YRVFE, YDSTA=>YDGEOMETRY%YRSTA,  &
 & YDLAP=>YDGEOMETRY%YRLAP, &
 & YDCSGEOM=>YDGEOMETRY%YRCSGEOM, YDCSGEOM_NB=>YDGEOMETRY%YRCSGEOM_NB,  &
 & YDSPGEOM=>YDGEOMETRY%YSPGEOM)
ASSOCIATE(NJKT4=>YDECUMF%NJKT4, NJKT5=>YDECUMF%NJKT5, &
 & NLOENG=>YDGEM%NLOENG, NGPTOT=>YDGEM%NGPTOT, &
 & NDGLG=>YDDIM%NDGLG, NGPBLKS=>YDDIM%NGPBLKS, &
 & NPROMA=>YDDIM%NPROMA, &
 & NFLEVG=>YDDIMV%NFLEVG, &
 & MYLATS=>YDMP%MYLATS, &
 & RLONDEP=>YDECUCONVCA%RLONDEP, RLATDEP=>YDECUCONVCA%RLATDEP, &
 & SLLONG=>YDECUCONVCA%SLLONG, SLLAT=>YDECUCONVCA%SLLAT, &
 & SLDLONG=>YDECUCONVCA%SLDLONG, SLDLAT=>YDECUCONVCA%SLDLAT, &
 & SLDDLAT=>YDECUCONVCA%SLDDLAT)
ZRADDEG=180._JPRB/RPI

IJKT=(NJKT4+NJKT5)/2    !get level in the middle of 850 and 500 hPa

ALLOCATE(ZLONDEPC(NPROMA,NFLEVG,NGPBLKS),ZLATDEPC(NPROMA,NFLEVG,NGPBLKS))
ALLOCATE(ZLDLAT(YDSL%NASLB1),ZLDLONG(YDSL%NASLB1))

ZLONDEPC=RLONDEP
ZLATDEPC=RLATDEP

ZLDLAT=SLDLAT
ZLDLONG=SLDLONG

!calculate latitude boundaries
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JLAT)
DO JLAT=1,NDGLG-1
  ZLATBOUND(JLAT)=0.5_JPRB*(YDCSGLEG%RLATIG(JLAT)+YDCSGLEG%RLATIG(JLAT+1))
ENDDO
!$OMP END PARALLEL DO

! work out departure GP in SL-halo array from lon/lat of departure point and do advection
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) &
!$OMP& PRIVATE(JK,ICEND,IBL,IOFF,JI,IGP,JYSS,JXSS,JSS,ZLATDEP,ZLONDEP)&
!$OMP& PRIVATE(ILATDEPG,ILATDEP,ZDLONG,ILONG,ILON,IGPDEP,ZDLONGDEP,ZDLATDEP,IXSS,IYSS)
DO JK=1,NGPTOT,NPROMA                           !
  ICEND=MIN(NPROMA,NGPTOT-JK+1)                 !
  IBL=(JK-1)/NPROMA+1                           ! loop over large scale gridpoints
  IOFF=JK                                       !

  DO JYSS=1,NIJH                          ! loop over small scale gridpoints
    DO JXSS=1,NIJH
      JSS=(JYSS-1)*NIJH+JXSS

!cdir nodep
      DO JI=1,ICEND                                 !
        IGP=IOFF+JI-1

        ZLATDEP=RLATDEP(JI,IJKT,IBL)+(-JYSS+0.5)*SLDLAT(YDSL%NSLCORE(IGP))/REAL(NIJH,JPRB)+SLDDLAT(YDSL%NSLCORE(IGP))
        ZLONDEP=RLONDEP(JI,IJKT,IBL)+(JXSS-0.5*NIJH-0.5)*SLDLONG(YDSL%NSLCORE(IGP))/REAL(NIJH,JPRB)

        IF (ZLONDEP >= 2*RPI) ZLONDEP=ZLONDEP-2*RPI
        IF (ZLONDEP < 0) ZLONDEP=ZLONDEP+2*RPI
        IF (ZLATDEP > 0.5*RPI) ZLATDEP=RPI-ZLATDEP
        IF (ZLATDEP < -0.5*RPI) ZLATDEP=-RPI-ZLATDEP

        !find global lat-index of dep point
        ILATDEPG=1
        DO WHILE (ZLATDEP < ZLATBOUND(ILATDEPG) .AND. ILATDEPG < NDGLG)!can this be replaced by something similar 
                                                                       !to ilon calc below?
          ILATDEPG=ILATDEPG+1
        ENDDO
        !find local lat-index of dep point
        ILATDEP=ILATDEPG-MYLATS(1)+1
        ZDLONG=2*RPI/REAL(NLOENG(ILATDEPG),JPRB)         !calculate delta longitude between gp on lat
        ILONG=NINT(ZLONDEP/ZDLONG)
        IF (YDSL%NSLSTA(ILATDEP) > ILONG+1) ILONG=ILONG+NLOENG(ILATDEPG)
        ILON=ILONG-YDSL%NSLSTA(ILATDEP)+2
        IF (ILON > YDSL%NSLONL(ILATDEP)) ILON=ILON-NLOENG(ILATDEPG)         !???? correct ??????

        IGPDEP=YDSL%NSLOFF(ILATDEP)+ILON

        IDEPLS(IGP,JSS)=IGPDEP

        ZDLONGDEP=ZLONDEP-(SLLONG(IGPDEP)-0.5*SLDLONG(IGPDEP))
        ZDLATDEP=(SLLAT(IGPDEP)+SLDDLAT(IGPDEP))-ZLATDEP

        IF (ZDLONGDEP >= 2*RPI) ZDLONGDEP=ZDLONGDEP-2*RPI

        IF (ZDLONGDEP>=SLDLONG(IGPDEP) .OR. ZDLONGDEP<0 .OR. ZDLATDEP>=SLDLAT(IGPDEP) .OR. ZDLATDEP<0) THEN
           WRITE(NULERR,'(A,4I6,6F9.3)') "dddd INCONSISTENCY: ",IGP,YDSL%NSLCORE(IGP),JSS,IGPDEP,ZLONDEP*ZRADDEG,&
                                 & ZLATDEP*ZRADDEG,ZDLONGDEP*ZRADDEG,ZDLATDEP*ZRADDEG,&
                                 & SLDLONG(IGPDEP)*ZRADDEG,SLDLAT(IGPDEP)*ZRADDEG
           WRITE(NULERR,'(6I6,4F9.3)') ILONG,ILON,  NLOENG(ILATDEPG),IGPDEP,YDSL%NSLONL(ILATDEP),YDSL%NSLSTA(ILATDEP),&
                                 & RLONDEP(JI,IJKT,IBL)*ZRADDEG,&
                                 & SLDLONG(YDSL%NSLCORE(IGP))*ZRADDEG,SLLONG(IGPDEP)*ZRADDEG,SLDLONG(IGPDEP)*ZRADDEG
        ENDIF

        IXSS=1
        DO WHILE (ZDLONGDEP > IXSS*ZDLONG/REAL(NIJH,JPRB))
          IXSS=IXSS+1
        ENDDO
        IYSS=1
        DO WHILE (ZDLATDEP > IYSS*SLDLAT(IGPDEP)/REAL(NIJH,JPRB))
          IYSS=IYSS+1
        ENDDO
        IDEPSS(IGP,JSS)=(IYSS-1)*NIJH+IXSS

        ! do the advection (arrays IDEPSS and IDEPLS could be saved)
        KFERT(IGP,JSS)=PFERTIN(IDEPLS(IGP,JSS),IDEPSS(IGP,JSS))
        KCELL(IGP,JSS)=PCELLIN(IDEPLS(IGP,JSS),IDEPSS(IGP,JSS))

      ENDDO ! JI

    ENDDO ! JXSS
  ENDDO ! JYSS
ENDDO ! JK
!$OMP END PARALLEL DO

! get the updated halo
LLINC=.FALSE.
IFIXSFLD(:)=0
IFLDSLB1=2*NIJH*NIJH
ALLOCATE(ZPB1(YDSL%NASLB1,IFLDSLB1))
ZPB1=0

!fill buffer
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JK,JLS,JSS)
DO JK=1,NGPTOT,NPROMA
  DO JSS=1,NIJH*NIJH
    DO JLS=JK,JK+MIN(NPROMA,NGPTOT-JK+1)-1
      ZPB1(YDSL%NSLCORE(JLS),JSS)=KCELL(JLS,JSS)
      ZPB1(YDSL%NSLCORE(JLS),JSS+NIJH*NIJH)=KFERT(JLS,JSS)
    ENDDO
  ENDDO
ENDDO
!$OMP END PARALLEL DO

!get halo
CALL SLCOMM(YDSL,IDUMARR,IFLDSLB1,LLINC,0,ZPB1)

!extract from buffer
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JK,JLS,JSS)
DO JK=1,YDSL%NASLB1,NPROMA
  DO JSS=1,NIJH*NIJH
    DO JLS=JK,JK+MIN(NPROMA,YDSL%NASLB1-JK+1)-1
      PCELLIN(JLS,JSS)=ZPB1(JLS,JSS)
      PFERTIN(JLS,JSS)=ZPB1(JLS,JSS+NIJH*NIJH)
    ENDDO
  ENDDO
ENDDO
!$OMP END PARALLEL DO
DEALLOCATE(ZPB1)
DEALLOCATE(ZLONDEPC,ZLATDEPC)
DEALLOCATE(ZLDLAT,ZLDLONG)

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('YOE_CUCONVCA:ADVECT_CA_SOPH',1,ZHOOK_HANDLE)
END SUBROUTINE ADVECT_CA_SOPH

!     ------------------------------------------------------
!     binary 2D output for debuging
!     ------------------------------------------------------
SUBROUTINE WRITE_FIELD(CDFILENAME, PFIELD, KPOINTS, KLEVS)

USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN  , ONLY : NULOUT

IMPLICIT NONE
INTEGER(KIND=JPIM),INTENT(IN)      :: KPOINTS, KLEVS
CHARACTER(LEN=*)  ,INTENT(IN)      :: CDFILENAME
REAL(KIND=JPRB)   ,INTENT(IN)      :: PFIELD(KPOINTS,KLEVS)
REAL(KIND=JPHOOK)                    :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('YOE_CUCONVCA:WRITE_FIELD',0,ZHOOK_HANDLE)
WRITE(NULOUT,*) 'Writing file ',CDFILENAME

OPEN(1142, FILE=CDFILENAME, FORM="unformatted", STATUS="REPLACE") 

WRITE(1142) KPOINTS
WRITE(1142) KLEVS
WRITE(1142) PFIELD

CLOSE(1142)

IF (LHOOK) CALL DR_HOOK('YOE_CUCONVCA:WRITE_FIELD',1,ZHOOK_HANDLE)
END SUBROUTINE WRITE_FIELD


END MODULE YOE_CUCONVCA
