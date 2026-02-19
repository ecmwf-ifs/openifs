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

MODULE YOMDIMO

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!**----------------------------------------------------------------
!**  DIMENSIONS OF ARRAYS RELATED TO OBSERVATIONS.
!**  DIMENSIONS DES TABLES RELIEES AUX OBSERVATIONS
!*
!*  OBSERVATIONS ARRAYS DIMENSIONS
!*  ------------------------------
!*    NOBTOT  : ACTUAL NUMBER OF OBSERVATIONS IN ROBSAR
!*    NOBTOTG : IDEM FOR SUM FOR ALL PROCESSORS
!*    NOBSLN  : ACTUAL LENGTH OF ROBSAR
!*    NACTIM  : ACTUAL NUMBER OF TIME SLOTS
!*
!*  OBSERVATIONS PROCESSING
!*  -----------------------
!*    NMXSOBT : MAXIMUM NUMBER OF SUB-OBSERVATIONS TYPES
!*    NMXAREA : MAXIMUM NUMBER OF AREAS
!*    NMXOBS  : MAXIMUM NUMBER OF OBSERVATIONS IN ROBSAR
!*    NMXSET  : MAXIMUM NUMBER OF OBSERVATIONS SET 
!*    NMXLEN  : MAXIMUM LENGTH OF ONE OBSERVATION SET
!*    NMXLEN  : MAXIMUM AMOUNT OF MEMORY USED FOR ONE OBSSET
!*    L_SPAN_KSETS_ACROSS_TSLOTS : .TRUE.  - span KSETs across timeslots (ecset.F90)
!*                                 .FALSE. - do not span KSETs across timeslots (ecset.F90)
!*    L_TASKOB_USE_SORTED_TIMES   : .TRUE. use dynamic (OpenMP) load balancing in TASKOB   (default)
!*    L_TASKOBAD_USE_SORTED_TIMES : .TRUE. use dynamic (OpenMP) load balancing in TASKOBAD (default)
!*    L_TASKOBTL_USE_SORTED_TIMES : .TRUE. use dynamic (OpenMP) load balancing in TASKOBTL (default)
!*    L_SORT_OBS_SETS_BY_TIMING   : .TRUE. use dynamic (OpenMP) load balancing in OOPS TL/AD obs operator (default)
!*    NOBTOV  : LOCAL NUMBER OF TOVS OBSERVATIONS 
!*    NOBTOVG : GLOBAL NUMBER OF TOVS OBSERVATIONS
!*    NOBTOVBP: NUMBER OF TOVS OBSERVATIONS, BY PROCESSOR
!*    NOBSCA  : LOCAL NUMBER OF SCATT. DATA 
!*    NOBNTV  : ACTUAL NUMBER OF NON-(TOVS/SCATT) OBSERVATIONS IN ROBSAR
!*
!*    NCSTO   : number of constant fields at obs points (CANARI and ECMWF)
!*
!*    NMXTCH  : NUMBER OF TOVS CHANNELS
!*    NMXTCHU : MAX NUMBER OF TOVS CHANNELS TO USE IN VAR ANALYSIS.

!*
!*  POSITIONING IN THE ARRAY OF MOD VAR AT OBS. POINTS
!*  --------------------------------------------------
!*    MULTI LEVEL FIELDS FOR NON-(TOVS/SCATT) OBS. 
!*      POSITION   POINTEE      DESCRIPTION
!*         0        GOMMV     GENERIC
!*         1        GOMU      U CONTRAVARIANT COMPONENT OF THE WIND
!*         2        GOMV      V      "           "               "
!*         3        GOMT      TEMPERATURE
!*         4        GOMQ      VAPOUR PHASE OF WATER
!*         5        GOMLW     CLOUD LIQUID WATER
!*         6        GOMCL     CLOUD COVER
!*
!*    SINGLE LEVEL FIELDS FOR NON-(TOVS/SCATT) OBS.
!*      POSITION   POINTEE      DESCRIPTION
!*         0        GOSMV     GENERIC
!*         1        GOSP      SURFACE PRESSURE
!*         2        GOSTS     SURFACE TEMPERATURE
!*         3        GOSWS     SOIL MOISTURE
!*         4        GOSSN     SNOW DEPTH
!*         5        GOSZ0     ROUGHNESS LENGTH
!*         6        GOSWL     VEG. SKIN RESERVOIR
!*
!*    MULTI LEVEL FIELDS FOR TOVS OBS.
!*      POSITION   POINTEE      DESCRIPTION
!*         0        GSMMV     GENERIC
!*         1        GSMT      TEMPERATURE
!*         2        GSMQ      VAPOUR PHASE OF WATER
!*         3        GSMLW     CLOUD LIQUID WATER
!*         4        GSMCL     CLOUD COVER
!*
!*    SINGLE LEVEL FIELDS FOR TOVS OBS.
!*      POSITION   POINTEE      DESCRIPTION
!*         0        GSSMV     GENERIC
!*         1        GSSP      SURFACE PRESSURE
!*         2        GSSTS     SURFACE TEMPERATURE
!*         3        GSSWS     SOIL MOISTURE
!*         4        GSSSN     SNOW DEPTH
!*         5        GSSZ0     ROUGHNESS LENGTH
!*         6        GSSWL     VEG. SKIN RESERVOIR
!*
!*    FIELDS FOR SCATTEROMETER OBS.
!*      POSITION   POINTEE      DESCRIPTION
!*         0        GSCMV     GENERIC
!*         1        GSCU      U WIND AT LOWEST MOD. LEV.
!*         2        GSCV      V WIND AT LOWEST MOD. LEV.
!*         3        GSCT      TEMPE  AT LOWEST MOD. LEV.
!*         4        GSCQ      SPECIF HUMID AT LOWEST MOD. LEV.
!*         5        GSCPS     SURFACE PRESSURE
!*         6        GSCTS     SURFACE TEMPERATURE
!*
!*
!*  DIMENSIONS DES TABLEAUX D'OBSERVATIONS
!*  --------------------------------------
!*    NOBTOT  : NOMBRE D'OBSERVATIONS EFFECTIVEMENT DANS ROBSAR
!*    NOBSLN  : LONGUEUR EFFECTIVE DU TABLEAU ROBSAR
!*    NACTIM  : NOMBRE DE FOURCHETTES TEMPORELLES
!*
!*  TRAITEMENT DES OBSERVATIONS
!*  ---------------------------
!*    NMXSOBT : NOMBRE MAXIMUM DE TYPES DE "SUB-OBSERVATIONS"
!*    NMXAREA : NOMBRE MAXIMUM DE ZONES
!*    NMXOBS  : NOMBRE MAXIMUM D'OBSERVATIONS DANS ROBSAR
!*    NMXSET  : NOMBRE MAXIMUM DE PAQUETS D'OBSERVATIONS (SAUF TOVS)
!*    NMXLEN  : LONGUEUR MAXIMALE D'UN PAQUET D'OBSERVATIONS (SAUF TOVS)
!*    NOBTOV  : NOMBRE EFFECTIF DE TOVS DANS ROBSAR
!*    NOBTOVG : NOMBRE EFFECTIF DE TOVS DANS ROBSAR GLOBALE
!*    NOBSCA  : NOMBRE EFFECTIF DE DONNEES SCATT. DANS ROBSAR
!*    NOBNTV  : NOMBRE EFFECTIF D'OBSERVATIONS NON TOVS DANS ROBSAR
!*
!*    NMXTCH  : NBRE DE CANAUX TOVS
!*    NMXTCHU : NBRE MAX DE CANAUX TOVS UTILISES DANS L'ANALYSE VAR.
!*-----------------------------------------------------------------

!  Position    Pointee         Description
!     0        ECSTO           generic
!     1        ECALBE          albedo
!     2        ECCORI          CORIOLIS
!     3        ECEMIS          emissivity
!     4        ECLSMA          land/sea mask
!     5        ECOROG          orography
!     9        ECARG           percentage of clay within the soil
!    10        ECSAB           percentage of sand within the soil
!    11        ECHV            resistance to evapotranspiration
!    12        ECZ0H           thermal roughness length *g

INTEGER(KIND=JPIM) :: NOBTOT
INTEGER(KIND=JPIM) :: NOBSLN
INTEGER(KIND=JPIM) :: NOBTOTG
INTEGER(KIND=JPIM) :: NACTIM
INTEGER(KIND=JPIM) :: NMXSOBT
INTEGER(KIND=JPIM) :: NMXAREA
INTEGER(KIND=JPIM) :: NMXSET
INTEGER(KIND=JPIM) :: NMXLEN
INTEGER(KIND=JPIM) :: NOBTOV
INTEGER(KIND=JPIM) :: NOBTOVG
INTEGER(KIND=JPIM), ALLOCATABLE :: NOBTOVBP(:)
INTEGER(KIND=JPIM) :: NOBSCA
INTEGER(KIND=JPIM) :: NOBNTV
INTEGER(KIND=JPIM) :: NOBNTVG
INTEGER(KIND=JPIM), ALLOCATABLE :: NOBNTVBP(:)
INTEGER(KIND=JPIM) :: NMXTCH
INTEGER(KIND=JPIM) :: NMXTCHU

LOGICAL :: L_SPAN_KSETS_ACROSS_TSLOTS
LOGICAL :: L_TASKOB_USE_SORTED_TIMES   = .TRUE.
LOGICAL :: L_TASKOBAD_USE_SORTED_TIMES = .TRUE.
LOGICAL :: L_TASKOBTL_USE_SORTED_TIMES = .TRUE.
LOGICAL :: L_SORT_OBS_SETS_BY_TIMING   = .TRUE.  ! Used in OOPS TL/AD observation operator
!     ------------------------------------------------------------------

END MODULE YOMDIMO
