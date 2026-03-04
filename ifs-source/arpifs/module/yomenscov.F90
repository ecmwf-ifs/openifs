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

MODULE YOMENSCOV
  USE PARKIND1 , ONLY : JPIM, JPRB
  IMPLICIT NONE
  SAVE
  PUBLIC
  !
  !     Author.
  !      -------
  !      Yann Michel, METEO-FRANCE/CNRM/GMAP, 12th June 2018
  !
  !* -----------!
  !* GRIDPOINT *!
  !* -----------!
  !* A structure for RF-based separable correlation model
  !* ------------------------------------------------------------
  !* 3d variable
  TYPE TYPE_GPLOC
     !* vertical
     INTEGER(KIND=JPIM)                 :: NPASSESV
     INTEGER(KIND=JPIM)                 :: NORDERV
     REAL(KIND=JPRB),      ALLOCATABLE  :: ALPHAV(:)
     REAL(KIND=JPRB),      ALLOCATABLE  :: MATGV(:,:)
     REAL(KIND=JPRB),      ALLOCATABLE  :: INVDEFV(:)
     REAL(KIND=JPRB),      ALLOCATABLE  :: DEFV(:)
     !* horizontal
     INTEGER(KIND=JPIM)                 :: NPASSESH
     INTEGER(KIND=JPIM)                 :: NORDERH
     REAL(KIND=JPRB),      ALLOCATABLE  :: ALPHAH(:,:)
     REAL(KIND=JPRB),      ALLOCATABLE  :: MATGH(:,:,:)
     !* normalization
     REAL(KIND=JPRB),      ALLOCATABLE  :: NORM(:)
  END TYPE TYPE_GPLOC
  !* ------------------------------------------------------------
  !* The localization itself
  !* ------------------------------------------------------------
  TYPE(TYPE_GPLOC) :: GP_LOC

  !* ------------------------------------------------------------
  !* A structure for ens data
  !* ------------------------------------------------------------
  TYPE TYPE_ENSDATA
     REAL(KIND=JPRB), POINTER :: GP3ENS(:,:,:,:,:)       => NULL()
     REAL(KIND=JPRB), POINTER :: GP2ENS(:,:,:,:)         => NULL()
  END TYPE TYPE_ENSDATA
  !* ------------------------------------------------------------
  !* The data itself
  !* ------------------------------------------------------------
  TYPE(TYPE_ENSDATA) :: ENS
  
  !* ------------------------------------------------------------
  !* Flags to debug/test
  !* ------------------------------------------------------------
  LOGICAL :: LAPPLY_LOC_VERT, LAPPLY_LOC_HORZ, LENSPERT2ONE
  !* ------------------------------------------------------------
  !* change of spatial resolution
  LOGICAL :: LDUAL_RES
  INTEGER(KIND=JPIM) :: NSMAX_ENS, NMSMAX_ENS
  !* ------------------------------------------------------------
  !* Flags to the localization method (spectral/gridpoint, u/v or vor/div)
  !* ------------------------------------------------------------
  LOGICAL :: LGPLOC, LSPLOC
  LOGICAL :: LUVLOC
  !* ------------------------------------------------------------
  !* Length-scales and weights
  !* ------------------------------------------------------------
  REAL(KIND=JPRB) :: HORZ_LOC
  REAL(KIND=JPRB) :: VERT_LOC
  REAL(KIND=JPRB) :: ENS_WEIGHT
  REAL(KIND=JPRB) :: STC_WEIGHT
  !* ------------------------------------------------------------
  !* Ensemble sizes
  !* ------------------------------------------------------------
  INTEGER(KIND=JPIM) :: NENS
  INTEGER(KIND=JPIM) :: NLAGGED
  INTEGER(KIND=JPIM) :: ILAG ! for openfainfo
  !* ----------!
  !* SPECTRAL *!
  !* ----------!
  !* ------------------------------------------------------------
  !* Structure for LOCALIZATION METHODS
  !* ------------------------------------------------------------
  TYPE TYPE_SPLOC
     !* METHOD 1
     REAL(KIND=JPRB),      ALLOCATABLE  :: HLOCOR(:,:)
     REAL(KIND=JPRB),      ALLOCATABLE  :: EVALS(:)
     REAL(KIND=JPRB),      ALLOCATABLE  :: EVECS(:,:)
     REAL(KIND=JPRB),      ALLOCATABLE  :: VLOCOR(:,:)
     REAL(KIND=JPRB),      ALLOCATABLE  :: NORM(:)
     REAL(KIND=JPRB),      ALLOCATABLE  :: NORM_SQRT(:)
  END TYPE TYPE_SPLOC
  !* ------------------------------------------------------------
  !* The localization itself
  !* ------------------------------------------------------------
  TYPE(TYPE_SPLOC) :: SP_LOC
  
END MODULE YOMENSCOV
