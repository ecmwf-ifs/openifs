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
MODULE YOEVDF

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------
!*    ** *YOEVDF* CONTAINS CONSTANTS NEEDED BY *VDF....*
!     ------------------------------------------------------------------
TYPE TVDF
REAL(KIND=JPRB) :: RLAM
REAL(KIND=JPRB) :: RVDIFTS
LOGICAL         :: LWDS 
REAL(KIND=JPRB) :: REPS1WDS 
REAL(KIND=JPRB) :: REPS2WDS 
REAL(KIND=JPRB) :: RETAWDS 
REAL(KIND=JPRB) :: RTOFDALPHA
REAL(KIND=JPRB) :: REISTHSC
INTEGER(KIND=JPIM) :: NSUBST
!---------------------------------------------------------------------
CONTAINS
  PROCEDURE, PASS :: PRINT => PRINT_CONFIGURATION 
END TYPE TVDF
!*     *YOEVDF* CONTAINS CONSTANTS NEEDED BY *VDF....*
!     FOR THE COMPUTATION OF VERTICAL DIFFUSION

!     A.C.M. BELJAARS      E.C.M.W.F.    14/12/89

!     OBUKHOV-L UPDATE     ACMB          26/03/90.
!     LWDS-upate           A.Beljaars    Jan-2014   

!     NAME        TYPE     DESCRIPTION
!     ----        ----     -----------

!     *RLAM*      REAL     *ASYMPTOTIC MIXING LENGTH FOR MOMENTUM
!     *RVDIFTS*   REAL     *FACTOR FOR TIME STEP WEIGHTING IN *VDF....*
!     *LWDS*      LOGICAL  .T. for Wood/Diamantakis/Staniforth scheme      
!     *REPS1WDS*  REAL     Epsilon1 in WDS       
!     *REPS2WDS*  REAL     Epsilon2 in WDS         
!     *RETAWDS*   REAL     Eta in WDS         
!     *REISTHSC   REAL     Threshold for Inversion strength (K) for Stratocumulus
!     *NSUBST*    INTEGER  Number of substeps in VDF           
!     ------------------------------------------------------------------

 !---------------------------------------------------------------------

CONTAINS 
  
SUBROUTINE PRINT_CONFIGURATION(SELF, KDEPTH, KOUTNO)
IMPLICIT NONE
CLASS(TVDF), INTENT(IN) :: SELF
INTEGER    , INTENT(IN) :: KDEPTH
INTEGER    , INTENT(IN) :: KOUTNO

INTEGER :: IDEPTHLOC

IDEPTHLOC = KDEPTH+2

WRITE(KOUTNO,*) REPEAT(' ',KDEPTH)    // 'model%yrml_phy_g%yrvdf : '
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RLAM = ', SELF%RLAM
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RVDIFTS = ', SELF%RVDIFTS
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LWDS = ', SELF%LWDS
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'REPS1WDS = ', SELF%REPS1WDS
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'REPS2WDS = ', SELF%REPS2WDS
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RETAWDS = ', SELF%RETAWDS
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NSUBST = ', SELF%NSUBST
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'REISTHSC = ', SELF%REISTHSC


END SUBROUTINE PRINT_CONFIGURATION

END MODULE YOEVDF
