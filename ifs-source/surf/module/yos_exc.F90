MODULE YOS_EXC
 
USE PARKIND1  ,ONLY : JPIM     ,JPRB

! (C) Copyright 2005- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------
!*    ** *YOS_EXC* CONTAINS CONSTANTS NEEDED BY *V....*
!     ------------------------------------------------------------------

TYPE :: TEXC
LOGICAL :: LELWDD           ! TRUE when longwave downward derivative is used for skin temperature
LOGICAL :: LELWTL           ! TRUE when longwave net radiation is tiled to the surface
LOGICAL :: LEOCWA           ! TRUE if WARM OCEAN LAYER PARAMETRIZATION active
LOGICAL :: LEOCCO           ! TRUE if COOL OCEAN SKIN PARAMETRIZATION active
LOGICAL :: LWCOU            ! TRUE if coupled to wave model
LOGICAL :: LWCOU2W          ! TRUE if coupled to wave model, there are feedbacks onto the atmosphere (momentum, warm layer scheme or ocean TKE scheme)
LOGICAL :: LWCOUHMF         ! TRUE if coupled to wave model, there is a direct feedback onto the ocean heat and moisture flux
LOGICAL :: LEOCSA           ! TRUE if SALINTY EFFECT ON SATURATION AT OCEAN SURFACE active
LOGICAL :: LEOCLA           ! TRUE if LANGMUIR CIRCULATION EFFIECT active
LOGICAL :: LSCMEC           ! TRUE if ECMWF Single Column Model
LOGICAL :: LROUGH           ! TRUE if externally specified roughness lengths
REAL(KIND=JPRB) :: REXTZ0M  ! roughness length for momentum (if LROUGH)
REAL(KIND=JPRB) :: REXTZ0H  ! roughness length for heat and moisture (if LROUGH)
REAL(KIND=JPRB) :: RKAP     ! VONKARMAN CONSTANT
REAL(KIND=JPRB) :: REPDU2   ! MINIMUM VELOCITY DIFFERENCE IN RI-NUMBER
REAL(KIND=JPRB) :: RPARZI   ! ANSATZ FOR PBL-H IN W* COMPUTATION
REAL(KIND=JPRB) :: RZ0ICE   ! ROUGHNESS OVER SEA ICE
REAL(KIND=JPRB) :: REPUST   ! MINIMUM FRICTION VELOCITY (SECURITY PARAMETER)
REAL(KIND=JPRB) :: RNU      ! SMOOTH SURFACE CONSTANT KINEMATIC AIR DENSITY
REAL(KIND=JPRB) :: RNUM     ! SMOOTH SURFACE CONSTANT IN Z0M=RNUM/u* (RNUM a faction of RNU)
REAL(KIND=JPRB) :: RNUH     ! SMOOTH SURFACE CONSTANT IN Z0H=RNUH/u* (RNUH a faction of RNU)
REAL(KIND=JPRB) :: RNUQ     ! SMOOTH SURFACE CONSTANT IN Z0Q=RNUQ/u* (RNUQ a faction of RNU)
REAL(KIND=JPRB) :: RCHAR    ! CHARNOCK CONSTANT
END TYPE TEXC

                            !  (simplified physics only)

!*     *YOS_EXC* CONTAINS CONSTANTS NEEDED BY *V....*
!     FOR THE COMPUTATION OF SURFACE DIFFUSION EXCHANGE COEFFICIENTS

!     A.C.M. BELJAARS      E.C.M.W.F.    14/12/89
!     OBUKHOV-L UPDATE     ACMB          26/03/90.
!     surf externalisation P. Viterbo    09/06/2005

END MODULE YOS_EXC
