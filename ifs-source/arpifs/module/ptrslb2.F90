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

MODULE PTRSLB2

USE PARKIND1  ,ONLY : JPIM

IMPLICIT NONE

SAVE

!     NFLDSLB2     - Number of fields in semi-lagrangian buffer 2
!     MSLBUF2      - whole SLB2-pointer.
!     MSLB2DBBC1   - SLB2-pointer for intermediate quantity ((p/mRT)*Jacobian
!                    part of Laplacian) in Lth level for diagnostic BBC in NH.
!     MSLB2DPHI1   - SLB2-pointer for (p/m Rd T)_L for diagnostic BBC in NH.
!     MSLB2USI     - SLB2-pointer for t SI quantity in U-wind eqn.
!     MSLB2VSI     - SLB2-pointer for t SI quantity in V-wind eqn.
!     MSLB2TSI     - SLB2-pointer for t SI quantity in T eqn.
!     MSLB2PDSI    - SLB2-pointer for t SI quantity in Pcha eqn.
!     MSLB2VDSI    - SLB2-pointer for t SI quantity in dcha eqn.
!     MSLB2SPSI    - SLB2-pointer for t SI quantity in continuity eqn.
!     MSLB2VVEL    - SLB2-pointer for vertical velocity.
!     MSLB2URL     - SLB2-pointer for gp U-wind used in the SL traj research. 
!     MSLB2VRL     - SLB2-pointer for gp V-wind used in the SL traj research.
!     MSLB2WRL     - SLB2-pointer for gp "etadot"-wind used in the SL traj research.
!     MSLB2ZRL     - SLB2-pointer for gp W-wind used in LSLDP_XYZ cartesian option
!     MSLB2URL5    - SLB2-pointer for gp U-wind used in the SL traj res (TL).
!     MSLB2VRL5    - SLB2-pointer for gp V-wind used in the SL traj res (TL).
!     MSLB2WRL5    - SLB2-pointer for gp "etadot" used in the SL traj res (TL).
!     MSLB2ZRL5    - SLB2-pointer for gp W-wind used in LSLDP_XYZ cartesian option (TL)
!     MSLB2USI5    - SLB2-pointer for t SI quantity in U-wind eqn
!                    (same as MSLB2USI but for trajectory).
!     MSLB2VSI5    - SLB2-pointer for t SI quantity in V-wind eqn
!                    (same as MSLB2VSI but for trajectory).
!     MSLB2U15     - SLB2-pointer for grid-point quantity in U-wind eqn
!                    (same as MSLB2U1 but for trajectory).
!     MSLB2V15     - SLB2-pointer for grid-point quantity in V-wind eqn
!                    (same as MSLB2V1 but for trajectory).
!     MSLB2T15     - SLB2-pointer for grid-point quantity in T eqn
!                    (same as MSLB2T1 but for trajectory).
!     MSLB2Q15     - SLB2-pointer for grid-point quantity in Q eqn
!                    (same as MSLB2Q1 but for trajectory).
!     MSLB2KAPPA   - SLB2-pointer for the KAPPA function controlling SLHD 
!     MSLB2KAPPAT  - SLB2-pointer for the KAPPA function controlling SLHD on T 
!                    (valid at t0)
!     MSLB2KAPPAM  - SLB2-pointer for horizontal momentum echange coef.
!                     used in 3D turbulence (valid at t+dt in F)
!     MSLB2KAPPAH  - SLB2-pointer for horizontal heat (scalar) echange coef.
!                     used in 3D turbulence (valid at t+dt in F)
!     MSLB2KAPPA5  - trajectory of the KAPPA function (MSLB2KAPPA)
!     MSLB2KAPPAT5 - trajectory of the KAPPA function for T (MSLB2KAPPAT)
!     MSLB2GWF     - SLB2-pointer to store full level "gw" at t or t-dt.
!     MSLB2GDW     - SLB2-pointer to store full level "g dw" at t or t-dt.
!     MSLB2GWS     - SLB2-pointer to store "g w_surf" at t or t-dt.
!     MSLB2STDDISU - SLB2-pointer for the zonal STDDIS coef. for COMAD 
!     MSLB2STDDISV - SLB2-pointer for the meridional STDDIS coef. for COMAD 
!     MSLB2STDDISW - SLB2-pointer for the vertical STDDIS coef. for COMAD 
!     MSLB2STDDISU5- trajectory of the STDDISU function
!     MSLB2STDDISV5- trajectory of the STDDISV function
!     MSLB2STDDISW5- trajectory of the STDDISW function

TYPE :: TPTRSLB2
INTEGER(KIND=JPIM) :: NFLDSLB2
INTEGER(KIND=JPIM) :: MSLBUF2
INTEGER(KIND=JPIM) :: MSLB2DBBC1
INTEGER(KIND=JPIM) :: MSLB2DPHI1
INTEGER(KIND=JPIM) :: MSLB2USI
INTEGER(KIND=JPIM) :: MSLB2VSI
INTEGER(KIND=JPIM) :: MSLB2TSI
INTEGER(KIND=JPIM) :: MSLB2PDSI
INTEGER(KIND=JPIM) :: MSLB2VDSI
INTEGER(KIND=JPIM) :: MSLB2SPSI
INTEGER(KIND=JPIM) :: MSLB2VVEL
INTEGER(KIND=JPIM) :: MSLB2URL
INTEGER(KIND=JPIM) :: MSLB2VRL
INTEGER(KIND=JPIM) :: MSLB2WRL
INTEGER(KIND=JPIM) :: MSLB2ZRL
INTEGER(KIND=JPIM) :: MSLB2URL5
INTEGER(KIND=JPIM) :: MSLB2VRL5
INTEGER(KIND=JPIM) :: MSLB2WRL5
INTEGER(KIND=JPIM) :: MSLB2ZRL5
INTEGER(KIND=JPIM) :: MSLB2USI5
INTEGER(KIND=JPIM) :: MSLB2VSI5
INTEGER(KIND=JPIM) :: MSLB2U15
INTEGER(KIND=JPIM) :: MSLB2V15
INTEGER(KIND=JPIM) :: MSLB2T15
INTEGER(KIND=JPIM) :: MSLB2Q15
INTEGER(KIND=JPIM) :: MSLB2KAPPA
INTEGER(KIND=JPIM) :: MSLB2KAPPAT
INTEGER(KIND=JPIM) :: MSLB2KAPPAM
INTEGER(KIND=JPIM) :: MSLB2KAPPAH
INTEGER(KIND=JPIM) :: MSLB2KAPPA5
INTEGER(KIND=JPIM) :: MSLB2KAPPAT5
INTEGER(KIND=JPIM) :: MSLB2GWF
INTEGER(KIND=JPIM) :: MSLB2GDW
INTEGER(KIND=JPIM) :: MSLB2GWS
INTEGER(KIND=JPIM) :: MSLB2STDDISU
INTEGER(KIND=JPIM) :: MSLB2STDDISV
INTEGER(KIND=JPIM) :: MSLB2STDDISW
INTEGER(KIND=JPIM) :: MSLB2STDDISU5
INTEGER(KIND=JPIM) :: MSLB2STDDISV5
INTEGER(KIND=JPIM) :: MSLB2STDDISW5
!--------------------------------------------------------------
CONTAINS
  
  PROCEDURE, PASS :: PRINT => PRINT_CONFIGURATION 

END TYPE TPTRSLB2
!==============================================================

!!TYPE(TPTRSLB2), POINTER :: YRPTRSLB2 => NULL()
CONTAINS 
  
SUBROUTINE PRINT_CONFIGURATION(SELF, KDEPTH, KOUTNO)
  IMPLICIT NONE
  CLASS(TPTRSLB2), INTENT(IN) :: SELF
  INTEGER        , INTENT(IN) :: KDEPTH
  INTEGER        , INTENT(IN) :: KOUTNO

  INTEGER :: IDEPTHLOC

  IDEPTHLOC = KDEPTH + 2
  
  WRITE(KOUTNO,*) REPEAT(' ',KDEPTH   ) // 'model%yrml_gconf%ptrslb2 : '
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NFLDSLB2 = ', SELF%NFLDSLB2
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLBUF2 = ', SELF%MSLBUF2
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB2DBBC1 = ', SELF%MSLB2DBBC1
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB2DPHI1 = ', SELF%MSLB2DPHI1
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB2USI = ', SELF%MSLB2USI
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB2VSI = ', SELF%MSLB2VSI
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB2TSI = ', SELF%MSLB2TSI
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB2PDSI = ', SELF%MSLB2PDSI
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB2VDSI = ', SELF%MSLB2VDSI
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB2SPSI = ', SELF%MSLB2SPSI
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB2VVEL = ', SELF%MSLB2VVEL
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB2URL = ', SELF%MSLB2URL
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB2VRL = ', SELF%MSLB2VRL
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB2WRL = ', SELF%MSLB2WRL
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB2ZRL = ', SELF%MSLB2ZRL
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB2URL5 = ', SELF%MSLB2URL5
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB2VRL5 = ', SELF%MSLB2VRL5
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB2WRL5 = ', SELF%MSLB2WRL5
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB2ZRL5 = ', SELF%MSLB2ZRL5
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB2USI5 = ', SELF%MSLB2USI5
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB2VSI5 = ', SELF%MSLB2VSI5
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB2U15 = ', SELF%MSLB2U15
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB2V15 = ', SELF%MSLB2V15
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB2T15 = ', SELF%MSLB2T15
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB2Q15 = ', SELF%MSLB2Q15
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB2KAPPA = ', SELF%MSLB2KAPPA
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB2KAPPAT = ', SELF%MSLB2KAPPAT
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB2KAPPAM = ', SELF%MSLB2KAPPAM
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB2KAPPAH = ', SELF%MSLB2KAPPAH
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB2KAPPA5 = ', SELF%MSLB2KAPPA5
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB2KAPPAT5 = ', SELF%MSLB2KAPPAT5
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB2GWF = ', SELF%MSLB2GWF
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB2GDW = ', SELF%MSLB2GDW
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB2GWS = ', SELF%MSLB2GWS
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB2STDDISU = ', SELF%MSLB2STDDISU
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB2STDDISV = ', SELF%MSLB2STDDISV
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB2STDDISW = ', SELF%MSLB2STDDISW

END SUBROUTINE PRINT_CONFIGURATION

END MODULE PTRSLB2
