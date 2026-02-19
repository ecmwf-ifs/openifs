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

MODULE YOMRCOEF

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!   -----------------------------------------------------------------

!*   Logical switches for writing and reading the radiation coefficients
!    in file or in core

!    LRCOEF                       : switch for write out and read
!                                   the radiation coefficients in file

!*   Logical switch connected with the simplified radiation

!    LTLADDIA    : switch for using diabatic adjoint and then reading
!                  rad.coef. at each time step when LRAYSP
!    LGLOBRAD    : switch to compute the global mean for the thermal
!                  radiation matrix

!*   Number of arrays for radiation coefficients

!    NG3SR                        : number of 3D fields (NFLEVG)
!    NGMTR                        : number of matrices (0:NFLEVG,0:NFLEVG)
!    NLATWR                       : number of latitude to be written
!    NLATRD                       : number of latitude to be read

!*   Buffer for radiation coefficients

!    NLENGSRB                       : length of buffer for 3D solar radiation
!                                     coefficients and correction for thermal
!                                     radiation
!    NLENGTRB                       : length of buffer for matrix of thermal
!                                     radiation coefficients
!    SOLRAD(NLENGSRB)               : buffer for 3D fields
!    THERRAD(NSLBR/NGPBLKS/1,NLENGTRB): buffer for matrix
!                                     for 3D fields

!*   Buffer for simplified thermal radiation

!    TRWEIGHT(NFLEVG+1,NFLEVG+1) : buffer for sum of mean weights
!                                  for the whole domain
!    TRMATSUM(NFLEVG+1,NFLEVG+1) : buffer for sum of thermal radiation
!                                  matrices for the whole domain

TYPE :: TRCOEF
LOGICAL :: LRCOEF
LOGICAL :: LTLADDIA
LOGICAL :: LGLOBRAD

INTEGER(KIND=JPIM) :: NG3SR
INTEGER(KIND=JPIM) :: NGMTR
INTEGER(KIND=JPIM) :: NLATWR
INTEGER(KIND=JPIM) :: NLATRD

INTEGER(KIND=JPIM) :: NLENGSRB
INTEGER(KIND=JPIM) :: NLENGTRB

REAL(KIND=JPRB),ALLOCATABLE:: SOLRAD(:)
REAL(KIND=JPRB),ALLOCATABLE:: THERRAD(:,:)

REAL(KIND=JPRB),ALLOCATABLE:: TRWEIGHT(:,:)
REAL(KIND=JPRB),ALLOCATABLE:: TRMATSUM(:,:)
!----------------------------------------------------------------------------
CONTAINS
  PROCEDURE, PASS :: PRINT => PRINT_CONFIGURATION 
END TYPE TRCOEF
!============================================================================

!!TYPE(TRCOEF), POINTER :: YRRCOEF => NULL()


!   ----------------------------------------------------------------

CONTAINS

SUBROUTINE PRINT_CONFIGURATION(SELF, KDEPTH, KOUTNO)
IMPLICIT NONE
CLASS(TRCOEF), INTENT(IN) :: SELF
INTEGER      , INTENT(IN) :: KDEPTH
INTEGER      , INTENT(IN) :: KOUTNO

INTEGER :: IDEPTHLOC

IDEPTHLOC = KDEPTH+2

WRITE(KOUTNO,*) REPEAT(' ',KDEPTH   ) // 'model%yrml_phy_rad%yrrcoef : '
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LRCOEF = ', SELF%LRCOEF
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LTLADDIA = ', SELF%LTLADDIA
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LGLOBRAD = ', SELF%LGLOBRAD
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NG3SR = ', SELF%NG3SR
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NGMTR = ', SELF%NGMTR
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NLATWR = ', SELF%NLATWR
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NLATRD = ', SELF%NLATRD
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NLENGSRB = ', SELF%NLENGSRB
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NLENGTRB = ', SELF%NLENGTRB
IF (ALLOCATED(SELF%SOLRAD)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SOLRAD ALLOCATED OF SHAPE ',     &
 &        SHAPE(SELF%SOLRAD), ' SUM = ', SUM(SELF%SOLRAD)
IF (ALLOCATED(SELF%THERRAD)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'THERRAD ALLOCATED OF SHAPE ',   &
 &        SHAPE(SELF%THERRAD), ' SUM = ', SUM(SELF%THERRAD)
IF (ALLOCATED(SELF%TRWEIGHT)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'TRWEIGHT ALLOCATED OF SHAPE ', &
 &        SHAPE(SELF%TRWEIGHT), ' SUM = ', SUM(SELF%TRWEIGHT)
IF (ALLOCATED(SELF%TRMATSUM)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'TRMATSUM ALLOCATED OF SHAPE ', &
 &        SHAPE(SELF%TRMATSUM), ' SUM = ', SUM(SELF%TRMATSUM)

END SUBROUTINE PRINT_CONFIGURATION

END MODULE YOMRCOEF
