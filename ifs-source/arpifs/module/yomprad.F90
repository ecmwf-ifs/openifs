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

MODULE YOMPRAD

USE PARKIND1 , ONLY : JPIM, JPRB, JPRD

IMPLICIT NONE

SAVE

!----------------------------------------------------------------
!*    new radiation grid data structures
!----------------------------------------------------------------

! NRESOL_ID : transformation package resolution identifier
! NSMAX     : spectral truncation
! NASM0     : address in a spectral array of (m, n=m)
! MYMS      : actual wave numbers handled by this processor
! NDGLG     : number of latitude rows
! NDLON     : length of a row of latitude near equator
! NDGSAL    : Local version of NDGSA , always 1
! NDGENL    : Number of latitude rows for which this process has grid points
! NDGSAH    : 1-NROWIDE
! NDGENH    : NDGENL+NROWIDE
! NDGSAG    : Global version of NDGSA
! NDGENG    : Global version of NDGEN
! NDSUR1    : over dimensioning of NDLON for technical reasons (at least 2)
! NDLSUR    : NDLON+NDSUR1
! NDGSUR    : number of additional rows at each pole for semi-lagrangian
!           : interpolation
! MYFRSTACTLAT : first actual lat on this PE in grid-point space,
!              : it is nfrstlat(my_region_ns)
! MYLSTACTLAT  : last actual lat on this PE in grid-point space,
!              : it is nlstlat(my_region_ns)
! NGPTOT    : number of local grid points
! NGPTOTG   : total number of grid points
! NRGRI     : number of grid points on a latitude row
! NLOENG    : as above but extended at poles for SL interpolation
! GELAM     : radiation grid geographic longitude "lambda"
! GELAT     : radiation grid geographic latitude "theta"
! GESLO     : sine of geographic longitude "sin(lambda)"
! GECLO     : cosine of geographic longitude "cos(lambda)
! GEMU      : sine of geographic latitude "sin(theta)"
! RMU       : mu              sin(theta)
! RSQM2     : SQRT(R1MU2)     cos(theta)
! RLATIG    : arcsin(mu)      theta  GLOBAL VIEW
! RLATI     : arcsin(mu)      theta
! RIPI      : bi-cubic interpolation coefficients

TYPE RADIATION_GRID_STRUCT
INTEGER(KIND=JPIM) :: NRESOL_ID, NGPTOT, NGPTOTG, &
 & NGPTOTMX, NSPEC2, NSMAX, &
 & NPTRFLOFF, NUMP, NDLON, &
 & NDGSAL, NDGENL, NDGSAH, NDGENH, &
 & NDGLG, NDGSAG, NDGENG, NDLSUR, &
 & NFRSTLOFF, NDSUR1,NDGSUR, &
 & MYFRSTACTLAT, MYLSTACTLAT  
INTEGER(KIND=JPIM), ALLOCATABLE, DIMENSION(:) :: NRGRI, NLOENG, NPTRFRSTLAT, NFRSTLAT, &
 & NLSTLAT, MYMS, NASM0  
INTEGER(KIND=JPIM), ALLOCATABLE, DIMENSION(:,:):: NSTA, NONL
REAL(KIND=JPRB),    ALLOCATABLE, DIMENSION(:)  :: GELAM, GELAT, GECLO, GESLO, &
 & RMU, RSQM2, RLATIG, RLATI
REAL(KIND=JPRD),    ALLOCATABLE, DIMENSION(:)  :: GEMU
REAL(KIND=JPRD),    ALLOCATABLE, DIMENSION(:,:):: RIPI
!----------------------------------------------------------------------------
CONTAINS
  PROCEDURE, PASS :: PRINT => PRINT_CONFIGURATION 
END TYPE RADIATION_GRID_STRUCT
!============================================================================

CONTAINS

SUBROUTINE PRINT_CONFIGURATION(SELF, KDEPTH, KOUTNO)
IMPLICIT NONE
CLASS(RADIATION_GRID_STRUCT), INTENT(IN) :: SELF
INTEGER                     , INTENT(IN) :: KDEPTH
INTEGER                     , INTENT(IN) :: KOUTNO

INTEGER :: IDEPTHLOC

IDEPTHLOC = KDEPTH+2

WRITE(KOUTNO,*) REPEAT(' ',KDEPTH   ) // 'model%yrml_phy_rad%radgrid : '
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NRESOL_ID = ', SELF%NRESOL_ID
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NGPTOT = ', SELF%NGPTOT
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NGPTOTG = ', SELF%NGPTOTG
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NGPTOTMX = ', SELF%NGPTOTMX
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NSPEC2 = ', SELF%NSPEC2
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NSMAX = ', SELF%NSMAX
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NPTRFLOFF = ', SELF%NPTRFLOFF
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NUMP = ', SELF%NUMP
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NDLON = ', SELF%NDLON
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NDGSAL = ', SELF%NDGSAL
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NDGENL = ', SELF%NDGENL
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NDGSAH = ', SELF%NDGSAH
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NDGENH = ', SELF%NDGENH
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NDGLG = ', SELF%NDGLG
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NDGSAG = ', SELF%NDGSAG
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NDGENG = ', SELF%NDGENG
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NDLSUR = ', SELF%NDLSUR
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NFRSTLOFF = ', SELF%NFRSTLOFF
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NDSUR1 = ', SELF%NDSUR1
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NDGSUR = ', SELF%NDGSUR
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MYFRSTACTLAT = ', SELF%MYFRSTACTLAT
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MYLSTACTLAT = ', SELF%MYLSTACTLAT
IF (ALLOCATED(SELF%NRGRI)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NRGRI ALLOCATED OF SHAPE ', SHAPE(SELF%NRGRI), &
 &       ' SUM = ', SUM(SELF%NRGRI)
IF (ALLOCATED(SELF%NLOENG)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NLOENG ALLOCATED OF SHAPE ', SHAPE(SELF%NLOENG), &
 &       ' SUM = ', SUM(SELF%NLOENG)
IF (ALLOCATED(SELF%NPTRFRSTLAT)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NPTRFRSTLAT ALLOCATED OF SHAPE ', &
SHAPE(SELF%NPTRFRSTLAT), ' SUM = ', SUM(SELF%NPTRFRSTLAT)
IF (ALLOCATED(SELF%NFRSTLAT)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NFRSTLAT ALLOCATED OF SHAPE ', &
SHAPE(SELF%NFRSTLAT), ' SUM = ', SUM(SELF%NFRSTLAT)
IF (ALLOCATED(SELF%NLSTLAT)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NLSTLAT ALLOCATED OF SHAPE ', SHAPE(SELF%NLSTLAT), &
 &       ' SUM = ', SUM(SELF%NLSTLAT)
IF (ALLOCATED(SELF%MYMS)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MYMS ALLOCATED OF SHAPE ', SHAPE(SELF%MYMS), &
 &       ' SUM = ', SUM(SELF%MYMS)
IF (ALLOCATED(SELF%NASM0)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NASM0 ALLOCATED OF SHAPE ', SHAPE(SELF%NASM0), &
 &       ' SUM = ', SUM(SELF%NASM0)
IF (ALLOCATED(SELF%NSTA)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NSTA ALLOCATED OF SHAPE ', SHAPE(SELF%NSTA), &
 &       ' SUM = ', SUM(SELF%NSTA)
IF (ALLOCATED(SELF%NONL)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NONL ALLOCATED OF SHAPE ', SHAPE(SELF%NONL), &
 &       ' SUM = ', SUM(SELF%NONL)
IF (ALLOCATED(SELF%GELAM)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'GELAM ALLOCATED OF SHAPE ', SHAPE(SELF%GELAM), &
 &       ' SUM = ', SUM(SELF%GELAM)
IF (ALLOCATED(SELF%GELAT)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'GELAT ALLOCATED OF SHAPE ', SHAPE(SELF%GELAT), &
 &       ' SUM = ', SUM(SELF%GELAT)
IF (ALLOCATED(SELF%GECLO)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'GECLO ALLOCATED OF SHAPE ', SHAPE(SELF%GECLO), &
 &       ' SUM = ', SUM(SELF%GECLO)
IF (ALLOCATED(SELF%GESLO)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'GESLO ALLOCATED OF SHAPE ', SHAPE(SELF%GESLO), &
 &       ' SUM = ', SUM(SELF%GESLO)
IF (ALLOCATED(SELF%RMU)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RMU ALLOCATED OF SHAPE ', SHAPE(SELF%RMU), &
 &       ' SUM = ', SUM(SELF%RMU)
IF (ALLOCATED(SELF%RSQM2)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RSQM2 ALLOCATED OF SHAPE ', SHAPE(SELF%RSQM2), &
 &       ' SUM = ', SUM(SELF%RSQM2)
IF (ALLOCATED(SELF%RLATIG)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RLATIG ALLOCATED OF SHAPE ', SHAPE(SELF%RLATIG), &
 &       ' SUM = ', SUM(SELF%RLATIG)
IF (ALLOCATED(SELF%RLATI)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RLATI ALLOCATED OF SHAPE ', SHAPE(SELF%RLATI), &
 &       ' SUM = ', SUM(SELF%RLATI)
IF (ALLOCATED(SELF%GEMU)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'GEMU ALLOCATED OF SHAPE ', SHAPE(SELF%GEMU), &
 &       ' SUM = ', SUM(SELF%GEMU)
IF (ALLOCATED(SELF%RIPI)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RIPI ALLOCATED OF SHAPE ', SHAPE(SELF%RIPI), &
 &       ' SUM = ', SUM(SELF%RIPI)

END SUBROUTINE PRINT_CONFIGURATION

END MODULE YOMPRAD
