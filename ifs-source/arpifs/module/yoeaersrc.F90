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

MODULE YOEAERSRC

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------
!*    ** *YOEAERSRC* - CONTROL OPTIONS FOR AEROSOLS' SOURCE
!     ------------------------------------------------------------------

TYPE :: TEAERSRC
INTEGER(KIND=JPIM) :: NDDUST, NSSALT, NAERWND, NDMSO, NPIST,NSO4SCHEME

LOGICAL :: LEPAERO, LAEROMIN, LOCNDMS

REAL(KIND=JPRB) :: RGELAV, RGEMUV, RDGLAV, RDGMUV
REAL(KIND=JPRB) :: RCLONV, RSLONV, RDCLONV, RDSLONV
REAL(KIND=JPRB) :: RLATVOL,RLONVOL
REAL(KIND=JPRB) :: RFCTDU, RFCTSS, RFCTDUR, RFCTSSR
REAL(KIND=JPRB) :: RAERDUB, RDMSMIN
REAL(KIND=JPRB) :: RCODECA, RSIDECA, RCOVSRA, RSIVSRA

REAL(KIND=JPRB) :: RSSFLX(9)
REAL(KIND=JPRB) :: RDDUSRC(9)
!---------------------------------------------------------------------
CONTAINS
  PROCEDURE, PASS :: PRINT => PRINT_CONFIGURATION 
END TYPE TEAERSRC
!=====================================================================

TYPE(TEAERSRC), POINTER :: YREAERSRC => NULL()

!     ------------------------------------------------------------------
! RSSFLX     : sea salt flux for (9)3-size bins (in mg m-2 s-1) for a 1 m s-1
!              wind speed, at 10 m height and 80% RH 
! Rlat/lonVOL: LAT/LON of a possible volcanic eruption
! RxxPHIL, PHOB: fractions of BC and OM hydrophilic or hydrophobic
! NDDUST, NSSALT : various formulations for dust and sea-salt emissions (test only)
! NDMSO      : 0= no oceanic DMS, 1= parametrized, 2=from source file
! NPIST      : Index for treatment of piston velocity
! RAERDUB    : base value for dust emission
! NAERWND    : 0 = 10m-wind for both DU and SS
!            : 1 = gust for SS, 10m-wind for DU
!            : 2 = gust for DU, 10m-wind for SS
!            : 3 = gust for both DU and SS
! RFCTSS, -DU: factors depending on choice of NAERWND
!     ------------------------------------------------------------------

CONTAINS

SUBROUTINE PRINT_CONFIGURATION(SELF, KDEPTH, KOUTNO)
IMPLICIT NONE
CLASS(TEAERSRC), INTENT(IN) :: SELF
INTEGER        , INTENT(IN) :: KDEPTH
INTEGER        , INTENT(IN) :: KOUTNO

INTEGER :: IDEPTHLOC

IDEPTHLOC = KDEPTH+2

WRITE(KOUTNO,*) REPEAT(' ',KDEPTH   ) // 'model%yrml_phy_aer%yreaersrc : '
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NDDUST = ', SELF%NDDUST
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NSSALT = ', SELF%NSSALT
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NAERWND = ', SELF%NAERWND
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NDMSO = ', SELF%NDMSO
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NPIST = ', SELF%NPIST
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NSO4SCHEME = ', SELF%NSO4SCHEME
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LEPAERO = ', SELF%LEPAERO
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LAEROMIN = ', SELF%LAEROMIN
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LOCNDMS = ', SELF%LOCNDMS
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RGELAV = ', SELF%RGELAV
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RGEMUV = ', SELF%RGEMUV
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RDGLAV = ', SELF%RDGLAV
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RDGMUV = ', SELF%RDGMUV
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RCLONV = ', SELF%RCLONV
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RSLONV = ', SELF%RSLONV
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RDCLONV = ', SELF%RDCLONV
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RDSLONV = ', SELF%RDSLONV
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RLATVOL = ', SELF%RLATVOL
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RLONVOL = ', SELF%RLONVOL
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RFCTDU = ', SELF%RFCTDU
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RFCTSS = ', SELF%RFCTSS
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RFCTDUR = ', SELF%RFCTDUR
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RFCTSSR = ', SELF%RFCTSSR
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RAERDUB = ', SELF%RAERDUB
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RDMSMIN = ', SELF%RDMSMIN
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RCODECA = ', SELF%RCODECA
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RSIDECA = ', SELF%RSIDECA
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RCOVSRA = ', SELF%RCOVSRA
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RSIVSRA = ', SELF%RSIVSRA
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RSSFLX SUM = ', SUM(SELF%RSSFLX)
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RDDUSRC SUM = ', SUM(SELF%RDDUSRC)


END SUBROUTINE PRINT_CONFIGURATION

END MODULE YOEAERSRC

