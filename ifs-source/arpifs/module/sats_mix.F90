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

MODULE SATS_MIX

! Purpose :
! -------
!    To initialise variables from the namelist NAMSATS.

! Interface :
! ---------
!    Empty.

! External :
! --------
!    POSNAM

! Method :
! ------
!    See Documentation.

! Reference :
! ---------
!    Object Oriented manners
 

! Author :
! ------
!    Ryad El Khatib *METEO-FRANCE* inspired by Yannick Tremolet *ECMWF*
!    Original : 23-Jul-2009

! Modifications :
! -------------
!  23/11/11 C. Lupu:  Add Met-9 radiance calculations to be performed in hretr
!  K. Yessad (July 2014): Move some variables.
!  19-Nove-2014  C. Lupu        Updates to sensors and platform lists from rttov_consts
!  M.Hamrud (Feb 2015) : LRTCALC_HRETR to true for MWHS (need Jacobian peak calc.)
!  J. Letertre-Danczak (Oct 2015) : Activate AHI
!  R. El Khatib 12-Aug-2016 NJPPF in namelist
!  C.Lupu (Feb 2017) : RTTOV 12.1 updates
!  C. Burrows (Apr 2018) : Activate ABI
!  F. Duruisseau 24-Sept-2018: Add BAYRAD namelist
!  C. Burrows (Jul 2019) : Add GIIRS
!  R. Eresmaa (Apr 2020) : Activate screening for HIRAS and IKFS2
!-----------------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN   , ONLY : NULOUT, NULNAM
USE YOMMP0   , ONLY : NPRINTLEV
USE RTTOV_CONST  , ONLY : NINST, INST_ID_HIRS, INST_ID_MSU, INST_ID_VTPR1, INST_ID_VTPR2, INST_ID_AIRS, &
 & INST_ID_IASI, INST_ID_MVIRI, INST_ID_SEVIRI, INST_ID_GOESIM, INST_ID_MTSATIM, INST_ID_CRIS, &
 & INST_ID_MWHS, INST_ID_IRAS, INST_ID_ABI, INST_ID_AHI, INST_ID_GIIRS, INST_ID_HIRAS, INST_ID_IKFS2
USE YOMSATS  , ONLY : NJPPF

IMPLICIT NONE
SAVE

PUBLIC
PRIVATE MISSING_INDICATOR


TYPE SATID_TABLE_T
  INTEGER(KIND=JPIM)  :: IBUFRID       ! WMO satellite ID
  INTEGER(KIND=JPIM)  :: IRTSERIES     ! Platform/series used by RTTOV
  INTEGER(KIND=JPIM)  :: IRTID         ! Satellite ID used by RTTOV
END TYPE SATID_TABLE_T

TYPE SUBT_TABLE_T
  INTEGER(KIND=JPIM)  :: ICODETYPE     
  INTEGER(KIND=JPIM)  :: IRTSERIES     ! Platform/series used by RTTOV
  INTEGER(KIND=JPIM)  :: IRTSUBTYPE    ! Subtype, for potential use
END TYPE SUBT_TABLE_T

INTEGER(KIND=JPIM), PARAMETER :: JPMXSATIDS=120
INTEGER(KIND=JPIM), PARAMETER :: JPMXSUBTS=50

!     MEMIS_SW     : switch for new ice /land classif method for emis
INTEGER(KIND=JPIM) :: MEMIS_SW=1

! Private missing indicator
INTEGER(KIND=JPIM) :: MISSING_INDICATOR=2147483647

LOGICAL :: LGHGRTTOV

! Define for which sensors the following tasks are performed in hretr:
LOGICAL :: LRTCALC_HRETR     (0:NINST-1) = .FALSE. ! radiance calculations
LOGICAL :: LCLD_RTCALC_SCREEN(0:NINST-1) = .FALSE. ! cloudy radiance calculations
LOGICAL :: LCLD_RTCALC_ASSIM (0:NINST-1) = .FALSE. ! cloudy assimilation

! Activate Lambertian treatment by sensor (clear-sky MW instruments only)
LOGICAL :: L_LAMBERTIAN(0:NINST-1) = .FALSE.

! Define table/namelists for BAYRAD
!

TYPE INV_ERROR_MOD_PTS
  INTEGER(KIND=JPIM)        :: CHAN = -1
  REAL(KIND=JPRB)           :: E_CLR = 1._JPRB
  REAL(KIND=JPRB)           :: E_SAT = 10._JPRB
  REAL(KIND=JPRB)           :: T_CLR = 230._JPRB
  REAL(KIND=JPRB)           :: T_SAT = 220._JPRB
  INTEGER(KIND=JPIM)        :: N = -1
END TYPE INV_ERROR_MOD_PTS


TYPE RH_ERROR_MOD_PTS
  REAL(KIND=JPRB),DIMENSION(2)    :: L_HIGH=(/-10._JPRB, 1._JPRB/)
  REAL(KIND=JPRB),DIMENSION(2)    :: L_LOW= (/ -2._JPRB, 1._JPRB/)
  REAL(KIND=JPRB),DIMENSION(2)    :: R_LOW= (/  2._JPRB, 1._JPRB/)
  REAL(KIND=JPRB),DIMENSION(2)    :: R_HIGH=(/ 10._JPRB, 1._JPRB/)
END TYPE RH_ERROR_MOD_PTS


TYPE OBS_ERR_TYPE
  LOGICAL                         :: LMODEL =  .FALSE.
  REAL(KIND=JPRB),DIMENSION(50)   :: LIST   =  -1._JPRB  ! Observation Error
END TYPE OBS_ERR_TYPE


TYPE SAT_BAY_RAD_TABLE_T
  LOGICAL                          :: LBAY =  .FALSE.           ! Active bay rad inv.
  INTEGER(KIND=JPIM)               :: KPROF=  1                 ! Nb profiles to perform
  INTEGER(KIND=JPIM)               :: NPLEV=  0                 ! Nb pressure levels for T Q
  INTEGER(KIND=JPIM),DIMENSION(90) :: BPLEV= -1                 ! Pressure levels
  INTEGER(KIND=JPIM)               :: NCHAN=  0                 ! Nb channels for T Q
  INTEGER(KIND=JPIM),DIMENSION(50) :: BCHAN= -1                 ! Channels selection for bay rad inv
  TYPE(OBS_ERR_TYPE)               :: VERR(999)
  REAL(KIND=JPRB)                  :: RETAC=  1._JPRB           ! Max ratio ABS(OBS-PSEUDO_OBS)/ABS(OBS-GUESS_SCATT)
  LOGICAL                          :: LINVMOD  = .FALSE.        ! Inversion Error Model
  LOGICAL                          :: LSCATTQC = .FALSE.        ! Active QC based on scattering simulations
  REAL(KIND=JPRB)                  :: RMIND = 0._JPRB
  REAL(KIND=JPRB)                  :: RFIND = 0._JPRB
  TYPE(INV_ERROR_MOD_PTS)          :: INV_ERR_M(50)
  TYPE(RH_ERROR_MOD_PTS)           :: RH_ERR_M(50)
  INTEGER(KIND=JPIM)               :: NSCATT_PARTICLES = 0      ! 0: Default way, >=1 : enable multi-particle simulation
  CHARACTER(LEN=256)               :: ISCATT_NAMES(50) = 'NONE'
END TYPE SAT_BAY_RAD_TABLE_T

TYPE(SAT_BAY_RAD_TABLE_T) :: BAYRAD_TABLE(0:JPMXSATIDS)

!--------------------------------


!     LCO2_DIAG_AIRS: switch for CO2-slicing on AIRS data
LOGICAL :: LCO2_DIAG_AIRS = .FALSE.
!     LCO2_DIAG_IASI: switch for CO2-slicing on IASI data
LOGICAL :: LCO2_DIAG_IASI = .FALSE.

! Expected number of fixed pressure levels for RTTOV if
! IFS-interpolation is used
INTEGER(KIND=JPIM) :: NLSAT = 44

! Switch whether a fatal error in RTTOV should trigger an abort
LOGICAL   :: LABORT_RTTOV_FAILURE = .FALSE.

LOGICAL   :: LNEW_HIRSCLD = .FALSE.

TYPE(SATID_TABLE_T) :: YLSATID_TABLE(JPMXSATIDS)

TYPE(SUBT_TABLE_T) :: YLSUBT_TABLE(JPMXSUBTS)

! Cloud sink variable: closest proximity of cloud top to surface pressure [hPa]
REAL(KIND=JPRB) :: CTOP_SFC_OFFSET = 20.0_JPRB

! Switch wether to read all the coefficient file or just a subset of channels
! from a file thtat contains all the channels (hyperspectral sounders only)
LOGICAL :: LPARTIAL_COEF_FILES = .FALSE.

!-----------------------------------------------------------------------------

CONTAINS

SUBROUTINE SUSATS(KMISS)

! KMISS : User-defined missing indicator for YLSATID_TABLE & YLSUBT_TABLE

INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KMISS

INTEGER(KIND=JPIM) :: IMISS, J
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "namsats.nam.h"

#include "posnam.intfb.h"

IF (LHOOK) CALL DR_HOOK('SATS_MIX:SUSATS',0,ZHOOK_HANDLE)

IF (PRESENT(KMISS)) THEN
  IMISS=KMISS
ELSE
  IMISS=MISSING_INDICATOR
ENDIF
WRITE(NULOUT, '('' SUSATS MISSING INDICATOR USED : '',I102)') IMISS

YLSATID_TABLE(:)%IBUFRID   = IMISS
YLSATID_TABLE(:)%IRTSERIES = IMISS
YLSATID_TABLE(:)%IRTID     = IMISS

YLSUBT_TABLE(:)%IRTSERIES  = IMISS
YLSUBT_TABLE(:)%ICODETYPE  = IMISS
YLSUBT_TABLE(:)%IRTSUBTYPE = IMISS

! Define sensors for which the radiances are computed in hretr:
LRTCALC_HRETR (INST_ID_HIRS) = .TRUE.
LRTCALC_HRETR (INST_ID_MSU) = .TRUE.
LRTCALC_HRETR (INST_ID_VTPR1) = .TRUE.
LRTCALC_HRETR (INST_ID_VTPR2) = .TRUE.
LRTCALC_HRETR (INST_ID_AIRS) = .TRUE.
LRTCALC_HRETR (INST_ID_IASI) = .TRUE.
LRTCALC_HRETR (INST_ID_MVIRI) = .TRUE.
LRTCALC_HRETR (INST_ID_SEVIRI) = .TRUE.
LRTCALC_HRETR (INST_ID_GOESIM) = .TRUE.
LRTCALC_HRETR (INST_ID_MTSATIM) = .TRUE.
LRTCALC_HRETR (INST_ID_CRIS) = .TRUE.
LRTCALC_HRETR (INST_ID_MWHS) = .TRUE.
LRTCALC_HRETR (INST_ID_IRAS) = .TRUE.
LRTCALC_HRETR (INST_ID_ABI) = .TRUE.
LRTCALC_HRETR (INST_ID_AHI) = .TRUE.
LRTCALC_HRETR (INST_ID_GIIRS) = .TRUE.
LRTCALC_HRETR (INST_ID_IKFS2) = .TRUE.
LRTCALC_HRETR (INST_ID_HIRAS) = .TRUE.

#ifdef NECSX
  NJPPF = 256
#else
! Smaller value uses less memory, larger value limits allocation/deallocation cycles
  NJPPF = 8
#endif
CALL POSNAM(NULNAM,'NAMSATS')
READ(NULNAM,NAMSATS)

WRITE(NULOUT, '('' MEMIS_SW = '',I2, &
 & '' LABORT_RTTOV_FAILURE = '',L2,'' LNEW_HIRSCLD = '',L2 )') &
 & MEMIS_SW, LABORT_RTTOV_FAILURE,LNEW_HIRSCLD
WRITE(NULOUT, '('' LRTCALC_HRETR = '',32L2)') LRTCALC_HRETR
WRITE(NULOUT, '('' LCLD_RTCALC_SCREEN = '',32L2)') LCLD_RTCALC_SCREEN
WRITE(NULOUT, '('' LCLD_RTCALC_ASSIM = '',32L2)') LCLD_RTCALC_ASSIM
WRITE(NULOUT, '('' LCO2_DIAG_AIRS ='',L2)') LCO2_DIAG_AIRS
WRITE(NULOUT, '('' LCO2_DIAG_IASI ='',L2)') LCO2_DIAG_IASI
IF (NPRINTLEV > 1) THEN
  WRITE(NULOUT, '('' YLSATID_TABLE(:)%IBUFRID  YLSATID_TABLE(:)%IRTSERIES &
   & YLSATID_TABLE(:)%IRTID '')')
  DO J=1, JPMXSATIDS
    WRITE(NULOUT, '(3I10)') YLSATID_TABLE(J)%IBUFRID,&
     & YLSATID_TABLE(J)%IRTSERIES,YLSATID_TABLE(J)%IRTID
  ENDDO
  WRITE(NULOUT, '('' YLSUBT_TABLE(:)%IRTSERIES  YLSUBT_TABLE(:)%ICODETYPE &
   & YLSUBT_TABLE(:)%IRTSUBTYPE '')')
  DO J=1, JPMXSUBTS
    WRITE(NULOUT, '(3I10)') YLSUBT_TABLE(J)%IRTSERIES,&
     & YLSUBT_TABLE(J)%ICODETYPE,YLSUBT_TABLE(J)%IRTSUBTYPE
  ENDDO
ENDIF
WRITE(NULOUT, '('' CTOP_SFC_OFFSET [hPa] = '',F7.2)') CTOP_SFC_OFFSET
WRITE(NULOUT, '('' LPARTIAL_COEF_FILES = '',L2)') LPARTIAL_COEF_FILES
WRITE(NULOUT, '('' NJPPF = '',I4)') NJPPF

IF (LHOOK) CALL DR_HOOK('SATS_MIX:SUSATS',1,ZHOOK_HANDLE)

END SUBROUTINE SUSATS

END MODULE SATS_MIX
