! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      MODULE YOWGRIBHD

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      IMPLICIT NONE

!*    ** *GRIB HEADERS:  COMMON (DEFAULT) GRIB HEADERS.

      INTEGER(KIND=JWIM), PARAMETER :: IMDLGRBID_G=106 !! see below the rule on how to select IMDLGRBID_G
      INTEGER(KIND=JWIM), PARAMETER :: IMDLGRBID_M=206
      INTEGER(KIND=JWIM) :: NDATE_TIME_WINDOW_END
      INTEGER(KIND=JWIM) :: NWINOFF
      INTEGER(KIND=JWIM) :: NGRIB_VERSION   ! grib version for output

      REAL(KIND=JWRB), PARAMETER :: PPMISS=-3.0_JWRB
      REAL(KIND=JWRB), PARAMETER :: PPEPS=1.0E-10_JWRB
      REAL(KIND=JWRB), PARAMETER :: PPREC=0.0_JWRB
      REAL(KIND=JWRB), PARAMETER :: PPRESOL=0.002_JWRB
      REAL(KIND=JWRB)            :: PPMIN_RESET=0.0_JWRB

      INTEGER(KIND=JWIM) :: NTENCODE

      INTEGER(KIND=JWIM), PARAMETER :: NGRBRESI=16
      INTEGER(KIND=JWIM), PARAMETER :: NGRBRESS=9

      CHARACTER(LEN=1), PARAMETER :: HOPERI='M' 
      CHARACTER(LEN=1), PARAMETER :: HOPERS='C' 

      LOGICAL :: LGRHDIFS
      LOGICAL :: LNEWLVTP
      LOGICAL :: LPADPOLES=.TRUE.
      LOGICAL :: LL_GRID_SIMPLE_MATRIX=.TRUE.

!*    VARIABLE.   TYPE.     PURPOSE.
!     ---------   -------   --------
!     IMDLGRBID_G INTEGER   GLOBAL MODEL IDENTIFICATION FOR GRIB CODING
!                           IT CAN ALSO BE MODIFIED IN THE INPUT NAMELIST.
!     IMDLGRBID_M INTEGER   LAW MODEL IDENTIFICATION FOR GRIB CODING
!                           IT CAN ALSO BE MODIFIED IN THE INPUT NAMELIST.
! 
! The generating process identification numbers (model numbers) for GRIB
! headers are allocated in pre-defined ranges for ECMWF GRIB coded
! fields. The field in the GRIB code for this number is 1 octet and with
! the value 255 indicating 'missing value', only numbers 1-254 are
! available for use.

! The atmospheric model allocated numbers are in the range 121 to 203.
! The Global wave model allocated numbers are in the range 104 to 120.
! The LAW model allocated numbers are in the range 204 to 220.
 
! The numbers 221 to 254 stay reserved for the moment.

! This pre-allocation was introduced to enable some Member States, which
! use the model number to identify products, to write their software in
! such a way that model number changes did not cause them problems with
! hard-coded model numbers. 

!     PPMISS      REAL      ALL SPECTRAL VALUES LESS OR EQUAL PPMISS ARE
!                           REPLACED BY THE MISSING DATA INDICATOR
!     PPEPS       REAL      SMALL NUMBER USED IN SPECTRAL PACKING OF 251
!     PPREC       REAL      REFERENCE VALUE FOR SPECTRAL PACKING OF 251
!     PPRESOL     REAL      MAXIMUN RESOLUTION POSSIBLE WHEN ENCODING 
!                           SPECTRA (PARAMETER 251).
!     PPMIN_RESET REAL      CAN BE USED TO SET THE MINIMUM OF PPMIN
!                           IN WGRIBOUT TO A LOWER VALUE.
!     NTENCODE    INTEGER   TOTAL NUMBER OF GRID POINTS FOR ENCODING
!     NGRBRESI    INTEGER   NUMBER OF BITS USED TO ENCODE INTEGRATED
!                           PARAMETERS
!     NGRBRESS    INTEGER   NUMBER OF BITS USED TO ENCODE SPECTRA
!     HOPERI      CHARACTER GRIB ENCODING ACTION FOR INTEGRATED FIELDS.
!     HOPERS      CHARACTER GRIB ENCODING ACTION FOR SPECTRA. 
!     LGRHDIFS    LOGICAL   IF TRUE THEN GRIB HEADER WILL USE INFORMATION 
!                           AS PROVIDED BY THE IFS.
!     LNEWLVTP    LOGICAL   IF TRUE THE NEW LEVTYPE DEFINITION WILL BE USED.
!     LPADPOLES   LOGICAL   TRUE IF POLES ARE PADDED WHEN SAVIND TO GRIB.
!     LL_GRID_SIMPLE_MATRIX IF TRUE THEN THE 2D SPECTRA WILL USE THE LEGACY grid_simple_matrix
!                           TO ENCODE THE 2D SPECTRA in GRIB1. THIS SHOULD BE PHASED OUT as soon as feasible!

 
! ----------------------------------------------------------------------
      END MODULE YOWGRIBHD
