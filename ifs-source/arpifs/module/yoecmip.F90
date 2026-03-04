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

MODULE YOECMIP

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

!     ------------------------------------------------------------------
!*     *YOECMIP* PARAMETERS FOR CMIP FORCINGS (after YOECMIP5)
!     ------------------------------------------------------------------
!     Shuting Yang  EC-Earth 9/02/2010 (after YOEOZOC)     
!     Hans Hersbach  9/08/2010 , to be used for AMIP runs as well
!     R. Senan/C. Roberts (Feb 2017) Moved YOECMIP5 to generic CMIP module 
!                                    YOECMIP with support for CMIP6 forcings.
!     O.Marsden  Jan 2017  Define a type containing (time-variable) arrays

!     NAME        TYPE   PURPOSE
!     ----        ----   -------
!     NO3CMIP     INT    Controls time-varying CMIP ozone forcing. Specified in NAERAD.
!                          0=no time-varying ozone.
!                          5=Read CMIP5 time-varying ozone
!                          6=Read CMIP6 time-varying ozone
!     NGHGCMIP    INT    Controls time-varying CMIP ghg forcing. Specified in NAERAD.
!                          0=Use default IFS GHG data.
!                          5=Read CMIP5 GHG data.
!                          6=Read CMIP6 GHG data.
!     NRCP        INT    RCP scenario used to extend CMIP5 and CMIP6 historical GHG.
!                        Specified in NAERAD.
!                          1=RCP 3-PD
!                          2=RCP 4.5
!                          3=RCP 6.0
!                          4=RCP 8.5
!     NCMIPFIXYR  INT    Positive values used to fix CMIP forcings at specified year.
!                        Specified in NAERAD.
!     GHGDATADIR  CHAR   Directory containing CMIP ghg data file. Specified in NAERAD.     
!     CSOLARDATA  CHAR   File containing TSI data used with NHINCSOL=4. Specified in NAERAD. 
!     CO3DATADIR  CHAR   Directory containing CMIP ozone data file. Specified in NAERAD.    
!     CO3DATAFIL  CHAR   Name of CMIP ozone data file.
!     NLON1_CMIP5 INT    Number of longitudes in CMIP5 ozone data 
!     NLAT1_CMIP5 INT    Number of latitudes in CMIP5 ozone data 
!     NLV1_CMIP5  INT    Number of levels in CMIP5 ozone data      
!     NLON1_CMIP6 INT    Number of longitudes in CMIP6 ozone data  
!     NLAT1_CMIP6 INT    Number of latitudes in CMIP6 ozone data   
!     NLV1_CMIP6  INT    Number of levels in CMIP6 ozone data      
!     NMONTH1     INT    Number of months in each ozone data file
!     NCURRYR     INT    Year of the ozone data currently stored in ZOZCL

!     ----------------------------------------------------------------

INTEGER(KIND=JPIM),PARAMETER :: NLON1_CMIP5=72
INTEGER(KIND=JPIM),PARAMETER :: NLAT1_CMIP5=37
INTEGER(KIND=JPIM),PARAMETER :: NLV1_CMIP5=24
INTEGER(KIND=JPIM),PARAMETER :: NLON1_CMIP6=144
INTEGER(KIND=JPIM),PARAMETER :: NLAT1_CMIP6=96
INTEGER(KIND=JPIM),PARAMETER :: NLV1_CMIP6=66
INTEGER(KIND=JPIM),PARAMETER :: NMONTH1=14


TYPE :: TECMIP

!! non-time-dependent stuff could be in a separate type, is it worth it?
!!TYPE :: TECMIP_CONF
INTEGER(KIND=JPIM) :: NRCP
INTEGER(KIND=JPIM) :: NO3CMIP
INTEGER(KIND=JPIM) :: NGHGCMIP
INTEGER(KIND=JPIM) :: NCMIPFIXYR
INTEGER(KIND=JPIM) :: NCURRYR=-1

CHARACTER (LEN = 260) ::  GHGDATADIR
CHARACTER (LEN = 260) ::  CSOLARDATA 
CHARACTER (LEN = 260) ::  CO3DATADIR
CHARACTER (LEN =  80) ::  CO3DATAFIL
!!END TYPE TECMIP_CONF

!!TYPE :: TECMIP
REAL(KIND=JPRB), ALLOCATABLE :: ZOZCL(:,:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: RSINC1(:)
REAL(KIND=JPRB), ALLOCATABLE :: RPROC1(:)
REAL(KIND=JPRB), ALLOCATABLE :: RLATCLI(:)
REAL(KIND=JPRB), ALLOCATABLE :: RLONCLI(:)
REAL(KIND=JPRB), ALLOCATABLE :: ROZT1(:,:,:)
END TYPE TECMIP

END MODULE YOECMIP























