! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      MODULE YOWWAMI

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      IMPLICIT NONE

!     **  *WAMINFO*                                               

      CHARACTER(LEN=14)  :: CBEGDT 
      CHARACTER(LEN=14)  :: CENDDT 
      CHARACTER(LEN=14)  :: CBPLTDT 
      CHARACTER(LEN=14)  :: CEPLTDT 

      CHARACTER(LEN=14)  :: CLSPDT 
      CHARACTER(LEN=14)  :: CRSTDT 
      INTEGER(KIND=JWIM) :: IANALPD 
      INTEGER(KIND=JWIM) :: IFOREPD 
      INTEGER(KIND=JWIM) :: IDELWIN 
      INTEGER(KIND=JWIM) :: IASSIM 
      INTEGER(KIND=JWIM) :: NFCST 
      INTEGER(KIND=JWIM) :: ISTAT(3) 

!     VARIABLE   TYPE     PURPOSE                                  
!     --------   ----     -------                                  
!     IBEGDT    INTEGER   BEGIN DATE OF RUN (YYMMDDHHMM)           
!     IENDDT    INTEGER   END DATE OF RUN (YYMMDDHHMM)             
!     IANALPD   INTEGER   ANALYSIS PERIOD IN SECONDS
!     IFOREPD   INTEGER   FORECAST PERIOD IN SECONDS
!     IDELWIN   INTEGER   WIND TIME STEP IN SECONDS
!     CDATER    CHAR*14   DATE FOR OUTPUT OF BOTH RESTART FILES
!                         IF NOT SET IT WILL BY DEFAULT TO CENDDT 
!     CDATES    CHAR*14   LAST DATE FOR OUTPUT OF RESTART FILES
!                         IF NOT SET IT WILL BY DEFAULT TO CENDDT 
!     CBPLTDT   CHAR*14   BEGIN DATE OF ANALYSIS (YYYYMMDDHHMM)        
!     CEPLTDT   CHAR*14   END DATE OF ANALYSIS                     
!     IASSIM    INTEGER   0 = NO DATA ASSIMILATION                 
!                         1 = DATA ASSIMILATION                    
!                             RECOVERY FILES BLSPXXXX, ETC)        
!     NFCST     INTEGER   0 = FORECAST RECOVERY FILES FETCHED      
!                              1 = ANALYSIS RECOVERY FILES FETCHED      
! ----------------------------------------------------------------------

      END MODULE YOWWAMI
