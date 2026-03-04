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

SUBROUTINE SUMETRIC(YDLAP,YDDIM,KSMAX,PMET,PMETDER,PMETKE)

!**** *SUMETRIC* - Compute metric for spectral norms computation.

!     Purpose.
!     --------
!        To compute the metric used for spectral norms computation, 
!        and to initialise the zonal truncation.

!**   Interface.
!     ----------
!        *CALL* *SUMETRIC(KINDPP)

!        Explicit arguments :
!        --------------------
!          KSMAX   : zonal truncation
!          PMET    : metric for basic fields
!          PMETDER : metric for derivatives
!          PMETKE  : metric for kinetic energy

!        Implicit arguments :  See modules above
!        --------------------

!     Method.
!     -------

!     Externals.
!     ----------
!      None.

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Ryad El Khatib *METEO-FRANCE*

!     Modifications.
!     --------
!        Original : 99-02-08 from SPNORM
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!     ------------------------------------------------------------------

USE YOMLAP   , ONLY : TLAP
USE YOMDIM   , ONLY : TDIM
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK

IMPLICIT NONE
TYPE(TLAP)        ,INTENT(IN)    :: YDLAP
TYPE(TDIM)        ,INTENT(IN)    :: YDDIM
INTEGER(KIND=JPIM),INTENT(OUT)   :: KSMAX 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PMET(0:YDDIM%NSMAX) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PMETDER(0:YDDIM%NSMAX) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PMETKE(0:YDDIM%NSMAX) 
INTEGER(KIND=JPIM) :: JM
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

!*       1.    ZONAL TRUNCATION
!              ----------------

IF (LHOOK) CALL DR_HOOK('SUMETRIC',0,ZHOOK_HANDLE)
ASSOCIATE(NSMAX=>YDDIM%NSMAX, &
 & RLAPDI=>YDLAP%RLAPDI, RLAPIN=>YDLAP%RLAPIN)
KSMAX=NSMAX

!*       2.    METRIC
!              ------

DO JM=0,NSMAX
  PMET   (JM) = 1.0_JPRB
  PMETDER(JM) = -RLAPDI(JM)
  PMETKE (JM) = -RLAPIN(JM)
ENDDO

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUMETRIC',1,ZHOOK_HANDLE)
END SUBROUTINE SUMETRIC
