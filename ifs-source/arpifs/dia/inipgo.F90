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

SUBROUTINE INIPGO(YDDIMV,YDMDDH,KBLOCKS,KNBARP)

!**** *INIPGO*  - DDH output : inititialize pseudo-grib coding

!     Purpose.
!     --------
!     Compute space needed for pseudo grib coding of DDH.
!     Allocate space for real and integer blocks for subsequent
!     pseudo grib coding and output of DDH data.
!     Reset pointers and counters used for coding.

!**   Interface.
!     ----------
!        *CALL* *INIPGO(..)*

!        Explicit arguments :     KBLOCKS - number of blocks : global=1
!        --------------------                                  zonal =NDHKD
!                                                              masks =NDHNOM
!                                 KNBARP  - number of records

!        Implicit arguments :      None.
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.   SUALPGO - memory manager allocation
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 93-03-01
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!     ------------------------------------------------------------------

USE YOMDIMV  , ONLY : TDIMV
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMMDDH  , ONLY : TMDDH
USE YOMPGO   , ONLY : NPTCH, NPTIN, NPTRE, NPTIND,&
 & NLENNAM, NLENRBL2, NLENREA, NLENIND, NLENINT, NLENCHA  

IMPLICIT NONE

TYPE(TDIMV)       , INTENT(IN) :: YDDIMV
TYPE(TMDDH)       , INTENT(INOUT):: YDMDDH
INTEGER(KIND=JPIM), INTENT(IN) :: KBLOCKS 
INTEGER(KIND=JPIM), INTENT(IN) :: KNBARP 
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "sualpgo.intfb.h"

!     ------------------------------------------------------------------

!*       1.  COMPUTE DIMENSIONS NEEDED - ALLOCATE SPACE.
!            -------------------------------------------

IF (LHOOK) CALL DR_HOOK('INIPGO',0,ZHOOK_HANDLE)
ASSOCIATE(NFLEVG=>YDDIMV%NFLEVG, &
 & NDHCSSU=>YDMDDH%NDHCSSU, NDHCVSU=>YDMDDH%NDHCVSU, NDHVS=>YDMDDH%NDHVS, &
 & NDHVV=>YDMDDH%NDHVV)
NLENREA=KBLOCKS*((NFLEVG+1)*(NDHCVSU+NDHVV)+(NDHCSSU+NDHVS))+1000
NLENIND=KNBARP
NLENINT=NLENIND+1000
NLENNAM=12

CALL SUALPGO

!     ------------------------------------------------------------------

!*       2.  RESET POINTERS.
!            ---------------

NPTCH=1
NPTIN=1
NPTRE=1
NPTIND=0
NLENRBL2=0

!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('INIPGO',1,ZHOOK_HANDLE)
END SUBROUTINE INIPGO
