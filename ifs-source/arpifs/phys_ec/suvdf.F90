! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction
SUBROUTINE SUVDF(YDVDF)

!     ------------------------------------------------------------------

!**   *SUVDF* IS THE SET-UP ROUTINE FOR COMMON BLOCK *YOEVDF*

!     A.C.M. BELJAARS         E.C.M.W.F.       2/11/89

!     PURPOSE
!     -------

!          THIS ROUTINE INITIALIZES THE CONSTANTS IN COMMON BLOCK
!     *YOEVDF*

!     INTERFACE.
!     ----------

!     CALL *SUVDF* FROM *SUPHEC*

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------

!     NONE.

!     REFERENCE.
!     ----------

!     MODIFICATIONS
!     -------------
!     J.-J. MORCRETTE         E.C.M.W.F.      91/07/14
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        A. Beljaars      Jan-2014 Inclusion of WDS numerics
!        P. Bechtold      Dec-2014 Inclusion of namelist call
!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPRB,JPIM
USE YOMHOOK   , ONLY : LHOOK,   DR_HOOK, JPHOOK

!USE YOEVDF   , ONLY : RLAM, RVDIFTS, &
!                    &LWDS, REPS1WDS, REPS2WDS ,RETAWDS, NSUBST
USE YOEVDF    , ONLY : TVDF
USE YOMDYNCORE, ONLY: RPLRG
USE YOMLUN    , ONLY : NULOUT, NULNAM

IMPLICIT NONE
REAL(KIND=JPRB),POINTER :: RLAM
REAL(KIND=JPRB),POINTER :: RVDIFTS
REAL(KIND=JPRB),POINTER :: RTOFDALPHA
REAL(KIND=JPRB),POINTER :: REISTHSC
LOGICAL,POINTER         :: LWDS 
INTEGER(KIND=JPIM),POINTER :: NSUBST

TYPE(TVDF), INTENT(INOUT),TARGET :: YDVDF
REAL(KIND=JPRB) :: ZC1,ZC2,ZP
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "namvdf.nam.h"
#include "posnam.intfb.h"

!     ------------------------------------------------------------------

!*         1.     SET FIRST SET OF CONSTANTS
!                 --------------------------

IF (LHOOK) CALL DR_HOOK('SUVDF',0,ZHOOK_HANDLE)
ASSOCIATE(REPS1WDS=>YDVDF%REPS1WDS,REPS2WDS=>YDVDF%REPS2WDS,&
 & RETAWDS=>YDVDF%RETAWDS)
RLAM => YDVDF%RLAM
RVDIFTS => YDVDF%RVDIFTS
LWDS => YDVDF%LWDS
NSUBST => YDVDF%NSUBST
REISTHSC=>YDVDF%REISTHSC

RLAM   =150._JPRB/RPLRG

RTOFDALPHA => YDVDF%RTOFDALPHA
!RTOFDALPHA=35. ! Omit _JPRB for bit-identicality
RTOFDALPHA=24._JPRB

!     ------------------------------------------------------------------

!*         2.      SET OTHER CONSTANTS
!                  -------------------

REISTHSC=5.0_JPRB ! Threshold of Inversion strength for Stratocumulus

!     ------------------------------------------------------------------

!*         3.      NUMERICS CONSTANTS
!                  -------------------
LWDS=.FALSE.       ! .T. for Wood, Diamantakis and Staniforth
RVDIFTS=1.5_JPRB   ! Overimplicitness constant in case LWDS=.F.

ZP=1.0_JPRB        ! P-constant controls WDS time integration
ZC1=SQRT(2._JPRB)
ZC2=1._JPRB+1._JPRB/ZC1
REPS1WDS=ZC2*(ZP+1._JPRB/ZC1+SQRT(ZP*(ZC1-1._JPRB)+0.5_JPRB))
REPS2WDS=ZC2*(ZP+1._JPRB/ZC1-SQRT(ZP*(ZC1-1._JPRB)+0.5_JPRB))
RETAWDS =ZC2*(1._JPRB+ZP)


! Set number of substeps for vertical diffusion
NSUBST=1

CALL POSNAM(NULNAM,'NAMVDF')
READ(NULNAM,NAMVDF)

IF (LWDS) THEN     
  NSUBST=1
ENDIF

WRITE(NULOUT,*)'SUVDF: LWDS=',LWDS,' NSUBST=',NSUBST,' RVDIFTS=',RVDIFTS,' RLAM=',RLAM,&
      &' RTOFDALPHA=',RTOFDALPHA

!     ------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUVDF',1,ZHOOK_HANDLE)
END SUBROUTINE SUVDF
