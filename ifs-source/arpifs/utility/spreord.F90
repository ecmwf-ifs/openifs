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

SUBROUTINE SPREORD(YDDIM,YDEDIM,YDLEP,KFLEV,PSPFILE,PSPBUF,LD_FILE_TO_MODEL)

!**** *SPREORD*  - Reorder spectral data from old to new (grib-like) ordering
!                  or the reverse.

!     Purpose.
!     --------
!           Reorder spectral data structure from file ordering to model ordering 
!           or vice-versa.
!           SM : model ordering is the global spectrum ordering used in the model.

!**   Interface.
!     ----------
!        *CALL* *SPREORD

!        Explicit arguments :
!        --------------------
!        PSPFILE          : Spectral array ready to be used for file
!        PSPBUF           : Spectral (buffer) array ready to be used in model
!        KFLEV            : Number of fields
!        LD_FILE_TO_MODEL : .TRUE. to order from file to model data structure

!        Implicit arguments :
!        --------------------
!        None

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!     None

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      Lars Isaksen *ECMWF*
!      Original : 95-06-10

!     Modifications.
!     --------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      R. El Khatib : 30-Mar-2010 : Optimisation (reverse arrays ordering)
!     ------------------------------------------------------------------

USE YOMDIM   , ONLY : TDIM
USE YEMDIM   , ONLY : TEDIM
USE YEMLAP   , ONLY : TLEP
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMCT0   , ONLY : LELAM

IMPLICIT NONE

TYPE(TDIM)        ,INTENT(IN)    :: YDDIM
TYPE(TEDIM)       ,INTENT(IN)    :: YDEDIM
TYPE(TLEP)        ,INTENT(IN)    :: YDLEP
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLEV 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPFILE(KFLEV,YDDIM%NSEFRE) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPBUF(KFLEV,YDDIM%NSPEC2G) 
LOGICAL           ,INTENT(IN)    :: LD_FILE_TO_MODEL 
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "espareord.intfb.h"
#include "spareord.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SPREORD',0,ZHOOK_HANDLE)
ASSOCIATE(NSEFRE=>YDDIM%NSEFRE, NSPEC2G=>YDDIM%NSPEC2G)
!     ------------------------------------------------------------------

IF (LELAM) THEN
  CALL ESPAREORD(YDDIM,YDEDIM,YDLEP,KFLEV,PSPFILE,PSPBUF,LD_FILE_TO_MODEL)
ELSE
  CALL SPAREORD(YDDIM,KFLEV,PSPFILE,PSPBUF,LD_FILE_TO_MODEL)
ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SPREORD',1,ZHOOK_HANDLE)
END SUBROUTINE SPREORD
