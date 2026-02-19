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

SUBROUTINE SUGRXFU(YDGEOMETRY,YDXFU,YDRIP,YDML_LBC,KFILE)

!**** *SUGRXFU*  - Initialize the instantaneous fluxes from *FA*

!     Purpose.
!     --------
!           Initialize the instantaneous fluxes fields of the model from FA.

!**   Interface.
!     ----------
!        *CALL* *SUGRXFU(....)*

!        Explicit arguments : KFILE - indicator for which file is to be read
!        ------------------   

!        Implicit arguments :
!        ------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS
!        Note de travail ARPEGE Nr 12 et 17

!     Author.
!     -------
!      R. El Khatib  *METEO-FRANCE*
!      Original : 97-12-05

!     Modifications.
!     --------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      G.Mozdzynski (Feb 2011): OOPS cleaning, use of derived type TGSGEOM
!      T. Wilhelmsson (Sept 2014) Geometry and setup refactoring.
!      P. Marguinaud  (Oct 2014)  Use IOXFU_MOD
!     ------------------------------------------------------------------

USE YOMRIP       , ONLY : TRIP
USE YOMXFU       , ONLY : TXFU
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK
USE IOXFU_MOD    , ONLY : NIOXFUACT_READ,&
                    &     IOXFU_COUNT,&
                    &     IOXFU_SELECTD
USE IOGRID_MOD   , ONLY : NIOGRIDCT_READ,&
                     &    IOGRID_SELECTF
USE IOFLDPTR_MOD , ONLY : IOFLDPTR
USE YEMLBC_MODEL   , ONLY : TELBC_MODEL

IMPLICIT NONE

TYPE (GEOMETRY)    , INTENT(IN)    :: YDGEOMETRY
TYPE (TRIP)        , INTENT(INOUT) :: YDRIP
TYPE (TXFU)        , INTENT(INOUT) :: YDXFU
TYPE (TELBC_MODEL)     , INTENT (IN)   :: YDML_LBC 
INTEGER (KIND=JPIM), INTENT (IN)   :: KFILE 

#include "rdfa2gp.intfb.h"

TYPE (IOFLDPTR),  ALLOCATABLE :: YLFLDGT (:)
REAL (KIND=JPRB), ALLOCATABLE :: ZBUFL (:,:)
INTEGER (KIND=JPIM) :: IFILE
INTEGER (KIND=JPIM) :: IFNUM

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('SUGRXFU',0,ZHOOK_HANDLE)


IF (KFILE == 7) THEN
  IFILE = 40
ELSE
  IFILE= KFILE
ENDIF

IFNUM = 0
CALL IOXFU_COUNT(YDGEOMETRY,YDXFU,NIOXFUACT_READ,IFNUM)

IF (IFNUM > 0) THEN

  ALLOCATE (YLFLDGT (IFNUM), ZBUFL (YDGEOMETRY%YRGEM%NGPTOT, IFNUM))

  CALL IOXFU_SELECTD(YDGEOMETRY, YDXFU, NIOXFUACT_READ, YLFLDGT)

  CALL RDFA2GP (YDGEOMETRY, YDRIP, IFNUM, ZBUFL, YLFLDGT%YFLDDSC, &
              & KFILE=IFILE, YDML_LBC=YDML_LBC)

  CALL IOGRID_SELECTF (YDGEOMETRY, NIOGRIDCT_READ, ZBUFL, YLFLDGT)
  
  DEALLOCATE (YLFLDGT, ZBUFL)

ENDIF 


IF (LHOOK) CALL DR_HOOK('SUGRXFU',1,ZHOOK_HANDLE)

END SUBROUTINE SUGRXFU

