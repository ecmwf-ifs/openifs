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

SUBROUTINE GPINO3CH(YDDIMV,YDOZO,YDDPHY,KSTART,KPROF,KPROMA,KGPLAT,PKOZO)

!**** *GPINO3CH* - Ozone chemistry.
 
!     Purpose.  
!     --------
!      Memory transfer: fills PKOZO with the content of buffer TOZ2DG.

!**   Interface.
!     ----------
!        *CALL* *GPINO3CH(...)

!        Explicit arguments :
!        --------------------

!        INPUT:
!         KSTART       : start of work
!         KPROF        : depth of work
!         KPROMA       : horizontal dimension
!         KGPLAT       : DM-global number of the latitude of point jrof=KSTART

!        OUTPUT:
!         PKOZO        : fields for photochemistery of ozone

!        Implicit arguments :
!        --------------------
  
!     Method.
!     -------
!        See documentation

!     Externals.  None. 
!     ----------
!      Called by EC_PHYS.
    
!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        ????

! Modifications.
! --------------
!   Original : ????
!   K. Yessad (Sep 2008): add missing comments, cleanings.
!   F. Vana : Cleaner 2D->3D array conversion.
! End Modifications
!------------------------------------------------------------------

USE YOMDIMV  , ONLY : TDIMV
USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMDPHY  , ONLY : TDPHY
USE YOMOZO   , ONLY : TOZO

! -------------------------------------------------------------------------

IMPLICIT NONE

TYPE(TDIMV)       ,INTENT(IN)    :: YDDIMV
TYPE(TDPHY)       ,INTENT(IN)    :: YDDPHY
TYPE(TOZO)        ,INTENT(IN)    :: YDOZO
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTART 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA
INTEGER(KIND=JPIM),INTENT(IN)    :: KGPLAT(KPROMA)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PKOZO(KPROMA,YDDIMV%NFLEVG,YDDPHY%NVCLIS) 

! -------------------------------------------------------------------------

INTEGER(KIND=JPIM) ::  IGL, JJ, JVAR, J1, J2
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! -------------------------------------------------------------------------

#include "abor1.intfb.h"

! -------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GPINO3CH',0,ZHOOK_HANDLE)
ASSOCIATE(NFLEVG=>YDDIMV%NFLEVG, &
 & NTOZ2D=>YDDPHY%NTOZ2D, NTOZ3D=>YDDPHY%NTOZ3D, NVCLIS=>YDDPHY%NVCLIS, &
 & TOZ2DG=>YDOZO%TOZ2DG)
! -------------------------------------------------------------------------

IF(NTOZ3D > 0) THEN
  CALL ABOR1('GPINIO3CH : NTOZ3D > 0 BROKEN WITH MESSAGE PASSING')
ENDIF

DO JVAR=1,NTOZ2D*NVCLIS*NFLEVG
  J1=MOD(JVAR-1,NFLEVG) +1
  J2=   (JVAR-1)/NFLEVG +1
  DO JJ=KSTART,KPROF
    IGL=KGPLAT(JJ)
    PKOZO(JJ,J1,J2)=TOZ2DG(JVAR,IGL)
  ENDDO
ENDDO

! -------------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('GPINO3CH',1,ZHOOK_HANDLE)
END SUBROUTINE GPINO3CH
