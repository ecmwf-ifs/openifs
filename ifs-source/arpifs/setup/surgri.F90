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

SUBROUTINE SURGRI(YDGEOMETRY,KULNAM,KULOUT,KDGL,KDGSA,KDGEN,KLOEN)

!**** *SURGRI*  - Routine to read or compute KLOEN from Namelist

!     Purpose.
!     --------
!        Reads from namelist NRGRI,which corresponds to NLOEN without
!        extra-rows, this interface is necessary because
!        static dimensioning seems necessary for namelists ; compute also
!        the reduced grid according to the allowed rate of aliasing given
!        in namelist

!**   Interface.
!     ----------

!     *CALL* SURGRI(...)

!        Explicit arguments :
!        --------------------
!     * INPUT:
!     KULNAM:  LOGICAL UNIT FOR NAMELIST
!     KULOUT:  LOGICAL UNIT FOR OUTPUT
!     KDGL  :  NUMBER OF GAUSSIAN LATITUDES (as NDGLG of YOMDIM)
!     KDGSA :  as NDGSAG of YOMDIM
!     KDGEN :  as NDGENG of YOMDIM

!     * OUTPUT:
!     KLOEN :  as NLOENG of YOMGEM

!        Implicit arguments :
!        --------------------

!     Method.
!     -------

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department DOcumentation of the IFS

!     Author.
!     -------
!      Philippe Courtier  *ECMWF*
!      Original : 92-03-18

!     Modifications.
!     --------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE PARDIM   , ONLY : JPMXGL
USE YOMRGRI  , ONLY : NRGRI

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY), INTENT(INOUT)    :: YDGEOMETRY
INTEGER(KIND=JPIM),INTENT(IN)    :: KDGSA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDGEN 
INTEGER(KIND=JPIM),INTENT(IN)    :: KULNAM 
INTEGER(KIND=JPIM),INTENT(IN)    :: KULOUT 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDGL 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KLOEN(KDGSA:KDGEN) 

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: ILOTOT, J, JGL
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "posnam.intfb.h"

#include "namrgri.nam.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SURGRI',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP)
ASSOCIATE(NDGNH=>YDDIM%NDGNH)
!     ------------------------------------------------------------------

!*         1. READ NAMELIST
!          ----------------

IF(KDGL >= JPMXGL)THEN
  CALL ABOR1('SURGRI:KDGL >= JPMXGL')
ENDIF

CALL POSNAM(KULNAM,'NAMRGRI')
READ(KULNAM,NAMRGRI)

!     ------------------------------------------------------------------

!*         2. MODIFY REDUCED GRID ACCORDING TO NAMELIST
!          --------------------------------------------

WRITE(KULOUT,FMT='('' CALLING SURGRI'')')
DO J=1,(KDGL+1)/2
  KLOEN(J)=NRGRI(J)
ENDDO
DO J=1,KDGL/2
  KLOEN(KDGL+1-J)=NRGRI(J)
ENDDO

IF(KDGSA <= 0)THEN
  KLOEN(0)=KLOEN(1)
  KLOEN(KDGL+1)=KLOEN(KDGL)
ENDIF

WRITE(KULOUT,FMT='('' (JGL,NRGRI) '')')
WRITE(KULOUT,FMT='(8(1X,''('',I4,I5,'')''))')&
 & (JGL,NRGRI(JGL),JGL=1,NDGNH)  

ILOTOT=0
DO JGL=1,KDGL
  ILOTOT=ILOTOT+NRGRI(JGL)
ENDDO
WRITE(KULOUT,FMT='('' NUMBER OF POINTS ON REDUCED GRID = '',I8)') ILOTOT
!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SURGRI',1,ZHOOK_HANDLE)
END SUBROUTINE SURGRI
