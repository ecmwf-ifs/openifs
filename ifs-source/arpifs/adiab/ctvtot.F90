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

SUBROUTINE CTVTOT(YDGEOMETRY,YGFL,KSTART,KPROF,PT,PGFL)

!**** *CTVTOT* - Computes T from 'TV' (TV defined as RT/Rd)

!     Purpose.
!     --------
!           COMPUTES T FROM 'TV'

!**   Interface.
!     ----------
!        *CALL* *CTVTOT(...)

!        Explicit arguments :
!        --------------------
!                              PT(NPROMA,NFLEVG)    - output: T
!                                                  - input: TV
!                              PQ(NPROMA,NFLEVG)    - Q

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Clive Temperton  *ECMWF*

!     Modifications.
!     --------------
!        Original : 94-05-16
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        Modified 04-02-06: Y. Seity : new arguments to GPRCP for AROME
!        M.Hamrud  15-Jan-2006  Revised GPRCP
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCST       , ONLY : RD
USE YOM_YGFL     , ONLY : TYPE_GFLD

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TYPE_GFLD)   ,INTENT(IN)    :: YGFL
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTART 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PT(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGFL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YGFL%NDIM) 

!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: ZR(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)

INTEGER(KIND=JPIM) :: IEND, IST, JLEV, JROF
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "gprcp.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CTVTOT',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP)
ASSOCIATE(NDIM=>YGFL%NDIM, &
 & NPROMA=>YDDIM%NPROMA, &
 & NFLEVG=>YDDIMV%NFLEVG)
!     ------------------------------------------------------------------

!*       1.    COMPUTES T FROM 'TV'
!              --------------------

IST =KSTART
IEND=KPROF

CALL GPRCP(NPROMA,IST,IEND,NFLEVG,PR=ZR,PGFL=PGFL)

DO JLEV = 1, NFLEVG
  DO JROF = IST, IEND
    PT(JROF,JLEV)=RD*PT(JROF,JLEV)/ZR(JROF,JLEV)
  ENDDO
ENDDO

!-----------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CTVTOT',1,ZHOOK_HANDLE)
END SUBROUTINE CTVTOT
