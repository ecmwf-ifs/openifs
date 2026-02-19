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

SUBROUTINE LATTE_STDDIS(YDGEOMETRY,YDDYNA,YDGMV,YDPTRSLB2,KST,KPROF,PDT,PGMV,PB2)

!**** *LATTE_STDDIS* input for corrected SL scheme COMAD.
!                    Computation of the stretching/shrinking functions
!                    in the fluid (based on the derivative of the velocity,
!                    cf definition of the 1D divergence/deformation)
!                    Functions are stored in PB2 to be used later under CALL_SL.

!     Purpose.
!     --------

!**   Interface.
!     ----------
!        *CALL* *LATTE_STDDIS(..)

!        Explicit arguments :
!        --------------------

!        INPUT:
!          KST      - first element of work.
!          KPROF    - depth of work.
!          PDT      - time step 
!          PGMV     - variables at time t, only wind derivatives will be used in PGMV

!        INPUT/OUTPUT:
!          PB2      - "SLBUF2" buffer.

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        Arpege documentation about semi-Lagrangian scheme.

!     Author.
!     -------
!        Original (S.Malardel): november 2013

!     Modifications.
!     --------------
!     F. Vana   26-Nov-2021   Reducing the content to a wrapper (do we really need this extra layer?)

!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMDYNA      , ONLY : TDYNA
USE YOMGMV       , ONLY : TGMV
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK
USE PTRSLB2      , ONLY : TPTRSLB2

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TDYNA)       ,INTENT(IN)    :: YDDYNA
TYPE(TGMV)        ,INTENT(INOUT) :: YDGMV
TYPE(TPTRSLB2)    ,INTENT(IN)    :: YDPTRSLB2
INTEGER(KIND=JPIM),INTENT(IN)    :: KST 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDT 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGMV(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDGMV%NDIMGMV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PB2(YDGEOMETRY%YRDIM%NPROMA,YDPTRSLB2%NFLDSLB2)
!     -------------------------------------------------------

REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

!     -------------------------------------------------------

#include "gp_stddis.intfb.h"

!     -------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('LATTE_STDDIS',0,ZHOOK_HANDLE)
ASSOCIATE(NPROMA=>YDGEOMETRY%YRDIM%NPROMA, NFLEVG=>YDGEOMETRY%YRDIMV%NFLEVG, &
 & YT0=>YDGMV%YT0, MSLB2STDDISW=>YDPTRSLB2%MSLB2STDDISW, &
 & MSLB2STDDISU=>YDPTRSLB2%MSLB2STDDISU, MSLB2STDDISV=>YDPTRSLB2%MSLB2STDDISV)
!     -------------------------------------------------------

!*       1.   Stretching functions computation

CALL GP_STDDIS(YDDYNA,NPROMA,KST,KPROF,NFLEVG,PDT,&
 & PGMV(1,1,YT0%MUL),PGMV(1,1,YT0%MDIV),&
 & PB2(1,MSLB2STDDISU),PB2(1,MSLB2STDDISV),PB2(1,MSLB2STDDISW))

!    --------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('LATTE_STDDIS',1,ZHOOK_HANDLE)
END SUBROUTINE LATTE_STDDIS
