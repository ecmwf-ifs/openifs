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

SUBROUTINE SUNH_VERTFESPLINE(YDGEOMETRY)

!**** *SUNH_VERTFESPLINE*  - Setup in Non-Hydrostatic model for
!                            VERTical Finite Element scheme with SPLINEs
!                            of general order

!     Purpose.
!     --------
!       Initialisation routine for the vertical finite element scheme
!       with general order B-splines.

!**   Interface.
!     ----------
!       *CALL* SUNH_VERTFESPLINE

!     Explicit arguments :
!     --------------------
!       None

!     Implicit arguments :
!     --------------------

!     Method.
!     -------

!     Externals.
!     ----------
!       see below

!     Reference.
!     ----------
!       ALADIN-NH documentation.

!     Author.
!     -------
!       Petra Smolikova (CHMI) after suvertfeb.
!       Original : 2017-09

!     Modifications.
!     --------------
!       J. Vivoda and P. Smolikova (Sep 2017): new options for VFE-NH
!       P. Smolikova (Sep 2020): VFE pruning.
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMDIMV  , ONLY : TDIMV
USE PARKIND1 , ONLY : JPRB ,JPIM
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN   , ONLY : NULOUT

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY), INTENT(INOUT) :: YDGEOMETRY
INTEGER(KIND=JPIM) :: IT_IN(0:YDGEOMETRY%YRDIMV%NFLEVG+1)
INTEGER(KIND=JPIM) :: IT_OUT(0:YDGEOMETRY%YRDIMV%NFLEVG+1)
INTEGER(KIND=JPIM) :: ITBC_IN (2), IBBC_IN (2)
INTEGER(KIND=JPIM) :: ITBC_OUT(2), IBBC_OUT(2)
REAL(KIND=JPRB)    :: ZETAF(0:YDGEOMETRY%YRDIMV%NFLEVG+1)
REAL(KIND=JPRB)    :: ZETAF_IN(0:YDGEOMETRY%YRDIMV%NFLEVG+1)

! for EIGSOL
REAL(KIND=JPRB) :: ZM(YDGEOMETRY%YRDIMV%NFLEVG+2, YDGEOMETRY%YRDIMV%NFLEVG+2)
REAL(KIND=JPRB) :: ZFR(YDGEOMETRY%YRDIMV%NFLEVG+2)
REAL(KIND=JPRB) :: ZFI(YDGEOMETRY%YRDIMV%NFLEVG+2)
REAL(KIND=JPRB) :: ZFN(YDGEOMETRY%YRDIMV%NFLEVG+2)
REAL(KIND=JPRB) :: ZMO(YDGEOMETRY%YRDIMV%NFLEVG+2,YDGEOMETRY%YRDIMV%NFLEVG+2)
REAL(KIND=JPRB) :: ZWO(YDGEOMETRY%YRDIMV%NFLEVG+3)
INTEGER(KIND=JPIM) :: IWO(YDGEOMETRY%YRDIMV%NFLEVG+3), I

INTEGER(KIND=JPIM) :: IER
REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

LOGICAL            :: LLIO, LLBC_MATERIAL_SURFACES

CHARACTER(LEN=15)  :: CLTAG

!     ------------------------------------------------------------------

#include "suvertfeb.intfb.h"
#include "eigsol.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUNH_VERTFESPLINE',0,ZHOOK_HANDLE)
ASSOCIATE(YDVETA=>YDGEOMETRY%YRVETA,YDVFE=>YDGEOMETRY%YRVFE,YDDIMV=>YDGEOMETRY%YRDIMV,YDCVER=>YDGEOMETRY%YRCVER)
ASSOCIATE(NFLEVG=>YDGEOMETRY%YRDIMV%NFLEVG)

!     ------------------------------------------------------------------

! this is switch for consistency of input assumptions and output values
! for example: when we assume zero derivative at the top
! we then assume also zero value in output of derivative operator
LLIO    = .TRUE.

LLBC_MATERIAL_SURFACES = .TRUE.

!  1. Preliminary calculations:
!  ----------------------------

ZETAF(0       ) = 0.0_JPRB
ZETAF(1:NFLEVG) = YDVETA%VFE_ETAF(1:NFLEVG)
ZETAF(NFLEVG+1) = 1.0_JPRB


ZETAF_IN(1:NFLEVG  ) = ZETAF(1:NFLEVG  )
IF(LLBC_MATERIAL_SURFACES)THEN
  ZETAF_IN(0         ) = 0.0_JPRB
  ZETAF_IN(  NFLEVG+1) = 1.0_JPRB
ELSE
  ZETAF_IN(0         ) = ZETAF(1)
  ZETAF_IN(  NFLEVG+1) = ZETAF(NFLEVG)
ENDIF

!*  2. Integral operators
!   ---------------------

IF(YDCVER%LVFE_ECMWF)THEN

  WRITE(NULOUT,*) "SUVERTFE RINTE"
  !-------------------------------
  CLTAG = "RINTE"
  IT_IN       = 0
  IT_OUT      = 0
  ITBC_IN (1) = 0 ; ITBC_IN (2) = 0
  IBBC_IN (1) = 0 ; IBBC_IN (2) = 0
  ITBC_OUT(1) = 1 ; ITBC_OUT(2) = 0
  IBBC_OUT(1) = 0 ; IBBC_OUT(2) = 0

  CALL SUVERTFEB(YDVETA,YDDIMV,YDVFE,YDCVER,CLTAG, &
   & .FALSE.,-1,NFLEVG,ZETAF(1:NFLEVG),IT_IN(1:NFLEVG), &
   & NFLEVG+1,ZETAF(1:NFLEVG+1),IT_OUT(1:NFLEVG+1), &
   & ITBC_IN,IBBC_IN,ITBC_OUT,IBBC_OUT,YDVFE%RINTE)
ELSE
  WRITE(NULOUT,*) "SUVERTFE RINTBF11"
  !----------------------------------
  CLTAG = "RINTBF11"

  IT_IN              = 0
  IT_IN (0         ) = 1
  IT_IN (  NFLEVG+1) = 1
  IT_OUT             = 0
  ITBC_IN (1) = 0 ; ITBC_IN (2) = 0
  IBBC_IN (1) = 0 ; IBBC_IN (2) = 0
  ITBC_OUT(1) = 1 ; ITBC_OUT(2) = 0
  IBBC_OUT(1) = 0 ; IBBC_OUT(2) = 0

  CALL SUVERTFEB(YDVETA,YDDIMV,YDVFE,YDCVER,CLTAG, &
   & .FALSE.,-1,NFLEVG+2,ZETAF_IN,IT_IN, &
   & NFLEVG+1,ZETAF(1:NFLEVG+1),IT_OUT(1:NFLEVG+1), &
   & ITBC_IN,IBBC_IN,ITBC_OUT,IBBC_OUT,YDVFE%RINTBF11)

  ZM = 0.0_JPRB
  ZM(1:NFLEVG+1,1:NFLEVG+2) = YDVFE%RINTBF11(1:NFLEVG+1,1:NFLEVG+2)

  CALL EIGSOL(NFLEVG+2,NFLEVG+2,ZM,ZFR,ZFI,0,ZMO,IWO,ZWO,IER)

  DO I = 1, NFLEVG + 2
    ZFN(I)=SQRT(ZFR(I)*ZFR(I)+ZFI(I)*ZFI(I))
    WRITE(NULOUT,'("RINTBF11 EIGENVALUES :: ",I5,3(1X,F12.6))') I, ZFN(I), ZFR(I), ZFI(I)
    CALL FLUSH(NULOUT)
  ENDDO
ENDIF

!*  3. First order derivative operators.
!   ------------------------------------

IF(YDCVER%LVFE_ECMWF)THEN
  WRITE(NULOUT,*) "SUVERTFE RDERI"
  ! ------------------------------
  CLTAG = "RDERI"

  IT_IN       = 0
  IT_OUT      = 0
  ITBC_IN (1) = 0 ; ITBC_IN (2) = 0
  IBBC_IN (1) = 0 ; IBBC_IN (2) = 0
  ITBC_OUT(1) = 1 ; ITBC_OUT(2) = 0
  IBBC_OUT(1) = 1 ; IBBC_OUT(2) = 0

  CALL SUVERTFEB(YDVETA,YDDIMV,YDVFE,YDCVER,CLTAG, &
   & .FALSE.,1,NFLEVG,ZETAF(1:NFLEVG),IT_IN(1:NFLEVG), &
   & NFLEVG,ZETAF(1:NFLEVG),IT_OUT(1:NFLEVG), &
   & ITBC_IN,IBBC_IN,ITBC_OUT,IBBC_OUT,YDVFE%RDERI)
ELSE
  WRITE(NULOUT,*) "SUVERTFE RDERBF11"
  ! -----------------------------------
  CLTAG = "RDERBF11"

  IT_IN              = 0
  IT_IN (0         ) = 1
  IT_IN (  NFLEVG+1) = 1
  IT_OUT             = 0
  ITBC_IN (1) = 0 ; ITBC_IN (2) = 0
  IBBC_IN (1) = 0 ; IBBC_IN (2) = 0

  IF(LLIO)THEN
    ! ITBC_OUT(1) = 1 ; ITBC_OUT(2) = 0
    ! IBBC_OUT(1) = 1 ; IBBC_OUT(2) = 0
    ITBC_OUT(1) = 0 ; ITBC_OUT(2) = 0
    IBBC_OUT(1) = 0 ; IBBC_OUT(2) = 0
  ELSE
    ITBC_OUT(1) = 0 ; ITBC_OUT(2) = 0
    IBBC_OUT(1) = 0 ; IBBC_OUT(2) = 0
  ENDIF

  CALL SUVERTFEB(YDVETA,YDDIMV,YDVFE,YDCVER,CLTAG, &
   & .FALSE.,1,NFLEVG+2,ZETAF,IT_IN, &
   & NFLEVG,ZETAF(1:NFLEVG),IT_OUT(1:NFLEVG), &
   & ITBC_IN,IBBC_IN,ITBC_OUT,IBBC_OUT,YDVFE%RDERBF11)

  IF (YDGEOMETRY%LNONHYD_GEOM) THEN

    WRITE(NULOUT,*) "SUVERTFE RDERBF00"
    ! -----------------------------------
    CLTAG = "RDERBF00"

    IT_IN       = 0
    IT_OUT      = 0
    ITBC_IN (1) = 0 ; ITBC_IN (2) = 0
    IBBC_IN (1) = 0 ; IBBC_IN (2) = 0
    ITBC_OUT(1) = 0 ; ITBC_OUT(2) = 0
    IBBC_OUT(1) = 0 ; IBBC_OUT(2) = 0

    CALL SUVERTFEB(YDVETA,YDDIMV,YDVFE,YDCVER,CLTAG, &
     & .FALSE.,1,NFLEVG+2,ZETAF,IT_IN, &
     & NFLEVG,ZETAF(1:NFLEVG),IT_OUT(1:NFLEVG), &
     & ITBC_IN,IBBC_IN,ITBC_OUT,IBBC_OUT,YDVFE%RDERBF00)

    YDVFE%RDERB = YDVFE%RDERBF00

    WRITE(NULOUT,*) "SUVERTFE RDERBF01"
    ! -----------------------------------
    CLTAG = "RDERBF01"

    IT_IN              = 0
    IT_IN (  NFLEVG+1) = 1
    IT_OUT             = 0
    ITBC_IN (1) = 0 ; ITBC_IN (2) = 0
    IBBC_IN (1) = 0 ; IBBC_IN (2) = 0

    ITBC_OUT(1) = 0 ; ITBC_OUT(2) = 0
    IF(LLIO)THEN
      IBBC_OUT(1) = 1 ; IBBC_OUT(2) = 0
    ELSE
      IBBC_OUT(1) = 0 ; IBBC_OUT(2) = 0
    ENDIF

    CALL SUVERTFEB(YDVETA,YDDIMV,YDVFE,YDCVER,CLTAG, &
     & .FALSE.,1,NFLEVG+2,ZETAF,IT_IN, &
     & NFLEVG,ZETAF(1:NFLEVG),IT_OUT(1:NFLEVG), &
     & ITBC_IN,IBBC_IN,ITBC_OUT,IBBC_OUT,YDVFE%RDERBF01)

    WRITE(NULOUT,*) "SUVERTFE RDERBF10"
    ! -----------------------------------
    CLTAG = "RDERBF10"

    IT_IN              = 0
    IT_IN (0         ) = 1
    IT_OUT             = 0
    ITBC_IN (1) = 0 ; ITBC_IN (2) = 0
    IBBC_IN (1) = 0 ; IBBC_IN (2) = 0

    IF(LLIO)THEN
      ITBC_OUT(1) = 1 ; ITBC_OUT(2) = 0
    ELSE
      ITBC_OUT(1) = 0 ; ITBC_OUT(2) = 0
    ENDIF
    IBBC_OUT(1) = 0 ; IBBC_OUT(2) = 0

    CALL SUVERTFEB(YDVETA,YDDIMV,YDVFE,YDCVER,CLTAG, &
     & .FALSE.,1,NFLEVG+2,ZETAF,IT_IN, &
     & NFLEVG,ZETAF(1:NFLEVG),IT_OUT(1:NFLEVG), &
     & ITBC_IN,IBBC_IN,ITBC_OUT,IBBC_OUT,YDVFE%RDERBF10)
  ENDIF
ENDIF


!*  4. Second order derivative operator

!   -----------------------------------

IF(YDGEOMETRY%LNONHYD_GEOM)THEN

  WRITE(NULOUT,*) "SUVERTFE RDDERBF01"
! ------------------------------------
  CLTAG = "RDDERBF01"

  IT_IN              = 0
  IT_IN (  NFLEVG+1) = 1
  IT_OUT             = 0
  ITBC_IN (1) = 0 ; ITBC_IN (2) = 0
  IBBC_IN (1) = 0 ; IBBC_IN (2) = 0

  ITBC_OUT(1) = 0 ; ITBC_OUT(2) = 0
  IF(LLIO)THEN
    IBBC_OUT(1) = 1 ; IBBC_OUT(2) = 0
  ELSE
    IBBC_OUT(1) = 0 ; IBBC_OUT(2) = 0
  ENDIF

  CALL SUVERTFEB(YDVETA,YDDIMV,YDVFE,YDCVER,CLTAG, &
   & .FALSE.,2,NFLEVG+2,ZETAF,IT_IN, &
   & NFLEVG, ZETAF(1:NFLEVG), IT_OUT(1:NFLEVG), &
   & ITBC_IN,IBBC_IN,ITBC_OUT,IBBC_OUT,YDVFE%RDDERBF01)

ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUNH_VERTFESPLINE',1,ZHOOK_HANDLE)
END SUBROUTINE SUNH_VERTFESPLINE
