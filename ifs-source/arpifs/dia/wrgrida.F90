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

SUBROUTINE WRGRIDA(YDGEOMETRY,YDSURF,YDXFU,YDMODEL,YDFACTX,CDFIC)

!     Purpose.
!     --------
!        Write out the model surface fields to ARPEGE file

!**   Interface.
!     ----------
!        *CALL* *WRGRIDA(KNUMER)

!        Explicit arguments : KNUMER : logical unit number of the file
!        --------------------

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        Note de travail ARPEGE NR 17

!     Author.
!     -------
!        R. El Khatib *METEO-FRANCE*

!     Modifications.
!     --------------
!   Original 18-Mar-2010 (From WRMLPPA)
!            14-Jan-2011 Use IOGRIDA_MOD (P.Marguinaud)
!            11-Sep-2012 Use IOFLDDESC_MOD (P.Marguinaud)
!            10-Oct-2013 Use FACTX (P.Marguinaud)
!            O. Marsden: June 2015 CY42 YRGMV, YRGFL, YRSURF, YRGMV5, and YRGFL5 are now passed by argument
!     -------------------------------------------------------------------------

USE GEOMETRY_MOD       , ONLY : GEOMETRY
USE SURFACE_FIELDS_MIX , ONLY : TSURF
USE YOMXFU             , ONLY : TXFU
USE TYPE_MODEL         , ONLY : MODEL
USE PARKIND1           , ONLY : JPIM     ,JPRB
USE YOMHOOK            , ONLY : LHOOK,   DR_HOOK, JPHOOK

USE IOGRIDA_MOD        , ONLY : IOGRIDA_COUNT,&
                             & IOGRIDA_SELECTD,&
                             & NIOGRIDACT_WRITE
USE IOGRID_MOD         , ONLY : IOGRID_SELECTF,&
                             & NIOGRIDCT_WRITE
USE YOMTAG             , ONLY : MTAG_MFIO_WRGRIDA

USE IOFLDPTR_MOD       , ONLY : IOFLDPTR

USE FACTX_MOD          , ONLY : FACTX

IMPLICIT NONE

TYPE (GEOMETRY), INTENT (IN)    :: YDGEOMETRY
TYPE(TSURF)    , INTENT (INOUT) :: YDSURF
TYPE(TXFU)     , INTENT (INOUT) :: YDXFU
TYPE(MODEL)    , INTENT (INOUT) :: YDMODEL
TYPE (FACTX)   , INTENT (INOUT) :: YDFACTX
CHARACTER(LEN=*),INTENT (IN)    :: CDFIC

#include "wrgp2fa.intfb.h"

REAL (KIND=JPRB), ALLOCATABLE :: ZREAL (:,:)
TYPE (IOFLDPTR),  ALLOCATABLE :: YLFLDPT (:) 
INTEGER (KIND=JPIM)           :: IFNUM

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('WRGRIDA',0,ZHOOK_HANDLE)

IFNUM = 0
CALL IOGRIDA_COUNT(YDGEOMETRY,YDSURF,YDMODEL%YRML_AOC%YRMCC,YDMODEL%YRML_PHY_G%YRDPHY,YDMODEL%YRML_PHY_MF%YRPHY, &
 & YDMODEL%YRML_PHY_MF%YRPHYDS,NIOGRIDACT_WRITE,IFNUM)

IF (IFNUM > 0) THEN

  ALLOCATE (YLFLDPT (IFNUM))

  CALL IOGRIDA_SELECTD(YDGEOMETRY, YDSURF, YDMODEL%YRML_AOC%YRMCC, YDMODEL%YRML_PHY_G%YRDPHY, YDMODEL%YRML_PHY_MF%YRPHY, &
 & YDMODEL%YRML_PHY_MF%YRPHYDS, NIOGRIDACT_WRITE, YLFLDPT)
 
  ALLOCATE (ZREAL (YDGEOMETRY%YRGEM%NGPTOT, IFNUM))

  CALL IOGRID_SELECTF (YDGEOMETRY, NIOGRIDCT_WRITE, ZREAL, YLFLDPT)

  CALL WRGP2FA(YDGEOMETRY,YDMODEL%YRML_GCONF%YRRIP,IFNUM,ZREAL,YLFLDPT%YFLDDSC,YDFACTX, &
 & CDFIC,KTAG=MTAG_MFIO_WRGRIDA)

  DEALLOCATE (ZREAL)

  DEALLOCATE (YLFLDPT)

ENDIF

IF (LHOOK) CALL DR_HOOK('WRGRIDA',1,ZHOOK_HANDLE)

END SUBROUTINE WRGRIDA

