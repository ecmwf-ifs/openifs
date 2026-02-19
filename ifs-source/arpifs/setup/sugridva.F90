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

SUBROUTINE SUGRIDVA(YDGEOMETRY,YDSURF,YDRIP,YDML_LBC)

USE YOMRIP             , ONLY : TRIP
USE GEOMETRY_MOD       , ONLY : GEOMETRY
USE SURFACE_FIELDS_MIX , ONLY : TSURF
USE PARKIND1           , ONLY : JPIM, JPRB
USE YOMHOOK            , ONLY : LHOOK, DR_HOOK, JPHOOK

USE IOGRIDVA_MOD       , ONLY : NIOGRIDVACT_READ,&
 &                              IOGRIDVA_COUNT,&
 &                              IOGRIDVA_SELECTD,&
 &                              IOGRIDVA_SELECTF
USE YEMLBC_MODEL         , ONLY : TELBC_MODEL
USE IOFLDDESC_MOD      , ONLY : IOFLDDESC


!**** *SUGRIDVA*  - Initialize the clim. gridpoint fields from *FA*

!     Purpose.
!     --------
!           Initialize the clim. gridpoint fields from *FA*

!**   Interface.
!     ----------
!        *CALL* *SUGRIDVA(...)*

!        Explicit arguments :
!        ------------------

!        Implicit arguments :
!        --------------------
!        See 'USE MODULE' above.

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!        OPENFA
!        RDGPFA
!        DISGRID
!        FAIRME
!        Is called by SUGRIDF

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS
!        Note de travail ARPEGE Nr 12 et 17

!     Author.
!     -------
!        RYAD EL KHATIB *METEO-FRANCE*
!        from SUGRCLIADM/SUGRIDVA 

!     Modifications.
!     --------------
!        ORIGINAL : 97-12-03
!        Y. Bouteloup : 02-09-11 xxO3ABC => Climatological ozone profiles
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        P. Marquet   : 02-07-17 add VCLIA for aerosol files
!        M.Hamrud      01-Jul-2006 Revised surface fields
!        R. El Khatib  26-Nov-2008  Cleaning
!        R. El Khatib : 23-Apr-2010 use disgrid_mod instead of disgrid
!        P. Marguinaud: 11-Sep-2012 Refactor using IOGRIDVA_MOD, IOFLDDESC_MOD and RDFA2GP
!     ------------------------------------------------------------------

IMPLICIT NONE

#include "rdfa2gp.intfb.h"

TYPE(GEOMETRY) , INTENT(IN)    :: YDGEOMETRY
TYPE(TSURF)    , INTENT(INOUT) :: YDSURF
TYPE(TRIP)     , INTENT(INOUT) :: YDRIP
TYPE(TELBC_MODEL)  , INTENT(IN)    :: YDML_LBC

TYPE(IOFLDDESC)  , ALLOCATABLE   :: YLFLDSC (:) 
REAL(KIND=JPRB)  , ALLOCATABLE   :: ZBUFL (:,:)

INTEGER (KIND=JPIM) :: IFNUM

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
!*       1.    PREPARATIONS.
!              ------------

IF (LHOOK) CALL DR_HOOK('SUGRIDVA',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP)
ASSOCIATE(NGPTOT=>YDGEM%NGPTOT, NGPTOTG=>YDGEM%NGPTOTG)
IFNUM = 0 
CALL IOGRIDVA_COUNT(YDGEOMETRY,YDSURF,NIOGRIDVACT_READ, IFNUM)

IF (IFNUM > 0) THEN

  ALLOCATE (YLFLDSC (IFNUM), ZBUFL (NGPTOT, IFNUM))

  CALL IOGRIDVA_SELECTD(YDGEOMETRY,YDSURF,NIOGRIDVACT_READ, YLFLDSC)

! read & distribute in local buffers
  CALL RDFA2GP (YDGEOMETRY, YDRIP, IFNUM, ZBUFL, YLFLDSC, KFILE=19_JPIM, YDML_LBC=YDML_LBC)
    
  

! fill model variables
  CALL IOGRIDVA_SELECTF(YDGEOMETRY,YDSURF,NIOGRIDVACT_READ, ZBUFL, YLFLDSC)
  
  DEALLOCATE (ZBUFL, YLFLDSC)

ENDIF 

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUGRIDVA',1,ZHOOK_HANDLE)
END SUBROUTINE SUGRIDVA
