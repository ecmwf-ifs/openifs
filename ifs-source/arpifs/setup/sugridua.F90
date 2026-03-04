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

SUBROUTINE SUGRIDUA(YDGEOMETRY,YDGFL,YDML_GCONF,YDML_LBC,KFILE)

!**** *SUGRIDUA*  - Initialize grid-point fields (GFL) from FA file

!     Purpose.
!     --------
!          To initialize grid-point fields (GFL) from FA file ; if not required in input
!          the fields are set to zero.

!     Interface.
!     ----------
!        *CALL* *SUGRIDUA

!        Explicit arguments :  
!        ------------------

!        KFILE : an indicator for which  file is to be read         

!        Implicit arguments :
!        ------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!        ABOR1

!     Reference.
!     ----------
!        Revised data flow in IFS/Arpege (By Mats Hamrud)
!        Arpege/Aladin files package (J. Clochard, R. El Khatib, D. Paradis)

!     Author.
!     -------
!      Eric Bazile *Meteo-France*
!      Original : 95-01-12

!     Modifications.
!     --------------
!      P.Smolikova 02-09-30 : Interface to INI2SCAN2M for d4 in NH
!      R. El Khatib 03-08-06 gfl.
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      Y. Seity 04-11-16 Initialisation of TKE for AROME
!      Y. Bouteloup 06-April-2006 : No read of GFL in coupling file when NCOUPLING/=1
!      R. El Khatib 10-Aug-2011 NIOLEVG management
!      P. Marguinaud 11-Sep-2012 Refactor using IOGRIDUA_MOD, IOFLDDESC_MOD and RDFA2GP
!      P. Marguinaud 10-Oct-2014 Reading coupling files using the IO server
!      O. Marsden: June 2015 CY42 YRGMV, YRGFL, YRSURF, YRGMV5, and YRGFL5 are now passed by argument
! End Modifications
!------------------------------------------------------------------------------

USE MODEL_GENERAL_CONF_MOD , ONLY : MODEL_GENERAL_CONF_TYPE
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMGFL       , ONLY : TGFL
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK
USE IOGRIDUA_MOD , ONLY : IOGRIDUA_COUNT,&
                        & IOGRIDUA_SELECTD,&
                        & NIOGRIDUACT_READ
USE IOGRID_MOD   , ONLY : IOGRID_SELECTF,&
                        & NIOGRIDCT_READ

USE IOFLDPTR_MOD , ONLY : IOFLDPTR
USE ERLBC_MOD    , ONLY : YSUGACTX
USE YEMLBC_MODEL   , ONLY : TELBC_MODEL

!------------------------------------------------------------------------------

IMPLICIT NONE

TYPE (GEOMETRY)    , INTENT (IN)    :: YDGEOMETRY
TYPE (TGFL)        , INTENT (INOUT) :: YDGFL
TYPE(MODEL_GENERAL_CONF_TYPE),INTENT(INOUT):: YDML_GCONF
TYPE (TELBC_MODEL)     , INTENT (IN)    :: YDML_LBC 
INTEGER (KIND=JPIM), INTENT (IN)    :: KFILE 
!------------------------------------------------------------------------------
TYPE (IOFLDPTR),  ALLOCATABLE :: YLFLDGT (:) 
REAL (KIND=JPRB), ALLOCATABLE :: ZGPBUFL (:,:) ! data fields to be read from file

INTEGER (KIND=JPIM) :: IFNUM
INTEGER (KIND=JPIM) :: IREP

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "rdfa2gp.intfb.h"
#include "sugridua_fixup.intfb.h"
#include "sugridua_map_part2.intfb.h"

IF (LHOOK) CALL DR_HOOK('SUGRIDUA',0,ZHOOK_HANDLE)

CALL SUGRIDUA_MAP_PART2(YDGEOMETRY,YDGFL,YDML_GCONF%YGFL,IREP,YSUGACTX)

IF (IREP == 0) THEN
  IF (LHOOK) CALL DR_HOOK('SUGRIDUA',1,ZHOOK_HANDLE)
  RETURN
ENDIF

IFNUM = 0
CALL IOGRIDUA_COUNT(YDGEOMETRY,YDGFL,YDML_GCONF%YGFL,NIOGRIDUACT_READ,IFNUM,KFILE=KFILE)

IF (IFNUM > 0) THEN

  ALLOCATE (YLFLDGT (IFNUM), ZGPBUFL (YDGEOMETRY%YRGEM%NGPTOT, IFNUM))

  CALL IOGRIDUA_SELECTD(YDGEOMETRY, YDGFL, YDML_GCONF%YGFL, NIOGRIDUACT_READ, YLFLDGT, KFILE=KFILE)
  
  CALL RDFA2GP (YDGEOMETRY, YDML_GCONF%YRRIP, IFNUM, ZGPBUFL, &
              & YLFLDGT%YFLDDSC, KFILE=KFILE, YDML_LBC=YDML_LBC)

  CALL IOGRID_SELECTF (YDGEOMETRY, NIOGRIDCT_READ, ZGPBUFL, YLFLDGT)

  DEALLOCATE (ZGPBUFL, YLFLDGT)

ENDIF

CALL SUGRIDUA_FIXUP(YDGEOMETRY,YDGFL,YDML_GCONF%YGFL,KFILE)


IF (LHOOK) CALL DR_HOOK('SUGRIDUA',1,ZHOOK_HANDLE)

END SUBROUTINE SUGRIDUA

