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

SUBROUTINE WRGRIDUA(YDGEOMETRY,YDGFL,YDXFU,YDML_GCONF,YDFACTX,CDFIC)

!**** *WRGRIDUA*  - Write out grid-point fields (GFL) to FA file

!     Purpose.
!     --------
!        To write out grid-point fields (GFL) to FA file           

!**   Interface.
!     ----------
!        *CALL* *WRGRIDUA*

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

!     Reference.
!     ----------
!        Revised data flow in IFS/Arpege (By Mats Hamrud)
!        Arpege/Aladin files package (J. Clochard, R. El Khatib, D. Paradis)

!     Author.
!     -------
!        *Meteo-France*

!     Modifications.
!     --------------
!        Original : 01-03-21 S. Ivatek-Sahdan
!        02-09-30 P. Smolikova - interface to INI2SCAN2M for d4 in NH
!        R. El Khatib 03-08-06 gfl.
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        P.Marguinaud  14-Jan-2011 Use WRGRIDUA_MOD
!        P.Marguinaud  11-Sep-2012 Use IOFLDDESC_MOD
!        P.Marguinaud  10-Oct-2013 Use FACTX
!        O. Marsden: June 2015 CY42 YRGMV, YRGFL, YRSURF, YRGMV5, and YRGFL5 are now passed by argument
!------------------------------------------------------------------------------

USE MODEL_GENERAL_CONF_MOD , ONLY : MODEL_GENERAL_CONF_TYPE
USE YOMXFU       , ONLY : TXFU
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMGFL  , ONLY : TGFL
USE PARKIND1, ONLY : JPIM, JPRB
USE YOMHOOK, ONLY : LHOOK, DR_HOOK, JPHOOK
USE IOGRIDUA_MOD, ONLY : IOGRIDUA_COUNT,&
                       & IOGRIDUA_SELECTD,&
                       & NIOGRIDUACT_WRITE
USE IOGRID_MOD,   ONLY : IOGRID_SELECTF,&
                       & NIOGRIDCT_WRITE

USE IOFLDPTR_MOD, ONLY : IOFLDPTR

USE YOMTAG, ONLY : MTAG_MFIO_WRGRIDUA

USE FACTX_MOD, ONLY : FACTX

IMPLICIT NONE

TYPE (GEOMETRY), INTENT (IN)    :: YDGEOMETRY
TYPE (TGFL) ,    INTENT (INOUT) :: YDGFL
TYPE(MODEL_GENERAL_CONF_TYPE),INTENT(INOUT):: YDML_GCONF
TYPE(TXFU)      ,INTENT(INOUT)  :: YDXFU
TYPE (FACTX),    INTENT (INOUT) :: YDFACTX
CHARACTER(LEN=*),INTENT(IN)     :: CDFIC

#include "wrgp2fa.intfb.h"


REAL (KIND=JPRB), ALLOCATABLE :: ZREAL (:,:)   ! data fields to be written to file
TYPE (IOFLDPTR),  ALLOCATABLE :: YLFLDGT (:) 
INTEGER (KIND=JPIM)           :: IFNUM

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('WRGRIDUA',0,ZHOOK_HANDLE)

IFNUM = 0
CALL IOGRIDUA_COUNT(YDGEOMETRY,YDGFL,YDML_GCONF%YGFL,NIOGRIDUACT_WRITE,IFNUM)

IF (IFNUM > 0) THEN

  ALLOCATE (YLFLDGT (IFNUM))

  CALL IOGRIDUA_SELECTD(YDGEOMETRY, YDGFL, YDML_GCONF%YGFL, NIOGRIDUACT_WRITE, YLFLDGT)

  ALLOCATE (ZREAL (YDGEOMETRY%YRGEM%NGPTOT, IFNUM))

  CALL IOGRID_SELECTF (YDGEOMETRY, NIOGRIDCT_WRITE, ZREAL, YLFLDGT)

  CALL WRGP2FA(YDGEOMETRY,YDML_GCONF%YRRIP,IFNUM,ZREAL,YLFLDGT%YFLDDSC,YDFACTX,CDFIC,KTAG=MTAG_MFIO_WRGRIDUA)

  DEALLOCATE (ZREAL)
  DEALLOCATE (YLFLDGT)

ENDIF

IF (LHOOK) CALL DR_HOOK('WRGRIDUA',1,ZHOOK_HANDLE)

END SUBROUTINE WRGRIDUA

