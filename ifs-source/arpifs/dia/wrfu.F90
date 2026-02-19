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

SUBROUTINE WRFU(YDGEOMETRY,YDCFU,YDXFU,YDRIP,YDFACTX,CDFIC)

!**** * WRFU  * - ROUTINE TO EXTRACT AND WRITE CUMULATED FLUXES

!     PURPOSE.
!     --------
!        WRITES YOMGFU ON DISK

!**   INTERFACE.
!     ----------
!        *CALL* *WRFU*

!     EXPLICIT ARGUMENTS :
!     --------------------
!         KNUMER : logical unit number of the file

!     IMPLICIT ARGUMENTS :
!     --------------------
!        COMMON  YOMGFU

!     METHOD.
!     -------
!     EXTERNALS.
!     ----------

!         FAIENC    (SEE MANIPULATION DE FICHIER ARPEGE)
!         EXTGPFT
!         GRNORMA

!     REFERENCE.
!     ----------
!        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS

!     AUTHOR.
!     -------
!      PHILIPPE COURTIER
!      ORIGINAL : 91-05-23

!     MODIFICATIONS.
!     --------------
!      Modified : 01-11 by S. Ivatek-Sahdan - bf - problem in ALADIN extension
!                        zone, all fluxes from NGPTOT_CAP+1 to NGPTOT => 0
!      R. El Khatib : 02-01-17 ALADIN E-zone treatment moved to cnt4.
!      R. El Khatib : 02-21-20 Cleaning
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      M.Hamrud      10-Jan-2004 CY28R1 Cleaning
!      P. Marguinaud: 14-Jan-2011 Use WRFU_MOD
!      P. Marguinaud: 11-Sep-2012 Use IOFLDDESC_MOD + MTAG_WRFU tag
!      P. Marguinaud: 10-Oct-2013 Use FACTX
!      O. Marsden: June 2015 CY42 YRGMV, YRGFL, YRSURF, YRGMV5, and YRGFL5 are now passed by argument
!        -----------------------------------------------------------

USE YOMCFU       , ONLY : TCFU
USE YOMRIP       , ONLY : TRIP
USE YOMXFU       , ONLY : TXFU
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE IOFU_MOD  ,ONLY : IOFU_COUNT,&
                    & IOFU_SELECTD,&
                    & NIOFUACT_WRITE
USE IOGRID_MOD,ONLY : IOGRID_SELECTF,&
                    & NIOGRIDCT_WRITE
USE YOMTAG    ,ONLY : MTAG_MFIO_WRFU

USE IOFLDPTR_MOD, ONLY : IOFLDPTR

USE FACTX_MOD, ONLY : FACTX

IMPLICIT NONE

TYPE (GEOMETRY), INTENT (IN)    :: YDGEOMETRY
TYPE(TCFU)      ,INTENT(INOUT)  :: YDCFU
TYPE(TRIP)      ,INTENT(INOUT)  :: YDRIP
TYPE(TXFU)      ,INTENT(INOUT)  :: YDXFU
TYPE (FACTX),    INTENT (INOUT) :: YDFACTX
CHARACTER(LEN=*),INTENT(IN)     :: CDFIC

#include "wrgp2fa.intfb.h"

REAL (KIND=JPRB), ALLOCATABLE :: ZREAL(:,:)
TYPE (IOFLDPTR),  ALLOCATABLE :: YLFLDPT (:)
INTEGER (KIND=JPIM)           :: IFNUM

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('WRFU',0,ZHOOK_HANDLE)

IFNUM = 0
CALL IOFU_COUNT(YDGEOMETRY,YDCFU,NIOFUACT_WRITE,IFNUM)

IF (IFNUM > 0) THEN

  ALLOCATE (YLFLDPT (IFNUM))

  CALL IOFU_SELECTD(YDGEOMETRY,YDCFU, NIOFUACT_WRITE, YLFLDPT)
 
  ALLOCATE (ZREAL (YDGEOMETRY%YRGEM%NGPTOT, IFNUM))
  CALL IOGRID_SELECTF (YDGEOMETRY, NIOGRIDCT_WRITE, ZREAL, YLFLDPT)
  CALL WRGP2FA(YDGEOMETRY,YDRIP,IFNUM,ZREAL,YLFLDPT%YFLDDSC,YDFACTX,CDFIC,KTAG=MTAG_MFIO_WRFU)

  DEALLOCATE (ZREAL)
  DEALLOCATE (YLFLDPT)

ENDIF


IF (LHOOK) CALL DR_HOOK('WRFU',1,ZHOOK_HANDLE)

END SUBROUTINE WRFU

