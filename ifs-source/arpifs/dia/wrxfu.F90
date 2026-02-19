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

SUBROUTINE WRXFU(YDGEOMETRY,YDXFU,YDRIP,YDFACTX,CDFIC)

USE YOMRIP       , ONLY : TRIP
USE YOMXFU       , ONLY : TXFU
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1, ONLY : JPIM, JPRB
USE YOMHOOK, ONLY : LHOOK, DR_HOOK, JPHOOK

USE IOXFU_MOD,  ONLY : IOXFU_COUNT,&
                     & IOXFU_SELECTD,&
                     & NIOXFUACT_WRITE
USE IOGRID_MOD, ONLY : IOGRID_SELECTF,&
                     & NIOGRIDCT_WRITE

USE IOFLDPTR_MOD, ONLY : IOFLDPTR

USE YOMTAG, ONLY : MTAG_MFIO_WRXFU

USE FACTX_MOD, ONLY : FACTX

!**** * WRXFU  * - ROUTINE TO EXTRAXCT AND WRITE INSTANTANEOUS FLUXES

!     PURPOSE.
!     --------
!        WRITES YOMYXFU ON DISK

!**   INTERFACE.
!     ----------
!        *CALL* *WRXFU*

!     EXPLICIT ARGUMENTS :
!     --------------------
!         KNUMER : logical unit number of the file

!     IMPLICIT ARGUMENTS :
!     --------------------
!        COMMON  YOMYXFU

!     METHOD.
!     -------
!     EXTERNALS.
!     ----------

!         FAIENC    (SEE MANIPULATION DE FICHIER ARPEGE)
!         EXTGPFT

!     REFERENCE.
!     ----------
!        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS

!     AUTHOR.
!     -------
!      O.CLARY
!      ORIGINAL : 06/92

!     MODIFICATIONS.
!     --------------
!      Modified : 01-11 by S. Ivatek-Sahdan - bf - problem in ALADIN extension
!                      zone, all fluxes from NGPTOT_CAP+1 to NGPTOT => 0
!      R. El Khatib : 02-01-17 ALADIN E-zone treatment moved to cnt4.
!      R. El Khatib : 02-21-20 Cleaning
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      M.Hamrud      10-Jan-2004 CY28R1 Cleaning
!      P. Marguinaud 14-Jan-2011 Use WRXFU_MOD
!      P. Marguinaud 11-Sep-2012 Use IOFLDDESC_MOD
!      P. Marguinaud 10-Oct-2013 Use FACTX
!        -----------------------------------------------------------

IMPLICIT NONE

TYPE (GEOMETRY), INTENT (IN)    :: YDGEOMETRY
TYPE(TRIP)      ,INTENT(INOUT)  :: YDRIP
TYPE(TXFU)      ,INTENT(INOUT)  :: YDXFU
TYPE (FACTX),    INTENT (INOUT) :: YDFACTX
CHARACTER(LEN=*),INTENT(IN)     :: CDFIC

#include "wrgp2fa.intfb.h"

REAL (KIND=JPRB), ALLOCATABLE :: ZREAL(:,:)
TYPE (IOFLDPTR),  ALLOCATABLE :: YLFLDGT (:)
INTEGER (KIND=JPIM)           :: IFNUM

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('WRXFU',0,ZHOOK_HANDLE)

IFNUM = 0
CALL IOXFU_COUNT(YDGEOMETRY,YDXFU,NIOXFUACT_WRITE,IFNUM)

IF (IFNUM > 0) THEN

  ALLOCATE (YLFLDGT (IFNUM))

  CALL IOXFU_SELECTD(YDGEOMETRY, YDXFU, NIOXFUACT_WRITE, YLFLDGT)
 
  ALLOCATE (ZREAL (YDGEOMETRY%YRGEM%NGPTOT, IFNUM))

  CALL IOGRID_SELECTF (YDGEOMETRY, NIOGRIDCT_WRITE, ZREAL, YLFLDGT)

  CALL WRGP2FA(YDGEOMETRY,YDRIP,IFNUM,ZREAL,YLFLDGT%YFLDDSC,YDFACTX,CDFIC,KTAG=MTAG_MFIO_WRXFU)

  DEALLOCATE (ZREAL)

  DEALLOCATE (YLFLDGT)

ENDIF

IF (LHOOK) CALL DR_HOOK('WRXFU',1,ZHOOK_HANDLE)

END SUBROUTINE WRXFU

