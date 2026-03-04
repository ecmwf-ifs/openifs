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

SUBROUTINE WRSPECA_FP(YDRIP,YDFPS,YDGEOMETRY)

!**** *WRSPECA_FP*  - Send spectral fields to FP server

!     Author. 
!     ------- 
!      Philippe Marguinaud *METEO-FRANCE*
!      Original : 15-04-2016


USE YOMRIP   , ONLY : TRIP
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK

USE YOMFP_SERV,   ONLY : FP_SERV
USE YOMTAG,       ONLY : MTAG_MFIO_WRSPECSPA
USE GEOMETRY_MOD, ONLY : GEOMETRY
USE YOMMP0,       ONLY : MYSETV

USE IOFLDPTR_MOD, ONLY : IOFLDPTR

USE IOSPECSPA_MOD, ONLY : NIOSPECSPACT_WRITE, &
                        & IOSPECSPA_COUNT,    &
                        & IOSPECSPA_SELECTD

IMPLICIT NONE

TYPE(TRIP)      ,INTENT(INOUT)  :: YDRIP
TYPE (FP_SERV),  INTENT (INOUT) :: YDFPS
TYPE (GEOMETRY), INTENT (IN)    :: YDGEOMETRY

#include "wrfld_fp.intfb.h"

INTEGER(KIND=JPIM) :: IFLDSPG

TYPE (IOFLDPTR), ALLOCATABLE :: YLFLDSTG (:)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('WRSPECA_FP',0,ZHOOK_HANDLE)

IFLDSPG = 0
CALL IOSPECSPA_COUNT (YDGEOMETRY, NIOSPECSPACT_WRITE, IFLDSPG)

ALLOCATE (YLFLDSTG (IFLDSPG))

CALL IOSPECSPA_SELECTD (YDGEOMETRY, NIOSPECSPACT_WRITE, YLFLDSTG)

CALL WRFLD_FP(YDFPS,YDGEOMETRY,YDRIP%TSTEP,MTAG_MFIO_WRSPECSPA,MYSETV,YLFLDSTG,YDFPS%YIOS_SP_CL,YDFPS%YSP_DISTR)

DEALLOCATE (YLFLDSTG)

IF (LHOOK) CALL DR_HOOK('WRSPECA_FP',1,ZHOOK_HANDLE)

END SUBROUTINE WRSPECA_FP

