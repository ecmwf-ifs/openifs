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

SUBROUTINE WRGP2FA_COMPRESS_MT (YDFA, KUNIT, YDFLDSC, PGPG, PGVALCO, YDCPDSC)

!**** *WRGP2FA_COMPRESS_MT*  - Compress global grid-point fields;
!                              This routine is called both by WRGP2FA (IO server not enabled),
!                              and by the IO server itfself

!     Author. 
!     ------- 
!      Philippe Marguinaud *METEO-FRANCE*
!      Original : 01-01-2011

!      Modifications 
!      P.Marguinaud : 11-09-2012 : Refactor using IOFLDDESC_MOD
!      P.Marguinaud : 10-10-2013 : Use IOCPTDESC_MOD
!      P.Marguinaud : 10-10-2014 : Use FACONO
!      P.Marguinaud : 04-10-2016 : Port to single precision


USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK
USE PARKIND1     , ONLY : JPIM, JPRB
USE FA_MOD       , ONLY : FA_COM, JPPRCM
USE IOFLDDESC_MOD, ONLY : IOFLDDESC
USE IOCPTDESC_MOD, ONLY : IOCPTDESC

IMPLICIT NONE

TYPE (FA_COM),       INTENT (INOUT)         :: YDFA
INTEGER (KIND=JPIM), INTENT (IN)            :: KUNIT
TYPE (IOFLDDESC),    INTENT (IN)            :: YDFLDSC 
REAL (KIND=JPRB),    INTENT (IN), TARGET    :: PGPG (:)
REAL (KIND=JPRB),    INTENT (OUT)           :: PGVALCO (:)
TYPE (IOCPTDESC),    INTENT (INOUT)         :: YDCPDSC

#include "facono_mt.h"

INTEGER (KIND=JPIM) :: IREP, INGRIB, INBPDG, INBCSP, ISTRON, IPUILA, IDMOPL, ILEV
INTEGER (KIND=JPIM) :: ILONGD
INTEGER (KIND=JPIM), PARAMETER :: INGRIB0 = 0

REAL (KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK ('WRGP2FA_COMPRESS_MT',0,ZHOOK_HANDLE)

IF (YDFLDSC%LIOLV) THEN
  ILEV = YDFLDSC%IOLEV
ELSE
  ILEV = YDFLDSC%ILEVG
ENDIF

CALL FAVEUR_MT (YDFA, IREP, KUNIT, INGRIB, INBPDG, INBCSP, ISTRON, IPUILA, IDMOPL) 

IF ((YDFLDSC%NGRIBL <= 3) .AND. (YDFLDSC%JBITS == 64)) THEN
  CALL FAGOTE_MT (YDFA, IREP, KUNIT, INGRIB0, INBPDG, INBCSP, ISTRON, IPUILA, IDMOPL)
ELSE
  CALL FAGOTE_MT (YDFA, IREP, KUNIT, YDFLDSC%NGRIBL, YDFLDSC%JBITS, YDFLDSC%JBITS, ISTRON, IPUILA, IDMOPL)
ENDIF

ILONGD = SIZE (PGVALCO) / JPPRCM
CALL FACONO_MT (YDFA, IREP, KUNIT, YDFLDSC%CPREF, ILEV, YDFLDSC%CSUFF,         &
              & PGPG, .FALSE., YDCPDSC%CNOMA, YDCPDSC%ILNOMA, PGVALCO, ILONGD, &
              & LDUNDF=YDFLDSC%LUNDF, PUNDF=YDFLDSC%XUNDF)
YDCPDSC%ILONGD = ILONGD

CALL FAGOTE_MT (YDFA, IREP, KUNIT, INGRIB, INBPDG, INBCSP, ISTRON, IPUILA, IDMOPL)

IF (LHOOK) CALL DR_HOOK ('WRGP2FA_COMPRESS_MT',1,ZHOOK_HANDLE)

END SUBROUTINE WRGP2FA_COMPRESS_MT

