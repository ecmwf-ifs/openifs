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

SUBROUTINE SUSPECA_MAP_PART2(YDGEOMETRY,YGFL,YDDYN,YDDYNA,KREP,YDSUSPCTX,YDSPEC)

!**** *SUSPECA_MAP_PART2*  - Wait for spectral fields sent by the IO server

!     Author. 
!     ------- 
!      Philippe Marguinaud *METEO FRANCE* 
!      Original : 01-10-2014

USE YOMDYN       , ONLY : TDYN
USE YOM_YGFL     , ONLY : TYPE_GFLD
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPRB, JPIM
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK

USE YOMMP0       , ONLY : MYSETV
USE YOMDYNA      , ONLY : TDYNA
USE IOSPECA_MOD  , ONLY : SUSPECA_MAP_CTX,&
 &                        IOSPECA_SELECTF,&
 &                        IOSPECA_FINISH,&
 &                        NSPECACT_READ
USE YOMIO_SERV   , ONLY : IO_SERV_C001
USE SPECTRAL_FIELDS_MOD

IMPLICIT NONE

TYPE (GEOMETRY),          INTENT (IN)    :: YDGEOMETRY
TYPE(TDYN),               INTENT (IN)    :: YDDYN
TYPE(TDYNA),              INTENT(IN)     :: YDDYNA
TYPE(TYPE_GFLD),          INTENT (IN)    :: YGFL
INTEGER (KIND=JPIM),      INTENT (OUT)   :: KREP
TYPE (SUSPECA_MAP_CTX),   INTENT (INOUT) :: YDSUSPCTX
TYPE (SPECTRAL_FIELD),    INTENT (INOUT) :: YDSPEC

#include "io_serv_map_recv_part2.intfb.h"
#include "suspeca_fixup.intfb.h"

INTEGER (KIND=JPIM), PARAMETER :: IFILE = 9_JPIM
INTEGER (KIND=JPIM) :: IFNUML
REAL (KIND=JPRB), POINTER :: ZSPBUFL (:,:)

REAL (KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK ('SUSPECA_MAP_PART2',0,ZHOOK_HANDLE)

KREP = 0

CALL IO_SERV_MAP_RECV_PART2 (IO_SERV_C001, YDSUSPCTX%YIOSMPP, ZSPBUFL)

IF (.NOT. ASSOCIATED (ZSPBUFL)) THEN
  KREP = -1
  IF (LHOOK) CALL DR_HOOK ('SUSPECA_MAP_PART2',1,ZHOOK_HANDLE)
  RETURN
ENDIF

IFNUML = SIZE (ZSPBUFL, 2)

CALL IOSPECA_SELECTF(YDGEOMETRY,YGFL,YDDYN,NSPECACT_READ, YDSUSPCTX%YIOSP, IFNUML, YDSUSPCTX%IVSETOFF (MYSETV),&
                    & 0, YDSUSPCTX%YFLDSC, ZSPBUFL, YDSPEC)

DEALLOCATE (ZSPBUFL)


! Compute Vorticity, Divergence & mean wind from geographical wind

CALL IOSPECA_FINISH(YDGEOMETRY,YDDYNA,NSPECACT_READ, YDSUSPCTX%YIOSP, YDSPEC)

! Release descriptors

DEALLOCATE (YDSUSPCTX%IVSETOFF, YDSUSPCTX%YFLDSC)

CALL SUSPECA_FIXUP(YDGEOMETRY,YGFL,IFILE,YDSPEC%GFL,YDSPEC%SP)

IF (LHOOK) CALL DR_HOOK ('SUSPECA_MAP_PART2',1,ZHOOK_HANDLE)

END SUBROUTINE SUSPECA_MAP_PART2

