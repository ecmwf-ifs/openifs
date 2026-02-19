! (C) Copyright 1989- Meteo-France.

SUBROUTINE SUFPV(YDGEOMETRY,YDNAMFPV)

!**** *SUFPV*  - INITIALIZE OUTPUT VERTIVAL GEOMETRY PARAMETRES (Full-POS)

!     PURPOSE.
!     --------
!        SETS DEFAULT VALUES, THEN READS NAMELIST NAMFPG AND CHECKS IT.

!**   INTERFACE. *CALL* *SUFPV*
!     ----------

!        EXPLICIT ARGUMENTS :
!        ------------------

!        IMPLICIT ARGUMENTS
!        ------------------
!        PARDIM, YOMDIM
!        PARFPOS
!        YOMFPG, NAMFPG

!     EXTERNALS.
!     ----------

!     AUTHOR.    RYAD EL KHATIB *METEO-FRANCE*
!     -------

!     MODIFICATIONS.
!     --------------
!      Original : 02-Oct-2017 from sufpg
!-----------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1,  ONLY : JPIM     ,JPRB
USE PARFPOS ,  ONLY : JPOSDOM, JPOSGL, JPOSLE
USE YOMHOOK,   ONLY : LHOOK,   DR_HOOK, JPHOOK
USE YOMLUN,    ONLY : NULOUT   ,NULNAM   ,NULERR
USE YOMFPG,    ONLY : TNAMFPV
USE YOMVERT,   ONLY : VP00
USE YOMCT0,    ONLY : LECMWF

!-----------------------------------------------------------------------

IMPLICIT NONE

!-----------------------------------------------------------------------

TYPE(GEOMETRY), INTENT(IN) :: YDGEOMETRY
TYPE(TNAMFPV), TARGET, INTENT(OUT) :: YDNAMFPV

REAL(KIND=JPRB), POINTER :: FPVALH(:)
REAL(KIND=JPRB), POINTER :: FPVBH(:)
INTEGER(KIND=JPIM), POINTER :: NFPLEV

! Sorry Doctor - pseudo-namelist variables
INTEGER(KIND=JPIM) :: NFPMAX(JPOSDOM)
INTEGER(KIND=JPIM) :: NFPFFTW(JPOSDOM)
INTEGER(KIND=JPIM) :: NFPFLT(JPOSDOM)
INTEGER(KIND=JPIM) :: NFPRGRI(JPOSGL)
INTEGER(KIND=JPIM) :: NFPHTYP
INTEGER(KIND=JPIM) :: NFPTTYP
INTEGER(KIND=JPIM) :: NMFPMAX
REAL(KIND=JPRB) :: FPMUCEN
REAL(KIND=JPRB) :: FPLOCEN
REAL(KIND=JPRB) :: FPSTRET
REAL(KIND=JPRB) :: FPLON0
REAL(KIND=JPRB) :: FPLAT0
REAL(KIND=JPRB) :: FPRPK
REAL(KIND=JPRB) :: FPNLGINC
LOGICAL :: LFPMAP
LOGICAL :: LFPMRT
INTEGER(KIND=JPIM) :: NFPDISTRIB

INTEGER(KIND=JPIM) :: JLEV

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!-----------------------------------------------------------------------

#include "sudefo_vv1.intfb.h"
#include "posnam.intfb.h"

!-----------------------------------------------------------------------

#include "namfpg.nam.h"

!-----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUFPV',0,ZHOOK_HANDLE)
ASSOCIATE(NFLEVG=>YDGEOMETRY%YRDIMV%NFLEVG, YDVAB=>YDGEOMETRY%YRVAB)

!-----------------------------------------------------------------------

WRITE(NULOUT,'(''- Set up F-post processing, vertical geometry -'')')

!*       0. POINTERS TO THE NAMELIST STRUCTURE
!           ----------------------------------

FPVALH(0:)=>YDNAMFPV%FPVALH(0:)
FPVBH(0:)=>YDNAMFPV%FPVBH(0:)
NFPLEV=>YDNAMFPV%NFPLEV

!*       1. SET DEFAULT VALUES
!           -------------------

NFPLEV = NFLEVG
DO JLEV = 0,NFLEVG
  FPVALH(JLEV) = YDVAB%VAH(JLEV)
  FPVBH (JLEV) = YDVAB%VBH(JLEV)
ENDDO


!*       2. READ NAMELIST
!           -------------

!     Read everything
CALL POSNAM(NULNAM,'NAMFPG')
READ(NULNAM,NAMFPG)

IF (NFPLEV > JPOSLE) THEN
  NFPLEV=JPOSLE
  WRITE(NULERR,*)  ' NFPLEV TRUNCATED TO JPOSLE=',JPOSLE
ENDIF

IF (NFPLEV /= NFLEVG) THEN
!       Set a default value for FPVALH,FPVBH if possible
  CALL SUDEFO_VV1(LECMWF,NFPLEV,FPVALH,FPVBH,VP00)
!       Posible modification of FPVALH,FPVBH from namelist
  CALL POSNAM(NULNAM,'NAMFPG')
  READ(NULNAM,NAMFPG)
ENDIF

!-----------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUFPV',1,ZHOOK_HANDLE)
END SUBROUTINE SUFPV
