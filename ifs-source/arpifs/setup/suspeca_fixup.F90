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

SUBROUTINE SUSPECA_FIXUP(YDGEOMETRY,YGFL,KFILE,PSPGFL,PSPSP)

!**** *SUSPECA_FIXUP*  - Massage spectral fields after reading them

!     Author. 
!     ------- 
!      Philippe Marguinaud *METEO FRANCE* 
!      Original : 01-10-2014
!      R. El Khatib abort if Ps is not Ln(Ps) instead of converting it


USE GEOMETRY_MOD, ONLY : GEOMETRY
USE PARKIND1    , ONLY : JPIM, JPRB
USE YOMHOOK     , ONLY : LHOOK, DR_HOOK, JPHOOK

USE YOM_YGFL    , ONLY : TYPE_GFLD
USE YOMCT0      , ONLY : N3DINI
USE YOMMP0      , ONLY : MYSETW, MYPROC
IMPLICIT NONE

TYPE (GEOMETRY),      INTENT (IN)    :: YDGEOMETRY
TYPE (TYPE_GFLD),     INTENT (IN)    :: YGFL
INTEGER (KIND=JPIM),  INTENT (IN)    :: KFILE
REAL (KIND=JPRB),     INTENT (OUT)   :: PSPGFL(:,:,:) 
REAL (KIND=JPRB),     INTENT (INOUT) :: PSPSP(:) 

#include "abor1.intfb.h"
#include "set2pe.intfb.h"

INTEGER (KIND=JPIM) :: JGFL, JLEV, JSP
INTEGER (KIND=JPIM) :: I3DINI
INTEGER(KIND=JPIM) :: ISETW0, IPRSPW0 ! Set-W processor and processor owning wavenumber 0

REAL (KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('SUSPECA_FIXUP',0,ZHOOK_HANDLE)

I3DINI=N3DINI

  ! special treatment of Spectral GFL not read in initial file
DO JGFL=1, YGFL%NUMFLDS
  IF (YGFL%YCOMP(JGFL)%LSP) THEN
    IF ( YGFL%YCOMP(JGFL)%NREQIN== -1) THEN
      DO JLEV=1,YDGEOMETRY%YRDIMV%NFLEVL
        DO JSP=2,YDGEOMETRY%YRDIM%NSPEC2
          PSPGFL(JLEV,JSP,YGFL%YCOMP(JGFL)%MPSP) = 0._JPRB
        ENDDO  
        IF (MYSETW == 1) THEN
          PSPGFL(JLEV,1,YGFL%YCOMP(JGFL)%MPSP)=YGFL%YCOMP(JGFL)%REFVALI
        ELSE
          PSPGFL(JLEV,1,YGFL%YCOMP(JGFL)%MPSP)=0._JPRB
        ENDIF
      ENDDO  
    ENDIF  
  ENDIF
ENDDO

DO JGFL=1, YGFL%NUMFLDS
  IF (YGFL%YCOMP(JGFL)%LSP) THEN
    IF ( YGFL%YCOMP(JGFL)%NREQIN== 0) THEN
      DO JLEV=1,YDGEOMETRY%YRDIMV%NFLEVL
        DO JSP=1,YDGEOMETRY%YRDIM%NSPEC2
          PSPGFL(JLEV,JSP,YGFL%YCOMP(JGFL)%MPSP) = 0._JPRB
        ENDDO 
      ENDDO  
    ENDIF  
  ENDIF
ENDDO

  ! Ensure the pronostic variable is Ln(Ps)

! unless KFILE=3 (file contains increments)
IF (KFILE /= 3) THEN
  IF (I3DINI == 0) THEN
    ISETW0=1
    CALL SET2PE(IPRSPW0,0,0,ISETW0,YDGEOMETRY%YRMP%NBSETSP)
    IF (MYPROC == IPRSPW0) THEN
      IF (PSPSP(1) > 5.E+4_JPRB) THEN
        CALL ABOR1('SUSPECA_FIXUP: THIS IS APPARENTLY Ps, NOT LN(Ps) !')
      ENDIF
    ENDIF
  ENDIF
ENDIF

IF (LHOOK) CALL DR_HOOK('SUSPECA_FIXUP',1,ZHOOK_HANDLE)

END SUBROUTINE SUSPECA_FIXUP
