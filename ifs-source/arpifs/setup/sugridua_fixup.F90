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

SUBROUTINE SUGRIDUA_FIXUP(YDGEOMETRY,YDGFL,YGFL,KFILE)

!**** *SUGRIDUA_FIXUP*  - Massage upper air fields after reading them

!     Author. 
!     ------- 
!      Philippe Marguinaud *METEO FRANCE* 
!      Original : 01-10-2014

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMGFL   , ONLY : TGFL
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOM_YGFL, ONLY : TYPE_GFLD

IMPLICIT NONE

TYPE(GEOMETRY),     INTENT(IN)    :: YDGEOMETRY
TYPE(TGFL),         INTENT(INOUT) :: YDGFL
TYPE(TYPE_GFLD)    ,INTENT(INOUT) :: YGFL
INTEGER(KIND=JPIM), INTENT(IN)    :: KFILE 
INTEGER(KIND=JPIM) :: JGFL, JLEV, JKGLO, IBL, IST, IEND, JROF
LOGICAL             :: LLREAD

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('SUGRIDUA_FIXUP',0,ZHOOK_HANDLE)

DO JKGLO = 1, YDGEOMETRY%YRGEM%NGPTOT, YDGEOMETRY%YRDIM%NPROMA

  IBL=(JKGLO-1)/YDGEOMETRY%YRDIM%NPROMA+1
  IST=1
  IEND=MIN(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRGEM%NGPTOT-JKGLO+1)

  DO JGFL=1,YGFL%NUMFLDS

    LLREAD = (YGFL%YCOMP(JGFL)%NREQIN==1).AND.((YGFL%YCOMP(JGFL)%NCOUPLING==1).OR.(KFILE/=9))

    IF (YGFL%YCOMP(JGFL)%LGP) THEN

      SELECT CASE (YGFL%YCOMP(JGFL)%NREQIN)

        CASE (1)

          IF (.NOT. LLREAD) THEN
            DO JLEV = 1, YDGEOMETRY%YRDIMV%NFLEVG
              DO JROF = IST, IEND
                YDGFL%GFL(JROF,JLEV,YGFL%YCOMP(JGFL)%MP,IBL)=0.0_JPRB
              ENDDO
            ENDDO
          ENDIF

        CASE (-1)

          DO JLEV = 1, YDGEOMETRY%YRDIMV%NFLEVG
            DO JROF = IST, IEND
              YDGFL%GFL(JROF,JLEV,YGFL%YCOMP(JGFL)%MP,IBL) = YGFL%YCOMP(JGFL)%REFVALI
            ENDDO
          ENDDO

        CASE DEFAULT

          DO JLEV = 1, YDGEOMETRY%YRDIMV%NFLEVG
            DO JROF = IST, IEND
              YDGFL%GFL(JROF,JLEV,YGFL%YCOMP(JGFL)%MP,IBL) = 0.0_JPRB
            ENDDO
          ENDDO

      END SELECT

    ENDIF

  ENDDO

ENDDO

IF (LHOOK) CALL DR_HOOK('SUGRIDUA_FIXUP',1,ZHOOK_HANDLE)

END SUBROUTINE SUGRIDUA_FIXUP

