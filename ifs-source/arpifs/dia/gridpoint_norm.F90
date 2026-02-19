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

SUBROUTINE GRIDPOINT_NORM(YDGEOMETRY,YDGRID,CDGREP,LDLEVS)

!     Purpose
!     -------
!       Computes norm of distributed grid point fields.

!     Author
!     ------
!       Y. Tremolet

!     Modifications
!     -------------
!       16-Jan-04: Original
!       G.Mozdzynski  19-Sept-2008 gpnorm optimisations
!     ------------------------------------------------------------------

USE GEOMETRY_MOD        , ONLY : GEOMETRY
USE PARKIND1            , ONLY : JPIM, JPRB
USE YOMHOOK             , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN              , ONLY : NULOUT
USE YOMMP0              , ONLY : MYPROC
USE GRIDPOINT_FIELDS_MIX, ONLY: GRIDPOINT_FIELD

IMPLICIT NONE

TYPE(GEOMETRY)       ,INTENT(IN)    :: YDGEOMETRY
TYPE(GRIDPOINT_FIELD),INTENT(IN)    :: YDGRID
CHARACTER(LEN=*)     ,INTENT(IN)    :: CDGREP
LOGICAL, OPTIONAL, INTENT(IN) :: LDLEVS

REAL(KIND=JPRB), ALLOCATABLE :: ZGP(:,:,:), ZNORMS(:,:)
REAL(KIND=JPRB) :: ZAVE, ZMIN, ZMAX
INTEGER(KIND=JPIM) :: JB,JF,JL,JJ,JN,JOFF,JLEN,IFLDS
LOGICAL :: LLEVS
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! ------------------------------------------------------------------
#include "gpnorm1.intfb.h"
! ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GRIDPOINT_NORM',0,ZHOOK_HANDLE)

IFLDS=YDGRID%NG3D*YDGRID%NFLEVG+YDGRID%NG2D
LLEVS=.TRUE.
IF (PRESENT(LDLEVS)) LLEVS=LDLEVS

IF (IFLDS>0) THEN
  ALLOCATE(ZGP(YDGRID%NPROMA,IFLDS,YDGRID%NGPBLKS))
  ALLOCATE(ZNORMS(3,IFLDS))

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JB,JOFF,JLEN,JN,JF,JL,JJ)
  DO JB=1,YDGRID%NGPBLKS
    JOFF=(JB-1)*YDGRID%NPROMA
    JLEN=MIN(YDGRID%NPROMA,YDGRID%NGPTOT-JOFF)
    JN=0
    DO JF=1,YDGRID%NG3D
      DO JL=1,YDGRID%NFLEVG
        JN=JN+1
        DO JJ=1,JLEN
          ZGP(JJ,JN,JB)=YDGRID%GP3D(JJ,JL,JF,JB)
        ENDDO
      ENDDO
    ENDDO
    DO JF=1,YDGRID%NG2D
      JN=JN+1
      DO JJ=1,JLEN
        ZGP(JJ,JN,JB)=YDGRID%GP2D(JJ,JF,JB)
      ENDDO
    ENDDO
  ENDDO
!$OMP END PARALLEL DO

  CALL GPNORM1(YDGEOMETRY,ZGP,IFLDS,LLEVS,PNORMS=ZNORMS)

  IF (MYPROC==1) THEN
    WRITE(NULOUT,*)CDGREP,' Grid Point Norms'
    IF (LLEVS) THEN
      WRITE(NULOUT,'(5(1X,A))') ' FIELD','LEV',&
        & '       AVERAGE       ','       MINIMUM       ','       MAXIMUM       '
      JN=0
      DO JF=1,YDGRID%NG3D
        DO JL=1,YDGRID%NFLEVG
          JN=JN+1
          WRITE(NULOUT,999)YDGRID%NGRIB3(JF),JL,ZNORMS(:,JN)
        ENDDO
      ENDDO
      JL=YDGRID%NFLEVG
      DO JF=1,YDGRID%NG2D
        JN=JN+1
        WRITE(NULOUT,999)YDGRID%NGRIB2(JF),JL,ZNORMS(:,JN)
      ENDDO
999   FORMAT((1X,I6),(1X,I3),3(1X,ES21.14))
    ELSE
      WRITE(NULOUT,'(4(1X,A))') ' FIELD',&
        & '       AVERAGE       ','       MINIMUM       ','       MAXIMUM       '
      JN=0
      DO JF=1,YDGRID%NG3D
        ZAVE=SUM(ZNORMS(1,JN+1:JN+YDGRID%NFLEVG))/YDGRID%NFLEVG
        ZMIN=MINVAL(ZNORMS(2,JN+1:JN+YDGRID%NFLEVG))
        ZMAX=MAXVAL(ZNORMS(3,JN+1:JN+YDGRID%NFLEVG))
        JN=JN+YDGRID%NFLEVG
        WRITE(NULOUT,888)YDGRID%NGRIB3(JF),ZAVE,ZMIN,ZMAX
      ENDDO
      DO JF=1,YDGRID%NG2D
        JN=JN+1
        WRITE(NULOUT,888)YDGRID%NGRIB2(JF),ZNORMS(:,JN)
      ENDDO
888   FORMAT(1X,I6,3(1X,ES21.14))
    ENDIF
  ENDIF

  DEALLOCATE(ZGP,ZNORMS)
ENDIF

IF (LHOOK) CALL DR_HOOK('GRIDPOINT_NORM',1,ZHOOK_HANDLE)
END SUBROUTINE GRIDPOINT_NORM
