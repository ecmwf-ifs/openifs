! (C) Copyright 2000- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE SETUP_GEOM_MOD
CONTAINS
SUBROUTINE SETUP_GEOM

USE PARKIND1  ,ONLY : JPRD, JPIM

USE TPM_GEN         ,ONLY : NOUT, NPRINTLEV
USE TPM_DIM         ,ONLY : R
USE TPM_DISTR       ,ONLY : D
USE TPM_FIELDS      ,ONLY : F
USE TPM_GEOMETRY    ,ONLY : G
!

IMPLICIT NONE

REAL(KIND=JPRD) :: ZSQM2(R%NDGL)
INTEGER(KIND=JPIM) :: IDGLU(0:R%NSMAX,R%NDGNH)
INTEGER(KIND=JPIM) :: JGL,JM,NSMAXLIN

LOGICAL    :: LLP1,LLP2

!     ------------------------------------------------------------------

IF(.NOT.D%LGRIDONLY) THEN

  LLP1 = NPRINTLEV>0
  LLP2 = NPRINTLEV>1

  IF(LLP1) WRITE(NOUT,*) '=== ENTER ROUTINE SETUP_GEOM ==='

  ALLOCATE (G%NMEN(R%NDGL))
  IF(LLP2)WRITE(NOUT,9) 'G%NMEN   ',SIZE(G%NMEN   ),SHAPE(G%NMEN   )

  NSMAXLIN = R%NDGL-1
  IF (R%NSMAX>=NSMAXLIN .OR. .NOT. G%LREDUCED_GRID) THEN
    ! linear or full grid
    DO JGL=1,R%NDGL
      G%NMEN(JGL) = MIN(R%NSMAX,(G%NLOEN(JGL)-1)/2)
    ENDDO
  ELSEIF (R%NSMAX>=R%NDGL*2/3-1) THEN
    ! quadratic grid
    ZSQM2(:) = 3*(NSMAXLIN-R%NSMAX)/R%NDGL*F%R1MU2(:)
    G%NMEN(1) = MIN(R%NSMAX,INT(REAL(G%NLOEN(1)-1,JPRD)/(2.0_JPRD+ZSQM2(1))))
    DO JGL=2,R%NDGNH
      G%NMEN(JGL) = MIN(R%NSMAX,MAX(G%NMEN(JGL-1),&
       &INT(REAL(G%NLOEN(JGL)-1,JPRD)/(2.0_JPRD+ ZSQM2(JGL)))))
    ENDDO
    !       * SOUTHERN HEMISPHERE :
    G%NMEN(R%NDGL) = MIN(R%NSMAX,INT(REAL(G%NLOEN(R%NDGL)-1,JPRD)/(2.0_JPRD+ZSQM2(R%NDGL))))
    DO JGL=R%NDGL-1, R%NDGNH+1, -1
      G%NMEN(JGL) = MIN(R%NSMAX,MAX(G%NMEN(JGL+1),&
       &INT(REAL(G%NLOEN(JGL)-1,JPRD)/(2.0_JPRD+ ZSQM2(JGL)))))
    ENDDO
  ELSE
    ! cubic grid
    ZSQM2(:) = F%R1MU2(:)
    G%NMEN(1) = MIN(R%NSMAX,INT(REAL(G%NLOEN(1)-1,JPRD)/(2.0_JPRD+ZSQM2(1)))-1)
    DO JGL=2,R%NDGNH
      G%NMEN(JGL) = MIN(R%NSMAX,MAX(G%NMEN(JGL-1),&
       &INT(REAL(G%NLOEN(JGL)-1,JPRD)/(2.0_JPRD+ ZSQM2(JGL)))-1))
    ENDDO
    !       * SOUTHERN HEMISPHERE :
    G%NMEN(R%NDGL) = MIN(R%NSMAX,INT(REAL(G%NLOEN(R%NDGL)-1,JPRD)/(2.0_JPRD+ZSQM2(R%NDGL)))-1)
    DO JGL=R%NDGL-1, R%NDGNH+1, -1
      G%NMEN(JGL) = MIN(R%NSMAX,MAX(G%NMEN(JGL+1),&
       &INT(REAL(G%NLOEN(JGL)-1,JPRD)/(2.0_JPRD+ ZSQM2(JGL)))-1))
    ENDDO
  ENDIF
  IF(LLP1) THEN
    WRITE(NOUT,FMT='('' (JGL,G%NLOEN,G%NMEN) '')')
    WRITE(NOUT,FMT='(8(1X,''('',I4,I4,I4,'')''))')&
     &(JGL,G%NLOEN(JGL),G%NMEN(JGL),JGL=1,R%NDGL)
  ENDIF
  ALLOCATE(G%NDGLU(0:R%NSMAX))
  IF(LLP2)WRITE(NOUT,9) 'G%NDGLU   ',SIZE(G%NDGLU   ),SHAPE(G%NDGLU   )
  IDGLU(:,:) = 0
  G%NDGLU(:) = 0
  DO JGL=1,R%NDGNH
    DO JM=0,G%NMEN(JGL)
      IDGLU(JM,JGL) = 1
    ENDDO
  ENDDO
  DO JM=0,R%NSMAX
    DO JGL=1,R%NDGNH
      G%NDGLU(JM) = G%NDGLU(JM)+IDGLU(JM,JGL)
    ENDDO
  ENDDO
  IF(LLP1) THEN
      WRITE(NOUT,FMT='('' (JM,G%NDGLU) '')')
    WRITE(NOUT,FMT='(10(1X,''('',I4,I4,'')''))')&
     &(JM,G%NDGLU(JM),JM=0,R%NSMAX)
  ENDIF

ENDIF

!     ------------------------------------------------------------------
9 FORMAT(1X,'ARRAY ',A10,' ALLOCATED ',8I8)

END SUBROUTINE SETUP_GEOM
END MODULE SETUP_GEOM_MOD
