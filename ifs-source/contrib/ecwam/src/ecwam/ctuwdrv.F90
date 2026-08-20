! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE CTUWDRV (DELPRO, MSTART, MEND,                  &
 &                  IJS, IJL, NINF, NSUP,                  &
 &                  BLK2GLO,                               &
 &                  WLATM1, WCORM1, DP,                    &
 &                  CGROUP_EXT, OMOSNH2KD_EXT,             &
 &                  COSPHM1_EXT, DEPTH_EXT, U_EXT, V_EXT )


! ----------------------------------------------------------------------

!**** *CTUWDRV*

!*    PURPOSE.
!     --------

!     CALLS CTUW FOR FREQUENCIES FR(MSTART) TO FR(MIND)

! -------------------------------------------------------------------

USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
USE YOWDRVTYPE  , ONLY : WVGRIDGLO

USE YOWCURR  , ONLY : LLCFLCUROFF
USE YOWGRID  , ONLY : NPROMA_WAM
USE YOWMPP   , ONLY : IRANK
USE YOWPARAM , ONLY : NANG     ,NFRE_RED
USE YOWSTAT  , ONLY : IREFRA
USE YOWTEST  , ONLY : IU06

USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
USE EC_LUN   , ONLY : NULERR
      
! ----------------------------------------------------------------------
      IMPLICIT NONE

#include "abort1.intfb.h"
#include "ctuw.intfb.h"

REAL(KIND=JWRB), INTENT(IN) :: DELPRO ! ADVECTION TIME STEP
INTEGER(KIND=JWIM), INTENT(IN) :: MSTART, MEND ! FREQUENCY START AND END INDEXES OVER WHICH THE WEIGHT WILL BE COMPUTED
INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL   ! GRID POINTS WITHIN A BLOCK
INTEGER(KIND=JWIM), INTENT(IN) :: NINF, NSUP ! GRID POINT WITH HALO EXTEND NINF:NSUP+1 
TYPE(WVGRIDGLO), INTENT(IN) :: BLK2GLO  ! BLOCK TO GRID TRANSFORMATION
REAL(KIND=JWRB), DIMENSION(NINF:NSUP,2), INTENT(IN) :: WLATM1  ! 1 - WLAT
REAL(KIND=JWRB), DIMENSION(NINF:NSUP,4), INTENT(IN) :: WCORM1  ! 1 - WCOR
REAL(KIND=JWRB), DIMENSION(NINF:NSUP,2), INTENT(IN) :: DP      ! COS PHI FACTOR
REAL(KIND=JWRB), DIMENSION(NINF:NSUP+1, NFRE_RED), INTENT(IN) :: CGROUP_EXT  ! GROUP VELOCITY
REAL(KIND=JWRB), DIMENSION(NINF:NSUP+1, NFRE_RED), INTENT(IN) :: OMOSNH2KD_EXT ! OMEGA / SINH(2KD)
REAL(KIND=JWRB), DIMENSION(NINF:NSUP+1), INTENT(IN) :: COSPHM1_EXT ! 1/COSPH
REAL(KIND=JWRB), DIMENSION(NINF:NSUP+1), INTENT(IN) :: DEPTH_EXT ! WATER DEPTH
REAL(KIND=JWRB), DIMENSION(NINF:NSUP+1), INTENT(IN) :: U_EXT ! U-COMPONENT OF SURFACE CURRENT
REAL(KIND=JWRB), DIMENSION(NINF:NSUP+1), INTENT(IN) :: V_EXT ! V-COMPONENT OF SURFACE CURRENT


INTEGER(KIND=JWIM) :: IJ, ICALL
INTEGER(KIND=JWIM) :: JKGLO, KIJS, KIJL, NPROMA, MTHREADS
!$    INTEGER,EXTERNAL :: OMP_GET_MAX_THREADS

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

LOGICAL :: LL2NDCALL
LOGICAL, DIMENSION(IJS:IJL) :: LCFLFAIL

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CTUWDRV',0,ZHOOK_HANDLE)

!! NPROMA=NPROMA_WAM
   MTHREADS=1
!$ MTHREADS=OMP_GET_MAX_THREADS()
   NPROMA=(IJL-IJS+1)/MTHREADS + 1

!$OMP   PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO, KIJS, KIJL, ICALL, IJ, LL2NDCALL)
DO JKGLO = IJS, IJL, NPROMA
  KIJS=JKGLO
  KIJL=MIN(KIJS+NPROMA-1,IJL)
  ICALL = 1
  CALL CTUW(DELPRO, MSTART, MEND,                          &
&           KIJS, KIJL, NINF, NSUP, LCFLFAIL(KIJS), ICALL, &
&           BLK2GLO,                                       &
&           WLATM1, WCORM1, DP,                            &
&           CGROUP_EXT, OMOSNH2KD_EXT,                     &
&           COSPHM1_EXT, DEPTH_EXT, U_EXT, V_EXT )


! WHEN SURFACE CURRENTS ARE USED AND LLCFLCUROFF IS TRUE
! THEN TRY TO SATISFY THE CFL CONDITION WITHOUT THE CURRENTS
! IF IT WAS VIOLATED IN THE FIRST PLACE
  IF (LLCFLCUROFF .AND. (IREFRA == 2 .OR. IREFRA == 3)) THEN
    LL2NDCALL=.FALSE.
    DO IJ=KIJS,KIJL
       IF (LCFLFAIL(IJ)) THEN
         LL2NDCALL=.TRUE.
         EXIT
       ENDIF
    ENDDO
    IF (LL2NDCALL) THEN
       ICALL = 2
       CALL CTUW(DELPRO, MSTART, MEND,                          &
&                KIJS, KIJL, NINF, NSUP, LCFLFAIL(KIJS), ICALL, &
&                BLK2GLO,                                       &
&                WLATM1, WCORM1, DP,                            &
&                CGROUP_EXT, OMOSNH2KD_EXT,                     &
&                COSPHM1_EXT, DEPTH_EXT, U_EXT, V_EXT )
    ENDIF
  ENDIF
ENDDO
!$OMP  END PARALLEL DO

DO IJ=IJS,IJL
  IF (LCFLFAIL(IJ)) THEN
    CALL FLUSH (IU06)
    WRITE(NULERR,*) '!!! ********************************* !!'
    WRITE(NULERR,*) '!!! WAVE MODEL HAS ABORTED !!!'
    WRITE(NULERR,*) '!!! FOLLOWING CFL CRITERION VIOLATION !!'
    WRITE(NULERR,*) '!!! FOR DELPRO = ', DELPRO
    WRITE(NULERR,*) '!!! ON PE ',IRANK
    WRITE(NULERR,*) '!!! ********************************* !!'
    WRITE(IU06,*) '!!! ********************************* !!'
    WRITE(IU06,*) '!!! WAVE MODEL HAS ABORTED !!!'
    WRITE(IU06,*) '!!! FOLLOWING CFL CRITERION VIOLATION !!'
    WRITE(IU06,*) '!!! FOR DELPRO = ', DELPRO
    WRITE(IU06,*) '!!! ON PE ',IRANK
    WRITE(IU06,*) '!!! ********************************* !!'
    CALL ABORT1
  ENDIF
ENDDO


IF (LHOOK) CALL DR_HOOK('CTUWDRV',1,ZHOOK_HANDLE)

END SUBROUTINE CTUWDRV
