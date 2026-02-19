! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE MCHUNK

!****  *MCHUNK* 

! DECOMPOSES A IJS to IJL BLOCK INTO CHUNKS OF MAXIMUM NPROMA_WAM POINTS

! -------------------------------------------------------------------

USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

USE YOWGRID  , ONLY : IJS, IJL, NPROMA_WAM, NCHNK, KIJL4CHNK,     &
 &                    IJFROMCHNK, ICHNKFROMIJ, IPRMFROMIJ
USE YOWMPP   , ONLY : IRANK
USE YOWTEST  , ONLY : IU06

USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

!----------------------------------------------------------------------

IMPLICIT NONE

#include "abort1.intfb.h"

INTEGER(KIND=JWIM) :: IJ, IPRM, ICHNK

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('MCHUNK',0,ZHOOK_HANDLE)

! NUMBER OF CHUNKS OF MAXIMUM NPROMA_WAM POINTS:
! ---------------------------------------------
NCHNK =  (IJL-IJS+1) / NPROMA_WAM
IF( NCHNK*NPROMA_WAM <= (IJL-IJS) ) NCHNK = NCHNK + 1

! MAXIMUM GRID POINT INDEX FOR EACH CHUNCK:
! ----------------------------------------
IF (ALLOCATED(KIJL4CHNK)) DEALLOCATE(KIJL4CHNK)
ALLOCATE(KIJL4CHNK(NCHNK))
DO ICHNK = 1, NCHNK
  KIJL4CHNK(ICHNK) = MIN(NPROMA_WAM, IJL-IJS+1-(ICHNK-1)*NPROMA_WAM )
ENDDO

! EQUIVALENCE BLOCK <=> CHUNKS
! -----------------------------
IF (ALLOCATED(IJFROMCHNK)) DEALLOCATE(IJFROMCHNK)
ALLOCATE(IJFROMCHNK(NPROMA_WAM,NCHNK))
IF (ALLOCATED(ICHNKFROMIJ)) DEALLOCATE(ICHNKFROMIJ)
ALLOCATE(ICHNKFROMIJ(IJS:IJL))
IF (ALLOCATED(IPRMFROMIJ)) DEALLOCATE(IPRMFROMIJ)
ALLOCATE(IPRMFROMIJ(IJS:IJL))

DO ICHNK = 1, NCHNK
  DO IPRM = 1, KIJL4CHNK(ICHNK)
    IJ = IJS + IPRM -1 + (ICHNK-1)*NPROMA_WAM 
    IJFROMCHNK(IPRM, ICHNK) = IJ
    ICHNKFROMIJ(IJ) = ICHNK
    IPRMFROMIJ(IJ) = IPRM
  ENDDO
  IF (KIJL4CHNK(ICHNK) < NPROMA_WAM) IJFROMCHNK(IPRM+1:NPROMA_WAM, ICHNK) = 0
ENDDO

  IF(IJS /= IJFROMCHNK(1,1) .OR. IJL /= IJFROMCHNK(KIJL4CHNK(NCHNK), NCHNK) ) THEN
    WRITE(IU06,*)'* MCHUNK : SERIOUS ISSUE WITH THE MODEL DECOMPOSITION FOR THE LOCAL PTS *'
    WRITE(0,*)'*************************************************************************'
    WRITE(0,*)'* IRANK = ',IRANK
    WRITE(0,*)'* MCHUNK : SERIOUS ISSUE WITH THE MODEL DECOMPOSITION FOR THE LOCAL PTS *'
    WRITE(0,*)'* THE FOLLOWING TWO NUMBERS SHOULD BE EQUAL !!!'
    WRITE(0,*)'* IJS = ', IJS
    WRITE(0,*)'* IJFROMCHNK(1,1) = ',IJFROMCHNK(1,1)
    WRITE(0,*)'* AND OR THE FOLLOWING TWO NUMBERS SHOULD BE EQUAL !!!'
    WRITE(0,*)'* IJL = ', IJL
    WRITE(0,*)'* IJFROMCHNK(KIJL4CHNK(NCHNK), NCHNK) = ', IJFROMCHNK(KIJL4CHNK(NCHNK), NCHNK)
    WRITE(0,*)'*                                                                       *'
    WRITE(0,*)'*************************************************************************'
    CALL ABORT1
  ENDIF


IF (LHOOK) CALL DR_HOOK('MCHUNK',1,ZHOOK_HANDLE)

END SUBROUTINE MCHUNK
