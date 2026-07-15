! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

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
USE EC_LUN   , ONLY : NULERR

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
  DO IPRM = KIJL4CHNK(ICHNK)+1, NPROMA_WAM
    IJFROMCHNK(IPRM, ICHNK) = 0
  ENDDO
ENDDO

  IF(IJS /= IJFROMCHNK(1,1) .OR. IJL /= IJFROMCHNK(KIJL4CHNK(NCHNK), NCHNK) ) THEN
    WRITE(IU06,*)'* MCHUNK : SERIOUS ISSUE WITH THE MODEL DECOMPOSITION FOR THE LOCAL PTS *'
    WRITE(NULERR,*)'*************************************************************************'
    WRITE(NULERR,*)'* IRANK = ',IRANK
    WRITE(NULERR,*)'* MCHUNK : SERIOUS ISSUE WITH THE MODEL DECOMPOSITION FOR THE LOCAL PTS *'
    WRITE(NULERR,*)'* THE FOLLOWING TWO NUMBERS SHOULD BE EQUAL !!!'
    WRITE(NULERR,*)'* IJS = ', IJS
    WRITE(NULERR,*)'* IJFROMCHNK(1,1) = ',IJFROMCHNK(1,1)
    WRITE(NULERR,*)'* AND OR THE FOLLOWING TWO NUMBERS SHOULD BE EQUAL !!!'
    WRITE(NULERR,*)'* IJL = ', IJL
    WRITE(NULERR,*)'* IJFROMCHNK(KIJL4CHNK(NCHNK), NCHNK) = ', IJFROMCHNK(KIJL4CHNK(NCHNK), NCHNK)
    WRITE(NULERR,*)'*                                                                       *'
    WRITE(NULERR,*)'*************************************************************************'
    CALL ABORT1
  ENDIF


IF (LHOOK) CALL DR_HOOK('MCHUNK',1,ZHOOK_HANDLE)

END SUBROUTINE MCHUNK
