! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE RADOZC ( YDEOZOC, KIDIA , KFDIA , KLON , KLEV,&
 & KRINT , KDLON , KSHIFT,&
 & PAPRS , PGEMU,&
 & POZON                 )  

!***********************************************************************
! CAUTION: THIS ROUTINE WORKS ONLY ON A NON-ROTATED, UNSTRETCHED GRID
!***********************************************************************

!**** *RADOZC* - COMPUTES DISTRIBUTION OF OZONE FROM CLIMATOLOGY

!     PURPOSE.
!     --------

!**   INTERFACE.
!     ----------
!        CALL *RADOZC* FROM *RADINT*

!        EXPLICIT ARGUMENTS :
!        --------------------
!     ==== INPUTS ===
!     ==== OUTPUTS ===

!        IMPLICIT ARGUMENTS :   NONE
!        --------------------

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------

!          NONE

!     REFERENCE.
!     ----------

!        SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
!        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE "I.F.S"

!     AUTHOR.
!     -------
!      J.-J. MORCRETTE  E.C.M.W.F.    95/01/25

!     MODIFICATIONS.
!     --------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!-----------------------------------------------------------------------

USE PARKIND1  , ONLY : JPIM, JPRB
USE YOMHOOK   , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOEOZOC   , ONLY : TEOZOC
USE YOMDYNCORE, ONLY : LAPE

IMPLICIT NONE

TYPE(TEOZOC)      ,INTENT(IN)    :: YDEOZOC
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KRINT 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSHIFT 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPRS(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEMU(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: POZON(KDLON,KLEV) 

!     -----------------------------------------------------------------

REAL(KIND=JPRB) :: ZOZLT(0:60) , ZOZON(KDLON,KLEV+1)
REAL(KIND=JPRB) :: ZRRR(0:59)

INTEGER(KIND=JPIM) :: IL, INLA, JC, JK, JL, IMAXC
REAL(KIND=JPRB)    :: ZPMR, ZSILAT, ZSIN
REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE
LOGICAL            :: LLATINT

#include "abor1.intfb.h"

!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------

!*         1.     LATITUDE INDEX WITHIN OZONE CLIMATOLOGY
!                 ---------------------------------------

IF (LHOOK) CALL DR_HOOK('RADOZC',0,ZHOOK_HANDLE)

ZSIN=PGEMU(KIDIA)
INLA=0
ZSILAT=-9999._JPRB
IF( LAPE ) THEN
  IF(ZSIN <= YDEOZOC%RSINC(1)) THEN
    INLA=1
    LLATINT=.FALSE.
  ELSEIF(ZSIN >= YDEOZOC%RSINC(64)) THEN
    INLA=64
    LLATINT=.FALSE.
  ELSE
    DO JL=63,1,-1
      IF (ZSIN <= YDEOZOC%RSINC(JL+1).AND.ZSIN >= YDEOZOC%RSINC(JL)) THEN
        INLA=JL
      ENDIF
    ENDDO
    ZSILAT=(ZSIN-YDEOZOC%RSINC(INLA))/(YDEOZOC%RSINC(INLA+1)-YDEOZOC%RSINC(INLA))
    LLATINT=.TRUE.
  ENDIF
ELSE
  DO JL=18,1,-1
    IF (ZSIN <= YDEOZOC%RSINC(JL+1).AND.ZSIN >= YDEOZOC%RSINC(JL)) THEN
      INLA=JL
    ENDIF
  ENDDO
ENDIF
IF (INLA == 0) THEN
  CALL ABOR1(' Problem with lat. interpolation in radozc!')
ENDIF

!     ------------------------------------------------------------------

!*         2.     LATITUDE INTERPOLATED FIELD
!                 ---------------------------

IF( LAPE ) THEN
  IMAXC=60
  IF(.NOT.LLATINT) THEN
    DO JC=0,IMAXC
      ZOZLT(JC)=YDEOZOC%ROZT(INLA,JC)
    ENDDO
  ELSE
    DO JC=0,IMAXC
      ZOZLT(JC)=YDEOZOC%ROZT(INLA,JC)+ZSILAT*(YDEOZOC%ROZT(INLA+1,JC)-YDEOZOC%ROZT(INLA,JC))
    ENDDO
  ENDIF
ELSE
  IMAXC=35
  ZSILAT=(ZSIN-YDEOZOC%RSINC(INLA))/(YDEOZOC%RSINC(INLA+1)-YDEOZOC%RSINC(INLA))
  IF(INLA == 18.OR.INLA == 1) THEN
    DO JC=0,IMAXC
      ZOZLT(JC)=YDEOZOC%ROZT(INLA,JC)
    ENDDO
  ELSE
    DO JC=0,IMAXC
      ZOZLT(JC)=YDEOZOC%ROZT(INLA,JC)+ZSILAT*(YDEOZOC%ROZT(INLA+1,JC)-YDEOZOC%ROZT(INLA,JC))
    ENDDO
  ENDIF
ENDIF

!     ------------------------------------------------------------------

!*         3.     VERTICAL INTERPOLATION 
!                 ----------------------

DO JC=0,IMAXC-1
  ZRRR(JC)=(1.0_JPRB/(YDEOZOC%RPROC(JC)-YDEOZOC%RPROC(JC+1)))*(ZOZLT(JC)-ZOZLT(JC+1))
ENDDO

IL=KSHIFT
DO JL=KIDIA,KFDIA,KRINT
  IL=IL+1
  DO JC=0,IMAXC-1
    DO JK=1,KLEV+1
      ZPMR=PAPRS(JL,JK)
      IF(ZPMR >= YDEOZOC%RPROC(JC).AND.ZPMR < YDEOZOC%RPROC(JC+1)) &
       & ZOZON(IL,JK)=ZOZLT(JC+1)+(ZPMR-YDEOZOC%RPROC(JC+1))*ZRRR(JC)  
    ENDDO
  ENDDO
ENDDO

IL=KSHIFT
DO JL=KIDIA,KFDIA,KRINT
  IL=IL+1
  DO JK=1,KLEV+1
    ZPMR=PAPRS(JL,JK)
    ZPMR=PAPRS(JL,JK)
    IF(ZPMR >= YDEOZOC%RPROC(IMAXC)) ZOZON(IL,JK)=ZOZLT(IMAXC)
  ENDDO
ENDDO

! INTEGRATION IN THE VERTICAL:
IL=KSHIFT
DO JL=KIDIA,KFDIA,KRINT
  IL=IL+1
  DO JK=1,KLEV
    POZON(IL,JK)=(PAPRS(JL,JK+1)-PAPRS(JL,JK))&
     & *(ZOZON(IL,JK)+ZOZON(IL,JK+1))*0.5_JPRB  
  ENDDO
ENDDO

!     -----------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('RADOZC',1,ZHOOK_HANDLE)
END SUBROUTINE RADOZC
