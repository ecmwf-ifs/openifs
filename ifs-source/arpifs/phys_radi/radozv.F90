! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE RADOZV ( YDECMIP, KIDIA , KFDIA , KLON , KLEV,&
 & KRINT , KDLON , KSHIFT,&
 & PAPRS , PGELAM, PGEMU,&
 & POZON                 )  

!***********************************************************************
! CAUTION: THIS ROUTINE WORKS ONLY ON A NON-ROTATED, UNSTRETCHED GRID
!***********************************************************************

!**** *RADOZV* - COMPUTES VARIABLE OZONE DISTRIBUTION

!     PURPOSE.
!     --------

!**   INTERFACE.
!     ----------
!        CALL *RADOZV* FROM *RADINT*

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
!     Shuting Yang  EC-Earth 9/02/2010
!     (after RADOZC, J.-J. MORCRETTE  E.C.M.W.F. 95/01/25)

!     MODIFICATIONS.
!     --------------
!     C. Roberts/R. Senan 26/01/2017 Support for CMIP6 forcings

!-----------------------------------------------------------------------

USE YOMLUN   , ONLY : NULOUT
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOECMIP  , ONLY : TECMIP, &
  &                   NLON1_CMIP5, NLAT1_CMIP5, NLV1_CMIP5, &
  &                   NLON1_CMIP6, NLAT1_CMIP6, NLV1_CMIP6

IMPLICIT NONE

TYPE(TECMIP)      ,INTENT(IN)    :: YDECMIP
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KRINT 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSHIFT 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPRS(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGELAM(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEMU(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: POZON(KDLON,KLEV) 
!     -----------------------------------------------------------------

!*       0.1   ARGUMENTS.
!              ----------

!     ----------------------------------------------------------------- 

!*       0.2   LOCAL ARRAYS.
!              -------------

REAL(KIND=JPRB), ALLOCATABLE :: ZOZLT(:,:) ! ozone on external data vertical coordinate
REAL(KIND=JPRB), ALLOCATABLE :: ZRRR(:)    ! 1/dp for external data vertical coordinate 
REAL(KIND=JPRB) :: ZWLON(KDLON), ZOZON(KDLON,KLEV+1)

INTEGER(KIND=JPIM) :: IL, INLA, JC, JK, JL, JIR
INTEGER(KIND=JPIM) :: IINLO1(KDLON), IINLO2(KDLON)

REAL(KIND=JPRB) :: ZPMR, ZSILAT, ZSIN, ZDLONR, ZLON, ZOZLT1, ZOZLT2
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

INTEGER(KIND=JPIM) :: NLON1, NLAT1 ,NLV1

#include "abor1.intfb.h"

!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------


IF (YDECMIP%NO3CMIP == 6) THEN ! CMIP6
  NLON1=NLON1_CMIP6
  NLAT1=NLAT1_CMIP6
  NLV1=NLV1_CMIP6
ELSE IF (YDECMIP%NO3CMIP == 5) THEN ! CMIP5 
  NLON1=NLON1_CMIP5
  NLAT1=NLAT1_CMIP5
  NLV1=NLV1_CMIP5
ELSE
  WRITE(NULOUT,*)"NO3CMIP:",YDECMIP%NO3CMIP
  CALL ABOR1('RADOZV: Value of NO3CMIP not supported')
ENDIF

IF (.NOT. ALLOCATED(ZOZLT))   ALLOCATE(ZOZLT(KDLON,0:NLV1+1))
IF (.NOT. ALLOCATED(ZRRR))    ALLOCATE(ZRRR(0:NLV1))


!*         1.     LATITUDE INDEX WITHIN OZONE CLIMATOLOGY
!                 ---------------------------------------

IF (LHOOK) CALL DR_HOOK('RADOZV',0,ZHOOK_HANDLE)

ZSIN=PGEMU(KIDIA)
INLA=0
ZSILAT=-9999._JPRB
DO JL=NLAT1-1,1,-1
  IF (ZSIN <= YDECMIP%RSINC1(JL+1).AND.ZSIN >= YDECMIP%RSINC1(JL)) THEN
    INLA=JL
  ENDIF
ENDDO

IF (INLA == 0) THEN
  CALL ABOR1(' Problem with lat. interpolation in radozv!')
ENDIF
ZSILAT=(ZSIN-YDECMIP%RSINC1(INLA))/(YDECMIP%RSINC1(INLA+1)-YDECMIP%RSINC1(INLA))

!
!        1a.     LONGITUDE INDEX WITHIN OZONE CLIMATOLOGY
!    ------------------------------------------------------------------
ZDLONR=YDECMIP%RLONCLI(2)-YDECMIP%RLONCLI(1)
IINLO1(:)=0
IINLO2(:)=0

IL=KSHIFT
DO JL=KIDIA,KFDIA,KRINT
  IL=IL+1
  ZLON=PGELAM(JL)
  DO JIR=1,NLON1-1
    IF ( ZLON < YDECMIP%RLONCLI(JIR+1) .AND. ZLON >= YDECMIP%RLONCLI(JIR) ) THEN
      ZWLON(IL)=(ZLON-YDECMIP%RLONCLI(JIR))/ZDLONR
      IINLO1(IL)=JIR
      IINLO2(IL)=JIR+1
    ENDIF
  ENDDO
  IF ( ZLON >= YDECMIP%RLONCLI(NLON1) ) THEN
    ZWLON(IL)=(ZLON-YDECMIP%RLONCLI(NLON1))/ZDLONR
    IINLO1(IL)=NLON1
    IINLO2(IL)=1
  ENDIF
ENDDO

!EC-EARTH: The following check is incorrect. None of the indexes maybe zero!
!IF (MAXVAL(IINLO1(:)) == 0 .OR. MAXVAL(IINLO2(:)) == 0) THEN
IF (ANY(IINLO1(KSHIFT+1:IL) == 0) .OR. ANY(IINLO2(KSHIFT+1:IL) == 0)) THEN
  CALL ABOR1(' Problem with lon. interpolation in radozv!')
ENDIF

!     ------------------------------------------------------------------

!*         2.     LATITUDE/LONGITUDE INTERPOLATED FIELD
!                 -------------------------------------

IF(INLA == NLAT1 .OR. INLA == 1) THEN
  DO JC=0,NLV1+1
    IL=KSHIFT
    DO JL=KIDIA,KFDIA,KRINT
      IL=IL+1
      ZOZLT(IL,JC)=YDECMIP%ROZT1(IINLO1(IL),INLA,JC)+ZWLON(IL)* &
        &         (YDECMIP%ROZT1(IINLO2(IL),INLA,JC)-YDECMIP%ROZT1(IINLO1(IL),INLA,JC))
    ENDDO
  ENDDO
ELSE
  DO JC=0,NLV1+1
    IL=KSHIFT
    DO JL=KIDIA,KFDIA,KRINT
      IL=IL+1
      ZOZLT1=YDECMIP%ROZT1(IINLO1(IL),INLA,JC)+ZWLON(IL)* &
        &   (YDECMIP%ROZT1(IINLO2(IL),INLA,JC)-YDECMIP%ROZT1(IINLO1(IL),INLA,JC))
      ZOZLT2=YDECMIP%ROZT1(IINLO1(IL),INLA+1,JC)+ZWLON(IL)* &
        &   (YDECMIP%ROZT1(IINLO2(IL),INLA+1,JC)-YDECMIP%ROZT1(IINLO1(IL),INLA+1,JC))
      ZOZLT(IL,JC)=ZOZLT1+ZSILAT*(ZOZLT2-ZOZLT1)
    ENDDO
  ENDDO
ENDIF

!     ------------------------------------------------------------------

!*         3.     VERTICAL INTERPOLATION 
!                 ----------------------

DO JC=0,NLV1
  ZRRR(JC)=(1.0_JPRB/(YDECMIP%RPROC1(JC)-YDECMIP%RPROC1(JC+1)))
ENDDO

IL=KSHIFT
DO JL=KIDIA,KFDIA,KRINT
  IL=IL+1
  DO JC=0,NLV1
    DO JK=1,KLEV+1
      ZPMR=PAPRS(JL,JK)
      IF(ZPMR >= YDECMIP%RPROC1(JC).AND.ZPMR < YDECMIP%RPROC1(JC+1)) &
       & ZOZON(IL,JK)=ZOZLT(IL,JC+1)+(ZPMR-YDECMIP%RPROC1(JC+1))*ZRRR(JC)* &
       &             (ZOZLT(IL,JC)-ZOZLT(IL,JC+1))
    ENDDO
  ENDDO
ENDDO

IL=KSHIFT
DO JL=KIDIA,KFDIA,KRINT
  IL=IL+1
  DO JK=1,KLEV+1
    ZPMR=PAPRS(JL,JK)
    IF(ZPMR >= YDECMIP%RPROC1(NLV1+1)) ZOZON(IL,JK)=ZOZLT(IL,NLV1+1)
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

! deallocate local arrays
DEALLOCATE(ZOZLT)
DEALLOCATE(ZRRR)

!     -----------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('RADOZV',1,ZHOOK_HANDLE)
END SUBROUTINE RADOZV
