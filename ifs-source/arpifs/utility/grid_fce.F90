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

SUBROUTINE GRID_FCE(KLAT_IN,KLON_IN,PLAT_IN,KIN,&
                  & KLAT_OUT,KLON_OUT,PLAT_OUT,KOUT,&
                  & PFLD_IN,PFLD_OUT)
!    Purpose.
!    --------
!      Wrapper of SUHIFCE for horizontal interpolation of forecast errors

!    Explicit arguments:
!    -------------------

!    Input:
!      KLAT_IN  - number of latitude rows of the input grid
!      KLON_IN  - number of longitudes for each row of the input grid
!      PLAT_IN  - latitude (radians) of each row of the input grid
!      KIN      - size of input array
!      KLAT_OUT - number of latitude rows of the output grid
!      KLON_OUT - number of longitudes for each row of the output grid
!      PLAT_OUT - latitude (radians) of each row of the output grid
!      KOUT     - size of output array
!      KGPTOT   - size of output array on each proc
!      PFLD_IN  - array of grid values at input resolution

!    Output:
!      PFLD_OUT - array of interpolated values


!    Author.
!    -------
!      S. Massart

!    Modifications.
!    --------------
!      Original: 10/11/2021
!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM, JPRB, JPRD
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCST   , ONLY : RPI

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KLAT_IN
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON_IN(KLAT_IN)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLAT_IN(KLAT_IN)
INTEGER(KIND=JPIM),INTENT(IN)    :: KIN
INTEGER(KIND=JPIM),INTENT(IN)    :: KLAT_OUT
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON_OUT(KLAT_OUT)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLAT_OUT(KLAT_OUT)
INTEGER(KIND=JPIM),INTENT(IN)    :: KOUT
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFLD_IN(KIN)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFLD_OUT(KOUT)


INTEGER(KIND=JPIM) :: IROF, JGL, JLON
REAL(KIND=JPRB)    :: ZLON0E(KLAT_IN)
REAL(KIND=JPRB)    :: ZDLONE(KLAT_IN)
REAL(KIND=JPRB)    :: ZLONM(KOUT)
REAL(KIND=JPRB)    :: ZLATM(KOUT)
REAL(KIND=JPRB)    :: ZEGRID(KIN,1)
REAL(KIND=JPRB)    :: ZIGRID(KOUT,1)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"




!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('GRID_FCE',0,ZHOOK_HANDLE)

!     *    1. INPUT GRID AND FIELD

  ZLON0E(:) = 0.0_JPRB
  ZDLONE(1:KLAT_IN) = 2.0_JPRB * RPI / KLON_IN(1:KLAT_IN)
  ZEGRID(1:KIN,1) = PFLD_IN(1:KIN)

!     *    2. OUTPUT GRID

  IROF = 0
  DO JGL = 1,KLAT_OUT
    DO JLON = 1,KLON_OUT(JGL)
      IROF=IROF+1
      ZLONM(IROF) = REAL(RPI,JPRD)*REAL(2*(JLON-1),JPRD)/REAL(KLON_OUT(JGL),JPRD)
      ZLATM(IROF) = PLAT_OUT(JGL)
    ENDDO
  ENDDO
  IF (IROF /= KOUT) CALL ABOR1('GRID_FCE: PROBLEM DIMENSION')

!     *    3. INTERPOLATION

  CALL SUHIFCE(ZLON0E,ZDLONE,PLAT_IN,KLAT_IN,KLON_IN,&
             & ZLONM,ZLATM,KOUT,&
             & KIN,1,ZEGRID,ZIGRID)
  PFLD_OUT(1:KOUT) = ZIGRID(1:KOUT,1)


IF (LHOOK) CALL DR_HOOK('GRID_FCE',1,ZHOOK_HANDLE)
!     ------------------------------------------------------------------
END SUBROUTINE GRID_FCE
