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

SUBROUTINE SUVERTFE1(YDGEOMETRY)

!**** *SUVERTFE1*  - Setup VERTical Finite Element scheme based on
!                    linear (1) elements

!     Purpose.
!     --------
!           Compute matrix for the vertical integrations based on
!           linear finite elements

!**   Interface.
!     ----------

!     *CALL* SUVERTFE1
!        Explicit arguments :
!        --------------------
!        None

!        Implicit arguments :
!        --------------------

!     Method.
!     -------

!     Externals.
!     ----------
!        MINV_8

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Mariano Hortal ECMWF

!     Modifications.
!     --------------
!        Original : 2000-06
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        K. Yessad (Sep 2008): use VFE_ETAF
!        T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!        F. Vana     17-Dec-2015 Support for single precision
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1 , ONLY : JPIM, JPRB, JPRD
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN   , ONLY : NULOUT   
USE YOMMP0   , ONLY : MYPROC

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY), INTENT(INOUT) :: YDGEOMETRY
REAL(KIND=JPRD) :: ZWORK(2*(YDGEOMETRY%YRDIMV%NFLEVG+1))
REAL(KIND=JPRD) :: ZAMAT(YDGEOMETRY%YRDIMV%NFLEVG+1,YDGEOMETRY%YRDIMV%NFLEVG+1)
REAL(KIND=JPRD) :: ZAMATI(YDGEOMETRY%YRDIMV%NFLEVG+1,YDGEOMETRY%YRDIMV%NFLEVG+1)
REAL(KIND=JPRD) :: ZBMAT(YDGEOMETRY%YRDIMV%NFLEVG+1,YDGEOMETRY%YRDIMV%NFLEVG+1)
REAL(KIND=JPRD) :: ZPROD(YDGEOMETRY%YRDIMV%NFLEVG+1,YDGEOMETRY%YRDIMV%NFLEVG+1)
REAL(KIND=JPRD) :: ZDERI(YDGEOMETRY%YRDIMV%NFLEVG+1,YDGEOMETRY%YRDIMV%NFLEVG+1)
REAL(KIND=JPRD) :: ZD(0:YDGEOMETRY%YRDIMV%NFLEVG+1)

INTEGER(KIND=JPIM) :: I, JI, JJ, JK, JLEV

REAL(KIND=JPRD) :: ZD0, ZDET, ZEPS

LOGICAL :: LLPRINT
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "minv_8.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUVERTFE1',0,ZHOOK_HANDLE)
ASSOCIATE(YDVAB=>YDGEOMETRY%YRVAB,YDVETA=>YDGEOMETRY%YRVETA,YDVFE=>YDGEOMETRY%YRVFE,YDSTA=>YDGEOMETRY%YRSTA, &
 & YDLAP=>YDGEOMETRY%YRLAP,YDCSGLEG=>YDGEOMETRY%YRCSGLEG, &
 & YDCSGEOM=>YDGEOMETRY%YRCSGEOM,YDCSGEOM_NB=>YDGEOMETRY%YRCSGEOM_NB,YDGSGEOM=>YDGEOMETRY%YRGSGEOM, &
 & YDGSGEOM_NB=>YDGEOMETRY%YRGSGEOM_NB, YDSPGEOM=>YDGEOMETRY%YSPGEOM)
ASSOCIATE(NFLEVG=>YDGEOMETRY%YRDIMV%NFLEVG, &
  & YDDIMV=>YDGEOMETRY%YRDIMV)!     ------------------------------------------------------------------

LLPRINT=.FALSE.

! Scheme No 1: integral by NLEVG+1 finite elements with mass matrix
!    DO JLEV=0,NFLEVG
!      VFE_ETAH(JLEV)=REAL(JLEV,JPRD)/NFLEVG
!      IF(JLEV > 0) YRVETA%VFE_ETAF(JLEV)=_HALF_*(VFE_ETAH(JLEV)+VFE_ETAH(JLEV-1))
!    ENDDO
DO JI=1,NFLEVG+1
  IF(JI == 1) THEN
    ZD(JI)=YDVETA%VFE_ETAF(1)
  ELSEIF(JI == NFLEVG+1) THEN
    ZD(JI)=(1.0_JPRD-YDVETA%VFE_ETAF(NFLEVG))
  ELSE
    ZD(JI)=YDVETA%VFE_ETAF(JI)-YDVETA%VFE_ETAF(JI-1)
  ENDIF
ENDDO
ZD(0)=0.0_JPRD
! #### first index is row number, second is column number
ZAMAT(:,:)=0.0_JPRD
DO JJ=1,NFLEVG+1
  IF(JJ > 1) ZAMAT(JJ,JJ-1)=ZD(JJ)/6._JPRD
  IF(JJ < NFLEVG+1) ZAMAT(JJ,JJ)=(ZD(JJ)+ZD(JJ+1))/3._JPRD
  IF(JJ < NFLEVG+1) ZAMAT(JJ,JJ+1)=ZD(JJ+1)/6._JPRD
ENDDO
ZAMAT(NFLEVG+1,NFLEVG+1)=ZD(NFLEVG+1)/3._JPRD

ZBMAT(:,:)=0.0_JPRD
DO JI=1,NFLEVG+1
  DO JJ=1,NFLEVG+1
    IF(JI == JJ+2) THEN
      ZBMAT(JJ,JI)=ZD(JJ+1)**2/24._JPRD
    ELSEIF(JI == JJ+1) THEN
      ZBMAT(JJ,JI)=ZD(JJ)**2/8._JPRD + ZD(JJ)*ZD(JJ+1)/4._JPRD + &
       & ZD(JJ+1)**2/8._JPRD  
    ELSEIF(JI == JJ) THEN
      ZBMAT(JJ,JI)=ZD(JJ-1)*ZD(JJ)/4._JPRD + &
       & 5._JPRD*ZD(JJ)**2/24._JPRD + &
       & ZD(JJ+1)*(ZD(JJ-1)+ZD(JJ))/4._JPRD  
    ELSEIF(JI < JJ) THEN
      ZBMAT(JJ,JI)=(ZD(JI-1)+ZD(JI))/4._JPRD*(ZD(JJ)+ZD(JJ+1))
    ENDIF
  ENDDO
ENDDO
ZBMAT(1,1)=5._JPRD*ZD(1)**2/24._JPRD + ZD(1)*ZD(2)/4._JPRD
DO JJ=2,NFLEVG
  ZBMAT(JJ,1)=ZD(1)*(ZD(JJ)+ZD(JJ+1))/4._JPRD
ENDDO
ZBMAT(NFLEVG+1,1)=ZD(1)*ZD(NFLEVG+1)/4._JPRD
DO JI=2,NFLEVG
  ZBMAT(NFLEVG+1,JI)=(ZD(JI-1)+ZD(JI))*ZD(NFLEVG+1)/4._JPRD
ENDDO
ZBMAT(NFLEVG+1,NFLEVG+1)=ZD(NFLEVG)*ZD(NFLEVG+1)/4._JPRD + &
 & ZD(NFLEVG+1)**2/3._JPRD  
ZBMAT(NFLEVG,NFLEVG+1)=ZD(NFLEVG)**2/8._JPRD + &
 & ZD(NFLEVG)*ZD(NFLEVG+1)/4._JPRD + ZD(NFLEVG+1)**2/6._JPRD  
ZDET=0.0_JPRD
ZEPS=100.0_JPRD*TINY(1.0_JPRD)
ZAMATI(:,:)=ZAMAT(:,:)*NFLEVG
CALL MINV_8(ZAMATI,NFLEVG+1,NFLEVG+1,ZWORK,ZDET,ZEPS,0,1)
ZAMATI(:,:)=ZAMATI(:,:)*NFLEVG
IF(MYPROC == 1.AND.LLPRINT) THEN
  WRITE(NULOUT,*) ' DETERMINANT OF A*NFLEVG =',ZDET
  WRITE(NULOUT,*) ' MATRIX A*100 :'
  DO JLEV=1,NFLEVG+1
    WRITE(NULOUT,'(100F6.1)') (100._JPRD*ZAMAT(JLEV,I),I=1,NFLEVG+1)
  ENDDO
  WRITE(NULOUT,*) ' INVERSE OF A :'
  DO JLEV=1,NFLEVG+1
    WRITE(NULOUT,'(100F6.0)') (ZAMATI(JLEV,I),I=1,NFLEVG+1)
  ENDDO
  WRITE(NULOUT,*) ' MATRIX B*1E6 :'
  DO JLEV=1,NFLEVG+1
    WRITE(NULOUT,'(100F6.1)') (1.0E6_JPRD*ZBMAT(JLEV,I),I=1,NFLEVG+1)
  ENDDO
ENDIF
ZPROD(:,:)=0.0_JPRD
DO JI=1,NFLEVG+1
  DO JJ=1,NFLEVG+1
    DO JK=1,NFLEVG+1
      ZPROD(JJ,JI)=ZPROD(JJ,JI)+ZAMATI(JJ,JK)*ZBMAT(JK,JI)
    ENDDO
  ENDDO
ENDDO
DO JI=1,NFLEVG
  DO JJ=1,NFLEVG+1
    IF(JI == 1) THEN
      YDVFE%RINTE(JJ,JI)=ZPROD(JJ,JI)+ ZPROD(JJ,JI+1)
    ELSE
      YDVFE%RINTE(JJ,JI)=ZPROD(JJ,JI+1)
    ENDIF
  ENDDO
ENDDO

IF(MYPROC == 1) THEN
  IF(LLPRINT) THEN
    WRITE(NULOUT,*) ' MATRIX FOR INTEGRALS (*1E4):'
    DO JLEV=1,NFLEVG+1
      WRITE(NULOUT,'(100F6.1)') (1.0E4_JPRD*YDVFE%RINTE(JLEV,I),I=1,NFLEVG)
    ENDDO
  ENDIF
  ZDERI(:,:)= ZPROD(:,:)
  ZDET=0.0_JPRD
  ZEPS=100.0_JPRD*TINY(1.0_JPRD)
  CALL MINV_8(ZDERI,NFLEVG+1,NFLEVG+1,ZWORK,ZDET,ZEPS,0,1)
  IF(LLPRINT) THEN
    WRITE(NULOUT,*) ' DETERMINANT=',ZDET
    WRITE(NULOUT,*) ' INVERSE OF MATRIX FOR INTEGRALS:'
    DO JLEV=1,NFLEVG+1
      WRITE(NULOUT,'(100F6.1)') (ZDERI(JLEV,I),I=1,NFLEVG+1)
    ENDDO
  ENDIF
  WRITE(NULOUT,*) ' DERIVATIVES OF FUNCTION f(x)=1'
  DO JLEV=1,NFLEVG
    ZD0=0.0_JPRD
    DO JK=1,NFLEVG+1
      ZD0=ZD0+ZDERI(JLEV+1,JK)
    ENDDO
    WRITE(NULOUT,*) ' LEVEL ',JLEV,' DERIVATIVE=',ZD0
  ENDDO
ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUVERTFE1',1,ZHOOK_HANDLE)
END SUBROUTINE SUVERTFE1
