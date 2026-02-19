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

SUBROUTINE SUVERTFE3D(YDGEOMETRY)

!**** *SUVERTFE3D*  - Setup VERTical Finite Element Derivative based on
!                    cubic (3) elements

!     Purpose.
!     --------
!           Compute matrix for the vertical differentiation based on
!           cubic B-spline finite elements

!**   Interface.
!     ----------

!     *CALL* SUVERTFE3D
!        Explicit arguments :
!        --------------------
!        None

!        Implicit arguments :
!        --------------------

!     Method.
!     -------

!     Externals.
!     ----------
!        MINV

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Mariano Hortal (after SUVERTFE3) ECMWF

!     Modifications.
!     --------------
!        Original : 2004-10
!        K. Yessad (Sep 2008): use VFE_ETAH and VFE_ETAF
!        T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1 , ONLY : JPIM, JPRB, JPRD
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY), INTENT(INOUT) :: YDGEOMETRY
INTEGER(KIND=JPIM) :: IDIM, J, JK, I_K, JKC, JD, JR, JKR, J1, J2, J3, JC
REAL(KIND=JPRD)    :: ZDELTA(-3:YDGEOMETRY%YRDIMV%NFLEVG+4)
REAL(KIND=JPRD)    :: ZFAA(4,-1:YDGEOMETRY%YRDIMV%NFLEVG+2), ZFBB(4,-1:YDGEOMETRY%YRDIMV%NFLEVG+2),&
                     &ZFCC(4,-1:YDGEOMETRY%YRDIMV%NFLEVG+2), ZFDD(4,-1:YDGEOMETRY%YRDIMV%NFLEVG+2)
REAL(KIND=JPRD)    :: ZDAA(4,-1:YDGEOMETRY%YRDIMV%NFLEVG+2), ZDBB(4,-1:YDGEOMETRY%YRDIMV%NFLEVG+2),&
                     &ZDCC(4,-1:YDGEOMETRY%YRDIMV%NFLEVG+2), ZDDD(4,-1:YDGEOMETRY%YRDIMV%NFLEVG+2)
REAL(KIND=JPRD)    :: ZA(YDGEOMETRY%YRDIMV%NFLEVG+4,YDGEOMETRY%YRDIMV%NFLEVG+4), &
                     &ZAINV(YDGEOMETRY%YRDIMV%NFLEVG+4,YDGEOMETRY%YRDIMV%NFLEVG+4),&
                     &ZWORK(2*(YDGEOMETRY%YRDIMV%NFLEVG+3)), ZPROD(YDGEOMETRY%YRDIMV%NFLEVG+4,YDGEOMETRY%YRDIMV%NFLEVG+4),&
                     &ZB(YDGEOMETRY%YRDIMV%NFLEVG+4,YDGEOMETRY%YRDIMV%NFLEVG+4)
REAL(KIND=JPRD)    :: ZDERI(YDGEOMETRY%YRDIMV%NFLEVG+4,YDGEOMETRY%YRDIMV%NFLEVG+4), &
                     &ZSF(YDGEOMETRY%YRDIMV%NFLEVG+4,YDGEOMETRY%YRDIMV%NFLEVG+4),&
                     &ZSFINV(YDGEOMETRY%YRDIMV%NFLEVG+4,YDGEOMETRY%YRDIMV%NFLEVG+4), ZSD(YDGEOMETRY%YRDIMV%NFLEVG+4, &
                     &YDGEOMETRY%YRDIMV%NFLEVG+4)
REAL(KIND=JPRD)    :: ZEPS, ZALPHA(4,4), ZBETA(2,2), ZHH(4,4)
REAL(KIND=JPRD)    :: ZRHS(2)
REAL(KIND=JPRD)    :: ZD0 , ZD0H, ZBMAX, ZD1, ZD2, ZD3, ZD4, ZDETB, ZDEL, ZDEL2, ZDEL3
REAL(KIND=JPRD)    :: ZDDD1, ZDDD2, ZDDD3, ZDDD4, ZDNH, ZDET
REAL(KIND=JPRD)    :: ZAAA2, ZBBB2, ZCCC2, ZAAA3, ZBBB3, ZCCC3
REAL(KIND=JPRD)    :: FHH, FDD
REAL(KIND=JPRD)    :: XA,XB,XC,XD,YA,YB,YC,YD, XX
REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

#include "minv_8.h"

! indefiniteINTEGRAL(EJ*EI)dx
FHH(XX,XA,XB,XC,XD,YA,YB,YC,YD)=&
 &   XX*(  XA*YA&
 &   + XX*( (XA*YB+XB*YA)/2._JPRD&
 &   + XX*( (XA*YC+XB*YB+XC*YA)/3._JPRD&
 &   + XX*( (XA*YD+XB*YC+XC*YB+XD*YA)/4._JPRD&
 &   + XX*( (      XB*YD+XC*YC+XD*YB)/5._JPRD&
 &   + XX*( (            XC*YD+XD*YC)/6._JPRD&
 & + XX*                     XD*YD /7._JPRD))))))  

! indefiniteINTEGRAL(EJ*EI')dx                '
FDD(XX,XA,XB,XC,XD,YA,YB,YC,YD)=&
 &   XX*(     XA*YB&
 &   + XX*( (2._JPRD*XA*YC+XB*YB         )/2._JPRD&
 &   + XX*( (3._JPRD*XA*YD+2._JPRD*XB*YC+XC*YB)/3._JPRD&
 &   + XX*( (3._JPRD*XB*YD+2._JPRD*XC*YC+XD*YB)/4._JPRD&
 &   + XX*( (3._JPRD*XC*YD+2._JPRD*XD*YC      )/5._JPRD&
 & + XX*( (3._JPRD*XD*YD               )/6._JPRD))))))

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUVERTFE3D',0,ZHOOK_HANDLE)
ASSOCIATE(YDVAB=>YDGEOMETRY%YRVAB,YDVETA=>YDGEOMETRY%YRVETA,YDVFE=>YDGEOMETRY%YRVFE,YDSTA=>YDGEOMETRY%YRSTA, &
 & YDLAP=>YDGEOMETRY%YRLAP,YDCSGLEG=>YDGEOMETRY%YRCSGLEG, &
 & YDCSGEOM=>YDGEOMETRY%YRCSGEOM,YDCSGEOM_NB=>YDGEOMETRY%YRCSGEOM_NB,YDGSGEOM=>YDGEOMETRY%YRGSGEOM, &
 & YDGSGEOM_NB=>YDGEOMETRY%YRGSGEOM_NB, YDSPGEOM=>YDGEOMETRY%YSPGEOM)
ASSOCIATE(NFLEVG=>YDGEOMETRY%YRDIMV%NFLEVG, &
  & YDDIMV=>YDGEOMETRY%YRDIMV)!     ------------------------------------------------------------------

IDIM=NFLEVG+4
ZEPS=100.0_JPRD*TINY(1.0_JPRD)

! compute distance between nodes (full eta levels)

ZDELTA( 0)=2._JPRD*(YDVETA%VFE_ETAF(1)-YDVETA%VFE_ETAH(0))
ZD0=ZDELTA( 0)
ZD0H=0.5_JPRD*ZD0
ZDELTA(-3)=ZDELTA(0)
ZDELTA(-2)=ZDELTA(0)
ZDELTA(-1)=ZDELTA(0)
DO JK=1,NFLEVG-1
  ZDELTA(JK)=YDVETA%VFE_ETAF(JK+1)-YDVETA%VFE_ETAF(JK)
ENDDO
ZDELTA(NFLEVG)=2._JPRD*(YDVETA%VFE_ETAH(NFLEVG) -YDVETA%VFE_ETAF(NFLEVG))
ZDNH=ZDELTA(NFLEVG)*0.5_JPRD
ZDELTA(NFLEVG+1)=ZDELTA(NFLEVG)
ZDELTA(NFLEVG+2)=ZDELTA(NFLEVG)
ZDELTA(NFLEVG+3)=ZDELTA(NFLEVG)
ZDELTA(NFLEVG+4)=ZDELTA(NFLEVG)

! compute basis functions    

ZBMAX=4._JPRD/6._JPRD 
  
DO JK=-1,NFLEVG+2
  ZD1=ZDELTA(JK-2)
  ZD2=ZDELTA(JK-1)
  ZD3=ZDELTA(JK  )
  ZD4=ZDELTA(JK+1)
! matrix alpha
  ZALPHA(1,1)=0.0_JPRD     
  ZALPHA(1,2)=0.0_JPRD      
  ZALPHA(1,3)=ZD3**3      
  ZALPHA(1,4)=ZD4**3+3._JPRD*ZD3*ZD4**2+3._JPRD*ZD3**2*ZD4  
  ZALPHA(2,1)=ZD1**2+2._JPRD*ZD1*ZD2
  ZALPHA(2,2)=ZD2**2
  ZALPHA(2,3)=ZD3**2
  ZALPHA(2,4)=ZD4**2+2._JPRD*ZD3*ZD4
  ZALPHA(3,1)= ZD1
  ZALPHA(3,2)= ZD2
  ZALPHA(3,3)=-ZD3
  ZALPHA(3,4)=-ZD4
  ZALPHA(4,1)=ZD1**3+3._JPRD*ZD1**2*ZD2+3._JPRD*ZD1*ZD2**2
  ZALPHA(4,2)=ZD2**3
  ZALPHA(4,3)=0.0_JPRD
  ZALPHA(4,4)=0.0_JPRD

! matrix beta
  ZBETA(1,1)=ZALPHA(2,2)-ZALPHA(2,1)*ZALPHA(4,2)/ZALPHA(4,1)
  ZBETA(1,2)=ZALPHA(2,3)-ZALPHA(2,4)*ZALPHA(1,3)/ZALPHA(1,4)
  ZBETA(2,1)=ZALPHA(3,2)-ZALPHA(3,1)*ZALPHA(4,2)/ZALPHA(4,1)
  ZBETA(2,2)=ZALPHA(3,3)-ZALPHA(3,4)*ZALPHA(1,3)/ZALPHA(1,4)

! right-hand side vector
  ZRHS(1)=-ZBMAX*(ZALPHA(2,1)/ZALPHA(4,1)+ZALPHA(2,4)/ZALPHA(1,4))
  ZRHS(2)=-ZBMAX*(ZALPHA(3,1)/ZALPHA(4,1)+ZALPHA(3,4)/ZALPHA(1,4))

! determinant of beta 
  ZDETB=ZBETA(1,1)*ZBETA(2,2)-ZBETA(1,2)*ZBETA(2,1)

! coefficients 
  ZDDD2=(ZRHS(1)*ZBETA(2,2)-ZRHS(2)*ZBETA(1,2))/ZDETB
  ZDDD3=(ZRHS(2)*ZBETA(1,1)-ZRHS(1)*ZBETA(2,1))/ZDETB
  ZDDD1=(ZBMAX - ZALPHA(4,2)*ZDDD2)/ZALPHA(4,1)    
  ZDDD4=(ZBMAX - ZALPHA(1,3)*ZDDD3)/ZALPHA(1,4) 

  ZAAA2=   ZD1**3*ZDDD1
  ZBBB2=3._JPRD*ZD1**2*ZDDD1
  ZCCC2=3._JPRD*ZD1   *ZDDD1

  ZAAA3=   ZD4**3*ZDDD4
  ZBBB3=3._JPRD*ZD4**2*ZDDD4
  ZCCC3=3._JPRD*ZD4   *ZDDD4

! coefficients of explicit cubic polynomials in each interval
! (x=0 at the left of the corresponding interval)
  
  ZFAA(1,JK)=0.0_JPRD
  ZFBB(1,JK)=0.0_JPRD
  ZFCC(1,JK)=0.0_JPRD
  ZFDD(1,JK)=ZDDD1
 
  ZFAA(2,JK)=ZAAA2
  ZFBB(2,JK)=ZBBB2
  ZFCC(2,JK)=ZCCC2
  ZFDD(2,JK)=ZDDD2  
  
  ZDEL=ZDELTA(JK)
  ZDEL2=ZDEL*ZDEL
  ZDEL3=ZDEL2*ZDEL
  ZFAA(3,JK)=ZAAA3+ZBBB3*ZDEL+   ZCCC3*ZDEL2   +ZDDD3*ZDEL3
  ZFBB(3,JK)=     -ZBBB3     -2._JPRD*ZCCC3*ZDEL -3._JPRD*ZDDD3*ZDEL2
  ZFCC(3,JK)=                    ZCCC3      +3._JPRD*ZDDD3*ZDEL
  ZFDD(3,JK)=                                  -ZDDD3

  ZDEL=ZDELTA(JK+1)
  ZDEL2=ZDEL*ZDEL
  ZDEL3=ZDEL2*ZDEL
  ZFAA(4,JK)=    ZDDD4*ZDEL3
  ZFBB(4,JK)=-3._JPRD*ZDDD4*ZDEL2
  ZFCC(4,JK)= 3._JPRD*ZDDD4*ZDEL
  ZFDD(4,JK)=   -ZDDD4 
ENDDO

ZDAA(:,:)=ZFAA(:,:)
ZDBB(:,:)=ZFBB(:,:)
ZDCC(:,:)=ZFCC(:,:)
ZDDD(:,:)=ZFDD(:,:)

! COMPUTE MATRIX SF ( SF(I,J)=BASISFUNC_J(X_I) )  (from FE space of the function
!                                                 to be differenciated into physical space)
! ---------------------------------------------

ZSF(:,:)=0.0_JPRD
DO JR=1,IDIM
  DO JC=1,IDIM
    JKC=JC-2
    IF(JC == JR-1) ZSF(JR,JC)=ZFAA(4,JKC)
    IF(JC == JR  ) ZSF(JR,JC)=ZFAA(3,JKC)
    IF(JC == JR+1) ZSF(JR,JC)=ZFAA(2,JKC)
  ENDDO
ENDDO

! Instead of the value at node -1 and NFLEVG+2 we should use the derivative
!   at node 1 and at NFLEVG:

ZSF(1,1)=ZFBB(4, 0)
ZSF(1,2)=ZFBB(3, 1)
ZSF(1,3)=ZFBB(2, 2)
ZSF(IDIM,IDIM-2)=ZFBB(4,NFLEVG-1)
ZSF(IDIM,IDIM-1)=ZFBB(3,NFLEVG  )
ZSF(IDIM,IDIM  )=ZFBB(2,NFLEVG+1)


! invert SF:       (ZSFINV projects the function to be integrated
!                  from physical space onto FE space)
ZDET=0.0_JPRD
ZSFINV(:,:)=ZSF(:,:)*3._JPRD/2._JPRD
CALL MINV_8(ZSFINV,IDIM,IDIM,ZWORK,ZDET,ZEPS,0,1)
ZSFINV(:,:)=ZSFINV(:,:)*3._JPRD/2._JPRD

! COMPUTE MATRIX SD ( SD(I,J)=BASISDER_J(X_I) )  (To project the derivative from FE space
!                                                to physical space 
! -------------------------------------------

ZSD(:,:)=0.0_JPRD
DO JR=1,IDIM
  DO JC=1,IDIM
    JKC=JC-2
    IF(JC == JR-1) ZSD(JR,JC)=ZDAA(4,JKC)
    IF(JC == JR  ) ZSD(JR,JC)=ZDAA(3,JKC)
    IF(JC == JR+1) ZSD(JR,JC)=ZDAA(2,JKC)
  ENDDO
ENDDO



! COMPUTE MATRIX A
! ----------------

ZA(:,:)=0.0_JPRD

DO JR=1,IDIM                                 ! JR is the row number of matrix A
  JKR=JR-2                                   ! JKR is the index of the test function
  ZHH(:,:)=0.0_JPRD                          ! ZHH is an auxiliary 4x4 matrix
  DO J=1,4                                   ! J is the index of the interval belonging to function JKR
    JD=JKR-3+J                               ! JD is the absolute number of interval J
    IF(JD >= 0 .AND. JD <= NFLEVG) THEN      ! Integration is from node 0 to node NFLEVG+1
      DO I_K=1,4                             ! I_K
        JKC=JKR+J-I_K                        ! JKC is the I_Kth basis number covering interval J
        IF(JKC >= -1 .AND. JKC <= NFLEVG+2) THEN
          ZHH(J,I_K)=FHH(ZDELTA(JD) &
             & ,ZFAA(J,JKR),ZFBB(J,JKR),ZFCC(J,JKR),ZFDD(J,JKR)&
             & ,ZFAA(I_K,JKC),ZFBB(I_K,JKC),ZFCC(I_K,JKC),ZFDD(I_K,JKC))  
        ENDIF
      ENDDO
    ENDIF
  ENDDO
  ZA(JR,JR  )=ZHH(1,1)+ZHH(2,2)+ZHH(3,3)+ZHH(4,4)
  IF(JR+1 <= IDIM) &
   & ZA(JR,JR+1)=      ZHH(2,1)+ZHH(3,2)+ZHH(4,3)  
  IF(JR+2 <= IDIM) &
   & ZA(JR,JR+2)=               ZHH(3,1)+ZHH(4,2)  
  IF(JR+3 <= IDIM) &
   & ZA(JR,JR+3)=                        ZHH(4,1)
  IF(JR-1 >= 1)    &
   & ZA(JR,JR-1)=      ZHH(1,2)+ZHH(2,3)+ZHH(3,4)
  IF(JR-2 >= 1)    &
   & ZA(JR,JR-2)=               ZHH(1,3)+ZHH(2,4)
  IF(JR-3 >= 1)    &
   & ZA(JR,JR-3)=                        ZHH(1,4)     

ENDDO

! invert A: 
ZDET=0.0_JPRD
ZAINV(:,:)=ZA(:,:)*NFLEVG*7._JPRD/2._JPRD
CALL MINV_8(ZAINV,IDIM,IDIM,ZWORK,ZDET,ZEPS,0,1)
ZAINV(:,:)=ZAINV(:,:)*NFLEVG*7._JPRD/2._JPRD

DO J1=1,IDIM
  DO J2=1,IDIM
    ZPROD(J1,J2)=0.0_JPRD
      DO J3=1,IDIM
        ZPROD(J1,J2)=ZPROD(J1,J2)+ZAINV(J1,J3)*ZA(J3,J2)
      ENDDO
  ENDDO
ENDDO

! COMPUTE MATRIX B
! ----------------

ZB(:,:)=0.0_JPRD

DO JR=1,IDIM                                 ! JR is the row number of matrix B
  JKR=JR-2                                   ! JKR is the index of the test function
  ZHH(:,:)=0.0_JPRD                               ! ZHH is an auxiliary 4x4 matrix
  DO J=1,4                                   ! J is the index of the interval belonging to function JKR
    JD=JKR-3+J                               ! JD is the absolute number of interval J
    IF(JD >= 0 .AND. JD <= NFLEVG) THEN      ! Integration is from node 0 to node NFLEVG+1
      DO I_K=1,4                             ! I_K 
        JKC=JKR+J-I_K                        ! JKC is the I_Kth basis number covering interval J
        IF(JKC >= -1 .AND. JKC <= NFLEVG+2) THEN
          ZHH(J,I_K)=FDD(ZDELTA(JD) &
             & ,ZFAA(J,JKR),ZFBB(J,JKR),ZFCC(J,JKR),ZFDD(J,JKR)&
             & ,ZFAA(I_K,JKC),ZFBB(I_K,JKC),ZFCC(I_K,JKC),ZFDD(I_K,JKC))  
        ENDIF
      ENDDO
    ENDIF
  ENDDO
  ZB(JR,JR  )=ZHH(1,1)+ZHH(2,2)+ZHH(3,3)+ZHH(4,4)
  IF(JR+1 <= IDIM) &
   & ZB(JR,JR+1)=      ZHH(2,1)+ZHH(3,2)+ZHH(4,3)  
  IF(JR+2 <= IDIM) &
   & ZB(JR,JR+2)=               ZHH(3,1)+ZHH(4,2)  
  IF(JR+3 <= IDIM) &
   & ZB(JR,JR+3)=                        ZHH(4,1)
  IF(JR-1 >= 1)    &
   & ZB(JR,JR-1)=      ZHH(1,2)+ZHH(2,3)+ZHH(3,4)
  IF(JR-2 >= 1)    &
   & ZB(JR,JR-2)=               ZHH(1,3)+ZHH(2,4)
  IF(JR-3 >= 1)    &
   & ZB(JR,JR-3)=                        ZHH(1,4)     

ENDDO

! COMPUTE PRODUCT SD*A^-1*B*SF^-1=ZDERI
! ------------------------------------

ZDERI(:,:)=0.0_JPRD
DO JC=1,IDIM
  DO JR=1,IDIM
    DO JK=1,IDIM
      ZDERI(JR,JC)=ZDERI(JR,JC)+ZB(JR,JK)*ZSFINV(JK,JC)
    ENDDO
  ENDDO
ENDDO

ZPROD(:,:)=0.0_JPRD
DO JC=1,IDIM
  DO JR=1,IDIM
    DO JK=1,IDIM
      ZPROD(JR,JC)=ZPROD(JR,JC)+ZAINV(JR,JK)*ZDERI(JK,JC)
    ENDDO
  ENDDO
ENDDO

ZDERI(:,:)=0.0_JPRD
DO JC=1,IDIM
  DO JR=1,IDIM
    DO JK=1,IDIM
      ZDERI(JR,JC)=ZDERI(JR,JC)+ZSD(JR,JK)*ZPROD(JK,JC)
    ENDDO
  ENDDO
ENDDO

! Boundary values chosen: linear extrapolation to nodes 0 and NFLEVG+1
!                         derivative at nodes 1 and NFLEVG
DO JR=1,NFLEVG
  DO JC=3,NFLEVG-2
    YDVFE%RDERI(JR,JC)=ZDERI(JR+2,JC+2)
  ENDDO
  YDVFE%RDERI(JR,1)=ZDERI(JR+2,3)+2._JPRD* ZDERI(JR+2,2)-&
                    &ZDERI(JR+2,1)/(YDVETA%VFE_ETAF(2)-YDVETA%VFE_ETAF(1))
  YDVFE%RDERI(JR,NFLEVG)=ZDERI(JR+2,NFLEVG+2)+2._JPRD* ZDERI(JR+2,NFLEVG+3)-&
                    &ZDERI(JR+2,NFLEVG+4)/(YDVETA%VFE_ETAF(NFLEVG-1)-YDVETA%VFE_ETAF(NFLEVG))
  YDVFE%RDERI(JR,2)= ZDERI(JR+2,4)-ZDERI(JR+2,2)+ZDERI(JR+2,1)/(YDVETA%VFE_ETAF(2)-YDVETA%VFE_ETAF(1))
  YDVFE%RDERI(JR,NFLEVG-1)= ZDERI(JR+2,NFLEVG+1)-ZDERI(JR+2,NFLEVG+3)+&
                    &ZDERI(JR+2,NFLEVG+4)/(YDVETA%VFE_ETAF(NFLEVG-1)-YDVETA%VFE_ETAF(NFLEVG))
ENDDO

!     ------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUVERTFE3D',1,ZHOOK_HANDLE)
END SUBROUTINE SUVERTFE3D
