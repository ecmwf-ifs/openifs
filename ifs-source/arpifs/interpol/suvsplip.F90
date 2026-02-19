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

SUBROUTINE SUVSPLIP(YDGEOMETRY,YDVSPLIP)

!**** *SUVSPLIP*  - Setup for vertical cubic spline interpolation
!                   (externalisable part of interpolator)

!     Purpose.
!     --------
!           Setup matrix for cubic B-spline interpolation in the vertical
!           (fills RVSPTRI ,RVSPC ,RFVV).
!           Intermediate quantities used for vertical interpolations.

!**   Interface.
!     ----------

!     *CALL* SUVSPLIP

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
!        Agathe Untch ECMWF

!     Modifications.
!     --------------
!        Original : 2000-06
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        K. Yessad (Jan 2011): cleanings.
!        T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMVSPLIP,     ONLY : TVSPLIP
USE PARKIND1 , ONLY : JPIM, JPRB, JPRD
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK

USE YOMLUN_IFSAUX, ONLY : NULOUT

! arp/ifs dependencies to be solved later.
USE YOMMP0   , ONLY : MYPROC

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY) , INTENT(IN)    :: YDGEOMETRY
TYPE(TVSPLIP),   INTENT(INOUT) :: YDVSPLIP

INTEGER(KIND=JPIM) :: IDIM,JR,JC,JK,JKC
REAL(KIND=JPRD) :: ZDELTA(-3:YDGEOMETRY%YRDIMV%NFLEVG+4),ZD1,ZD2,ZD3,ZD4
REAL(KIND=JPRD) :: ZDEL, ZDEL2, ZDEL3
REAL(KIND=JPRD) :: ZAAA2,ZBBB2,ZCCC2,ZDDD2,ZAAA3,ZBBB3,ZCCC3,ZDDD3
REAL(KIND=JPRD) :: ZDDD1,ZDDD4
REAL(KIND=JPRD) :: ZALPHA(4,4), ZBETA(2,2), ZRHS(2), ZBMAX, ZDETB
REAL(KIND=JPRD) :: ZSF(YDGEOMETRY%YRDIMV%NFLEVG+2,YDGEOMETRY%YRDIMV%NFLEVG+2), &
 & ZSFINV(YDGEOMETRY%YRDIMV%NFLEVG+2,YDGEOMETRY%YRDIMV%NFLEVG+2)
REAL(KIND=JPRD) :: ZWORK(2*(YDGEOMETRY%YRDIMV%NFLEVG+2)), ZEPS, ZDET

LOGICAL, PARAMETER :: LLPRINT = .FALSE.

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "minv_8.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUVSPLIP',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

ALLOCATE(YDVSPLIP%RVSPTRI (YDGEOMETRY%YRDIMV%NFLEVG,-1:1))
ALLOCATE(YDVSPLIP%RVSPC (10))
ALLOCATE(YDVSPLIP%RFVV (4,-1:YDGEOMETRY%YRDIMV%NFLEVG+2,4))

ASSOCIATE(YDVAB=>YDGEOMETRY%YRVAB,YDVETA=>YDGEOMETRY%YRVETA,YDVFE=>YDGEOMETRY%YRVFE,YDSTA=>YDGEOMETRY%YRSTA, &
 & YDLAP=>YDGEOMETRY%YRLAP,YDCSGLEG=>YDGEOMETRY%YRCSGLEG, &
 & YDCSGEOM=>YDGEOMETRY%YRCSGEOM,YDCSGEOM_NB=>YDGEOMETRY%YRCSGEOM_NB,YDGSGEOM=>YDGEOMETRY%YRGSGEOM, &
 & YDGSGEOM_NB=>YDGEOMETRY%YRGSGEOM_NB, YDSPGEOM=>YDGEOMETRY%YSPGEOM)
ASSOCIATE(NFLEVG=>YDGEOMETRY%YRDIMV%NFLEVG, &
 & RFVV=>YDVSPLIP%RFVV, &
 & RVSPC=>YDVSPLIP%RVSPC, RVSPTRI=>YDVSPLIP%RVSPTRI,YDDIMV=>YDGEOMETRY%YRDIMV)
IDIM=NFLEVG+2
ZEPS=100.0_JPRD*TINY(1.0_JPRD)

! compute distance between nodes (full eta levels)

ZDELTA( 0)=2._JPRD*(YDVETA%VETAF(1)-YDVETA%VETAH(0))
ZDELTA(-3)=ZDELTA( 0)
ZDELTA(-2)=ZDELTA( 0)
ZDELTA(-1)=ZDELTA( 0)
DO JK=1,NFLEVG-1
  ZDELTA(JK)=YDVETA%VETAF(JK+1)-YDVETA%VETAF(JK)
ENDDO
ZDELTA(NFLEVG  )=2._JPRD*(YDVETA%VETAH(NFLEVG) -YDVETA%VETAF(NFLEVG))
ZDELTA(NFLEVG+1)=ZDELTA(NFLEVG  )
ZDELTA(NFLEVG+2)=ZDELTA(NFLEVG  )
ZDELTA(NFLEVG+3)=ZDELTA(NFLEVG  )
ZDELTA(NFLEVG+4)=ZDELTA(NFLEVG  )

! compute basis functions    

RFVV(:,:,:)=0.0_JPRD

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

  ZAAA2=        ZD1**3*ZDDD1
  ZBBB2=3._JPRD*ZD1**2*ZDDD1
  ZCCC2=3._JPRD*ZD1   *ZDDD1

  ZAAA3=        ZD4**3*ZDDD4
  ZBBB3=3._JPRD*ZD4**2*ZDDD4
  ZCCC3=3._JPRD*ZD4   *ZDDD4

! coefficients of explicit cubic polynomials in each interval
  
  RFVV(1,JK,1)=0.0_JPRD
  RFVV(1,JK,2)=0.0_JPRD
  RFVV(1,JK,3)=0.0_JPRD
  RFVV(1,JK,4)=ZDDD1
 
  RFVV(2,JK,1)=ZAAA2
  RFVV(2,JK,2)=ZBBB2
  RFVV(2,JK,3)=ZCCC2
  RFVV(2,JK,4)=ZDDD2  
  
  ZDEL=ZDELTA(JK)
  ZDEL2=ZDEL*ZDEL
  ZDEL3=ZDEL2*ZDEL
  RFVV(3,JK,1)=ZAAA3+ZBBB3*ZDEL+        ZCCC3*ZDEL2        +ZDDD3*ZDEL3
  RFVV(3,JK,2)=     -ZBBB3     -2._JPRD*ZCCC3*ZDEL -3._JPRD*ZDDD3*ZDEL2
  RFVV(3,JK,3)=                         ZCCC3      +3._JPRD*ZDDD3*ZDEL
  RFVV(3,JK,4)=                                            -ZDDD3

  ZDEL=ZDELTA(JK+1)
  ZDEL2=ZDEL*ZDEL
  ZDEL3=ZDEL2*ZDEL
  RFVV(4,JK,1)=         ZDDD4*ZDEL3
  RFVV(4,JK,2)=-3._JPRD*ZDDD4*ZDEL2
  RFVV(4,JK,3)= 3._JPRD*ZDDD4*ZDEL
  RFVV(4,JK,4)=        -ZDDD4 

ENDDO !end loop over nodes

! COMPUTE MATRIX SF ( SF(I,J)=BASISFUNC_J(X_I) )
! ---------------------------------------------

ZSF(:,:)=0.0_JPRD
DO JR=2,IDIM-1
  DO JC=1,IDIM
    JKC=JC-1

    IF(JC == JR-1) ZSF(JR,JC)=RFVV(4,JKC,1)
    IF(JC == JR  ) ZSF(JR,JC)=RFVV(3,JKC,1)
    IF(JC == JR+1) ZSF(JR,JC)=RFVV(2,JKC,1)

  ENDDO
ENDDO

! derivatives

! f' 
!!ZSF(1,1)=RFVV(4,0,2)
!!ZSF(1,2)=RFVV(3,1,2)
!!ZSF(1,3)=RFVV(2,2,2)

! f''
ZSF(1,1)=2.0_JPRD*RFVV(4,0,3)
ZSF(1,2)=2.0_JPRD*RFVV(3,1,3)
ZSF(1,3)=2.0_JPRD*RFVV(2,2,3)

! f'
!!ZSF(IDIM,IDIM-2)=RFVV(4,NFLEVG-1,2)
!!ZSF(IDIM,IDIM-1)=RFVV(3,NFLEVG  ,2)
!!ZSF(IDIM,IDIM  )=RFVV(2,NFLEVG+1,2)

! f''
ZSF(IDIM,IDIM-2)=2.0_JPRD*RFVV(4,NFLEVG-1,3)
ZSF(IDIM,IDIM-1)=2.0_JPRD*RFVV(3,NFLEVG  ,3)
ZSF(IDIM,IDIM  )=2.0_JPRD*RFVV(2,NFLEVG+1,3)

IF(MYPROC == 1.AND.LLPRINT) THEN
  ! invert SF: 
  ZDET=0.0_JPRD
  ZSFINV(:,:)=ZSF(:,:)
  CALL MINV_8(ZSFINV,IDIM,IDIM,ZWORK,ZDET,ZEPS,0,1)

! print S and S^-1
  WRITE(NULOUT,*) ' DETERMINANT OF SF =',ZDET
  WRITE(NULOUT,*) ' MATRIX SF :'
  DO JR=1,IDIM
    WRITE(NULOUT,'(100F9.3)') (ZSF(JR,JC),JC=1,IDIM)
  ENDDO
  WRITE(NULOUT,*) ' INVERSE OF SF :'
  DO JR=1,IDIM
    WRITE(NULOUT,'(100F9.3)') (ZSFINV(JR,JC),JC=1,IDIM)
  ENDDO
ENDIF

! tridiagonal band of SF for tridiagonal solver

RVSPTRI(:,:)=0.0_JPRD

! diagonal
RVSPTRI(1,0)=ZSF(2,2)-ZSF(2,1)*ZSF(1,2)/ZSF(1,1)
DO JR=2,NFLEVG-1
  RVSPTRI(JR,0)=ZSF(JR+1,JR+1)
ENDDO
RVSPTRI(NFLEVG,0)=ZSF(NFLEVG+1,NFLEVG+1)&
 & -ZSF(NFLEVG+1,NFLEVG+2)*ZSF(NFLEVG+2,NFLEVG+1)/ZSF(NFLEVG+2,NFLEVG+2)  

! lower off-diagonal
DO JR=2,NFLEVG-1
  RVSPTRI(JR,-1)=ZSF(JR+1,JR)
ENDDO
RVSPTRI(NFLEVG,-1)=ZSF(NFLEVG+1,NFLEVG) &
 & -ZSF(NFLEVG+1,NFLEVG+2)*ZSF(NFLEVG+2,NFLEVG)/ZSF(NFLEVG+2,NFLEVG+2)  

! upper off-diagonal
RVSPTRI(1,1)= ZSF(2,3)-ZSF(2,1)*ZSF(1,3)/ZSF(1,1)
DO JR=2,NFLEVG-1
  RVSPTRI(JR,1)=ZSF(JR+1,JR+2)
ENDDO 

IF(MYPROC == 1.AND.LLPRINT) THEN
  WRITE(NULOUT,*) ' tridiagonal band of SF:'
  WRITE(NULOUT,'(100F9.3)') (RVSPTRI(JR,-1),JR=1,NFLEVG)
  WRITE(NULOUT,'(100F9.3)') (RVSPTRI(JR, 0),JR=1,NFLEVG)
  WRITE(NULOUT,'(100F9.3)') (RVSPTRI(JR,+1),JR=1,NFLEVG)
ENDIF

! constants needed in VSPLTRANS

RVSPC(1)=  1.0_JPRD/ZSF(1,1)
RVSPC(2)= -ZSF(1,2)/ZSF(1,1)
RVSPC(3)= -ZSF(1,3)/ZSF(1,1)
RVSPC(4)= -ZSF(2,1)/ZSF(1,1)
RVSPC(5)=  1.0_JPRD/ZSF(NFLEVG+2,NFLEVG+2)
RVSPC(6)= -ZSF(NFLEVG+2,NFLEVG  )/ZSF(NFLEVG+2,NFLEVG+2)
RVSPC(7)= -ZSF(NFLEVG+2,NFLEVG+1)/ZSF(NFLEVG+2,NFLEVG+2)
RVSPC(8)= -ZSF(NFLEVG+1,NFLEVG+2)/ZSF(NFLEVG+2,NFLEVG+2)
RVSPC(9)= 1.0_JPRD/(YDVETA%VETAF(2)-YDVETA%VETAF(1))
RVSPC(10)= 1.0_JPRD/(YDVETA%VETAF(NFLEVG)-YDVETA%VETAF(NFLEVG-1))

IF(MYPROC == 1.AND.LLPRINT) THEN
  WRITE(NULOUT,*) ' constants RVSPC: ',(RVSPC(JR),JR=1,10)
ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUVSPLIP',1,ZHOOK_HANDLE)
END SUBROUTINE SUVSPLIP
