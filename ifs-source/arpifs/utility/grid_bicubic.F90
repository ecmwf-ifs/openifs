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

SUBROUTINE GRID_BICUBIC(KLAT_IN,KLON_IN,PLAT_IN,KIN,&
                      & KLAT_OUT,KLON_OUT,PLAT_OUT,KOUT,&
                      & PFLD_IN,PFLD_OUT)  

!    Purpose.
!    --------
!      Bicubic interpolation of grid field from model resolution
!      to grid defined by input parameter. Assumes the whole
!      field is in memory.

!    Explicit arguments:
!    -------------------

!    Input:
!      KLAT_IN  - number of latitude rows of the input grid
!      KLON_IN  - number of longitudes for each row of the input grid
!      PLAT_IN  - latitude (radians) of each row of the input grid
!      KOUT     - size of input array
!      KLAT_OUT - number of latitude rows of the output grid
!      KLON_OUT - number of longitudes for each row of the output grid
!      PLAT_OUT - latitude (radians) of each row of the output grid
!      KOUT     - size of output array
!      PFLD_IN  - array of grid values at input resolution

!    Output:
!      PFLD_OUT - array of interpolated values

!    Author.
!    -------
!      Y.Tremolet

!    Modifications.
!    --------------
!      Original: 24/04/01
!        M.Hamrud      01-Oct-2003 CY28 Cleaning

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMCST, ONLY : RPI
USE YOMLUN, ONLY : NULOUT

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

INTEGER(KIND=JPIM) :: JLAT, JLON, ILATS(4), ILAT, II, JJ, JIN, J1, J2, J3, J4
REAL(KIND=JPRB) :: ZLAT, ZOUT, ZDLONOUT, ZDLONIN, ZXI, ZP0, ZP1, ZM1, ZM2, &
 & Z1,Z2,Z3,Z4,ZX1,ZX2,ZX3,ZX4,ZD1,ZD2,ZD3,ZD4,ZZ1,ZZ2,ZZ3,ZZ4  
REAL(KIND=JPRB), ALLOCATABLE :: ZTMP(:,:)
INTEGER(KIND=JPIM) :: IND(KLAT_IN)
REAL(KIND=JPRB) :: ZZ, ZMINH, ZMAXH, ZAVGH, ZMINL, ZMAXL, ZAVGL
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GRID_BICUBIC',0,ZHOOK_HANDLE)

CALL GSTATS(1944,0)
II=1
DO JLAT=1,KLAT_IN
  IND(JLAT)=II
  II=II+KLON_IN(JLAT)
ENDDO

II=KLON_OUT(1)
DO JLAT=2,KLAT_OUT
  II=MAX(II,KLON_OUT(JLAT))
ENDDO
ALLOCATE(ZTMP(II,4))

II=0
DO JLAT=1,KLAT_OUT
  ZLAT=PLAT_OUT(JLAT)
  ILAT=1
  DO WHILE (ILAT<KLAT_IN .AND. ZLAT<PLAT_IN(ILAT))
    ILAT=ILAT+1
  ENDDO

! Check the bounds
! ... assuming nobody runs the IFS with less than 4 latitude rows !!!
  ILAT=MIN(ILAT,KLAT_IN-1)
  ILAT=MAX(ILAT,3)

  ILATS(1)=ILAT-2
  ILATS(2)=ILAT-1
  ILATS(3)=ILAT
  ILATS(4)=ILAT+1

  ZDLONOUT=2.0*RPI/KLON_OUT(JLAT)
  DO JJ=1,4
    ZDLONIN=2.0*RPI/KLON_IN(ILATS(JJ))
    DO JLON=1,KLON_OUT(JLAT)
      ZOUT=(JLON-1)*ZDLONOUT
      JIN=INT(ZOUT/ZDLONIN)
      J1=MOD(JIN-1,KLON_IN(ILATS(JJ)))
      J2=MOD(JIN  ,KLON_IN(ILATS(JJ)))
      J3=MOD(JIN+1,KLON_IN(ILATS(JJ)))
      J4=MOD(JIN+2,KLON_IN(ILATS(JJ)))
      ZXI=(J2-1)*ZDLONIN
      ZP0=(ZOUT-ZXI)/ZDLONIN
      ZP1=ZP0+1.0
      ZM1=ZP0-1.0
      ZM2=ZP0-2.0
      ZTMP(JLON,JJ)=-ZP0*ZM1*ZM2*PFLD_IN(IND(ILATS(JJ))+J1)/6. &
       & +ZP1*ZM1*ZM2*PFLD_IN(IND(ILATS(JJ))+J2)/2. &
       & -ZP1*ZP0*ZM2*PFLD_IN(IND(ILATS(JJ))+J3)/2. &
       & +ZP1*ZP0*ZM1*PFLD_IN(IND(ILATS(JJ))+J4)/6.  
    ENDDO
  ENDDO

! zi is x_i, zxi is (x-x_i)

  Z1=PLAT_IN(ILATS(1))
  Z2=PLAT_IN(ILATS(2))
  Z3=PLAT_IN(ILATS(3))
  Z4=PLAT_IN(ILATS(4))

  ZX1=ZLAT-Z1
  ZX2=ZLAT-Z2
  ZX3=ZLAT-Z3
  ZX4=ZLAT-Z4

  ZD1=(Z1-Z2)*(Z1-Z3)*(Z1-Z4)
  ZD2=(Z2-Z1)*(Z2-Z3)*(Z2-Z4)
  ZD3=(Z3-Z1)*(Z3-Z2)*(Z3-Z4)
  ZD4=(Z4-Z1)*(Z4-Z2)*(Z4-Z3)

  ZZ1=ZX2*ZX3*ZX4/ZD1
  ZZ2=ZX1*ZX3*ZX4/ZD2
  ZZ3=ZX1*ZX2*ZX4/ZD3
  ZZ4=ZX1*ZX2*ZX3/ZD4

  DO JLON=1,KLON_OUT(JLAT)
    II=II+1
    PFLD_OUT(II)=ZZ1 * ZTMP(JLON,1) &
     & + ZZ2 * ZTMP(JLON,2) &
     & + ZZ3 * ZTMP(JLON,3) &
     & + ZZ4 * ZTMP(JLON,4)  
  ENDDO
ENDDO

IF (ALLOCATED(ZTMP)) DEALLOCATE(ZTMP)

! This check is too slow. Should be done only for
! fields which actually require it.

ZMINH=PFLD_IN(1)
ZMAXH=PFLD_IN(1)
ZAVGH=PFLD_IN(1)
DO II=2,KIN
  ZZ=PFLD_IN(II)
  ZMINH=MIN(ZZ,ZMINH)
  ZMAXH=MAX(ZZ,ZMAXH)
  ZAVGH=ZAVGH+ZZ
ENDDO
ZAVGH=ZAVGH/II

ZMINL=PFLD_OUT(1)
ZMAXL=PFLD_OUT(1)
ZAVGL=PFLD_OUT(1)
DO II=2,KOUT
  ZZ=PFLD_OUT(II)
  ZMINL=MIN(ZZ,ZMINL)
  ZMAXL=MAX(ZZ,ZMAXL)
  ZAVGL=ZAVGL+ZZ
ENDDO
ZAVGL=ZAVGL/II

IF (ZMINH>=0. .AND. ZMINL<=0.) THEN
  WRITE(NULOUT,999)'IN =',ZMINH,ZMAXH,ZAVGH
  WRITE(NULOUT,999)'OUT=',ZMINL,ZMAXL,ZAVGL

  DO II=1,KOUT
    PFLD_OUT(II)=MAX(PFLD_OUT(II),ZMINH)
  ENDDO
ENDIF
CALL GSTATS(1944,1)

999 FORMAT('GREPTRAJ MIN,MAX,AVG ',A4,3(2X,ES11.4))

IF (LHOOK) CALL DR_HOOK('GRID_BICUBIC',1,ZHOOK_HANDLE)
END SUBROUTINE GRID_BICUBIC
