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

SUBROUTINE GRID_BICONSERV (PVAH,PVBH, KLAT_IN,KLON_IN,PLAT_IN,KIN,&
                         & KLAT_OUT,KLON_OUT,PLAT_OUT,KOUT,&
                         & PFLD_IN,PFLD_OUT, &
                         & LDPAR,KLEVEL,PSP_IN,PSP_OUT)

!    Purpose.
!    --------
!      Conserving bicubic interpolation of grid field from model 
!      resolution to grid defined by input parameter. Assumes the 
!      whole field is in memory. Assumes the field is in density
!      units (upper air fields multiplied by the pressure thickness
!      dp of the layer, surface pressure as ps, not lnps).

!    Explicit arguments:
!    -------------------

!    Input:
!      PVAH/PVBH- A and B of vertical coordinates
!      KLAT_IN  - number of latitude rows of the input grid
!      KLON_IN  - number of longitudes for each row of the input grid
!      PLAT_IN  - latitude (radians) of each row of the input grid
!      KIN      - size of input array
!      KLAT_OUT - number of latitude rows of the output grid
!      KLON_OUT - number of longitudes for each row of the output grid
!      PLAT_OUT - latitude (radians) of each row of the output grid
!      KOUT     - size of output array
!      PFLD_IN  - array of grid values at input resolution

!    Output:
!      PFLD_OUT - array of interpolated values

!    Optional:
!      LDPAR    - parity for extra-polar calculations
!      KLEVEL   - vertical level
!      PSP_IN   - surface pressure at input  res. (for conserving 3D interp)
!      PSP_OUT  - surface pressure at output res. (for conserving 3D interp)

!    Author.
!    -------
!      E. Holm - modification of GRID_BICUBIC by Y.Tremolet

!    Modifications.
!    --------------
!      Original: 27/05/04
!      20050507 Elias Holm: Introduce order and limiter in call and
!                           generalize code for up to quintic interpolation
!        Y.Tremolet    26-Jan-2005 Added optional arguments
!      20080204 Elias Holm: Surface pressure on input/output and
!                           simplified
!      20081015 Elias Holm: Simplify and reorder calculations for speed
!      20130312 Elias Holm: Minor cleaning
!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM, JPRB, JPRD
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCST   , ONLY : RPI
USE ECSORT_MIX,ONLY : KEYSORT

IMPLICIT NONE

REAL (KIND=JPRB), INTENT(IN)     :: PVAH(0:)
REAL (KIND=JPRB), INTENT(IN)     :: PVBH(0:) 
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
LOGICAL, OPTIONAL, INTENT(IN)    :: LDPAR
INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KLEVEL
REAL(KIND=JPRB), OPTIONAL, INTENT(IN) :: PSP_IN(KIN)
REAL(KIND=JPRB), OPTIONAL, INTENT(IN) :: PSP_OUT(KOUT)

INTEGER(KIND=JPIM) :: IND_IN  (KLAT_IN)
INTEGER(KIND=JPIM) :: ILON_IN (-2:KLAT_IN+3)
REAL(KIND=JPRB)    :: ZLAT_IN (-2:KLAT_IN+3)
REAL(KIND=JPRB)    :: ZLON_IN (-2:KLAT_IN+3)
REAL(KIND=JPRB)    :: ZDLAT_IN(-2:KLAT_IN+3), ZDLATI_IN(-2:KLAT_IN+3)
REAL(KIND=JPRB)    :: ZDLON_IN(-2:KLAT_IN+3), ZDLONI_IN(-2:KLAT_IN+3)
REAL(KIND=JPRB)    :: ZHLAT_IN(-2:KLAT_IN+3)
REAL(KIND=JPRB)    :: ZHLON_IN(-2:KLAT_IN+3)
REAL(KIND=JPRB)    :: ZDHLAT_IN(-2:KLAT_IN+3), ZDHLATI_IN(-2:KLAT_IN+3)
REAL(KIND=JPRB)    :: ZG_IN   (-2:KLAT_IN+3)
INTEGER(KIND=JPIM) :: ILAT_IO  (0:KLAT_OUT)
REAL(KIND=JPRB)    :: ZLAT_OUT (0:KLAT_OUT), ZLON_OUT (1:KLAT_OUT)
REAL(KIND=JPRB)    :: ZDLAT_OUT(1:KLAT_OUT), ZDLON_OUT(1:KLAT_OUT)
REAL(KIND=JPRB)    :: ZHLAT_OUT(1:KLAT_OUT), ZHLON_OUT(1:KLAT_OUT)
REAL(KIND=JPRB)    :: ZG_OUT   (1:KLAT_OUT)

INTEGER(KIND=JPIM) :: J, JLAT, JLON, JNUM, JINX, JINY, JINYD, JOUTY,IRET
INTEGER(KIND=JPIM) :: JXMAX_IN, JXMAX_OUT, JYMAX_OUT,IOFF,IBL
INTEGER(KIND=JPIM) :: JBL, IBL_OUT, JINYSTA, JINYEND, IBLATSTA_OUT(1:KLAT_OUT)
INTEGER(KIND=JPIM) :: IBLATEND_OUT(1:KLAT_OUT), IBLON_OUT(1:KLAT_OUT),IOFFB(1:KLAT_OUT+1)
INTEGER(KIND=JPIM) :: INDEXB(KLAT_OUT)
REAL(KIND=JPRB)    :: ZSUM, ZOUT, ZX, ZY, ZY2, ZDHI
REAL(KIND=JPRB)    :: ZC(0:3), ZDC, ZDL, ZDLL, ZDR, ZDRR, ZCC2, ZCC3
REAL(KIND=JPRB)    :: ZMINMOD, ZCURVLIM, ZLIML, ZLIMR, ZLIM, ZPARITY, ZPRIMADD
REAL(KIND=JPRD)    :: ZFRAC, ZFR
REAL(KIND=JPRB), ALLOCATABLE :: ZPRIMX(:,:), ZPRIMY(:,:), ZPRIM(:,:), ZFLD_IN(:)
REAL(KIND=JPRB), ALLOCATABLE :: ZDIFX(:,:), ZDIFY(:,:), ZOUTPRE(:)
LOGICAL :: LLCONS3D
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
ZMINMOD(ZFR)=REAL((0.5_JPRD+SIGN(0.5_JPRD,ZFR))*MIN(1._JPRD,ABS(ZFR)),JPRB)
ZFRAC(ZDC,ZDR)=REAL(ZDC,JPRD)/SIGN(ABS(REAL(ZDR,JPRD))+1.E-30_JPRD,REAL(ZDR,JPRD))
ZCURVLIM(ZDL,ZDR)=SIGN(MIN(ABS(ZDL),ABS(ZDR)/3._JPRB),ZDL)

IF (LHOOK) CALL DR_HOOK('GRID_BICONSERV',0,ZHOOK_HANDLE)

!     *    0. CONSTANTS AND FUNCTIONS

ZCC2=1.0_JPRB/1536._JPRB
ZCC3=1.0_JPRB/15360._JPRB

ZPARITY=1.0_JPRB
LLCONS3D=PRESENT(PSP_IN) .AND. PRESENT(PSP_OUT) .AND. PRESENT(KLEVEL)
IF (PRESENT(LDPAR)) THEN
  IF (.NOT.LDPAR) ZPARITY=-1.0_JPRB
ENDIF

ALLOCATE(ZFLD_IN(KIN))
CALL GSTATS(1470,0)
IF (LLCONS3D) THEN
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(J)
  DO J=1,KIN
!!    ZFLD_IN(J)=PFLD_IN(J)*( (YRVAB%VAH(KLEVEL)-YRVAB%VAH(KLEVEL-1)) &
!!            &            +(YRVAB%VBH(KLEVEL)-YRVAB%VBH(KLEVEL-1))*PSP_IN(J) )
    ZFLD_IN(J)=PFLD_IN(J)*( (PVAH(KLEVEL)-PVAH(KLEVEL-1)) &
            &            +(PVBH(KLEVEL)-PVBH(KLEVEL-1))*PSP_IN(J) )
  ENDDO
!$OMP END PARALLEL DO
ELSE
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(J)
  DO J=1,KIN
    ZFLD_IN(J)=PFLD_IN(J)
  ENDDO
!$OMP END PARALLEL DO
ENDIF
CALL GSTATS(1470,1)

!     *    1. INPUT GRID AND FIELD

!*         1.1 Geometry of input grid
!              Assuming lat_(j+.5) = 0.5*(lat_(j) + lat_(j+1))

CALL GSTATS(1943,0)
ZLAT_IN(0)      = SIGN(RPI/2._JPRB,PLAT_IN(1))
ZLAT_IN(KLAT_IN)= SIGN(RPI/2._JPRB,PLAT_IN(KLAT_IN))
JNUM=1
JXMAX_IN=0
DO JLAT=1,KLAT_IN
  IF (JLAT<KLAT_IN) ZLAT_IN(JLAT)=0.5_JPRB*(PLAT_IN(JLAT)+PLAT_IN(JLAT+1))
  ZDLAT_IN(JLAT)=ABS(ZLAT_IN(JLAT-1)-ZLAT_IN(JLAT))
  ILON_IN(JLAT)=KLON_IN(JLAT)
  ZDLON_IN(JLAT)=2._JPRB*RPI/ILON_IN(JLAT)
  IND_IN(JLAT)=JNUM
  JNUM=JNUM+ILON_IN(JLAT)
  JXMAX_IN=MAX(JXMAX_IN,ILON_IN(JLAT))
!              Starting longitude (must be less or equal to zero)
  ZLON_IN(JLAT)=-ZDLON_IN(JLAT)/2._JPRB
!              Scale factors
  ZHLON_IN(JLAT) = (SIN(ZLAT_IN(JLAT-1))-SIN(ZLAT_IN(JLAT)))/ZDLAT_IN(JLAT)
  ZHLAT_IN(JLAT) = 1._JPRB
  ZG_IN(JLAT)=ZHLON_IN(JLAT)*ZHLAT_IN(JLAT)*ZDLON_IN(JLAT)*ZDLAT_IN(JLAT)
ENDDO

!          1.1.1 Extension over poles
! Starting longitude for extension over the poles for even ILON_IN is 0.,
! and for odd ILON_IN shift by half a gridcell to -RPI/ILON_IN.

DO JINY = 1, 3
  ILON_IN (1-JINY)       = ILON_IN (JINY)
  ZDLAT_IN(1-JINY)       = ZDLAT_IN(JINY)
  ZDLON_IN(1-JINY)       = ZDLON_IN(JINY)
  ZHLAT_IN(1-JINY)       = ZHLAT_IN(JINY)
  ZHLON_IN(1-JINY)       = ZHLON_IN(JINY)
  ZG_IN   (1-JINY)       = ZG_IN   (JINY)
  ZLAT_IN (1-JINY)       = 2._JPRB*ZLAT_IN(0)-ZLAT_IN(JINY-1)
  ZLON_IN(1-JINY)        = ZLON_IN(JINY)&
                          &+MOD(ILON_IN(1-JINY),2)*RPI/ILON_IN(1-JINY)
  ILON_IN (KLAT_IN+JINY) = ILON_IN (KLAT_IN+1-JINY)
  ZDLAT_IN(KLAT_IN+JINY) = ZDLAT_IN(KLAT_IN+1-JINY)
  ZDLON_IN(KLAT_IN+JINY) = ZDLON_IN(KLAT_IN+1-JINY)
  ZHLAT_IN(KLAT_IN+JINY) = ZHLAT_IN(KLAT_IN+1-JINY)
  ZHLON_IN(KLAT_IN+JINY) = ZHLON_IN(KLAT_IN+1-JINY)
  ZG_IN   (KLAT_IN+JINY) = ZG_IN   (KLAT_IN+1-JINY)
  ZLAT_IN (KLAT_IN+JINY) = 2._JPRB*ZLAT_IN(KLAT_IN)-ZLAT_IN(KLAT_IN-JINY)
  ZLON_IN (KLAT_IN+JINY) = ZLON_IN (KLAT_IN+1-JINY)&
                          &+MOD(ILON_IN(KLAT_IN+JINY),2)*RPI/ILON_IN(KLAT_IN+JINY)
ENDDO

!              1.1.2 Inverses
ZDLONI_IN(:)=1./ZDLON_IN(:)
ZDLATI_IN(:)=1./ZDLAT_IN(:)
ZDHLAT_IN(:)=ZDLAT_IN(:)*ZHLAT_IN(:)
ZDHLATI_IN(:)=1./ZDHLAT_IN(:)

!*         1.2 Primitive function of input field with respect to longitude
!              Valid at grid-cell longitude boundaries.
!              Make an extension zone of +/-3 grid-cells for interpolation.
!              The extension over the poles shifts fields by RPI.

!          1.2.1 Primitive function
ALLOCATE(ZPRIMX(-2:JXMAX_IN+3,-2:KLAT_IN+3))
ALLOCATE(ZDIFX (-2:JXMAX_IN+2,-2:KLAT_IN+3))
CALL GSTATS(1943,1)

CALL GSTATS(1470,0)
!$OMP PARALLEL DO SCHEDULE(STATIC,1) PRIVATE(JLAT,ZSUM,JLON)
DO JLAT=1,KLAT_IN
  ZSUM=0._JPRB
  DO JLON=1,ILON_IN(JLAT)
    ZSUM=ZSUM+ZFLD_IN(IND_IN(JLAT)-1+JLON)*ZG_IN(JLAT)
    ZPRIMX(JLON,JLAT)=ZSUM
  ENDDO
ENDDO
!$OMP END PARALLEL DO
CALL GSTATS(1470,1)

!          1.2.2 The extension over the poles shifts fields by RPI.

CALL GSTATS(1943,0)
DO JINY = 1, 3
  ZSUM=0._JPRB
  DO JLON=1,ILON_IN(1-JINY)
    JINX = 1 + MOD( ILON_IN(1-JINY)/2 + JLON - 1, ILON_IN(1-JINY) )
    ZSUM=ZSUM+ZPARITY*ZFLD_IN(IND_IN(JINY)-1+JINX)*ZG_IN(1-JINY)
    ZPRIMX(JLON,1-JINY) = ZSUM
  ENDDO
  ZSUM=0.
  DO JLON=1,ILON_IN(KLAT_IN+JINY)
    JINX = 1+MOD( ILON_IN(KLAT_IN+JINY)/2 + JLON - 1, ILON_IN(KLAT_IN+JINY) )
    ZSUM=ZSUM+ZPARITY*ZFLD_IN(IND_IN(KLAT_IN+1-JINY)-1+JINX)*ZG_IN(KLAT_IN+JINY)
    ZPRIMX(JLON,KLAT_IN+JINY) = ZSUM
  ENDDO
ENDDO

!          1.2.3 Extension in x-direction

DO JINX = 1, 3
  DO JLAT=-2,KLAT_IN+3
    ZPRIMX(1-JINX            ,JLAT) =-ZPRIMX(ILON_IN(JLAT),JLAT)&
                                    &+ZPRIMX(ILON_IN(JLAT)+1-JINX,JLAT)
    ZPRIMX(ILON_IN(JLAT)+JINX,JLAT) = ZPRIMX(ILON_IN(JLAT),JLAT)&
                                    &+ZPRIMX(JINX,JLAT)
  ENDDO
ENDDO
CALL GSTATS(1943,1)

!              1.2.5 Difference of primitive function
CALL GSTATS(1470,0)
!$OMP PARALLEL DO SCHEDULE(STATIC,1) PRIVATE(JLAT,JLON)
DO JLAT=-2,KLAT_IN+3
  DO JLON=-2,ILON_IN(JLAT)+2
    ZDIFX(JLON,JLAT)=ZPRIMX(JLON+1,JLAT)-ZPRIMX(JLON,JLAT)
  ENDDO
ENDDO
!$OMP END PARALLEL DO
CALL GSTATS(1470,1)

!     *    2. OUTPUT GRID AND FIELD

!*         2.1 Geometry of output grid
!              Output grid can be finer or coarser than input grid.
!              Therefore dimensioning [KLAT_IN/KLAT_OUT]+4+1
!              to allow for 4th order interpolation of primitive
!              function wrt latitude at top and bottom of grid-cell
!              boundaries (output latitudes).

!          2.1.1 Grid blocks with constant longitude spacing

CALL GSTATS(1943,0)
IBL_OUT = 1
IBLATSTA_OUT(1) = 1
IBLON_OUT(1) = KLON_OUT(1)
DO JLAT=2,KLAT_OUT
  IF((KLON_OUT(JLAT)/=KLON_OUT(JLAT-1)) .OR. (JLAT == KLAT_OUT/2))THEN
    IBL_OUT = IBL_OUT + 1
    IBLATEND_OUT(IBL_OUT-1) = JLAT-1
    IBLATSTA_OUT(IBL_OUT) = JLAT
    IBLON_OUT(IBL_OUT) = KLON_OUT(JLAT)
  ENDIF
ENDDO
IBLATEND_OUT(IBL_OUT)=KLAT_OUT

!          2.1.2 Latitudes and scale factors

ZLAT_OUT(0)       = SIGN(RPI/2._JPRB,PLAT_OUT(1))
ZLAT_OUT(KLAT_OUT)= SIGN(RPI/2._JPRB,PLAT_OUT(KLAT_OUT))
DO JLAT=1,KLAT_OUT
  IF (JLAT<KLAT_OUT) ZLAT_OUT(JLAT)=0.5_JPRB*(PLAT_OUT(JLAT)+PLAT_OUT(JLAT+1))
  ZDLAT_OUT(JLAT)=ABS(ZLAT_OUT(JLAT-1)-ZLAT_OUT(JLAT))
  ZDLON_OUT(JLAT)=2._JPRB*RPI/KLON_OUT(JLAT)
!              Starting longitude (must be less or equal to zero)
  ZLON_OUT(JLAT)=-ZDLON_OUT(JLAT)/2._JPRB
!              Scale factors
  ZHLON_OUT(JLAT) = ABS(SIN(ZLAT_OUT(JLAT-1))-SIN(ZLAT_OUT(JLAT)))/ZDLAT_OUT(JLAT)
  ZHLAT_OUT(JLAT) = 1._JPRB
  ZG_OUT(JLAT) = ZHLON_OUT(JLAT)*ZHLAT_OUT(JLAT)*ZDLON_OUT(JLAT)*ZDLAT_OUT(JLAT)
ENDDO

!          2.1.3 Closest input row above each output row

ILAT_IO(0)=0
DO JLAT=1,KLAT_OUT
  ILAT_IO(JLAT)=ILAT_IO(JLAT-1)
  DO WHILE (ZLAT_IN(ILAT_IO(JLAT)) > ZLAT_OUT(JLAT))
    ILAT_IO(JLAT)=ILAT_IO(JLAT)+1
  ENDDO
  ILAT_IO(JLAT)=ILAT_IO(JLAT)-1
ENDDO
JXMAX_OUT=0
JYMAX_OUT=0
IOFF=0
DO JBL=1, IBL_OUT
  JXMAX_OUT = MAX(JXMAX_OUT,IBLON_OUT(JBL))
  JYMAX_OUT = MAX(JYMAX_OUT,IBLATEND_OUT(JBL)-IBLATSTA_OUT(JBL)+1)
  IOFFB(JBL)=IOFF
  DO JOUTY = 1, IBLATEND_OUT(JBL)-IBLATSTA_OUT(JBL)+1
    IOFF=IOFF+IBLON_OUT(JBL)
  ENDDO
ENDDO
IOFFB(IBL_OUT+1)=IOFF
CALL KEYSORT(IRET,IOFFB(2:IBL_OUT+1),IBL_OUT,INDEX=INDEXB(1:IBL_OUT),&
 & INIT=.TRUE.,DESCENDING=.TRUE.)

CALL GSTATS(1943,1)
  
!*         2.2 Interpolation of primitive function

CALL GSTATS(1470,0)
JNUM=0
!$OMP PARALLEL &
!$OMP&PRIVATE(JBL,IBL,JINYSTA,JINYEND,JINYD,JLON,JINY,ZOUT,ZPRIMADD,JINX,ZX, &
!$OMP& ZDLL,ZDL,ZDC,ZDR,ZDRR,ZC,ZLIMR,ZLIML,ZLIM,ZDHI, &
!$OMP& JOUTY,JLAT,ZY,ZY2,JNUM, &
!$OMP& ZPRIM,ZDIFY,ZPRIMY,ZOUTPRE)
ALLOCATE(ZPRIM(0:JXMAX_OUT,-3:3+KLAT_IN))
ALLOCATE(ZDIFY(1:JXMAX_OUT,-3:3+KLAT_IN))
ALLOCATE(ZPRIMY(0:JXMAX_OUT,0:JYMAX_OUT))
ALLOCATE(ZOUTPRE(0:JXMAX_OUT))
!$OMP DO SCHEDULE(STATIC,1)
DO JBL=1, IBL_OUT
  IBL=INDEXB(JBL)

!          2.2.1 Interpolate PF of input field with respect to longitude
!                to output grid-cell longitude boundaries.
!                (output longitude / input latitude)

  JINYSTA=ILAT_IO(IBLATSTA_OUT(IBL)-1)
  JINYEND=ILAT_IO(IBLATEND_OUT(IBL))
  JINYD=JINYEND-JINYSTA
  DO JLON=0,IBLON_OUT(IBL)
    ZPRIM(JLON,-3)=0.
    ZOUTPRE(JLON)=ZLON_OUT(IBLATSTA_OUT(IBL))+JLON*ZDLON_OUT(IBLATSTA_OUT(IBL))
  ENDDO
  DO JINY=JINYSTA-2,JINYEND+3
    DO JLON=0,IBLON_OUT(IBL)
      ZOUT =  ZOUTPRE(JLON)-ZLON_IN(JINY)
      ZPRIMADD=0._JPRB
      IF (ZOUT<0._JPRB) THEN
        ZOUT    = ZOUT+ILON_IN(JINY)*ZDLON_IN(JINY)
        ZPRIMADD=-ZPRIMX(ILON_IN(JINY),JINY)
      ENDIF
      IF (ZOUT>ILON_IN(JINY)*ZDLON_IN(JINY)) THEN
        ZOUT    = ZOUT-ILON_IN(JINY)*ZDLON_IN(JINY)
        ZPRIMADD= ZPRIMX(ILON_IN(JINY),JINY)
      ENDIF
      JINX=INT(ZOUT*ZDLONI_IN(JINY))
      ZX=-1._JPRB+2._JPRB*(ZOUT-JINX*ZDLON_IN(JINY))*ZDLONI_IN(JINY)
!          2.2.1.1 Basic coefficients
      ZDLL=ZDIFX(JINX-2,JINY)
      ZDL =ZDIFX(JINX-1,JINY)
      ZDC =ZDIFX(JINX  ,JINY)
      ZDR =ZDIFX(JINX+1,JINY)
      ZDRR=ZDIFX(JINX+2,JINY)
      ZC(0) = (ZPRIMX(JINX+1,JINY)+ZPRIMX(JINX  ,JINY))*.5
      ZC(1) = ZDC*.5
!          2.2.1.2 Economization of the power series (5th to 3rd order)
      ZC(2)=(130.*(ZDR-ZDL)-17.*(ZDRR-ZDLL))*ZCC2
      ZC(3)=(452.*(ZDR-2.*ZDC+ZDL)-33.*(ZDRR-2.*ZDC+ZDLL))*ZCC3
!          2.2.1.3 Limiting
      ZC(3) = ZCURVLIM(ZC(3),ZC(2))
      ZLIMR = ZMINMOD( ZFRAC( (ZDR-ZDC)*.5,       2.*( ZC(2)+ZC(3) ) ) )
      ZLIML = ZMINMOD( ZFRAC( (ZDL-ZDC)*.5, ZLIMR*2.*(-ZC(2)+ZC(3) ) ) )
      ZLIM  = ZLIMR*ZLIML
      ZC(2) = ZLIM*ZC(2)
      ZC(3) = ZLIM*ZC(3)
!          2.2.1.5 Evaluation
      ZPRIM(JLON,JINY-JINYSTA)= ZC(0)+ZC(1)*ZX   &
                 & +(ZX*ZX-1.)*(ZC(2)+ZC(3)*ZX)   &
                 & +ZPRIMADD
    ENDDO
  ENDDO

!          2.2.2  Differentiate wrt longitude

  DO JINY=-2,3+JINYD
    ZDHI = 1./(ZDLON_OUT(IBLATSTA_OUT(IBL))*ZHLON_IN(JINYSTA+JINY))
    DO JLON=IBLON_OUT(IBL),1,-1
      ZPRIM(JLON,JINY) = (ZPRIM(JLON,JINY)-ZPRIM(JLON-1,JINY))*ZDHI
    ENDDO
  ENDDO

!          2.2.3 PF wrt latitude and scalefactor*cellarea (G=h_y, dA=dy ==> G*dA)
!                (output longitude / input latitude)

  DO JINY=-2,3+JINYD
    DO JLON=1, IBLON_OUT(IBL)
      ZPRIM(JLON,JINY)= ZPRIM(JLON,JINY-1)+ZPRIM(JLON,JINY)
    ENDDO
  ENDDO
  DO JINY=-2,2+JINYD
    DO JLON=1, IBLON_OUT(IBL)
      ZDIFY(JLON,JINY) = (ZPRIM(JLON,JINY+1)-ZPRIM(JLON,JINY))*ZDHLATI_IN(JINYSTA+JINY+1)
    ENDDO
  ENDDO

!          2.2.4 Interpolation PF to output latitude cell boundaries

  DO JOUTY = 0, IBLATEND_OUT(IBL)-IBLATSTA_OUT(IBL)+1
    JLAT =         JOUTY+IBLATSTA_OUT(IBL)-1
    JINY = ILAT_IO(JOUTY+IBLATSTA_OUT(IBL)-1)-JINYSTA
    ZY=-1.+2.*(-ZLAT_OUT(JLAT)+ZLAT_IN(JINYSTA+JINY))*ZDLATI_IN(JINYSTA+JINY+1)
    ZY2 = (ZY*ZY-1.)
    DO JLON=1, IBLON_OUT(IBL)
!          2.2.4.1 Basic coefficients
      ZDLL=ZDIFY(JLON,JINY-2)
      ZDL =ZDIFY(JLON,JINY-1)
      ZDC =ZDIFY(JLON,JINY  )
      ZDR =ZDIFY(JLON,JINY+1)
      ZDRR=ZDIFY(JLON,JINY+2)
      ZC(0) = (ZPRIM(JLON,JINY+1)+ZPRIM(JLON,JINY))*.5
      ZC(1) = ZDC*.5
!          2.2.1.2 Economization of the power series (5th to 3rd order)
      ZC(2)=(130.*(ZDR-ZDL)-17.*(ZDRR-ZDLL))*ZCC2
      ZC(3)=(452.*(ZDR-2.*ZDC+ZDL)-33.*(ZDRR-2.*ZDC+ZDLL))*ZCC3
!          2.2.1.3 Limiting
      ZC(3) = ZCURVLIM(ZC(3),ZC(2))
      ZLIMR = ZMINMOD( ZFRAC( (ZDR-ZDC)*.5,       2.*( ZC(2)+ZC(3) ) ) )
      ZLIML = ZMINMOD( ZFRAC( (ZDL-ZDC)*.5, ZLIMR*2.*(-ZC(2)+ZC(3) ) ) )
      ZLIM  = ZLIMR*ZLIML
      ZC(2) = ZLIM*ZC(2)
      ZC(3) = ZLIM*ZC(3)
!          2.2.5.4 Back to real coefficients
      ZC(1:3) = ZDHLAT_IN(JINYSTA+JINY+1)*ZC(1:3)
!          2.2.4.5 Evaluation
      ZPRIMY(JLON,JOUTY)= ZC(0)+ZC(1)*ZY+ZY2*(ZC(2)+ZC(3)*ZY)
    ENDDO
  ENDDO

!          2.2.5 Take difference of primitive function to get output field
!                (output longitude / output latitude)

  JNUM=IOFFB(IBL)
  DO JOUTY = 1, IBLATEND_OUT(IBL)-IBLATSTA_OUT(IBL)+1
    JLAT = JOUTY+IBLATSTA_OUT(IBL)-1
    ZDHI = 1./(ZDLAT_OUT(JLAT)*ZHLAT_OUT(JLAT))
    DO JLON=1, IBLON_OUT(IBL)
      JNUM=JNUM+1
      IF (LLCONS3D) THEN
!!        PFLD_OUT(JNUM) =  ((ZPRIMY(JLON,JOUTY)-ZPRIMY(JLON,JOUTY-1))*ZDHI) &
!!         & /( (YRVAB%VAH(KLEVEL)-YRVAB%VAH(KLEVEL-1)) &
!!         & +(YRVAB%VBH(KLEVEL)-YRVAB%VBH(KLEVEL-1))*PSP_OUT(JNUM) )
        PFLD_OUT(JNUM) =  ((ZPRIMY(JLON,JOUTY)-ZPRIMY(JLON,JOUTY-1))*ZDHI) &
         & /( (PVAH(KLEVEL)-PVAH(KLEVEL-1)) &
         & +(PVBH(KLEVEL)-PVBH(KLEVEL-1))*PSP_OUT(JNUM) )
      ELSE
        PFLD_OUT(JNUM) =  (ZPRIMY(JLON,JOUTY)-ZPRIMY(JLON,JOUTY-1))*ZDHI
      ENDIF
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
IF (ALLOCATED(ZPRIM ))  DEALLOCATE(ZPRIM )
IF (ALLOCATED(ZPRIMY))  DEALLOCATE(ZPRIMY)
IF (ALLOCATED(ZDIFY))   DEALLOCATE(ZDIFY)
IF (ALLOCATED(ZOUTPRE)) DEALLOCATE(ZOUTPRE)
!$OMP END PARALLEL
CALL GSTATS(1470,1)

IF (ALLOCATED(ZPRIMX)) DEALLOCATE(ZPRIMX)
IF (ALLOCATED(ZFLD_IN)) DEALLOCATE(ZFLD_IN)
IF (ALLOCATED(ZDIFX)) DEALLOCATE(ZDIFX)


IF (LHOOK) CALL DR_HOOK('GRID_BICONSERV',1,ZHOOK_HANDLE)

END SUBROUTINE GRID_BICONSERV
