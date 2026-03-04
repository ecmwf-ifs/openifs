! (C) Copyright 1989- Meteo-France.

SUBROUTINE SUFPG2(YDNAMFPINT,YDNAMFPSCI,KFPDOM,YDFPUSERGEO,YDGEOMETRY,KFPRGP,PCO_EZO,LDPOSHOR,PFPLA,PFPLO,PFPGM,PFPNORX,PFPNORY, &
 & KNUMD,PFPMS)

!**** *SUFPG2*  - INITIALIZE OUTPUT GEOMETRY ARRAYS (Full-POS)

!     PURPOSE.
!     --------
!        COMPUTES THE COORDINATES OF OUTPUT POINTS IN THE GEOMETRY OF THE
!        MODEL, THE OUTPUT MAP FACTOR AND THE RELATIVE ROTATION MATRIX.
!        WHEN AN OUTPUT POINT IS AT A GEOGRAPHICAL POLE,
!        THE EQUATION OF THE COMPASS IS DEGENERATING,
!        SO A DIRECT (AND EXACT) COMPUTATION IS NEEDED.

!**   INTERFACE. *CALL* *SUFPG2*
!     ----------

!        EXPLICIT ARGUMENTS :
!        ------------------
!         All dummy arguments are DM-global.

!         KFPDOM           : number of subdomains
!         KFPRGP           : number of points on the (complete) output grid(s)
!         PCO_EZO          : number of rows to shift the virtual extension zone
!                            from the limits toward the center of the input domain
!         LDPOSHOR         : .TRUE. if actual horizontal processing
!         PFPLA, PFPLO     : relative locations (latitudes/longitudes or x/y)
!                            of output points
!         PFPGM            : output map factor
!         PFPNORX, PFPNORY : relative rotation matrix
!         KNUMD            : Subdomain index
!         PFPMS            : Resolution of the output grids (mean resolution of the grid devided by the map factor), in meters.

!        IMPLICIT ARGUMENTS
!        ------------------
!        PARDIM, YOMDIM
!        PARFPOS
!        YOMFPC
!        YOMFPG, NAMFPG
!        YOMGEM
!        YOMCST

!     EXTERNALS. CORDON ECORDON TRARECA SUGENORD ELALO2XY
!     ----------

!     Called by SUFPG.

!     AUTHOR.    RYAD EL KHATIB *METEO-FRANCE*
!     -------
!      970304 ORIGINAL from previous SUFPG

!     MODIFICATIONS.
!     --------------
!      01-04-14 M. Janousek: Replacement of EGGX by EGGPACK functions
!      R. El Khatib : 01-08-07 Pruning options
!      O.Spaniel    : 03-04-15 cleaning
!      JD. Gril 08-02-2005 : Modifs for Mercator Rotated Tilted
!                            cleaning NSOTRP, test FPRPK
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      Y.Seity       19-Nov-2007 Initialisation of LAT1 LAT2 LON1 LON2 before calling EGGX_n
!      K. Yessad (Jan 2010): externalisation of group EGGX in XRD/IFSAUX
!      R. El Khatib : 31-Jul-2012 LFPNEWGRID instead of LDSAME
!      R. El Khatib 31-Aug-2012 new E-zone management
!      R. El Khatib 27-Sep-2013 Boyd periodization in Fullpos-2
!      R. El Khatib 04-Aug-2014 Pruning of the conf. 927/928
!      R. El Khatib 27-Jul-2016 minor bugfix for the case NFPBOYD ; recode LWIDER_DOM
!      R. El Khatib 11-Dec-2018 : fix virtual location of E-zone if spline method used for biperiodicization
!-----------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMLUN   , ONLY : NULOUT
USE YOMCT0   , ONLY : LELAM    ,LRPLANE
USE YOMFPC        , ONLY : TNAMFPINT, TNAMFPSCI, LTRACEFP
USE TYPE_FPUSERGEO, ONLY : TFPUSERGEO
USE YOMCST   , ONLY : RPI      ,RA
USE MPL_MODULE, ONLY: MPL_BROADCAST

IMPLICIT NONE

TYPE (TNAMFPINT),  INTENT(IN) :: YDNAMFPINT
TYPE (TNAMFPSCI),  INTENT(IN) :: YDNAMFPSCI
INTEGER(KIND=JPIM),INTENT(IN)    :: KFPDOM
TYPE (TFPUSERGEO), INTENT(IN)    :: YDFPUSERGEO(KFPDOM)
TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
INTEGER(KIND=JPIM),INTENT(IN)    :: KFPRGP
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCO_EZO
LOGICAL           ,INTENT(OUT)   :: LDPOSHOR
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFPLA(KFPRGP)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFPLO(KFPRGP)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFPGM(KFPRGP)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFPNORX(KFPRGP)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFPNORY(KFPRGP)
INTEGER(KIND=JPIM),INTENT(OUT)   :: KNUMD(KFPRGP)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFPMS(KFPRGP)

!-----------------------------------------------------------------------

!     ZGELAT, ZGELAM   : real latitudes and longitudes in radians

!     ZGESLAT          : sine of the real latitudes
!     ZGESLAM, ZGECLAM : sine and cosine of the real longitudes

!     ZSLAT            : sine of the relative latitudes
!     ZSLAM, ZCLAM     : sine and cosine of the relative longitudes

!     ZFPMU            : sine of latitudes on the "work" pp. sphere
!     ZFPSLO, ZFPCLO   : sine and cosine of longitudes on the "work" pp. sphere

!     ZFPNX2, ZFPNY2   : geographical compass, from the output geometry.
!     ZFPNX1, ZFPNY1   : geographical compass, from the input geometry.

REAL(KIND=JPRB), ALLOCATABLE :: ZGNORX(:), ZGNORY(:)
REAL(KIND=JPRB), ALLOCATABLE :: ZGELAMG(:,:), ZGELATG(:,:)
REAL(KIND=JPRB), POINTER :: ZPTRGELAM(:,:,:), ZPTRGELAT(:,:,:)
REAL(KIND=JPRB) :: ZGELAT(KFPRGP), ZGELAM(KFPRGP)
REAL(KIND=JPRB) :: ZGESLAT(KFPRGP), ZGESLAM(KFPRGP), ZGECLAM(KFPRGP)
REAL(KIND=JPRB), ALLOCATABLE :: ZSLAT(:), ZSLAM(:), ZCLAM(:)
REAL(KIND=JPRB), ALLOCATABLE :: ZFPMU(:), ZFPSLO(:), ZFPCLO(:)
REAL(KIND=JPRB), ALLOCATABLE :: ZFPNX1(:), ZFPNY1(:)
REAL(KIND=JPRB) :: ZFPNX2(KFPRGP), ZFPNY2(KFPRGP)

REAL(KIND=JPRB), ALLOCATABLE :: ZZLA(:), ZZLO(:)

INTEGER(KIND=JPIM) :: IERR, IPTR, ISOTRP, ISTTYP, ITRIV, IPTR2, IFPGUX, IFPLUX, IOFF, ISTART
INTEGER(KIND=JPIM) :: J, JDOM, JGL, JI, JIN, JLON, JOUT, JJ, IGIVO, IST1, ILON, ILAT
INTEGER(KIND=JPIM) :: IFPLSUR, IFPGSUR, ISW, INE, ILEN
INTEGER(KIND=JPIM) :: ITO(1), ISHIFT

LOGICAL :: LLSUBSET, LLNEWGRID

REAL(KIND=JPRB) :: Z180SPI, Z2C, Z2PI, ZC2M1, ZC2P1, ZLATE, &
 & ZCAFPCLO(1), ZCAFPLO, ZCAFPSLA(1), ZCAFPSLO(1), ZCLACEN,&
 & ZCLOCEN, ZDLA, ZDLO, ZDUM,&
 & ZEPS, ZFCLACEN, ZFCLOCEN, ZFSLACEN,&
 & ZFSLOCEN, ZGEFPCLO(1), ZGEFPSLA(1), ZGEFPSLO(1), ZGELOCEN,&
 & ZGEMUCEN, ZLONW, ZLONE, ZLATS, ZLATN, ZLON0, ZLAT0, ZLONC, ZLATC,&
 & ZDELX, ZDELY, ZLON, ZPIS180, ZPISNL,&
 & ZRELA, ZRELO, ZSIGN, ZSLACEN, ZSLOCEN, ZUS2C,&
 & ZLON1, ZLAT1, ZLONR, ZLATR, ZBETA, ZRPK, ZZ, ZFPMS
! ZFPMS : mean resolution of the output grids (squared average between zonal and meridian resolutions), in meters

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!-----------------------------------------------------------------------

#include "eggx_n.h"
#include "fpbipere.h"
#include "egath_grid.h"
#include "abor1.intfb.h"

#include "elalo2xy.intfb.h"
#include "sugenord.intfb.h"





!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUFPG2',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDGEM=>YDGEOMETRY%YRGEM, &
 & YDEGEO=>YDGEOMETRY%YREGEO)
ASSOCIATE(LFPML_STD=>YDNAMFPINT%LFPML_STD, LFPML_LAN=>YDNAMFPINT%LFPML_LAN, &
 & LFPML_SEA=>YDNAMFPINT%LFPML_SEA, &
 & NFPBOYD=>YDNAMFPSCI%NFPBOYD, &
 & NGPTOT=>YDGEM%NGPTOT, NGPTOTG=>YDGEM%NGPTOTG, NSTTYP=>YDGEM%NSTTYP, RMUCEN=>YDGEM%RMUCEN, RSTRET=>YDGEM%RSTRET, &
 & RLOCEN=>YDGEM%RLOCEN, &
 & NDLON=>YDDIM%NDLON, NDGUXG=>YDDIM%NDGUXG, NDLUXG=>YDDIM%NDLUXG, NRESOL=>YDDIM%NRESOL, &
 & NLAT=>YDFPUSERGEO%NLAT, NLON=>YDFPUSERGEO%NLON, NFPLUX=>YDFPUSERGEO%NFPLUX, &
 & NFPGUX=>YDFPUSERGEO%NFPGUX, RLONC=>YDFPUSERGEO%RLONC, RLATC=>YDFPUSERGEO%RLATC, &
 & RDELX=>YDFPUSERGEO%RDELX, RDELY=>YDFPUSERGEO%RDELY, NFPBWX=>YDFPUSERGEO%NFPBWX, &
 & NFPBWY=>YDFPUSERGEO%NFPBWY, CFPGRID=>YDFPUSERGEO%CFPGRID,&
 & NFPTTYP=>YDFPUSERGEO%NFPTTYP, FPMUCEN=>YDFPUSERGEO%FPMUCEN, &
 & FPLOCEN=>YDFPUSERGEO%FPLOCEN, LFPMAP=>YDFPUSERGEO%LFPMAP, FPSTRET=>YDFPUSERGEO%FPSTRET, &
 & FPLON0=>YDFPUSERGEO%FPLON0, FPLAT0=>YDFPUSERGEO%FPLAT0, &
 & LFPMRT=>YDFPUSERGEO%LFPMRT, LFPBIPER=>YDFPUSERGEO%LFPBIPER, &
 & ELON0=>YDEGEO%ELON0, ELON1=>YDEGEO%ELON1, ELON2=>YDEGEO%ELON2, ELONC=>YDEGEO%ELONC, &
 & ELAT0=>YDEGEO%ELAT0, ELAT1=>YDEGEO%ELAT1, ELAT2=>YDEGEO%ELAT2, ELATC=>YDEGEO%ELATC, &
 & EDELX=>YDEGEO%EDELX, EDELY=>YDEGEO%EDELY, LMRT=>YDEGEO%LMRT, LMAP=>YDEGEO%LMAP)

!*       0. PREPARATIONS
!           ------------

Z2PI=2.0_JPRB*RPI
ZPIS180=RPI/180._JPRB
Z180SPI=1.0_JPRB/ZPIS180

!*       1. GEOGRAPHICAL COORDINATES, OUTPUT MAP FACTOR AND GEOGRAPHICAL
!           COMPASS, SEEN FROM OUTPUT GEOMETRY
!           ------------------------------------------------------------

IF (CFPGRID(1) == 'LALON') THEN
  IPTR = 0
  DO J = 1, KFPDOM

    ZLATS=RLATC(J)-RDELY(J)*0.5_JPRB*REAL(NLAT(J)-1,JPRB)
    ZLONW=MOD(RLONC(J)-RDELX(J)*0.5_JPRB*REAL(NLON(J)-1,JPRB),360._JPRB)
    IF (ZLONW >= 360._JPRB) ZLONW=ZLONW-360._JPRB

    ZLATS = ZLATS * ZPIS180
    ZLONW = ZLONW * ZPIS180
    ZDELX = RDELX(J) * ZPIS180
    ZDELY = RDELY(J) * ZPIS180
    ZFPMS=(RPI*RA/180._JPRB)*SQRT(0.5_JPRB*(RDELX(J)*RDELX(J)+RDELY(J)*RDELY(J)))
    DO JGL = 1, NLAT(J)
      IOFF=NLON(J)*(JGL-1)
      DO JLON = 1, NLON(J)
        ZZ = ZLONW+MOD(REAL(JLON-1,JPRB)*ZDELX,Z2PI)
        IF (ZZ >= Z2PI) ZZ = ZZ-Z2PI
        ZGELAM(JLON+IOFF+IPTR) = ZZ
        ZGESLAM(JLON+IOFF+IPTR) = SIN(ZZ)
        ZGECLAM(JLON+IOFF+IPTR) = COS(ZZ)
        ZZ=ZLATS+REAL(JGL-1,JPRB)*ZDELY
        ZGELAT(JLON+IOFF+IPTR) = ZZ
        ZGESLAT(JLON+IOFF+IPTR) = SIN(ZZ)
        PFPGM (JLON+IOFF+IPTR) = 1.0_JPRB
        IF (LFPML_STD.OR.LFPML_LAN.OR.LFPML_SEA) THEN
          PFPMS (JLON+IOFF+IPTR) = ZFPMS
        ENDIF
        ZFPNX2(JLON+IOFF+IPTR) = 0.0_JPRB
        ZFPNY2(JLON+IOFF+IPTR) = 1.0_JPRB
      ENDDO
    ENDDO
    KNUMD(IPTR+1:IPTR+NLAT(J)*NLON(J))=J
    IPTR=IPTR+NLAT(J)*NLON(J)
    IF (LTRACEFP) THEN
      WRITE(UNIT=NULOUT,FMT='(''ZFPMS('',I2.2,'') = '',F10.2)') J, ZFPMS
    ENDIF
  ENDDO

ELSEIF (CFPGRID(1) == 'LELAM') THEN

  IPTR = 1
  ZLON0 = FPLON0(1) * ZPIS180
  ZLAT0 = FPLAT0(1) * ZPIS180
  ZLONR=0._JPRB
  ZLATR=0._JPRB
  ZLATS=0._JPRB
  ZLATN=0._JPRB
  ZLONW=0._JPRB
  ZLONE=0._JPRB
  ZBETA=0._JPRB
  IGIVO=0_JPIM
  ISOTRP=0_JPIM
  IST1=1
  JJ=-1
  IF (LFPMRT(1)) JJ=JJ-1
  DO J = 1, KFPDOM
    ZLONC = RLONC(J) * ZPIS180
    ZLATC = RLATC(J) * ZPIS180
    ZDELX = RDELX(J)
    ZDELY = RDELY(J)
!   Interpolation area :
    IF (NFPBOYD /= 0) THEN
      IFPLUX=NFPLUX(J)+2*NFPBWX(J)+2*(NLON(J)-NFPLUX(J))
      IFPGUX=NFPGUX(J)+2*NFPBWY(J)+2*(NLAT(J)-NFPGUX(J))
      IFPLSUR=IFPLUX
      IFPGSUR=IFPGUX
    ELSE
      IFPLUX=NFPLUX(J)
      IFPGUX=NFPGUX(J)
      IFPLSUR=NLON(J)
      IFPGSUR=NLAT(J)
    ENDIF
    WRITE(NULOUT,*) 'Call EGGX_N by SUFPG2'
!   Computation will return the locations of points on the interpolation area
    IF (LFPMAP(J)) THEN
      ZRPK=HUGE(1._JPRB)
      ! will be computed by eggx
    ELSE
      ZRPK=-9._JPRB
    ENDIF
    CALL EGGX_N(RPI,RA,JJ,ZLONR,ZLATR,ZBETA,ZLONW,ZLATS,&
     & ZLONE,ZLATN,ZLON0,ZLAT0,ZRPK,NULOUT,ISOTRP,IGIVO,&
     & ZGELAM(IPTR),ZGELAT(IPTR),PFPGM(IPTR),ZFPNX2(IPTR),&
     & ZFPNY2(IPTR),IST1,IFPLSUR,IST1,IFPGSUR,IST1,IFPLUX,&
     & IST1,IFPGUX,ZDELX,ZDELY,ZLONC,ZLATC)
    IF (NFPBOYD /= 0) THEN
      ISW=((NFPBWY(J)+NLAT(J)-NFPGUX(J))*IFPLUX)+NFPBWX(J)+NLON(J)-NFPLUX(J)+1
      INE=((NFPBWY(J)+NLAT(J))*IFPLUX)-NFPBWX(J)-(NLON(J)-NFPLUX(J))
    ELSE
      ISW=1
      INE=(NFPGUX(J)*NFPLUX(J))+(NFPGUX(J)-1)*(NLON(J)-NFPLUX(J))
    ENDIF
    WRITE(NULOUT,'('' SOUTH-WEST LATITUDE   = '',E21.15,'' NORTH-EAST LATITUDE   = '',E21.15)') ZGELAT(ISW), ZGELAT(INE)
    WRITE(NULOUT,'('' SOUTH-WEST LONGITUDE  = '',E21.15,'' NORTH-EAST LONGITUDE  = '',E21.15)') ZGELAM(ISW), ZGELAM(INE)
    WRITE(NULOUT,'('' SOUTH-WEST MAP FACTOR = '',E21.15,'' NORTH-EAST MAP FACTOR = '',E21.15)') PFPGM(ISW),  PFPGM(INE)
    WRITE(NULOUT,'('' SOUTH-WEST COMPASS/X  = '',E21.15,'' NORTH-EAST COMPASS/X  = '',E21.15)') ZFPNX2(ISW), ZFPNX2(INE)
    WRITE(NULOUT,'('' SOUTH-WEST COMPASS/Y  = '',E21.15,'' NORTH-EAST COMPASS/Y  = '',E21.15)') ZFPNY2(ISW), ZFPNY2(INE)
    ZLATN = ZLATN * Z180SPI
    ZLATS = ZLATS * Z180SPI
    ZLONE = ZLONE * Z180SPI
    ZLONW = ZLONW * Z180SPI
    IF (ZLONE < ZLONW) THEN
      ZLONW=ZLONW-360._JPRB
    ENDIF
    IF(LTRACEFP) THEN
      WRITE(UNIT=NULOUT,FMT='('' Angles RLATN, RLATS, RLONW, RLONE, are in degrees '')')
      WRITE(UNIT=NULOUT,FMT='(''  ( J    RLATN   RLATS   RLONW   RLONE  )'')')
      WRITE(UNIT=NULOUT,FMT='('' ('',I2,1X,4F8.2,'')'')') J,ZLATN,ZLATS,ZLONW,ZLONE
    ENDIF
    IF (NFPBOYD == 0) THEN
      ! Fill the extension zone with appropriate geographical values (spline method)
      ! Set the coordinates of E-zone points at the same location but
      ! as far as possible from the target C+I domain. Shift this location toward the inside by 1 row
      ! for numerical security.
      IF (LELAM.AND.LFPBIPER(J)) THEN
        ALLOCATE(ZGELAMG(NGPTOTG,1))
        ALLOCATE(ZGELATG(NGPTOTG,1))
        ITO(1)=1
        ZPTRGELAM(1:NGPTOT,1:1,1:1)=>YDGEOMETRY%YRGSGEOM_NB%GELAM(1:NGPTOT)
        ZPTRGELAT(1:NGPTOT,1:1,1:1)=>YDGEOMETRY%YRGSGEOM_NB%GELAT(1:NGPTOT)
        CALL EGATH_GRID(PGPG=ZGELAMG,KRESOL=NRESOL,KFGATHG=1,KTO=ITO,PGP=ZPTRGELAM)
        CALL EGATH_GRID(PGPG=ZGELATG,KRESOL=NRESOL,KFGATHG=1,KTO=ITO,PGP=ZPTRGELAT)
        ! PCO_EZO should not be helpful anymore. For security however, it is re-used if negative.
        ! Then its absolute value is the number of rows the E-zone is shifted toward the center
        ! of the input domain (this is exactly what I wanted to do !!!). REK
        IF (PCO_EZO < 0._JPRB) THEN
          ISHIFT=MAX(1,MIN(MIN(-NINT(PCO_EZO),NDLUXG/2),NDGUXG/2))
        ELSE
          ISHIFT=1
        ENDIF
        IF (ISHIFT > 1) THEN
          WRITE(NULOUT,'(''VIRTUAL LOCATION OF E-ZONE WILL BE SHIFTED BY '',I5,'' ROWS '')') ISHIFT
        ENDIF
        IF (ZLATC > ELATC) THEN
          IF (ZLONC > ELONC) THEN
            ! South-West :
            ZLATE=ZGELATG(1+ISHIFT,1)
            ZLONE=ZGELAMG(1+ISHIFT,1)
          ELSE
            ! South-East :
            ZLATE=ZGELATG(NDLON+NDLUXG-ISHIFT,1)
            ZLONE=ZGELAMG(NDLON+NDLUXG-ISHIFT,1)
          ENDIF
        ELSE
          IF (ZLONC > ELONC) THEN
            ! North-West :
            ZLATE=ZGELATG(NDLON*(NDGUXG-1-ISHIFT)+1+ISHIFT,1)
            ZLONE=ZGELAMG(NDLON*(NDGUXG-1-ISHIFT)+1+ISHIFT,1)
          ELSE
            ! North-East :
            ZLATE=ZGELATG(NDLON*(NDGUXG-1-ISHIFT)+NDLUXG-ISHIFT,1)
            ZLONE=ZGELAMG(NDLON*(NDGUXG-1-ISHIFT)+NDLUXG-ISHIFT,1)
          ENDIF
        ENDIF
        DEALLOCATE(ZGELAMG)
        DEALLOCATE(ZGELATG)
        CALL MPL_BROADCAST(ZLONE,KTAG=1,KROOT=1,CDSTRING='SUFPG2:')
        CALL MPL_BROADCAST(ZLATE,KTAG=1,KROOT=1,CDSTRING='SUFPG2:')
      ELSE
        ZLONE=MOD(ZLONC+RPI,2._JPRB*RPI)
        ZLATE=-ZLATC
      ENDIF
      WRITE(NULOUT,'(''VIRTUAL LOCATION OF E-ZONE IN DEGREES : LONGITUDE = '',E15.7,'' LATITUDE = '',E15.7)') &
       & ZLONE*Z180SPI,ZLATE*Z180SPI
      DO JGL=1,NLAT(J)
        ISTART=1+NFPLUX(J)*MIN(1,NFPGUX(J)/JGL)
        IOFF=(JGL-1)*NLON(J)
        DO JLON=ISTART,IFPLSUR
          ZGELAM(IOFF+JLON)=ZLONE
          ZGELAT(IOFF+JLON)=ZLATE
          ZFPNX2(IOFF+JLON)=0._JPRB
          ZFPNY2(IOFF+JLON)=1._JPRB
        ENDDO
      ENDDO
    ENDIF
    ILEN=IFPLSUR*IFPGSUR
    IPTR2=IPTR-1 + ILEN
    IF (.NOT.(LFPMAP(J) .OR. (LELAM .AND. LMAP))) THEN
      PFPGM(IPTR:IPTR2)=1._JPRB
    ELSEIF (NFPBOYD == 0) THEN
      CALL FPBIPERE(NFPLUX(J),NFPGUX(J),IFPLSUR,IFPGSUR,1,ILEN,PFPGM(IPTR:IPTR2),0,LDZON=.TRUE.)
    ENDIF
    ZFPMS=SQRT(0.5_JPRB*(RDELX(J)*RDELX(J)+RDELY(J)*RDELY(J)))
    IF (LFPML_STD.OR.LFPML_LAN.OR.LFPML_SEA) THEN
      DO JI = IPTR, IPTR+ILEN-1
        PFPMS(JI)=ZFPMS/PFPGM(JI)
      ENDDO
    ENDIF
    KNUMD(IPTR:IPTR+ILEN-1)=J
    IPTR = IPTR + ILEN
    IF (LTRACEFP) THEN
      WRITE(UNIT=NULOUT,FMT='(''ZFPMS('',I2.2,'') = '',F10.2)') J, ZFPMS
    ENDIF
  ENDDO
  DO J = 1, IPTR-1
    ZGELAM(J) = MOD(ZGELAM(J),Z2PI)
  ENDDO

  DO J = 1, KFPRGP
    ZGESLAT(J) = SIN(ZGELAT(J))
    ZGESLAM(J) = SIN(ZGELAM(J))
    ZGECLAM(J) = COS(ZGELAM(J))
  ENDDO

ELSEIF (CFPGRID(1) == 'GAUSS') THEN

  ALLOCATE(ZFPMU(KFPRGP))
  ALLOCATE(ZFPSLO(KFPRGP))
  ALLOCATE(ZFPCLO(KFPRGP))

  Z2C = 2.0_JPRB* FPSTRET(1)
  ZUS2C = 1.0_JPRB/ Z2C
  ZC2M1 = FPSTRET(1) * FPSTRET(1) - 1.0_JPRB
  ZC2P1 = FPSTRET(1) * FPSTRET(1) + 1.0_JPRB
  ZFPMS = RPI*RA/REAL(NLAT(1),JPRB)
  IPTR = 0
  DO JGL = 1,NLAT(1)
    ZPISNL = Z2PI / YDFPUSERGEO(1)%NFPRGRI(JGL)
    DO JLON = 1, YDFPUSERGEO(1)%NFPRGRI(JGL)
      IPTR = IPTR+1
      ZFPMU(IPTR)=YDFPUSERGEO(1)%FPMU(JGL)
      PFPGM(IPTR)=(ZC2P1 + ZC2M1*ZFPMU(IPTR)) * ZUS2C
      IF (LFPML_STD.OR.LFPML_LAN.OR.LFPML_SEA) THEN
        PFPMS(IPTR)=ZFPMS/PFPGM(IPTR)
      ENDIF
      ZLON=REAL(JLON-1,JPRB)
      ZFPSLO(IPTR) = SIN(ZPISNL * ZLON)
      ZFPCLO(IPTR) = COS(ZPISNL * ZLON)
    ENDDO
  ENDDO
  KNUMD(1:KFPRGP)=1 ! jdom actually

  ZFSLACEN = FPMUCEN(1)
  ZFCLACEN = SQRT(1.0_JPRB - ZFSLACEN*ZFSLACEN)
  ZFCLOCEN = COS(FPLOCEN(1))
  ZFSLOCEN = SIN(FPLOCEN(1))

  !OIFS
  ZGESLAT  = ZFPMU
  ZGESLAM  = ZFPSLO
  ZGECLAM  = ZFPCLO
  !OIFS






  CALL SUGENORD(NFPTTYP(1),KFPRGP,FPSTRET(1),FPMUCEN(1),FPLOCEN(1),&
   & ZFPMU,ZFPSLO,ZFPCLO,ZFPNX2,ZFPNY2)

  DO J = 1, KFPRGP

    ZGELAT(J) = ASIN(MAX(-1.0_JPRB,MIN(1.0_JPRB,ZGESLAT(J))))
    ZGELAM(J) = MOD( ACOS(MAX(-1.0_JPRB,MIN(1.0_JPRB,ZGECLAM(J))))&
     & *SIGN(1.0_JPRB,ZGESLAM(J))&
     & - Z2PI*MIN(0.0_JPRB,SIGN(1.0_JPRB,ZGESLAM(J))), Z2PI)
  ENDDO
  IF (LTRACEFP) THEN
    WRITE(UNIT=NULOUT,FMT='(''ZFPMS('',I2.2,'') = '',F10.2)') 1, ZFPMS
  ENDIF

ENDIF

!*       3. RELATIVE OUTPUT COORDINATES AND GEOGRAPHICAL COMPASS,
!           SEEN FROM INPUT GEOMETRY
!           -----------------------------------------------------

! Activate interpolations
LLNEWGRID = ANY(YDFPUSERGEO(:)%LFPCOORD)

!       3.1 LIMITED AREA CASE

IF (LELAM.AND.LLNEWGRID) THEN

  IF (LMAP.AND.LRPLANE) THEN

    ALLOCATE(ZGNORX(KFPRGP))
    ALLOCATE(ZGNORY(KFPRGP))

    ZLON0=ELON0*Z180SPI ; ZLAT0=ELAT0*Z180SPI
    ZLON1=ELON1*Z180SPI ; ZLAT1=ELAT1*Z180SPI
    ZLONC=ELONC*Z180SPI ; ZLATC=ELATC*Z180SPI
    WRITE(NULOUT,*) 'Call ELALO2XY by SUFPG2'

    CALL ELALO2XY(KFPRGP,ZLON0,ZLAT0,ZLON1,ZLAT1,ZLONC,ZLATC,LMRT,&
     & ZGELAM,ZGELAT,PFPLO,PFPLA,ZGNORX,ZGNORY)

    DO J = 1, KFPRGP
      PFPNORX(J)=-ZGNORX(J)*ZFPNY2(J)+ZGNORY(J)*ZFPNX2(J)
      PFPNORY(J)=ZGNORY(J)*ZFPNY2(J)+ZGNORX(J)*ZFPNX2(J)
    ENDDO

    DEALLOCATE(ZGNORX)
    DEALLOCATE(ZGNORY)

  ELSE

    ! LMAP=.F. and LRPLANE=.F. are treated identicaly
    IPTR=0
    DO JDOM = 1, KFPDOM
      IF (NFPBOYD /= 0) THEN
        IFPLUX=NFPLUX(JDOM)+2*NFPBWX(JDOM)+2*(NLON(JDOM)-NFPLUX(JDOM))
        IFPGUX=NFPGUX(JDOM)+2*NFPBWY(JDOM)+2*(NLAT(JDOM)-NFPGUX(JDOM))
      ELSE
        IFPLUX=NFPLUX(JDOM)
        IFPGUX=NFPGUX(JDOM)
      ENDIF
      DO JGL = 1, IFPGUX
        DO JLON = 1, IFPLUX
          IPTR=IPTR+1
          PFPLO(IPTR) = REAL(JLON-1,JPRB)*RDELX(JDOM)
          PFPLA(IPTR) = REAL(JGL-1,JPRB)*RDELY(JDOM)
        ENDDO
      ENDDO
    ENDDO
    PFPNORX(:)=ZFPNX2(:)
    PFPNORY(:)=ZFPNY2(:)

  ENDIF

!       3.2 GLOBAL CASE

ELSEIF (.NOT.LELAM.AND.LLNEWGRID) THEN
!       3.2.1 Relative output coordinates

  ALLOCATE(ZSLAT(KFPRGP))
  ALLOCATE(ZSLAM(KFPRGP))
  ALLOCATE(ZCLAM(KFPRGP))

  ZSLACEN =  RMUCEN
  ZCLACEN = SQRT(1.0_JPRB-RMUCEN*RMUCEN)
  ZCLOCEN = COS(RLOCEN)
  ZSLOCEN = SIN(RLOCEN)

  !OIFS
  ZSLAT = ZGESLAT
  ZSLAM = ZGESLAM
  ZCLAM = ZGECLAM
  !OIFS





!       3.2.2 Direction of North, seen from input geometry

  ALLOCATE(ZFPNX1(KFPRGP))
  ALLOCATE(ZFPNY1(KFPRGP))

  CALL SUGENORD(NSTTYP,KFPRGP,RSTRET,RMUCEN,RLOCEN,&
   & ZSLAT,ZSLAM,ZCLAM,ZFPNX1,ZFPNY1)

!       3.2.3 Interface to output arrays

  DO J = 1, KFPRGP
    PFPLA(J) = ASIN(MAX(-1.0_JPRB,MIN(1.0_JPRB,ZSLAT(J))))
    PFPLO(J) = ACOS(MAX(-1.0_JPRB,MIN(1.0_JPRB,ZCLAM(J))))*SIGN(1.0_JPRB,ZSLAM(J))&
     & - Z2PI*MIN(0.0_JPRB,SIGN(1.0_JPRB,ZSLAM(J)))
    PFPNORX(J)=-ZFPNX1(J)*ZFPNY2(J)+ZFPNY1(J)*ZFPNX2(J)
    PFPNORY(J)=ZFPNY1(J)*ZFPNY2(J)+ZFPNX1(J)*ZFPNX2(J)
  ENDDO

  DEALLOCATE(ZFPNX1)
  DEALLOCATE(ZFPNY1)

!       3.3 NO ACTUAL INTERPOLATIONS

ELSE

  IF (LELAM) THEN
    IPTR=0
    IF (NFPBOYD /= 0) THEN
      IFPLUX=NFPLUX(1)+2*NFPBWX(1)+2*(NLON(1)-NFPLUX(1))
      IFPGUX=NFPGUX(1)+2*NFPBWY(1)+2*(NLAT(1)-NFPGUX(1))
    ELSE
      IFPLUX=NFPLUX(1)
      IFPGUX=NFPGUX(1)
    ENDIF
    DO JGL = 1, IFPGUX
      DO JLON = 1, IFPLUX
        IPTR=IPTR+1
        PFPLO(IPTR) = REAL(JLON-1,JPRB)*EDELX
        PFPLA(IPTR) = REAL(JGL-1,JPRB)*EDELY
      ENDDO
    ENDDO
  ELSE
    IPTR=0
    DO JGL = 1,NLAT(1)
      ZPISNL = Z2PI / YDFPUSERGEO(1)%NFPRGRI(JGL)
      DO JLON = 1, YDFPUSERGEO(1)%NFPRGRI(JGL)
        IPTR = IPTR+1
        PFPLA(IPTR)=ASIN(MAX(-1.0_JPRB,MIN(1.0_JPRB,REAL(YDFPUSERGEO(1)%FPMU(JGL),JPRB))))
        ZLON=REAL(JLON-1,JPRB)
        PFPLO(IPTR) = ZPISNL * ZLON
      ENDDO
    ENDDO
  ENDIF

  DO J = 1, KFPRGP
    PFPNORX(J) = 0.0_JPRB
    PFPNORY(J) = 1.0_JPRB
  ENDDO

ENDIF

!*       4. OUTPUT POINTS AT GEOGRAPHICAL POLES
!           -----------------------------------

IF (LLNEWGRID) THEN

!       Here it is considered for the time being that LALON grid cannot
!       be a rotated grid
  ITRIV=0
  IERR=0
  ZEPS=EPSILON(1.0_JPRB)*10._JPRB
  IF (CFPGRID(1) == 'GAUSS'.AND.NFPTTYP(1) == 2) THEN
    ZGEMUCEN=FPMUCEN(1)
    ZGELOCEN=FPLOCEN(1)
  ELSE
    ZGEMUCEN=1._JPRB
    ZGELOCEN=0._JPRB
  ENDIF
  IF (.NOT.LELAM) THEN
    IF (CFPGRID(1) == 'LALON') THEN
      IF (NSTTYP == 1) THEN
!             Trivial : poles are the sames
        ITRIV=1
      ELSE
!             Direct computation : poles are distincts
        ITRIV=2
        ZGEMUCEN=1._JPRB
        ZGELOCEN=0._JPRB
      ENDIF
    ELSEIF (CFPGRID(1) == 'GAUSS') THEN
      IF (NSTTYP /= NFPTTYP(1)) THEN
!             Direct computation : poles are distincts
        ITRIV=2
      ELSEIF (NSTTYP == 1.AND.NFPTTYP(1) == 1) THEN
!             Trivial : poles are the sames
        ITRIV=1
      ELSE
        IF (ABS(FPMUCEN(1)-RMUCEN) <= ZEPS  .AND.&
           & ABS(FPLOCEN(1)-RLOCEN) <= ZEPS       ) THEN
!               Trivial : poles are the sames
          ITRIV=1
        ELSE
!               Direct computation : poles are distincts
          ITRIV=2
        ENDIF
      ENDIF
    ELSE
      IERR=1
    ENDIF
  ELSEIF (LELAM.AND..NOT.LRPLANE) THEN
    IF (CFPGRID(1) == 'LALON') THEN
!  NB: The new EGGX does not use NROTEQ (it is virtually 0)
!      Therefore the rotated lat-lon domains should be defined maybe
!      using ELON0,ELAT0 - not yet coded
      ITRIV=1
    ELSEIF (CFPGRID(1) == 'GAUSS') THEN
      IF (NFPTTYP(1) == 1) THEN
!             Trivial : poles are the sames
        ITRIV=1
      ELSEIF (NFPTTYP(1) == 2) THEN
!             Direct computation from EGGX formulations : poles are distincts
        ITRIV=3
      ELSE
!             Direct computation from EGGX formulations possible only if poles
!             are distincts
        CALL ABOR1(' SUFPG2: ITRIV == 3 NOT YET CODED')
      ENDIF
    ELSE
      IERR=1
    ENDIF
  ELSE
    IF (CFPGRID(1) == 'GAUSS'.AND.NFPTTYP(1) == 2) THEN
      IERR=1
    ENDIF
  ENDIF

  DO JI=1, KFPRGP
    IF (ABS(COS(ZGELAT(JI))) < ZEPS) THEN

      IF (IERR == 1) THEN
!             No way out !
        WRITE(NULOUT,*) 'AN OUTPUT POINT AT A GEOGRAPHICAL POLE !'
        WRITE(NULOUT,*) 'THIS MAKES WIND TO BE UNCOMPUTABLE '
        WRITE(NULOUT,*) 'CHANGE THE GEOMETRY TO AVOID SUCH A'
        WRITE(NULOUT,*) 'SINGULAR POINT'
        CALL ABOR1(' SUFPG2: AN OUTPUT POINT AT A GEOGRAPHICAL POLE')
      ELSE
        IF (ITRIV == 1) THEN
!                 Trivial case
          PFPNORX(JI) = 0.0_JPRB
          PFPNORY(JI) = 1.0_JPRB
        ELSEIF (ITRIV == 2) THEN
!               Exact computation at poles
          IF (ABS(ABS(SIN(ZGELAT(JI)))-ZGEMUCEN) < ZEPS) THEN
!                 The point is at a pole of its grid
            ZSIGN=SIGN(1.0_JPRB,SIN(ZGELAT(JI)))
            IF (CFPGRID(1) == 'LALON') THEN
              PFPNORX(JI) = - SIN(ZGELAM(JI))
              PFPNORY(JI) = - COS(ZGELAM(JI))*ZSIGN
            ELSE
!                   Lobatto grid
              PFPNORX(JI) = - ZFPSLO(JI)
              PFPNORY(JI) = - ZFPCLO(JI)*ZSIGN
            ENDIF
          ELSE
!                 general case
            ZGEFPSLA(1)=ZGEMUCEN
            ZGEFPSLO(1)=SIN(ZGELOCEN)
            ZGEFPCLO(1)=COS(ZGELOCEN)
            ZSLACEN =  RMUCEN
            ZCLACEN = SQRT(1.0_JPRB-RMUCEN*RMUCEN)
            ZCLOCEN = COS(RLOCEN)
            ZSLOCEN = SIN(RLOCEN)

            !OIFS
            ZCAFPSLA = ZGEFPSLA
            ZCAFPSLO = ZGEFPSLO
            ZCAFPCLO = ZGEFPCLO
            !OIFS





            ZCAFPLO= ACOS(MAX(-1.0_JPRB,MIN(1.0_JPRB,ZCAFPCLO(1))))&
             & *SIGN(1.0_JPRB,ZCAFPSLO(1)) &
             & - Z2PI*MIN(0.0_JPRB,SIGN(1.0_JPRB,ZCAFPSLO(1)))
            ISTTYP=3
            CALL SUGENORD(ISTTYP,1,ZDUM,ZCAFPSLA(1),ZCAFPLO,&
             & ZSLAT(JI),ZSLAM(JI),ZCLAM(JI),&
             & PFPNORX(JI),PFPNORY(JI))
          ENDIF
        ELSEIF (ITRIV == 3) THEN
!               Exact computation at poles
          CALL ABOR1(' SUFPG2: ITRIV == 3 NOT YET CODED')
        ENDIF
      ENDIF

    ENDIF
  ENDDO

ENDIF

IF (ALLOCATED(ZFPMU)) DEALLOCATE(ZFPMU)
IF (ALLOCATED(ZFPSLO)) DEALLOCATE(ZFPSLO)
IF (ALLOCATED(ZFPCLO)) DEALLOCATE(ZFPCLO)
IF (ALLOCATED(ZSLAT))  DEALLOCATE(ZSLAT)
IF (ALLOCATED(ZSLAM))  DEALLOCATE(ZSLAM)
IF (ALLOCATED(ZCLAM))  DEALLOCATE(ZCLAM)
IF (ALLOCATED(ZFPNX1))  DEALLOCATE(ZFPNX1)
IF (ALLOCATED(ZFPNY1))  DEALLOCATE(ZFPNY1)

!*       5. COMPOSE ROTATION AND OUTPUT MAP FACTOR IF GAUSSIAN GRID IN OUTPUT
!           -----------------------------------------------------------------

IF (CFPGRID(1) == 'GAUSS') THEN
  DO J = 1, KFPRGP
    PFPNORX(J)=PFPNORX(J)/PFPGM(J)
    PFPNORY(J)=PFPNORY(J)/PFPGM(J)
  ENDDO
ENDIF

!*       6. COMPUTE LDPOSHOR
!           ----------------

IF (LLNEWGRID) THEN
!       Check if the output grid is a sub-set of the input grid
  LLSUBSET=.TRUE.
  ZEPS=1.E-4_JPRB
  IPTR=0
  IF (LELAM) THEN
    ALLOCATE(ZZLA(NDGUXG*NDLUXG),ZZLO(NDGUXG*NDLUXG))
    ZRELA=1.0_JPRB/EDELY
    ZRELO=1.0_JPRB/EDELX
    DO JGL = 1, NDGUXG
      DO JLON = 1, NDLUXG
        IPTR=IPTR+1
        ZZLO(IPTR) = REAL(JLON-1,JPRB)*EDELX
        ZZLA(IPTR) = REAL(JGL-1,JPRB)*EDELY
      ENDDO
    ENDDO
    ILAT=NLAT(1)
    IF (NFPGUX(1)>0) ILAT=NFPGUX(1)
    ILON=NLON(1)
    IF (NFPLUX(1)>0) ILON=NFPLUX(1)
    DO JGL = 1, ILAT
      DO JLON = 1, ILON
        JOUT=JLON+(JGL-1)*NLON(1)
        ZDLA=2.0_JPRB*ZEPS
        ZDLO=2.0_JPRB*ZEPS
        DO JIN=1,IPTR
          ZDLA=MIN(ZDLA,(ABS(PFPLA(JOUT)-ZZLA(JIN)))*ZRELA)
          ZDLO=MIN(ZDLO,(ABS(PFPLO(JOUT)-ZZLO(JIN)))*ZRELO)
        ENDDO
        IF (ZDLA > ZEPS.OR.ZDLO > ZEPS) THEN
          LLSUBSET=.FALSE.
          EXIT
        ENDIF
      ENDDO
    ENDDO
    DEALLOCATE(ZZLA)
    DEALLOCATE(ZZLO)
  ELSE
    LLSUBSET=.FALSE.
    WRITE(NULOUT,*) 'THE OUTPUT GRID IS ASSUMED NOT TO BE A SUB-SET OF THE INPUT ONE'
  ENDIF
ELSE
  LLSUBSET=.FALSE.
ENDIF

IF (LLSUBSET) THEN
  WRITE(NULOUT,*) 'THE OUTPUT GRID IS A SUB-SET OF THE INPUT ONE'
ENDIF

LDPOSHOR=LLNEWGRID.AND..NOT.LLSUBSET

IF(LTRACEFP) THEN
  WRITE(UNIT=NULOUT,FMT='('' LFPNEWGRID = '',L2)') LLNEWGRID
ENDIF
!-----------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUFPG2',1,ZHOOK_HANDLE)
END SUBROUTINE SUFPG2
