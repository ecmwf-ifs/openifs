! (C) Copyright 1995- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction
SUBROUTINE SUFCDF

USE YOMFORC1S, ONLY : JPSTPFC  ,UFI      ,VFI      ,TFI      ,&
     &            QFI      ,PSFI     ,SRFFI    ,TRFFI    ,R30FI    ,&
     &            S30FI    ,R30FI_C  ,S30FI_C  ,DTIMFC   ,RTSTFC   ,&
     &            NSTPFC,DIMFORC, CO2FI
USE YOMRIP   , ONLY : RTIMTR
USE YOMCT01S , ONLY : NSTOP    ,NSTART
USE YOMCST   , ONLY : RTT ,RDAY ,RG ,RETV, RLVTT, RLSTT, RMD, RMCO2
USE YOMLUN1S , ONLY : NULOUT   ,NULNAM
USE YOMDYN1S , ONLY : TSTEP
USE YOMLOG1S , ONLY : LDBGS1   ,NDIMCDF
USE YOMGF1S  , ONLY : RALT     ,RZUV
USE YOMGC1S  , ONLY : LMASK
USE YOMDPHY  , ONLY : NPOI,    NLON      ,NLAT     ,NLALO
USE YOMGPD1S , ONLY : VFZ0F,   VFGEO
USE YOETHF   , ONLY : R2ES     ,R3LES    ,R3IES    ,R4LES    ,&
     &            R4IES    ,R5LES    ,R5IES    ,RVTMP2   ,R5ALVCP  ,&
     &            R5ALSCP  ,RALVDCP  ,RALSDCP  ,RTWAT    ,RTICE    ,&
     &            RTWAT_RTICE_R, RTICECU, RTWAT_RTICECU_R
USE PARKIND1  ,ONLY : JPIM     ,JPRB, JPRD
USE YOEPHY    ,ONLY : LEAIRCO2COUP
USE YOEPHY    ,ONLY : LEINTWIND
USE YOMHOOK   ,ONLY : DR_HOOK, JPHOOK, LHOOK
USE MPL_MODULE

#ifdef DOC

!**** *SUFCDF  * - ROUTINE TO SETUP THE ATMOSPHERIC FORCING DATA
!                  FROM NETCDF INPUT

!     PURPOSE.
!     --------
!        Initialize the common block YOMFORC

!**   INTERFACE.
!     ----------
!        *CALL* *SUFCDF*

!     EXPLICIT ARGUMENTS :  
!     --------------------

!     IMPLICIT ARGUMENTS :
!     --------------------
!        COMMON  YOMFORC
!        COMMON  YOMRIP
!        COMMON  YOMLUN1S
!        COMMON  YOMDYN1S
!        COMMON  YOMGC1S
!        COMMON  YOMLOG1S

!     METHOD.
!     -------
!     EXTERNALS.
!     ----------
!        IYMD2C

!     REFERENCE.
!     ----------
!        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS

!     AUTHOR.
!     -------
!        JEAN-FRANCOIS MAHFOUF AND PEDRO VITERBO  *ECMWF*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 95-03-13
!        BART VD HURK (KNMI): READING NETCDF INPUT
!        ANNA AGUSTI-PANAREDA: 2020-11-17 VARIABLE ATMOSPHERIC CO2 FORCING

#endif
IMPLICIT NONE
#include "rdfvar.intfb.h"
#include "minmax.intfb.h"

CHARACTER CHEADER*400
CHARACTER*100 CNAME

!* Netcdf interface
!
REAL*4,ALLOCATABLE :: ZREAL3(:,:)
REAL*4,ALLOCATABLE :: ZLATI(:)
REAL*4,ALLOCATABLE :: ZLONI(:)
INTEGER ISTART3(3),ICOUNT3(3),IDUMAR(10)
INTEGER ISTART2(2),ICOUNT2(2)
INTEGER NCID,IERR,NIDLON,NILON,NIDLAT,NILAT,NIDTIM,NITIM,NVARID,IDUM,&
     &       NVARS,NDIM,NGATTS,IRECDIM

!* Other variables
REAL(KIND=JPRB),ALLOCATABLE :: ZSLMI(:)
REAL(KIND=JPRB),ALLOCATABLE :: ZGEOPD(:)
REAL(KIND=JPRB),ALLOCATABLE :: ZREAL3D(:,:)
INTEGER(KIND=JPIM) :: IFHR,IFMIN,IFSS,IDATST,ISTHR,ISTMIN,ISTSS,&
     &       ITIMST,IRDST,JT,NMX,NMY,NHORI,IVAR,JL
REAL(KIND=JPRB) :: ZFORREF,ZFORST,ZFORLS,ZTSTRT,ZFCEN,&
     &       ZFAC,ZDZ,ZQS,ZCOR,ZRH,ZFRAC
LOGICAL LWIND2D,LCTPF,LSNOWF,LOINTP ! spatial interpolation (only when forcing = 2D)

INTEGER(KIND=JPIM) :: MYPROC, NPROC
REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

#include "fcttim.h"
#include "namforc1s.h"
#include "netcdf.inc"
#include "fcttre.h"


IF (LHOOK) CALL DR_HOOK('SUFCDF',0,ZHOOK_HANDLE)

MYPROC = MPL_MYRANK()
NPROC  = MPL_NPROC()
!---------------------------------------------------------------------

! CONTENTS OF NETCDF FILE
!  lat:units = "degrees_north";long_name = "latitude";
!  lon:units = "degrees_east";long_name = "longitude";
!  time:units = "seconds";long_name = "Seconds since 19790101:00.00";
!  Tair:units = "K";long_name = "Temperature";
!  CO2air:units = "kg/kg";long_name = "Atmospheric CO2 mixing ratio";
!  Qair:units = "kg/kg";long_name = "Specific humidity";
!  Wind_E:units = "m/s";long_name = "Wind speed u";
!  Wind_N:units = "m/s";long_name = "Wind speed v";
!  SWdown:units = "W/m2";long_name = "Downward shortwave radiation";
!  LWdown:units = "W/m2";long_name = "Downward longwave radiation";
!  Rainf:units = "kg/m2s";long_name = "Rainfall";
!  Snowf:units = "kg/m2s";long_name = "Snowfall";
!  PSurf:units = "Pa";long_name = "Pressure";
!  Ctpf:units = "", long_name =" convective precipitation fraction"

!---------------------------------------------------------------------

!     1. READ DATA PROPERTIES

LCTPF=.FALSE.
LWIND2D=.FALSE.
LSNOWF=.FALSE.

ZPHISTA=10.
ZUV=10.0_JPRB
ZDTFORC=3600.

IFYYYY=1979
IFMM=1
IFDD=1
IFTIM=0

NDIMFORC=NDIMCDF
LOADIAB=.FALSE.
LWIND2D=.FALSE.
LSNOWF=.FALSE.
LCTPF=.FALSE.
CFORCU   ='forcing'
CFORCV   ='forcing'
CFORCT   ='forcing'
CFORCQ   ='forcing'
CFORCC   ='forcing'
CFORCP   ='forcing'
CFORCG   ='forcing'
CFORCRAIN='forcing'
CFORCSNOW='forcing'
CFORCSW  ='forcing'
CFORCLW  ='forcing'

REWIND(NULNAM)
READ(NULNAM,NAMFORC)

!! Move to modules ! 
DTIMFC =ZDTFORC
RALT = ZPHISTA
RZUV = ZUV
DIMFORC=NDIMFORC


LOINTP=.FALSE.
WRITE(NULOUT,*)'RALT = ',RALT
WRITE(NULOUT,*)'RZUV = ',RZUV
WRITE(NULOUT,*)'DTIMFC  = ',DTIMFC
WRITE(NULOUT,*)'DIMFORC  = ',DIMFORC


!* read forcing variables one by one
!* check presence of specific fields:
!* 2D wind 
LWIND2D=.FALSE. 
NCID = NCOPN(CFORCU, NCNOWRIT, IERR)
IF( IERR == 0 ) THEN
  CALL NCINQ(NCID, IDUM, NVARS, IDUM, IDUM, IERR)
  DO IVAR=1,NVARS
    CALL NCVINQ(NCID,IVAR,CNAME,IDUM,IDUM,&
       &     IDUMAR,IDUM,IERR)
    IF(CNAME(1:6).EQ.'Wind_E')THEN
      LWIND2D=.TRUE.
    ENDIF
  ENDDO
ENDIF
CALL NCCLOS(NCID,IERR)

IF(LWIND2D)THEN

!* wind speed u
  CNAME='Wind_E'
  CALL RDFVAR(CFORCU,CNAME,UFI)
  CNAME='Wind_N'
  CALL RDFVAR(CFORCV,CNAME,VFI)
ELSE

!* total wind speed
  CNAME='Wind'
  CALL RDFVAR(CFORCU,CNAME,UFI)
  VFI=0.0_JPRB
ENDIF

!* interpolate wind to correct height

IF(ZUV.NE.RALT)THEN
!   DO JL=1,NPOI
!     ZFAC=LOG(RALT/VFZ0F(JL))/LOG(ZUV/VFZ0F(JL))
!     UFI(JL,1:NSTPFC)=ZFAC*UFI(JL,1:NSTPFC)
!     VFI(JL,1:NSTPFC)=ZFAC*VFI(JL,1:NSTPFC)
!   ENDDO
  IF (LEINTWIND) THEN
    WRITE(NULOUT,*)'SUFCDF: WIND DATA INTERPOLATED FROM ',ZUV,&
      &     ' TO ',RALT,' M'
    WRITE(NULOUT,*)'We are not sure if this is correct! '
    WRITE(NULOUT,*)'Check sufcdf and callpar1s'
  ELSE
    WRITE(NULOUT,*)'SUFCDF: WIND DATA INTERPOLATED FROM ',ZUV,&
      &     ' TO ',RALT,' M'
    WRITE(NULOUT,*)'Set LEINTWIND=.TRUE. in NAMPHY1S' 
    CALL ABOR1('Surfcdf:')
  ENDIF
ENDIF

!* temperature
CNAME='Tair'
CALL RDFVAR(CFORCT,CNAME,TFI)

!* specific humidity
CNAME='Qair'
CALL RDFVAR(CFORCQ,CNAME,QFI)

IF (LEAIRCO2COUP) THEN
  WRITE(NULOUT,*) 'sufcdf: LEAIRCO2COUP = ',LEAIRCO2COUP
  !* atmospheric CO2
  WRITE(NULOUT,*) 'sufcdf: reading CO2 forcing '
  CNAME='CO2air'
  CALL RDFVAR(CFORCC,CNAME,CO2FI)
ENDIF

!* surface pressure
CNAME='PSurf'
CALL RDFVAR(CFORCP,CNAME,PSFI)

!* ADIABATIC HEIGHT CORRECTION

IF(LOADIAB)THEN
! read surface geopotential        
  NCID = NCOPN(CFORCG, NCNOWRIT, IERR)
  ISTART2=1
  ICOUNT2(1)=NILON
  ICOUNT2(2)=NILAT
  NHORI=NILON*NILAT
  ALLOCATE (ZREAL3(NHORI,1))
  ALLOCATE (ZREAL3D(NLALO,1))
  ALLOCATE (ZSLMI(NHORI))
  ZSLMI=1.
  ALLOCATE (ZLATI(NILAT))
  ALLOCATE (ZLONI(NILON))
  NVARID = NCVID(NCID, 'lat', IERR)
  CALL NCVGT(NCID, NVARID, 1, NILAT, ZLATI, IERR)
  NVARID = NCVID(NCID, 'lon', IERR)
  CALL NCVGT(NCID, NVARID, 1, NILON, ZLONI, IERR)

  NVARID = NCVID(NCID, 'geopot', IERR)
  CALL NCVGT(NCID, NVARID, ISTART2,ICOUNT2, ZREAL3(1,1), IERR)
  CALL NCCLOS(NCID,IERR)
  IF(LOINTP)THEN
    !CALL INTPF(NILON,NILAT,NLALO,ZLONI,ZLATI,&
    ! &               ZREAL3(1,1),ZSLMI,ZREAL3D(1,1))
    WRITE(NULOUT,*)'LOINTP option not possivle in sufcdf'
    CALL ABOR1('SUFCDF:')
  ELSE
    ZREAL3D(:,1)=ZREAL3(1,1)
  ENDIF
  IF( MYPROC == 1 ) THEN
    CALL MINMAX('GEOPOT',ZREAL3D(:,1),NMX,NMY,LMASK,NULOUT)
  ENDIF
  ALLOCATE (ZGEOPD(NPOI))
  ZGEOPD=PACK(ZREAL3D(:,1),LMASK)
  DO JL=1,NPOI
    ZDZ=(VFGEO(JL)-ZGEOPD(JL))/RG
    DO JT=1,NSTPFC
      ZQS=FOEEW(TFI(JL,JT))/PSFI(JL,JT)
      ZCOR=1.0_JPRB/(1.0_JPRB-RETV*ZQS)
      ZQS=ZQS*ZCOR
      ZRH=MIN(1._JPRB,QFI(JL,JT)/ZQS)
      TFI(JL,JT)=TFI(JL,JT)-0.0065_JPRB*ZDZ
      PSFI(JL,JT)=PSFI(JL,JT)-10.0_JPRB*ZDZ
      ZQS=FOEEW(TFI(JL,JT))/PSFI(JL,JT)
      ZCOR=1.0_JPRB/(1.0_JPRB-RETV*ZQS)
      ZQS=ZQS*ZCOR
      QFI(JL,JT)=ZRH*ZQS
    ENDDO
  ENDDO
  DEALLOCATE(ZGEOPD,ZREAL3,ZREAL3D,ZSLMI,ZLONI,ZLATI)
ENDIF
   
!* shortwave radiation
CNAME='SWdown'
CALL RDFVAR(CFORCSW,CNAME,SRFFI)

!* longwave radiation
CNAME='LWdown'
CALL RDFVAR(CFORCLW,CNAME,TRFFI)

!* rainfall (all large scale)
CNAME='Rainf'
CALL RDFVAR(CFORCRAIN,CNAME,R30FI)

!* snowfall (all large scale)
!* check presence of snowfall precipitation 
LSNOWF=.FALSE.
NCID = NCOPN(CFORCSNOW, NCNOWRIT, IERR)
IF( IERR == 0 ) THEN
  CALL NCINQ(NCID, IDUM, NVARS, IDUM, IDUM, IERR)
  DO IVAR=1,NVARS
    CALL NCVINQ(NCID,IVAR,CNAME,IDUM,IDUM,&
      &     IDUMAR,IDUM,IERR)
    IF(CNAME(1:5).EQ.'Snowf')THEN
      LSNOWF=.TRUE.
    ENDIF
  ENDDO
ENDIF
CALL NCCLOS(NCID,IERR)

IF (LSNOWF) THEN
  CNAME='Snowf'
  CALL RDFVAR(CFORCSNOW,CNAME,S30FI)
ELSE
  WRITE(NULOUT,*) " SNOWF DEDUCED FROM TAIR AND RAINF "
  S30FI=0.0_JPRB
  DO JL=1,NPOI
  DO JT=1,NSTPFC
    IF(TFI(JL,JT).LT.RTT)THEN
      S30FI(JL,JT)=R30FI(JL,JT)
      R30FI(JL,JT)=0.0_JPRB
    ENDIF
  ENDDO
  ENDDO
ENDIF


!* check presence of convective precipitation fraction # should be in the rainf file 
LCTPF=.FALSE.
NCID = NCOPN(CFORCRAIN, NCNOWRIT, IERR)
IF( IERR == 0 ) THEN
  CALL NCINQ(NCID, IDUM, NVARS, IDUM, IDUM, IERR)
  DO IVAR=1,NVARS
    CALL NCVINQ(NCID,IVAR,CNAME,IDUM,IDUM,&
       &     IDUMAR,IDUM,IERR)
    IF(CNAME(1:4).EQ.'Ctpf')THEN
      LCTPF=.TRUE.
    ENDIF
  ENDDO
ENDIF
CALL NCCLOS(NCID,IERR)

IF (LCTPF) THEN
  CNAME='Ctpf'
  ALLOCATE(ZREAL3D(NPOI,JPSTPFC))
  CALL RDFVAR(CFORCRAIN,CNAME,ZREAL3D)
  IF (MINVAL(ZREAL3D).LT.0.0_JPRB) THEN
    WRITE(NULOUT,*) " CONVECTIVE FRACTION < 0 , ABORT",MINVAL(ZREAL3D)
!     CALL ABOR1('SUFCDF:')
  ENDIF
  IF (MAXVAL(ZREAL3D).GT.1.0_JPRB) THEN
    WRITE(NULOUT,*) " CONVECTIVE FRACTION > 1 , ABORT",MAXVAL(ZREAL3D)
!     CALL ABOR1('SUFCDF:')
  ENDIF
  WRITE(NULOUT,*) "Correcting Rainf/Snoww using Ctpf "
  DO JL=1,NPOI
    DO JT = 1, NSTPFC
      ZFRAC=MIN(1.0_JPRB,MAX(0._JPRB,ZREAL3D(JL,JT)))
      R30FI_C(JL,JT)=ZFRAC*R30FI(JL,JT)
      S30FI_C(JL,JT)=ZFRAC*S30FI(JL,JT)
      R30FI(JL,JT)=(1.0_JPRB - ZFRAC)*R30FI(JL,JT)
      S30FI(JL,JT)=(1.0_JPRB - ZFRAC)*S30FI(JL,JT)
    ENDDO
  ENDDO
ELSE
  WRITE(NULOUT,*) " CONVECTIVE TP NOT PRESENT, SET TO 0. "
  R30FI_C(:,1:NSTPFC)=0.0_JPRB
  S30FI_C(:,1:NSTPFC)=0.0_JPRB
ENDIF

!Convert atmospheric CO2 from ppm to kg/kg
IF (LEAIRCO2COUP) THEN
  DO JL=1,NPOI
  DO JT=1,NSTPFC
      CO2FI(JL,JT)=CO2FI(JL,JT)*RMCO2/(RMD*1000000._JPRB)
  ENDDO
  ENDDO
ENDIF



WRITE(NULOUT,*) " FORCING DATA READ FOR ",NSTPFC," FORCING STEPS"

IF (LHOOK) CALL DR_HOOK('SUFCDF',1,ZHOOK_HANDLE)

RETURN
END SUBROUTINE SUFCDF
