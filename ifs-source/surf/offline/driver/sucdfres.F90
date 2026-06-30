SUBROUTINE SUCDFRES

#ifdef DOC
! (C) Copyright 2000- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

!**** *SUCDFFL * - Routine to initialize NetCDF output

!     Purpose.
!     --------
!           Initialize NetCDF outputfile restartout.nc (restart output file)

!**   Interface.
!     ----------
!        *CALL* *SUCDFRES

!     Explicit arguments :
!     --------------------
!        none

!     Method.
!     -------

!     Externals.
!     ----------

!     Reference.
!     ----------

!     Author.
!     -------
!        Bart vd HURK (KNMI)

!     Modifications.
!     --------------
!        Original     : 2000-7-17
!       E. Dutra     2008, add FLAKE restart fields
!       E. Dutra     11/2009 add LAI restart
!       Emanuel Dutra, June 2014: update netcd inteface and add netcdf4 
!       Anna Agusti-Panareda, June 2021: Add C3/C4 photosynthetic pathway type (Ctype)

#endif
USE PARKIND1  ,ONLY : JPIM     ,JPRB,   JPRD
USE YOMHOOK   ,ONLY : LHOOK    ,DR_HOOK, JPHOOK

USE YOMLUN1S , ONLY : NPOSRES  ,NULOUT   ,CTUNITS,NCTYPE 
USE YOMDPHY  , ONLY : NCSS, NLAT     ,NLON      ,NTILES ,NVHILO,NCLIMDAT

USE YOMLOG1S , ONLY : CMODID   ,CVERID   ,NDIMCDF
USE YOMRIP   , ONLY : NINDAT   ,NSSSSS
USE YOEPHY   , ONLY : LEFLAKE  ,LECTESSEL,LESNML, LEC4MAP
USE NETCDF
USE NETCDF_UTILS, ONLY: NCERROR,TYPE_NCDF_VAR,INIT_NCDF_VAR,NC_DEF_VAR,&
                        DEF_NC_FILE,TYPE_NCDF_DIMSID
USE MPL_MODULE
IMPLICIT NONE

!* -- Local variables 
!* -- NETCDF RELATED
INTEGER(KIND=JPIM) :: NPOS
INTEGER(KIND=JPIM) :: VARID,IDIMID2(2),IDIM2,IDIMID3(3),IDIM3

!* -- OTHERS
TYPE(TYPE_NCDF_VAR)          :: YD_VARINFO
TYPE(TYPE_NCDF_DIMSID)       :: YD_DIMSID
INTEGER(KIND=JPIM)           :: NVARS2D,NVARS3D,IVAR
CHARACTER(LEN=20)            :: CVARS2D(40),CVARS3D(40)
INTEGER(KIND=JPIM)           :: MYPROC, NPROC
REAL(KIND=JPHOOK)              :: ZHOOK_HANDLE

ASSOCIATE( NLATID=>YD_DIMSID%NLATID, NLONID=>YD_DIMSID%NLONID, NLEVSID=>YD_DIMSID%NLEVSID, &
         & NTILID=>YD_DIMSID%NTILID, NVTID=>YD_DIMSID%NVTID  , NMONID=>YD_DIMSID%NMONID  , &
         & NTIMID=>YD_DIMSID%NTIMID, NLEVSNID=>YD_DIMSID%NLEVSNID )

IF (LHOOK) CALL DR_HOOK('SUCDFRES',0,ZHOOK_HANDLE)

MYPROC = MPL_MYRANK()
NPROC  = MPL_NPROC()

IF( MYPROC == 1 ) THEN

!* -- 1. create netcdf file

  CALL NCERROR( NF90_CREATE('restartout.nc',NCTYPE,NPOSRES),'Creating restartout.nc file')
  WRITE(NULOUT,*)'NETCDF-FILE restartout.nc OPENED ON UNIT ',NPOSRES
  NPOS=NPOSRES

!* -- 2. add header to the file and global attributes
  CALL DEF_NC_FILE(NPOS,YD_DIMSID)
 
!* =====================================
!*       2D Fields lat/lon (or lon only)
IDIMID2(1) = NLONID
IDIMID2(2) = NLATID
IF (NDIMCDF == 2 ) THEN
  IDIM2 = 2
ELSE
  IDIM2 = 1
ENDIF

  NVARS2D=31
  CVARS2D(1:NVARS2D)=(/'Mask     ','z0m      ','lz0h     ','landsea  ','geopot   ',&
                      'cvl      ','cvh      ','tvl      ','tvh      ','sotype   ',&
                      'sdor     ','sst      ','seaice   ','CanopInt ','SWE      ',&
                      'SnowT    ','SAlbedo  ','snowdens ','AvgSurfT ','TLICE    ',&
                      'TLMNW    ','TLWML    ','TLBOT    ','TLSF     ','HLICE    ',&
                      'HLML     ','LDEPTH   ','CLAKE    ','x        ','CLAKEF   ',&
                      'cu       '/)

  IF (LEC4MAP) THEN

    NVARS2D=NVARS2D+1
    CVARS2D(NVARS2D)='Ctype    '

  ENDIF

  DO IVAR=1,NVARS2D
    CALL INIT_NCDF_VAR(YD_VARINFO,TRIM(CVARS2D(IVAR)),IACCUR_IN=NF90_DOUBLE,&
                       CCOORD_IN="lat lon")
    CALL NC_DEF_VAR(NPOS,YD_VARINFO,IDIMID2(1:IDIM2),VARID)
  ENDDO

!* =====================================
!*       3D Fields - nlevs : lat/lon/level

  IF (NDIMCDF == 2 ) THEN
    IDIMID3(1) = NLONID
    IDIMID3(2) = NLATID
    IDIMID3(3) = NLEVSID
    IDIM3=3
  ELSE
    IDIMID3(1) = NLONID
    IDIMID3(2) = NLEVSID
    IDIM3=2
  ENDIF
  NVARS3D=3
  CVARS3D(1:NVARS3D)=(/'SoilTemp ','SoilMoist','icetemp  '/)
  DO IVAR=1,NVARS3D
    CALL INIT_NCDF_VAR(YD_VARINFO,TRIM(CVARS3D(IVAR)),IACCUR_IN=NF90_DOUBLE,&
                       CCOORD_IN="nlevs lat lon")
    CALL NC_DEF_VAR(NPOS,YD_VARINFO,IDIMID3(1:IDIM3),VARID)
  ENDDO

!* =====================================
!*       3D Fields - month : lat/lon/month

  IF (NDIMCDF == 2 ) THEN
    IDIMID3(1) = NLONID
    IDIMID3(2) = NLATID
    IDIMID3(3) = NMONID
    IDIM3=3
  ELSE
    IDIMID3(1) = NLONID
    IDIMID3(2) = NMONID
    IDIM3=2
  ENDIF
  NVARS3D=4
  CVARS3D(1:NVARS3D)=(/'Malbedo','Mlail  ','Mlaih  ','fwet   '/)
  DO IVAR=1,NVARS3D
    CALL INIT_NCDF_VAR(YD_VARINFO,TRIM(CVARS3D(IVAR)),IACCUR_IN=NF90_DOUBLE,&
                       CCOORD_IN="month lat lon")
    CALL NC_DEF_VAR(NPOS,YD_VARINFO,IDIMID3(1:IDIM3),VARID)
  ENDDO

!* =====================================
!*       3D Fields: veg. type : lat/lon/vegetation type fraction 

  IF (NDIMCDF == 2 ) THEN
    IDIMID3(1) = NLONID
    IDIMID3(2) = NLATID
    IDIMID3(3) = NVTID
    IDIM3=3
  ELSE
    IDIMID3(1) = NLONID
    IDIMID3(2) = NVTID
    IDIM3=2
  ENDIF
  NVARS3D=9
  CVARS3D(1:NVARS3D)=(/'r0vt    ','Anday   ','Anfm    ','Rstr    ',&
                        'Rstr2   ','Blast   ','laivt   ','Biomstr ',&
                        'Biomstr2'/)
  DO IVAR=1,NVARS3D
    CALL INIT_NCDF_VAR(YD_VARINFO,TRIM(CVARS3D(IVAR)),IACCUR_IN=NF90_DOUBLE,&
                       CCOORD_IN="vtype lat lon")
    CALL NC_DEF_VAR(NPOS,YD_VARINFO,IDIMID3(1:IDIM3),VARID)
  ENDDO

!* =====================================
!*       3D Fields: tiled : lat/lon/tile

  IF (NDIMCDF == 2 ) THEN
    IDIMID3(1) = NLONID
    IDIMID3(2) = NLATID
    IDIMID3(3) = NTILID
    IDIM3=3
  ELSE
    IDIMID3(1) = NLONID
    IDIMID3(2) = NTILID
    IDIM3=2
  ENDIF
  NVARS3D=5
  CVARS3D(1:NVARS3D)=(/'Qh    ','Evap  ','Momu  ','Momv  ','SkinT '/)
  DO IVAR=1,NVARS3D
    CALL INIT_NCDF_VAR(YD_VARINFO,TRIM(CVARS3D(IVAR)),IACCUR_IN=NF90_DOUBLE,&
                       CCOORD_IN="tile lat lon")
    CALL NC_DEF_VAR(NPOS,YD_VARINFO,IDIMID3(1:IDIM3),VARID)
  ENDDO

!* =====================================
!*       3D Fields: snow : lat/lon/nlevsn
  IF (LESNML) THEN
    IF (NDIMCDF == 2 ) THEN
      IDIMID3(1) = NLONID
      IDIMID3(2) = NLATID
      IDIMID3(3) = NLEVSNID
      IDIM3=3
    ELSE
      IDIMID3(1) = NLONID
      IDIMID3(2) = NLEVSNID
      IDIM3=2
    ENDIF
    NVARS3D=4
    CVARS3D(1:NVARS3D)=(/'SnowTML     ','SWEML       ','snowdensML  ',&
                       'slwML       '/)
    DO IVAR=1,NVARS3D
      CALL INIT_NCDF_VAR(YD_VARINFO,TRIM(CVARS3D(IVAR)),IACCUR_IN=NF90_DOUBLE,&
                        CCOORD_IN="nlevsn lat lon")
      CALL NC_DEF_VAR(NPOS,YD_VARINFO,IDIMID3(1:IDIM3),VARID)
    ENDDO
  ENDIF


!* -- END OF DEFINITIONS 

  CALL NCERROR( NF90_ENDDEF(NPOS) )

ENDIF

END ASSOCIATE

CALL MPL_BARRIER()

IF (LHOOK) CALL DR_HOOK('SUCDFRES',1,ZHOOK_HANDLE)

END SUBROUTINE SUCDFRES
