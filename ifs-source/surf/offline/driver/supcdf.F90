SUBROUTINE SUPCDF
USE PARKIND1  ,ONLY : JPIM     ,JPRB, JPRD
USE YOMHOOK   ,ONLY : LHOOK    ,DR_HOOK, JPHOOK

#ifdef DOC
! (C) Copyright 2000- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

!**** *SUPCDF * - Routine to initialize prognostic NetCDF output

!     Purpose.
!     --------
!           Initialize NetCDF outputfile o_gg.nc
!           Initialize NetCDF outputfile o_ocp.nc

!**   Interface.
!     ----------
!        *CALL* *SUPCDF

!        Explicit arguments :
!        --------------------
!        none

!     Method.
!     -------
!        Fields are stored in the model as NLON NLAT NLEVS TIME
!        and written to the NetCDF file in reverse order

!     Externals.
!     ----------

!     Reference.
!     ----------

!     Author.
!     -------
!        Bart vd HURK (KNMI)

!     Modifications.
!     --------------
!        Original     : 2000-7-14
!        Y. Takaya (ECMWF) add ocean mixed layer model part
!        Emanuel Dutra, May 2014: update netcd inteface and add netcdf4 
#endif

USE YOMLUN1S , ONLY : NPOSGG   ,NPOSOCP  ,NULOUT   ,NCTYPE
USE YOMDPHY  , ONLY : NCSS, NLAT     ,NLON, NCOM
USE YOMLOG1S , ONLY : NDIMCDF,LWRGG, LWROCP
USE YOEPHY   , ONLY : LEFLAKE,LESNML,LEURBAN
USE MPL_MODULE
USE NETCDF
USE NETCDF_UTILS, ONLY: NCERROR,TYPE_NCDF_VAR,INIT_NCDF_VAR,NC_DEF_VAR,&
                        DEF_NC_FILE,TYPE_NCDF_DIMSID
IMPLICIT NONE

!* -- Local variables 
!* -- NETCDF RELATED
INTEGER(KIND=JPIM) :: NPOS
INTEGER(KIND=JPIM) :: VARID,DIMIDS3(3),IDIMID3(3),IDIMID4(4)
INTEGER(KIND=JPIM) :: IDIM3,IDIM4

!* -- OTHERS
TYPE(TYPE_NCDF_VAR)          :: YD_VARINFO
TYPE(TYPE_NCDF_DIMSID)       :: YD_DIMSID

INTEGER(KIND=JPIM)           :: IVAR
INTEGER(KIND=JPIM),PARAMETER :: NVARSMAX=20 
INTEGER(KIND=JPIM)           :: NVARS3D,NVARS2D,NVARS3DSN
CHARACTER(LEN=20)            :: CVARS2D(NVARSMAX),CVARS3D(NVARSMAX),CVARS3DSN(NVARSMAX)
INTEGER(KIND=JPIM)           :: MYPROC, NPROC
REAL(KIND=JPHOOK)              :: ZHOOK_HANDLE


IF (LHOOK) CALL DR_HOOK('SUPCDF',0,ZHOOK_HANDLE)

MYPROC = MPL_MYRANK()
NPROC  = MPL_NPROC()

ASSOCIATE( NLATID=>YD_DIMSID%NLATID, NLONID=>YD_DIMSID%NLONID, NLEVSID=>YD_DIMSID%NLEVSID, &
         & NTILID=>YD_DIMSID%NTILID, NVTID=>YD_DIMSID%NVTID  , NMONID=>YD_DIMSID%NMONID  , &
         & NTIMID=>YD_DIMSID%NTIMID, NLEVSNID=>YD_DIMSID%NLEVSNID )

IF( MYPROC == 1 ) THEN
!* -- Initialize variables
NVARS2D=14
CVARS2D(1:NVARS2D)=(/'AvgSurfT    ','CanopInt    ','SWE         ',&
                     'SnowT       ','SAlbedo     ','snowdens    ',&
                     'TLICE       ','TLMNW       ','TLWML       ',&
                     'TLBOT       ','TLSF        ','HLICE       ',&
                     'HLML        ','slw         '/)
NVARS3D=3
CVARS3D(1:NVARS3D)=(/'SoilMoist   ','SoilTemp    ','icetemp     '/)

NVARS3DSN=4
CVARS3DSN(1:NVARS3DSN)=(/'SnowTML     ','SWEML       ','snowdensML  ',&
                     'slwML       '/)

!* -- 1. create netcdf file
  IF(LWRGG)THEN
    CALL NCERROR( NF90_CREATE('o_gg.nc',NCTYPE,NPOSGG),'Creating o_gg.nc file')
    WRITE(NULOUT,*)'NETCDF-FILE o_gg.nc OPENED ON UNIT ',NPOSGG
  ENDIF

  IF(LWRGG)THEN
    NPOS = NPOSGG

!* -- 2. add header to the file and global attributes
    CALL DEF_NC_FILE(NPOS,YD_DIMSID)

  !* -- general dimension definitions for 1d/2d grids 
    IF(NDIMCDF == 2)THEN
      IDIMID3(1) = NLONID
      IDIMID3(2) = NLATID
      IDIMID3(3) = NTIMID
      IDIM3=3
      IDIMID4(1) = NLONID
      IDIMID4(2) = NLATID
      IDIMID4(3) = NLEVSID
      IDIMID4(4) = NTIMID
      IDIM4=4
    ELSE
      IDIMID3(1) = NLONID
      IDIMID3(2) = NTIMID
      IDIM3=2
      IDIMID4(1) = NLONID
      IDIMID4(2) = NLEVSID
      IDIMID4(3) = NTIMID
      IDIM4=3
    ENDIF

  !* =====================================
  !*       2D variables 

    DO IVAR=1,NVARS2D
      CALL INIT_NCDF_VAR(YD_VARINFO,TRIM(CVARS2D(IVAR)),CCOORD_IN="time lat lon")
      CALL NC_DEF_VAR(NPOS,YD_VARINFO,IDIMID3(1:IDIM3),VARID)
    ENDDO

!   !* =====================================
!   !*       3D variables 
    DO IVAR=1,NVARS3D
      CALL INIT_NCDF_VAR(YD_VARINFO,TRIM(CVARS3D(IVAR)),CCOORD_IN="time nlevs lat lon")
      CALL NC_DEF_VAR(NPOS,YD_VARINFO,IDIMID4(1:IDIM4),VARID)
    ENDDO
  
!   !* =====================================
!   !*       3D variables snow
    IF (LESNML) THEN
  
    ! change DIMIDS4
      IF ( IDIM4 == 4 ) THEN
        IDIMID4(3) = NLEVSNID
      ELSE 
        IDIMID4(2) = NLEVSNID
      ENDIF
    
      DO IVAR=1,NVARS3DSN
        CALL INIT_NCDF_VAR(YD_VARINFO,TRIM(CVARS3DSN(IVAR)),CCOORD_IN="time nlevsn lat lon")
        CALL NC_DEF_VAR(NPOS,YD_VARINFO,IDIMID4(1:IDIM4),VARID)
      ENDDO 
    ENDIF

!* -- END OF DEFINITIONS 

    CALL NCERROR( NF90_ENDDEF(NPOS) )

  ENDIF

ENDIF

END ASSOCIATE

CALL MPL_BARRIER()

IF (LHOOK) CALL DR_HOOK('SUPCDF',1,ZHOOK_HANDLE)

END SUBROUTINE SUPCDF
