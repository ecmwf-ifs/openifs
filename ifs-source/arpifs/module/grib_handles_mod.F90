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

MODULE GRIB_HANDLES_MOD
! Purpose.
! --------
! Define type for keeping GRIB handles
! Keep the global GRIB samples in YGRIB_SAMPLES  to be ready to clone from

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 27 Sep 2017

!----------------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM,JPRB
!!USE TYPE_GEOMETRY,ONLY : GEOMETRY

TYPE TYPE_GRIB_HANDLES

!!TYPE(GEOMETRY),POINTER :: GEOM=>NULL()

INTEGER(KIND=JPIM) :: NGRIB_HANDLE_GG=-99
INTEGER(KIND=JPIM) :: NGRIB_HANDLE_GG2=-99
INTEGER(KIND=JPIM) :: NGRIB_HANDLE_GG_ML=-99
INTEGER(KIND=JPIM) :: NGRIB_HANDLE_SH=-99
INTEGER(KIND=JPIM) :: NGRIB_HANDLE_SH_ML=-99
INTEGER(KIND=JPIM) :: NGRIB_HANDLE_DIAG=-99

END TYPE TYPE_GRIB_HANDLES

TYPE(TYPE_GRIB_HANDLES) :: YGRIB_SAMPLES

CONTAINS

!============================================================================

SUBROUTINE READ_GRIB_SAMPLES

!**** *READ_GRIB_SAMPLES* - Routine to read GRIB API sample files

!     Purpose. Read sample files. These should never be modified,
!      just something to clone from
!     --------

!**   Interface.
!     ----------
!        *CALL* *READ_GRIB_SAMPLES*

!        Explicit arguments :  None.
!        --------------------

!        Implicit arguments :  
!        --------------------

!     Method.  
!     -------  

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 27 Sep 2017
!----------------------------------------------------------------------------

USE GRIB_API_INTERFACE, ONLY : IGRIB_NEW_FROM_SAMPLES, IGRIB_GET_MESSAGE_SIZE, &
 & IGRIB_GET_MESSAGE, IGRIB_NEW_FROM_MESSAGE, JPKSIZE_T, IGRIB_RELEASE, &
 & IGRIB_GET_VALUE, IGRIB_SET_VALUE
USE YOMHOOK           , ONLY : LHOOK, DR_HOOK, JPHOOK
USE MPL_MODULE        , ONLY : MPL_BROADCAST
USE YOMMP0            , ONLY : MYPROC
IMPLICIT NONE
INTEGER(KIND=JPIM), PARAMETER  :: JPGRIBHANDLES=6
INTEGER(KIND=JPIM) :: JGRIB, IOMASTER, ITAG, IGRIB2_TABLES_VERSION_LATEST
INTEGER(KIND=JPIM) :: IGRIB_HANDLE, IGRIB_HANDLEX, IGRIB_LENX, INTLEN, ILENINT
INTEGER(KIND=JPKSIZE_T)         :: IGRIB_LEN
INTEGER(KIND=JPIM), ALLOCATABLE :: IGRIBMSG(:)
CHARACTER(LEN=32) :: CLGRIB_SAMPL_FNAME(JPGRIBHANDLES)
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!----------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GRIB_HANDLES_MOD:READ_GRIB_SAMPLES',0,ZHOOK_HANDLE)

CLGRIB_SAMPL_FNAME(1) = 'gg_sfc_grib1' ! Edition 1
CLGRIB_SAMPL_FNAME(2) = 'gg_sfc_grib2' ! Edition 2
CLGRIB_SAMPL_FNAME(3) = 'gg_ml'        ! Edition 2 
CLGRIB_SAMPL_FNAME(4) = 'sh_sfc'       ! Edition 1
CLGRIB_SAMPL_FNAME(5) = 'sh_ml'        ! Edition 2
CLGRIB_SAMPL_FNAME(6) = 'diag'         ! Pseudo-GRIB
IOMASTER = 1
ITAG = 666
INTLEN = 4
DO JGRIB = 1, JPGRIBHANDLES
  IF(MYPROC == IOMASTER) THEN
    CALL IGRIB_NEW_FROM_SAMPLES(IGRIB_HANDLE,TRIM(CLGRIB_SAMPL_FNAME(JGRIB)))
    IF (JGRIB == 2.OR.JGRIB == 3.OR.JGRIB == 5) THEN
      ! Use latest tables version for the GRIB-2 samples
      CALL IGRIB_GET_VALUE(IGRIB_HANDLE,'tablesVersionLatest',IGRIB2_TABLES_VERSION_LATEST)
      CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'tablesVersion',IGRIB2_TABLES_VERSION_LATEST)
    ENDIF
    CALL IGRIB_GET_MESSAGE_SIZE(IGRIB_HANDLE,IGRIB_LEN)
    IGRIB_LENX = IGRIB_LEN
  ENDIF
  CALL MPL_BROADCAST(IGRIB_LENX, KTAG=ITAG,&
   & KROOT=IOMASTER,CDSTRING = 'GRIB_HANDLES_MOD:READ_GRIB_SAMPLES')
  ILENINT = (IGRIB_LENX+INTLEN-1)/INTLEN
  ALLOCATE(IGRIBMSG(ILENINT))
  IF(MYPROC == IOMASTER) THEN
    CALL IGRIB_GET_MESSAGE(IGRIB_HANDLE,IGRIBMSG)
    CALL IGRIB_RELEASE(IGRIB_HANDLE)
  ENDIF
  CALL MPL_BROADCAST(IGRIBMSG, KTAG=ITAG+1,&
   & KROOT=IOMASTER,CDSTRING = 'GRIB_HANDLES_MOD:READ_GRIB_SAMPLES 2')
  CALL IGRIB_NEW_FROM_MESSAGE(IGRIB_HANDLEX,IGRIBMSG)
  DEALLOCATE(IGRIBMSG)
  SELECT CASE (JGRIB)
  CASE(1)
    YGRIB_SAMPLES%NGRIB_HANDLE_GG    = IGRIB_HANDLEX
  CASE(2)
    YGRIB_SAMPLES%NGRIB_HANDLE_GG2   = IGRIB_HANDLEX
  CASE(3)
    YGRIB_SAMPLES%NGRIB_HANDLE_GG_ML = IGRIB_HANDLEX
  CASE(4)
    YGRIB_SAMPLES%NGRIB_HANDLE_SH    = IGRIB_HANDLEX
  CASE(5)
    YGRIB_SAMPLES%NGRIB_HANDLE_SH_ML = IGRIB_HANDLEX
  CASE(6)
    YGRIB_SAMPLES%NGRIB_HANDLE_DIAG  = IGRIB_HANDLEX
  END SELECT
ENDDO

IF (LHOOK) CALL DR_HOOK('GRIB_HANDLES_MOD:READ_GRIB_SAMPLES',1,ZHOOK_HANDLE)
!----------------------------------------------------------------------------

END SUBROUTINE READ_GRIB_SAMPLES
!============================================================================
SUBROUTINE RELEASE_GRIB_HANDLES(YDGBH)
USE GRIB_API_INTERFACE, ONLY : IGRIB_RELEASE
USE YOMHOOK           , ONLY : LHOOK, DR_HOOK, JPHOOK
IMPLICIT NONE
TYPE(TYPE_GRIB_HANDLES),INTENT(INOUT) :: YDGBH
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!----------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GRIB_HANDLES_MOD:RELEASE_GRIB_HANDLES',0,ZHOOK_HANDLE)

IF(YDGBH%NGRIB_HANDLE_GG    /= -99) CALL IGRIB_RELEASE(YDGBH%NGRIB_HANDLE_GG)
IF(YDGBH%NGRIB_HANDLE_GG2   /= -99) CALL IGRIB_RELEASE(YDGBH%NGRIB_HANDLE_GG2)
IF(YDGBH%NGRIB_HANDLE_GG_ML /= -99) CALL IGRIB_RELEASE(YDGBH%NGRIB_HANDLE_GG_ML)
IF(YDGBH%NGRIB_HANDLE_SH    /= -99) CALL IGRIB_RELEASE(YDGBH%NGRIB_HANDLE_SH)
IF(YDGBH%NGRIB_HANDLE_SH_ML /= -99) CALL IGRIB_RELEASE(YDGBH%NGRIB_HANDLE_SH_ML)
IF(YDGBH%NGRIB_HANDLE_DIAG  /= -99) CALL IGRIB_RELEASE(YDGBH%NGRIB_HANDLE_DIAG)

IF (LHOOK) CALL DR_HOOK('GRIB_HANDLES_MOD:RELEASE_GRIB_HANDLES',1,ZHOOK_HANDLE)
END SUBROUTINE RELEASE_GRIB_HANDLES
!============================================================================
END MODULE GRIB_HANDLES_MOD
