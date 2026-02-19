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

SUBROUTINE WRITE_SPEC_GRIB(YDGEOMETRY,YDRIP,CDFILE,YDSPEC,CDMODE,KSECLOC)

!     Purpose.
!     --------
!       Gather global versions of spectral arrays on one PE
!       and write it in grib file.

!     Author.
!     -------
!       Yannick Tremolet *ECMWF*

!     Modifications.
!     --------------
!       Original  : 29-01-01
!       21-05-02 Y.Tremolet: Use SPECTRAL_FIELD type.
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        M.Hamrud      10-Jan-2004 CY28R1 Cleaning
!        Y.Tremolet    23-Aug-2004 Fix problem for very low resolution
!        Y.Tremolet    04-Nov-2004 Added optional mode for appending to file
!        Y.Tremolet    10-Nov-2004 Option to force ECMWF local GRIB definition
!        Y.Tremolet    22-Nov-2004 Change SFC fields into ML 1 (for MARS).
!        Y.Tremolet    30-Mar-2005 PREGRBENC -> PRESET_GRIB_HEAD
!        M.Hamrud      27-Oct-2009 Use IOSTREAM
!        F.Bouttier    31-Jul-2012 LAM case still not working, but no longer fatal
! ------------------------------------------------------------------

USE YOMRIP             , ONLY : TRIP
USE GEOMETRY_MOD       , ONLY : GEOMETRY
USE PARKIND1           , ONLY : JPIM     ,JPRB
USE YOMHOOK            , ONLY : LHOOK,   DR_HOOK, JPHOOK
USE YOMCT0             , ONLY : LELAM
USE SPECTRAL_FIELDS_MOD, ONLY : ASSIGNMENT(=), SPECTRAL_FIELD
USE IOSTREAM_MIX       , ONLY : SETUP_IOSTREAM, SETUP_IOREQUEST, IO_PUT,&
 &                              CLOSE_IOSTREAM, TYPE_IOSTREAM , TYPE_IOREQUEST,&
 &                              CLOSE_IOREQUEST

IMPLICIT NONE

TYPE(GEOMETRY)      ,INTENT(IN)    :: YDGEOMETRY
TYPE(TRIP)          ,INTENT(IN)    :: YDRIP
CHARACTER(LEN=*)    ,INTENT(IN)    :: CDFILE
TYPE(SPECTRAL_FIELD),INTENT(IN)    :: YDSPEC
CHARACTER(LEN=1)  , OPTIONAL, INTENT(IN) :: CDMODE
INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KSECLOC(:)

INTEGER(KIND=JPIM) :: ILEVS3D(YDSPEC%NFLEVG),IBSET3D(YDSPEC%NFLEVG),JLEV,ILEVS2D(YDSPEC%NS2G)
INTEGER(KIND=JPIM) :: IOPROC
CHARACTER(LEN=1) :: CLMODE

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
TYPE(TYPE_IOSTREAM) :: YL_IOSTREAM
TYPE(TYPE_IOREQUEST) :: YL_IOREQUEST

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('WRITE_SPEC_GRIB',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP)
ASSOCIATE(NBSETLEV=>YDMP%NBSETLEV)
CLMODE="w"
IF (PRESENT(CDMODE)) CLMODE=CDMODE
IOPROC=1
DO JLEV=1,YDSPEC%NFLEVG
  ILEVS3D(JLEV)=JLEV
ENDDO
DO JLEV=1,YDSPEC%NFLEVG
  IBSET3D(JLEV) = NBSETLEV(ILEVS3D(JLEV))
ENDDO
ILEVS2D(:)=1

CALL SETUP_IOSTREAM(YL_IOSTREAM,'CIO',TRIM(CDFILE),CDMODE=CLMODE,&
 & KIOMASTER=IOPROC)
CALL SETUP_IOREQUEST(YL_IOREQUEST,'SPECTRAL_FIELDS',LDGRIB=.TRUE.,CDLEVTYPE='ML',&
 & KGRIB3D=YDSPEC%NGRIB3(1:YDSPEC%NS3D),KLEVS3D=ILEVS3D,KBSET3D=IBSET3D,&
 & KGRIB2D=YDSPEC%NGRIB2(1:YDSPEC%NS2G),KLEVS2D=ILEVS2D,&
 & KSMAX=YDSPEC%NSMAX,PTSTEP=YDRIP%TSTEP)
IF (LELAM) THEN    
  ! not implemented yet
ELSE
  CALL IO_PUT(YL_IOSTREAM,YL_IOREQUEST,PR2=YDSPEC%SP2D,PR3=YDSPEC%SP3D)
ENDIF
CALL CLOSE_IOREQUEST(YL_IOREQUEST)
CALL CLOSE_IOSTREAM(YL_IOSTREAM)
  

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('WRITE_SPEC_GRIB',1,ZHOOK_HANDLE)
END SUBROUTINE WRITE_SPEC_GRIB
