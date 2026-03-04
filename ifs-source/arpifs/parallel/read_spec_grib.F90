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

SUBROUTINE READ_SPEC_GRIB(YDMP,CDFILE,YDSPEC)

!     Purpose.
!     --------
!       Read spectral fields from grib file and distribute them.

!     Author.
!     -------
!       Yannick Tremolet *ECMWF*

!     Modifications.
!     --------------
!       Original  : 13-02-03
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        Y.Tremolet    08-Jul-2004 Truncate/pad fields if resolution differ
!        Y.Tremolet    11-Aug-2005 Do not assume order of fields in file and
!                                  distribute work in the same way as SUSPECG
!        M.Hamrud      27-Oct-2009 Use IOSTREAM
! ------------------------------------------------------------------

USE YOMMP              , ONLY : TMP
USE PARKIND1           , ONLY : JPIM, JPRB
USE YOMHOOK            , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN             , ONLY : NULOUT
USE SPECTRAL_FIELDS_MOD, ONLY : ASSIGNMENT(=), SPECTRAL_FIELD
USE IOSTREAM_MIX       , ONLY : SETUP_IOSTREAM, SETUP_IOREQUEST, IO_GET,&
 &                              CLOSE_IOSTREAM, TYPE_IOSTREAM , TYPE_IOREQUEST,&
 &                              CLOSE_IOREQUEST

IMPLICIT NONE

TYPE(TMP)           ,INTENT(IN)    :: YDMP
CHARACTER(LEN=*)    ,INTENT(IN)    :: CDFILE
TYPE(SPECTRAL_FIELD),INTENT(INOUT) :: YDSPEC

INTEGER(KIND=JPIM) :: ILEVS3D(YDSPEC%NFLEVG),IBSET3D(YDSPEC%NFLEVG),JLEV
INTEGER(KIND=JPIM) :: IOPROC
TYPE(TYPE_IOSTREAM) :: YL_IOSTREAM
TYPE(TYPE_IOREQUEST) :: YL_IOREQUEST

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

IF (LHOOK) CALL DR_HOOK('READ_SPEC_GRIB',0,ZHOOK_HANDLE)
ASSOCIATE(NBSETLEV=>YDMP%NBSETLEV)
WRITE(NULOUT,*)'READ_SPEC_GRIB: reading file ',CDFILE

IOPROC=1
DO JLEV=1,YDSPEC%NFLEVG
  ILEVS3D(JLEV)=JLEV
ENDDO
DO JLEV=1,YDSPEC%NFLEVG
  IBSET3D(JLEV) = NBSETLEV(ILEVS3D(JLEV))
ENDDO
CALL SETUP_IOSTREAM(YL_IOSTREAM,'CIO',TRIM(CDFILE),CDMODE='r',&
 & KIOMASTER=IOPROC)
CALL SETUP_IOREQUEST(YL_IOREQUEST,'SPECTRAL_FIELDS',LDGRIB=.TRUE.,CDLEVTYPE='ML',&
 & LDINTERP=.TRUE., &
 & KGRIB3D=YDSPEC%NGRIB3(1:YDSPEC%NS3D),KLEVS3D=ILEVS3D,KBSET3D=IBSET3D,&
 & KGRIB2D=YDSPEC%NGRIB2(1:YDSPEC%NS2G), &
 & KSMAX=YDSPEC%NSMAX)
CALL IO_GET(YL_IOSTREAM,YL_IOREQUEST,PR2=YDSPEC%SP2D,PR3=YDSPEC%SP3D)
CALL CLOSE_IOREQUEST(YL_IOREQUEST)
CALL CLOSE_IOSTREAM(YL_IOSTREAM)

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('READ_SPEC_GRIB',1,ZHOOK_HANDLE)
! ------------------------------------------------------------------
END SUBROUTINE READ_SPEC_GRIB
