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

SUBROUTINE SUOFNAME(KSTEP,PTSTEP,KDIGITS,CDPATH,CDLABEL,CDOFNAME,CDINC)

!**** *SUOFNAME*  - Setup an output filename

!     Purpose.
!     --------
!           To compute an output filename

!     Explicit arguments :
!     --------------------

!           KSTEP    : (input, optional) time step. If missing, extension is 'INIT'
!           KDIGITS  : (input optional) number of digits for timestamp if timestep >= 0
!           CDPATH   : (input optional) path prefix (ending with a "/")
!           CDLABEL  : (input) file label
!           CDOFNAME : (output) filename (including timestamp)
!           CDINC    : (output optional) timestamp

!     Externals.
!     ----------

!     Author.
!     -------
!        Ryad EL KHATIB *METEO-FRANCE*
!        Original : 23-Oct-2008 from inifaout and others

!     Modifications.
!     --------------
!        P.Marguinaud: 26-Apr-2012 : Add minutes in filename
!        H.Varella   : 25-July-2012: Introduction of the LWRIBVEC_FULL key
!        P.Marguinaud: 10-Oct-2013 : More minutes in filename
!        K. Yessad (July 2014): Move some variables.
!        R. El Khatib  10-Dec-2015   KSTEP optional, CDPATH optional
!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCT0   , ONLY : CNMEXP, LINFLAT
USE YOMVAR   , ONLY : MBGVEC, LWRIBVEC_FULL
USE YOMLUN   , ONLY : NULOUT
USE QACTEX   , ONLY : LAEINC

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN), OPTIONAL  :: KSTEP
REAL(KIND=JPRB)   ,INTENT(IN), OPTIONAL  :: PTSTEP
INTEGER(KIND=JPIM),INTENT(IN), OPTIONAL  :: KDIGITS
CHARACTER(LEN=*)  ,INTENT(IN), OPTIONAL  :: CDPATH
CHARACTER(LEN=*)  ,INTENT(IN)            :: CDLABEL
CHARACTER(LEN=*)  ,INTENT(OUT)           :: CDOFNAME
CHARACTER(LEN=*)  ,INTENT(OUT), OPTIONAL :: CDINC

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: IINC, ILENINC, ILENOFNAME, IDIGITS, IST
CHARACTER(LEN=5)   :: CLSTEP
CHARACTER(LEN=16)  :: CLINC
CHARACTER(LEN=5)   :: CLMEMBER

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "get_clinc.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUOFNAME',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

!*           COMPUTE TIME STAMP & OUTPUT FILE NAME
!            -------------------------------------

IDIGITS=4
IF (PRESENT(KDIGITS)) THEN
  IF (KDIGITS > 6) THEN
    CALL ABOR1('SUOFNAME:KDIGITS COULD NOT EXCEED 6')
  ELSE
    IDIGITS=MAX(1,KDIGITS)
  ENDIF
ENDIF

CDOFNAME=' '
IST=6-IDIGITS+1

IF (CDLABEL=='ICMCM') THEN

! coupled fields for COM method
  IF (PRESENT(KSTEP)) THEN
    IF (KSTEP < 0) THEN
      ILENOFNAME=LEN(CDLABEL)+4+LEN(CLSTEP)
      IF (ILENOFNAME > LEN(CDOFNAME)) CALL ABOR1('SUOFNAME:CDOFNAME TOO SHORT')
      WRITE(CLSTEP,'(I5.4)') KSTEP
      CDOFNAME=CDLABEL//CNMEXP(1:4)//CLSTEP
      IF (PRESENT(CDINC)) THEN
        ILENINC=LEN(CLSTEP)
        IF (ILENINC > LEN(CDINC)) CALL ABOR1('SUOFNAME:CDINC TOO SHORT')
        CDINC=CLSTEP
      ENDIF
    ELSE
      ILENOFNAME=LEN(CDLABEL)+5+IDIGITS
      IF (ILENOFNAME > LEN(CDOFNAME)) CALL ABOR1('SUOFNAME:CDOFNAME TOO SHORT')
      WRITE(CLINC,'(I6.6)') IINC
      CDOFNAME=CDLABEL//CNMEXP(1:4)//'+'//CLINC(IST:6)
      IF (PRESENT(CDINC)) THEN
        IF (IDIGITS > LEN(CDINC)) CALL ABOR1('SUOFNAME:CDINC TOO SHORT')
        CDINC=CLINC(IST:6)
      ENDIF
    ENDIF
  ELSE
    ILENOFNAME=LEN(CDLABEL)+8
    IF (ILENOFNAME > LEN(CDOFNAME)) CALL ABOR1('SUOFNAME:CDOFNAME TOO SHORT')
    CDOFNAME=CDLABEL//CNMEXP(1:4)//'INIT'
    IF (PRESENT(CDINC)) THEN
      CDINC='INIT'
    ENDIF
  ENDIF

ELSE

! 'standard' ICMSH case
  IF (PRESENT(KSTEP)) THEN
    IF (KSTEP < 0) THEN
      WRITE(CLSTEP,'(I5.4)') KSTEP
      IF (PRESENT(CDPATH)) THEN
        ILENOFNAME=LEN_TRIM(CDPATH)+LEN(CDLABEL)+4+LEN(CLSTEP)
        IF (ILENOFNAME > LEN(CDOFNAME)) CALL ABOR1('SUOFNAME:CDOFNAME TOO SHORT')
        IF (PRESENT(CDPATH)) THEN
          CDOFNAME=TRIM(CDPATH)//CDLABEL//CNMEXP(1:4)//CLSTEP
        ELSE
          CDOFNAME=CDLABEL//CNMEXP(1:4)//CLSTEP
        ENDIF
      ELSE
        ILENOFNAME=LEN(CDLABEL)+4+LEN(CLSTEP)
        IF (ILENOFNAME > LEN(CDOFNAME)) CALL ABOR1('SUOFNAME:CDOFNAME TOO SHORT')
        CDOFNAME=CDLABEL//CNMEXP(1:4)//CLSTEP
      ENDIF
      IF (PRESENT(CDINC)) THEN
        ILENINC=LEN(CLSTEP)
        IF (ILENINC > LEN(CDINC)) CALL ABOR1('SUOFNAME:CDINC TOO SHORT')
        CDINC=CLSTEP
      ENDIF
    ELSE
      IF (PRESENT(PTSTEP)) THEN
        CALL GET_CLINC (CLINC,KSTEP,IDIGITS,PTSTEP)
      ELSE
        CALL GET_CLINC (CLINC,KSTEP,IDIGITS)
      ENDIF
      IF (PRESENT(CDPATH)) THEN
        CDOFNAME=TRIM(CDPATH)//CDLABEL//CNMEXP(1:4)//'+'//TRIM(CLINC)
      ELSE
        CDOFNAME=CDLABEL//CNMEXP(1:4)//'+'//TRIM(CLINC)
      ENDIF
      IF (PRESENT (CDINC)) THEN
        IF (LEN_TRIM (CLINC) > LEN (CDINC)) CALL ABOR1('SUOFNAME:CDINC TOO SHORT')
        CDINC = TRIM (CLINC)
      ENDIF
    ENDIF
  ELSE
    ILENOFNAME=LEN_TRIM(CDPATH)+LEN(CDLABEL)+8
    IF (ILENOFNAME > LEN(CDOFNAME)) CALL ABOR1('SUOFNAME:CDOFNAME TOO SHORT')
    IF (PRESENT(CDPATH)) THEN
      CDOFNAME=TRIM(CDPATH)//CDLABEL//CNMEXP(1:4)//'INIT'
    ELSE
      CDOFNAME=CDLABEL//CNMEXP(1:4)//'INIT'
    ENDIF
    IF (PRESENT(CDINC)) THEN
      CDINC='INIT'
    ENDIF
  ENDIF

  IF ((LINFLAT.AND.LAEINC) .OR. LWRIBVEC_FULL) THEN
    WRITE(NULOUT,*)MBGVEC
    WRITE(CLMEMBER(1:5),'(''_m'',I3.3)')MBGVEC
    WRITE(NULOUT,*)'CLMEMBER :',CLMEMBER
    WRITE(NULOUT,*)'CDOFNAME AVANT MODIF :',TRIM(CDOFNAME)
    CDOFNAME=TRIM(CDOFNAME)//CLMEMBER
    WRITE(NULOUT,*)'CDOFNAME APRES MODIF :',TRIM(CDOFNAME)
  ENDIF

ENDIF

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUOFNAME',1,ZHOOK_HANDLE)
END SUBROUTINE SUOFNAME
