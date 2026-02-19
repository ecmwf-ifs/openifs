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

SUBROUTINE PRINT_GFP(CDTEXT,CDNAME,YD)

!**** *PRINT_GFP*  - Print GFP information (from SUAFN3)

!     Author. 
!     ------- 
!      Philippe Marguinaud *METEO FRANCE* 
!      Original : 01-10-2014
!      R. El Khatib 08-Dec-2015 Interoperability GRIB2 vs FA
!      R. El Khatib 09-Sep-2016 better interoperability GRIB2 vs FA

USE YOMHOOK, ONLY : LHOOK, DR_HOOK, JPHOOK
USE PARKIND1, ONLY : JPRB
USE YOMLUN, ONLY : NULOUT
USE YOMCT0, ONLY : LARPEGEF
USE TYPE_FPDSPHYS, ONLY : FPDSPHY

IMPLICIT NONE

! This tiny subroutine enables any sequence of the type.

CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: CDTEXT
CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: CDNAME
TYPE(FPDSPHY),    OPTIONAL, INTENT(IN) :: YD

CHARACTER(LEN=12) :: CLNAME
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('PRINT_GFP',0,ZHOOK_HANDLE)
IF (LARPEGEF) THEN
  IF (PRESENT(CDTEXT).AND.PRESENT(CDNAME).AND.PRESENT(YD)) THEN
    CLNAME=CDNAME//' '
    WRITE(NULOUT,FMT='(1X,A25,'' : '',A12,1X,A16,3X,I2,5X,I2,3X,I4,5X,A8,  &
     & 1X,I1,4X,I2,6X,L1,6X,L1)')                                          &
     & CDTEXT,CLNAME, YD%CLNAME, YD%IBITS, YD%INTER, YD%IORDR,             &
     & YD%CLPAIR, YD%IMASK, YD%IANO, YD%LLMON, YD%LLSRF  
  ELSE
    WRITE(UNIT=NULOUT,FMT='(3X,''Description'',12X,'' : '','' TYPE NAME  '', &
     & 1X,''    %CLNAME     '',1X,''%IBITS'',1X,''%INTER'',1X,''%IORDR'',1X, &
     & ''%CLPAIR'',1X,''%IMASK'',1X,''%IANO'',1X,''%LLMON'',1X,''%LLSRF'')')  
  ENDIF
ELSE
  IF (PRESENT(CDTEXT).AND.PRESENT(CDNAME).AND.PRESENT(YD)) THEN
    IF (YD%IGRIB > 0) THEN
      CLNAME=CDNAME//' '
      WRITE(NULOUT,FMT='(1X,A25,'' : '',A12,I6,4X,I2,5X,I2,3X,I4,5X,A8, &
       & 1X,I1,4X,I2,6X,L1,6X,L1)')                                     &
       & CDTEXT,CLNAME, YD%IGRIB, YD%IBITS, YD%INTER, YD%IORDR,         &
       & YD%CLPAIR, YD%IMASK, YD%IANO, YD%LLMON, YD%LLSRF  
    ENDIF
  ELSE
    WRITE(UNIT=NULOUT,FMT='(3X,''Description'',12X,'' : '','' TYPE NAME  '', &
     & 1X,''%IGRIB'',1X,''%IBITS'',1X,''%INTER'',1X,''%IORDR'',1X,           &
     & ''%CLPAIR'',1X,''%IMASK'',1X,''%IANO'',1X,''%LLMON'',1X,''%LLSRF'')')  
  ENDIF
ENDIF

IF (LHOOK) CALL DR_HOOK('PRINT_GFP',1,ZHOOK_HANDLE)

END SUBROUTINE PRINT_GFP
