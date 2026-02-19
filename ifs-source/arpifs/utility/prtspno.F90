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

SUBROUTINE PRTSPNO(KULOUT,KFLEV,KLEV,KNSMAX,CDVAR,PX)

!**** *PRTSPNO* - Print spectrum

!     Purpose.
!     --------
!        Print spectrum

!**   Interface.
!     ----------
!        *CALL* *PRTSPNO(KULOUT,KFLEV,KLEV,KNSMAX,CDVAR,PX)

!        Explicit arguments :
!        --------------------
!           KULOUT : Logical unit for outputs.                   (input)
!           KFLEV  : Dimensioning                                (input)
!           KLEV   : Number of levels to print                   (input)
!           KNSMAX : Truncation                                  (input)
!           CDVAR  : Name of the variable                        (input)
!           PX     : Array containing the spectrums              (input)

!        Implicit arguments : none.
!        --------------------

!     Method.
!     -------

!     Externals. None
!     ----------
!      Called by SPNORM.

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Philippe Courtier  *ECMWF*

!     Modifications.
!     --------------
!        Original : 92-12-23
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KFLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KNSMAX 
INTEGER(KIND=JPIM),INTENT(IN)    :: KULOUT 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
CHARACTER(LEN=4)  ,INTENT(IN)    :: CDVAR 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PX(KFLEV,0:KNSMAX) 
INTEGER(KIND=JPIM) :: I1, I2, J1, JLEV, JN
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

!*       1.    Prints.
!              -------

IF (LHOOK) CALL DR_HOOK('PRTSPNO',0,ZHOOK_HANDLE)
I1=(KNSMAX+1)/10
I2=MOD(KNSMAX+1,10)
DO JLEV=1,KLEV
  DO J1=1,I1
    WRITE(KULOUT,'(1X,A4,I3,10(2X,E10.4))')&
     & CDVAR,JLEV,(PX(JLEV,(J1-1)*10+JN-1),JN=1,10)  
  ENDDO
  WRITE(KULOUT,'(1X,A4,I3,10(2X,E10.4))')&
   & CDVAR,JLEV,(PX(JLEV,I1*10+JN-1),JN=1,I2)  
ENDDO

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('PRTSPNO',1,ZHOOK_HANDLE)
END SUBROUTINE PRTSPNO

