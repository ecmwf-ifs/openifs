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

SUBROUTINE SCORDO(KORDR,PFREQR,PFREQI,PWHOU,PWORK)

!**** *SCORDO* - Order frequencies and corresponding normal modes.

!     Purpose.   To sort the frequencies and corresponding normal modes
!     --------   in ascending order.

!**   Interface.
!     ----------
!        *CALL* *SCORDO(KORDR,PFREQR,PFREQI,PWHOU,PWORK)*

!        Explicit arguments :  KORDR - number of modes (input)
!        --------------------  PFREQR - real part of frequencies (I-O)
!                              PFREQI - imaginary part of frequencies (I-O)
!                              PWHOU  - normal modes (I-O)
!                              PWORK  - work array

!        Implicit arguments :  NONE.
!        --------------------

!     Method.
!     -------

!     Externals.  none
!     ----------  

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS
!      P.Courtier and J.F. Geleyn: A global model with varaible
!      resolution. Application to a shallow-water model.
!      Qua.J.Roi.Met.Soc. 1988.

!     Author.
!     -------
!        Mats Hamrud and Philippe Courtier  *ECMWF*

!     Modifications.
!     --------------
!        Original : 88-05-25
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KORDR 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFREQR(KORDR) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFREQI(KORDR) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PWHOU(KORDR,KORDR) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PWORK(KORDR,KORDR) 
REAL(KIND=JPRB) :: ZONER(KORDR),ZONEI(KORDR)
INTEGER(KIND=JPIM) :: INDEX(KORDR)

INTEGER(KIND=JPIM) :: J1, J2
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "xutrii.intfb.h"

!     ------------------------------------------------------------------

!*       1.    SORT FREQUENCIES.
!              ----------------

IF (LHOOK) CALL DR_HOOK('SCORDO',0,ZHOOK_HANDLE)
CALL XUTRII(KORDR,PFREQR,INDEX)

!     ------------------------------------------------------------------

!*       2.    REORDER FREQUENCIES AND NORMAL MODES.
!              -------------------------------------

DO J1=1,KORDR
  ZONER(J1)=PFREQR(INDEX(J1))
  ZONEI(J1)=PFREQI(INDEX(J1))
ENDDO

DO J1=1,KORDR
  PFREQR(J1)=ZONER(J1)
  PFREQI(J1)=ZONEI(J1)
ENDDO

DO J1=1,KORDR
  DO J2=1,KORDR
    PWORK(J2,J1)= PWHOU(J2,INDEX(J1))
  ENDDO
ENDDO

DO J1=1,KORDR
  DO J2=1,KORDR
    PWHOU(J2,J1)=PWORK(J2,J1)
  ENDDO
ENDDO

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SCORDO',1,ZHOOK_HANDLE)
END SUBROUTINE SCORDO

